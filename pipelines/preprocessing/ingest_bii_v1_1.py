#!/usr/bin/env python3
"""Copy Biodiversity Intactness v1.1 into GCS and register it as a public EE collection.

Each source object is streamed from source.coop to
``gs://landler-opendata-gee-assets/datasets/bii/v1_1/`` without local staging. The stream
recomputes the source's S3 ETag (the MD5 for single-part objects; the multipart MD5 of
per-part MD5s otherwise) together with its CRC32C and SHA-256. A matching ETag proves the
received bytes equal the source. GCS's server-side CRC32C must then equal the stream's
CRC32C. The source COGs are registered unchanged as
``projects/landler-open-data/assets/datasets/bii/v1_1/bii_<year>``.

All Google calls use the active ``gcloud`` identity (user login or service account).
Reruns are idempotent: verified objects and registrations are checked, not rewritten.
An interrupted upload leaves no object because of ``if_generation_match=0``.

Example::

    python pipelines/preprocessing/ingest_bii_v1_1.py run --work-dir DIR
    python pipelines/preprocessing/ingest_bii_v1_1.py verify --work-dir DIR
"""

from __future__ import annotations

import argparse
import base64
import csv
import fcntl
import hashlib
import http.client
import json
import math
import multiprocessing
import subprocess
import time
import urllib.error
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path

import ee
import google.auth.credentials
import google_crc32c
import numpy as np
import rasterio
from google.api_core.exceptions import NotFound
from google.cloud import storage
from googleapiclient.errors import HttpError
from rasterio.windows import Window

SOURCE_URL = "https://data.source.coop/vizzuality/biodiversity-intactness-100m-v1-1"
BUCKET = "landler-opendata-gee-assets"
BILLING_PROJECT = "landler-open-data"
EE_PROJECT = "landler-open-data"
PREFIX = "datasets/bii/v1_1"
COLLECTION = f"projects/{EE_PROJECT}/assets/{PREFIX}"
YEARS = tuple(range(2017, 2026))
DOCS = {
    "README.md": "text/markdown",
    "biodiversity-intactness-100m-white-paper-v1.1.pdf": "application/pdf",
}
BAND = "bii"
MIB = 1024 * 1024
READ_CHUNK = 8 * MIB
UPLOAD_CHUNK = 64 * MIB
SAMPLES_PER_YEAR = 20
ZERO_SAMPLES = 4
SAMPLE_OVERVIEW_LEVEL = 6
EE_UNMASK = -9999
# Cloudflare rejects urllib's default User-Agent with HTTP 403.
USER_AGENT = {
    "User-Agent": "eii-bii-ingest/1.0 (+https://github.com/landler-io/ecosystem-integrity-index)"
}

CREDITS = {
    "title": "Global 100m Projections of Biodiversity Intactness 2017-2025 (v1.1)",
    "citation": (
        "Francis Gassert, Joe Mazzarello, and Samantha Hyde. 2026. Global 100m Projections of "
        "Biodiversity Intactness for the years 2017-2025 (v1.1). "
        "https://source.coop/vizzuality/biodiversity-intactness-100m-v1-1"
    ),
    "producers": "Impact Observatory; Vizzuality",
    "license": "CC-BY-4.0",
    "license_url": "https://creativecommons.org/licenses/by/4.0/",
    "source_dataset_url": "https://source.coop/vizzuality/biodiversity-intactness-100m-v1-1",
    "methodology_url": f"{SOURCE_URL}/biodiversity-intactness-100m-white-paper-v1.1.pdf",
    "code_url": "https://github.com/Vizzuality/biodiversity-intactness-100m-update",
    "values": "0-1; 1 = most intact. Source declares no nodata; non-land pixels hold 0.",
    "dataset_version": "v1.1",
}


class GcloudCredentials(google.auth.credentials.Credentials):
    """Borrow short-lived tokens from whichever identity ``gcloud`` has active."""

    def refresh(self, request: object) -> None:  # noqa: ARG002
        self.token = subprocess.run(
            ["gcloud", "auth", "print-access-token"], check=True, capture_output=True, text=True
        ).stdout.strip()
        # google-auth compares naive UTC; gcloud tokens live 60 minutes.
        self.expiry = datetime.now(timezone.utc).replace(tzinfo=None) + timedelta(minutes=45)


def log(message: str) -> None:
    print(f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} {message}", flush=True)


def part_size(size: int, parts: int) -> int:
    """Infer the whole-MiB S3 multipart part size that yields ``parts`` parts."""
    for mib in range(1, 5121):
        if math.ceil(size / (mib * MIB)) == parts:
            return mib * MIB
    raise ValueError(f"No whole-MiB part size splits {size} bytes into {parts} parts")


class StreamHasher:
    """Recompute an S3 ETag plus CRC32C and SHA-256 over one sequential byte stream."""

    def __init__(self, size: int, etag: str) -> None:
        self.size, self.expected_etag = size, etag
        parts = etag.partition("-")[2]
        self.parts = int(parts) if parts else 0
        self.part_size = part_size(size, self.parts) if self.parts else size
        self.length = 0
        self._part, self._part_length, self._digests = hashlib.md5(), 0, []
        self._sha256, self._crc32c = hashlib.sha256(), google_crc32c.Checksum()

    def update(self, data: bytes) -> None:
        self.length += len(data)
        self._sha256.update(data)
        self._crc32c.update(data)
        view = memoryview(data)
        while view:
            take = view[: self.part_size - self._part_length]
            self._part.update(take)
            self._part_length += len(take)
            view = view[len(take) :]
            if self.parts and self._part_length == self.part_size:
                self._digests.append(self._part.digest())
                self._part, self._part_length = hashlib.md5(), 0

    def result(self) -> dict:
        if not self.parts:
            etag = self._part.hexdigest()
        else:
            digests = self._digests + ([self._part.digest()] if self._part_length else [])
            etag = f"{hashlib.md5(b''.join(digests)).hexdigest()}-{len(digests)}"
        return {
            "size": self.length,
            "etag": etag,
            "sha256": self._sha256.hexdigest(),
            "crc32c": base64.b64encode(self._crc32c.digest()).decode(),
        }

    def verified(self) -> dict:
        result = self.result()
        if result["size"] != self.size or result["etag"] != self.expected_etag:
            raise ValueError(f"Stream does not match source: {result} vs {self.expected_etag}")
        return result


def head(url: str) -> dict:
    with urllib.request.urlopen(
        urllib.request.Request(url, headers=USER_AGENT, method="HEAD"), timeout=60
    ) as r:
        return {
            "url": url,
            "size": int(r.headers["Content-Length"]),
            "etag": r.headers["ETag"].strip('"'),
            "last_modified": r.headers["Last-Modified"],
        }


def stream_source(source: dict, sink=None, max_failures: int = 30) -> dict:
    """Hash (and optionally forward) the source, resuming dropped reads with HTTP ranges."""
    hasher = StreamHasher(source["size"], source["etag"])
    response, failures, reported = None, 0, 0
    while hasher.length < source["size"]:
        try:
            if response is None:
                request = urllib.request.Request(
                    source["url"],
                    headers={
                        **USER_AGENT,
                        "Range": f"bytes={hasher.length}-",
                        "If-Match": f'"{source["etag"]}"',
                    },
                )
                response = urllib.request.urlopen(request, timeout=120)
                if response.status != 206 or response.headers["ETag"].strip('"') != source["etag"]:
                    raise ValueError(f"Unexpected range response {response.status}")
            chunk = response.read(READ_CHUNK)
            if not chunk:
                raise ConnectionError("Source closed the stream early")
        except urllib.error.HTTPError as error:
            if error.code < 500:
                raise
            response, failures = None, failures + 1
        except (OSError, http.client.HTTPException) as error:
            response, failures = None, failures + 1
            log(f"retry {failures} at byte {hasher.length} of {source['url']}: {error}")
        else:
            if sink is not None:
                sink(chunk)
            hasher.update(chunk)
            if hasher.length - reported >= 1024 * MIB * 4:
                reported = hasher.length
                log(f"{source['url']}: {hasher.length / source['size']:.1%}")
            continue
        if failures > max_failures:
            raise RuntimeError(f"Too many read failures: {source['url']}")
        time.sleep(min(60, 2**failures))
    if response is not None:
        response.close()
    return hasher.verified()


def object_record(source: dict, blob: storage.Blob, hashes: dict) -> dict:
    return {
        "source_key": source["url"].removeprefix(SOURCE_URL + "/"),
        "source_url": source["url"],
        "source_size": source["size"],
        "source_etag": source["etag"],
        "source_last_modified": source["last_modified"],
        "gcs_uri": f"gs://{blob.bucket.name}/{blob.name}",
        "gcs_generation": str(blob.generation),
        "gcs_crc32c": blob.crc32c,
        "gcs_size": blob.size,
        "sha256": hashes["sha256"],
    }


def copy_object(bucket: storage.Bucket, name: str, content_type: str, write: bool) -> dict:
    """Stream one source object into GCS once and prove the stored bytes match it."""
    source = head(f"{SOURCE_URL}/{name}")
    blob = bucket.blob(f"{PREFIX}/{'' if name.endswith('.tif') else 'docs/'}{name}")
    try:
        blob.reload()
    except NotFound:
        if not write:
            raise
        log(f"copying {source['url']} ({source['size']} bytes) to {blob.name}")
        blob.metadata = {"source_url": source["url"], "source_etag": source["etag"]}
        with blob.open(
            "wb",
            chunk_size=UPLOAD_CHUNK,
            ignore_flush=True,
            if_generation_match=0,
            content_type=content_type,
        ) as writer:
            hashes = stream_source(source, writer.write)
        blob.reload()
    else:
        metadata = blob.metadata or {}
        if metadata.get("source_etag") != source["etag"] or blob.size != source["size"]:
            raise ValueError(f"Existing object differs from source; refusing: {blob.name}")
        if metadata.get("sha256"):
            hashes = {"crc32c": metadata["crc32c"], "sha256": metadata["sha256"]}
        elif write:
            log(f"re-hashing source for unverified existing object {blob.name}")
            hashes = stream_source(source)
        else:
            raise ValueError(f"Object was never verified: {blob.name}")
    if blob.crc32c != hashes["crc32c"] or blob.size != source["size"]:
        raise ValueError(f"GCS bytes differ from the verified source stream: {blob.name}")
    if not (blob.metadata or {}).get("sha256"):
        blob.metadata = {
            **(blob.metadata or {}),
            "crc32c": hashes["crc32c"],
            "sha256": hashes["sha256"],
            "source_etag_verified": "true",
        }
        blob.patch(if_generation_match=blob.generation)
    return object_record(source, blob, hashes)


def inspect_cog(url: str) -> dict:
    """Read the grid and COG layout from the source, whose bytes equal the GCS copy."""
    with (
        rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_MAX_RETRY=5),
        rasterio.open(f"/vsicurl/{url}") as dataset,
    ):
        info = {
            "layout": dataset.tags(ns="IMAGE_STRUCTURE").get("LAYOUT"),
            "compression": dataset.tags(ns="IMAGE_STRUCTURE").get("COMPRESSION"),
            "block_shapes": [list(shape) for shape in dataset.block_shapes],
            "overviews": dataset.overviews(1),
            "count": dataset.count,
            "dtype": dataset.dtypes[0],
            "nodata": dataset.nodata,
            "crs": dataset.crs.to_string(),
            "transform": list(dataset.transform)[:6],
            "shape": [dataset.height, dataset.width],
        }
    if (
        info["layout"] != "COG"
        or info["count"] != 1
        or info["block_shapes"] != [[512, 512]]
        or not info["overviews"]
    ):
        raise ValueError(f"Source is not a single-band tiled COG with overviews: {info}")
    return info


def image_properties(year: int, obj: dict) -> dict:
    return {
        **CREDITS,
        "year": year,
        "band_units": "unitless intactness fraction",
        "source_url": obj["source_url"],
        "source_etag": obj["source_etag"],
        "source_size": obj["source_size"],
        "source_last_modified": obj["source_last_modified"],
        "gcs_uri": obj["gcs_uri"],
        "gcs_generation": obj["gcs_generation"],
        "gcs_crc32c": obj["gcs_crc32c"],
        "gcs_sha256": obj["sha256"],
    }


def collection_properties() -> dict:
    return {
        **CREDITS,
        "description": (
            f"{CREDITS['title']}. Annual images bii_<year>, band '{BAND}'. Unmodified source "
            f"COGs copied to gs://{BUCKET}/{PREFIX}/. {CREDITS['values']}"
        ),
        "years": f"{YEARS[0]}-{YEARS[-1]}",
        "gcs_prefix": f"gs://{BUCKET}/{PREFIX}/",
    }


def registration_manifest(year: int, obj: dict, grid: dict) -> dict:
    band = {"id": BAND, "tilesetId": "data", "tilesetBandIndex": 0}
    if grid["nodata"] is not None:
        band["missingData"] = {"values": [grid["nodata"]]}
    return {
        "name": f"{COLLECTION}/bii_{year}",
        "tilesets": [{"id": "data", "sources": [{"uris": [obj["gcs_uri"]]}]}],
        "bands": [band],
        "startTime": f"{year}-01-01T00:00:00Z",
        "endTime": f"{year + 1}-01-01T00:00:00Z",
        "properties": image_properties(year, obj),
    }


def registration_errors(asset: dict, manifest: dict, obj: dict, grid: dict) -> list[str]:
    """List every difference between a registered EE image and its manifest/source grid."""
    errors = []
    uris = [u for t in asset.get("tilesets", []) for s in t["sources"] for u in s["uris"]]
    if uris != [f"{obj['gcs_uri']}#{obj['gcs_generation']}"]:
        errors.append(f"uris {uris} do not pin generation {obj['gcs_generation']}")
    if [band["id"] for band in asset.get("bands", [])] != [BAND]:
        errors.append("band ids differ")
    for key, value in manifest["properties"].items():
        if asset.get("properties", {}).get(key) != value:
            errors.append(f"property {key} differs")
    for key in ("startTime", "endTime"):
        if asset.get(key) != manifest[key]:
            errors.append(f"{key} differs")
    for band in asset.get("bands", []):
        if band["dataType"].get("precision") != "FLOAT" or grid["dtype"] != "float32":
            errors.append("dtype differs")
        affine = band["grid"]["affineTransform"]
        transform = [
            affine.get(k, 0)
            for k in ("scaleX", "shearX", "translateX", "shearY", "scaleY", "translateY")
        ]
        if band["grid"].get("crsCode") != grid["crs"]:
            errors.append("crs differs")
        if not np.allclose(transform, grid["transform"], rtol=0, atol=1e-12):
            errors.append("transform differs")
        dimensions = band["grid"]["dimensions"]
        if [dimensions["height"], dimensions["width"]] != grid["shape"]:
            errors.append("dimensions differ")
    return errors


def asset_or_none(asset_id: str) -> dict | None:
    try:
        return ee.data.getAsset(asset_id)
    except ee.EEException as error:
        if isinstance(error.__context__, HttpError) and error.__context__.resp.status == 404:
            return None
        raise


def ensure_collection(write: bool) -> None:
    """Create missing parents, then set the collection's credit properties."""
    parts = COLLECTION.split("/")
    for length in range(4, len(parts) + 1):
        name = "/".join(parts[:length])
        kind = "IMAGE_COLLECTION" if length == len(parts) else "FOLDER"
        existing = asset_or_none(name)
        if existing is None and write:
            ee.data.createAsset({"type": kind}, name)
        elif existing is None or existing["type"] != kind:
            raise ValueError(f"Missing or wrong-type destination parent: {name}")
    if ee.data.getAsset(COLLECTION).get("properties", {}) != collection_properties():
        if not write:
            raise ValueError("Collection properties differ")
        ee.data.updateAsset(COLLECTION, {"properties": collection_properties()}, ["properties"])


def choose_pixels(
    overview: np.ndarray, shape: tuple[int, int], seed: int, count: int, zeros: int
) -> list[tuple[int, int]]:
    """Pick full-resolution pixels under mostly positive and some zero overview cells."""
    rng = np.random.default_rng(seed)
    scale = (shape[0] / overview.shape[0], shape[1] / overview.shape[1])
    chosen = []
    for cells, wanted in (
        (np.argwhere(overview > 0), count - zeros),
        (np.argwhere(overview == 0), zeros),
    ):
        for row, col in cells[rng.choice(len(cells), min(wanted, len(cells)), replace=False)]:
            chosen.append(
                (
                    min(int((row + rng.random()) * scale[0]), shape[0] - 1),
                    min(int((col + rng.random()) * scale[1]), shape[1] - 1),
                )
            )
    return sorted(chosen)


def sample_check(year: int, obj: dict, grid: dict) -> list[dict]:
    """Compare source pixels (HTTP range reads) with the registered EE image."""
    url = obj["source_url"]
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_MAX_RETRY=5):
        with rasterio.open(f"/vsicurl/{url}", overview_level=SAMPLE_OVERVIEW_LEVEL) as dataset:
            overview = dataset.read(1)
        pixels = choose_pixels(overview, tuple(grid["shape"]), year, SAMPLES_PER_YEAR, ZERO_SAMPLES)
        with rasterio.open(f"/vsicurl/{url}") as dataset:
            source = [float(dataset.read(1, window=Window(c, r, 1, 1))[0, 0]) for r, c in pixels]
    transform = rasterio.Affine(*grid["transform"])
    points = ee.FeatureCollection(
        [
            ee.Feature(
                ee.Geometry.Point(list(transform * (c + 0.5, r + 0.5)), grid["crs"]), {"i": i}
            )
            for i, (r, c) in enumerate(pixels)
        ]
    )
    sampled = (
        ee.Image(f"{COLLECTION}/bii_{year}")
        .unmask(EE_UNMASK, sameFootprint=False)
        .sampleRegions(collection=points, projection=ee.Projection(grid["crs"], grid["transform"]))
        .getInfo()["features"]
    )
    target = {f["properties"]["i"]: f["properties"][BAND] for f in sampled}
    rows = [
        {"row": r, "col": c, "source": source[i], "ee": target.get(i)}
        for i, (r, c) in enumerate(pixels)
    ]
    bad = [
        row
        for row in rows
        if row["ee"] is None or not math.isclose(row["source"], row["ee"], rel_tol=1e-7, abs_tol=0)
    ]
    if bad or len(rows) != SAMPLES_PER_YEAR:
        raise ValueError(f"bii_{year}: EE/source sample mismatch {bad or len(rows)}")
    return rows


def connect() -> storage.Bucket:
    credentials = GcloudCredentials()
    ee.Initialize(credentials=credentials, project=EE_PROJECT)
    client = storage.Client(project=BILLING_PROJECT, credentials=credentials)
    bucket = client.bucket(BUCKET, user_project=BILLING_PROJECT)
    bucket.reload()
    if (
        bucket.location.lower() != "us-central1"
        or bucket.storage_class != "STANDARD"
        or not bucket.requester_pays
    ):
        raise ValueError("Require Standard storage in us-central1 with Requester Pays")
    return bucket


def process_year(year: int, write: bool, work_dir: Path) -> dict:
    """Copy, inspect, register, and sample one year; safe to repeat."""
    bucket = connect()
    obj = copy_object(bucket, f"bii_{year}.tif", "image/tiff", write)
    grid = inspect_cog(obj["source_url"])
    manifest = registration_manifest(year, obj, grid)
    asset = asset_or_none(manifest["name"])
    if asset is None:
        if not write:
            raise ValueError(f"Not registered: {manifest['name']}")
        ee.data.startExternalImageIngestion(manifest, allow_overwrite=False)
        asset = ee.data.getAsset(manifest["name"])
    errors = registration_errors(asset, manifest, obj, grid)
    if errors:
        raise ValueError(f"{manifest['name']}: {errors}")
    record = {
        **obj,
        "year": year,
        "ee_id": manifest["name"],
        "cog": grid,
        "samples": sample_check(year, obj, grid),
    }
    (work_dir / "years").mkdir(exist_ok=True)
    (work_dir / "years" / f"bii_{year}.json").write_text(json.dumps(record, indent=2) + "\n")
    log(f"bii_{year} verified: gen {obj['gcs_generation']}, {SAMPLES_PER_YEAR} samples match")
    return record


def execute(work_dir: Path, years: list[int], write: bool) -> dict:
    bucket = connect()
    docs = [copy_object(bucket, name, kind, write) for name, kind in DOCS.items()]
    ensure_collection(write)
    records, failures = [], []
    with ProcessPoolExecutor(len(years), mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = {pool.submit(process_year, year, write, work_dir): year for year in years}
        for future in as_completed(futures):
            try:
                records.append(future.result())
            except Exception as error:  # Let the other copies finish; rerun resumes.
                failures.append(futures[future])
                log(f"bii_{futures[future]} FAILED: {error!r}")
    if failures:
        raise RuntimeError(f"Years failed; rerun to resume: {sorted(failures)}")
    acl = ee.data.getAssetAcl(COLLECTION)
    if not acl.get("all_users_can_read"):
        if not write:
            raise ValueError("Collection is not public")
        ee.data.setAssetAcl(COLLECTION, {**acl, "all_users_can_read": True})
    records.sort(key=lambda record: record["year"])
    inventory = {
        "verified_at": datetime.now(timezone.utc).isoformat(),
        "mode": "run" if write else "verify",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "software": {
            "earthengine": ee.__version__,
            "rasterio": rasterio.__version__,
            "gdal": rasterio.__gdal_version__,
        },
        "collection": COLLECTION,
        "collection_public": True,
        "collection_properties": ee.data.getAsset(COLLECTION).get("properties", {}),
        "docs": docs,
        "images": records,
    }
    (work_dir / "inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
    columns = ["source_key", "source_size", "source_etag", "gcs_uri", "gcs_generation"]
    columns += ["gcs_crc32c", "sha256", "ee_id"]
    with (work_dir / "inventory.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(docs + records)
    return inventory


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run", "verify"))
    parser.add_argument("--work-dir", type=Path, default=Path("bii_ingest"))
    parser.add_argument("--years", type=int, nargs="+", default=list(YEARS))
    args = parser.parse_args()
    if not set(args.years) <= set(YEARS):
        parser.error(f"years must be within {YEARS[0]}-{YEARS[-1]}")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    with (args.work_dir / ".lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        inventory = execute(args.work_dir, args.years, args.command == "run")
    log(f"{args.command} complete: {len(inventory['images'])} images in {COLLECTION}")


if __name__ == "__main__":
    main()
