#!/usr/bin/env python3
"""Ingest CHELSA V2.1 1981-2010 climatologies for the v2 potential-NPP predictors.

CHELSA provides the four bioclim predictors (bio1, bio4, bio12, bio15) and, with
``pet_penman_mean``, aridity (P / PET, computed in ``build_predictor_stack.py``). The sixth
CHELSA predictor, ``chelsa_npp`` (``CHELSA_npp_1981-2010_V.2.1.tif`` at ``BASE_URL``), is not
ingested here: it was uploaded to Earth Engine before this pipeline and exported to
``datasets/chelsa_npp/v1_0/`` in the bucket, where ``sync_predictor_sources.py`` reads it.
Licence: CC0 1.0 (EnviDat, CHELSA climatologies
doi:10.16904/envidat.228 and CHELSA-BIOCLIM+ doi:10.16904/envidat.332). Cite Karger et al.
(2017) Sci. Data 4, 170122 and Brun et al. (2022) Earth Syst. Sci. Data 14, 5573-5603.

Each source GeoTIFF (30 arcsec, uint16 with scale/offset tags, deflate stripes) is downloaded,
checked against the server's Content-Length, rewritten as a ZSTD COG with the same values,
dtype, grid and scale/offset, and uploaded without overwriting to
``gs://<bucket>/datasets/chelsa_<var>/v1_0/chelsa_<var>_1981_2010_v2-1.tif``. The COG is also
kept under ``<root>/chelsa_<var>/``, the layout ``sync_predictor_sources.py`` mirrors. A
manifest with source URL, checksums and GCS generation goes to ``data/ingest_records/chelsa``.

Notes: bio12 (kg m-2 yr-1, scale 0.1) saturates at 65535 = 6553.5 mm in the wettest places.
The bioclim layers have values over the ocean; pet_penman_mean has nodata 65535 there.

Example:
    python pipelines/preprocessing/ingest_chelsa_v2_1.py predictor_stack \\
        --download-dir sources/chelsa_v2_1
"""

from __future__ import annotations

import argparse
import base64
import contextlib
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import google_crc32c
import rasterio
import rasterio.shutil
from google.api_core.exceptions import PreconditionFailed
from google.cloud import storage

REPO = Path(__file__).resolve().parents[2]
RECORD = REPO / "data" / "ingest_records" / "chelsa"
BUCKET = "landler-opendata-gee-assets"
BILLING_PROJECT = "landler-open-data"
BASE_URL = "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio"
VARIABLES = ("bio1", "bio4", "bio12", "bio15", "pet_penman_mean")
VERSION = "v1_0"


def hashes(path: Path) -> dict[str, str]:
    sha, crc = hashlib.sha256(), google_crc32c.Checksum()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            sha.update(chunk)
            crc.update(chunk)
    return {"sha256": sha.hexdigest(), "crc32c": base64.b64encode(crc.digest()).decode()}


def download(url: str, target: Path) -> None:
    with urllib.request.urlopen(urllib.request.Request(url, method="HEAD")) as r:
        size = int(r.headers["Content-Length"])
    if target.exists() and target.stat().st_size == size:
        return
    part = target.with_name(target.name + ".part")
    urllib.request.urlretrieve(url, part)
    if part.stat().st_size != size:
        raise OSError(f"{url}: got {part.stat().st_size} bytes, expected {size}")
    part.rename(target)


def to_cog(source: Path, target: Path) -> None:
    if target.exists():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + ".part")
    rasterio.shutil.copy(
        source,
        part,
        driver="COG",
        compress="ZSTD",
        predictor="YES",
        blocksize=512,
        overview_resampling="average",
        bigtiff="IF_SAFER",
        num_threads="ALL_CPUS",
    )
    with rasterio.open(source) as a, rasterio.open(part) as b:
        assert (a.dtypes, a.shape, a.transform, a.scales, a.offsets, a.nodata) == (
            b.dtypes,
            b.shape,
            b.transform,
            b.scales,
            b.offsets,
            b.nodata,
        ), source
    part.rename(target)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="Local mirror root (predictor_stack)")
    parser.add_argument("--download-dir", type=Path, required=True)
    args = parser.parse_args()

    args.download_dir.mkdir(parents=True, exist_ok=True)
    bucket = storage.Client(project=BILLING_PROJECT).bucket(BUCKET, user_project=BILLING_PROJECT)
    records = []
    for var in VARIABLES:
        source_name = f"CHELSA_{var}_1981-2010_V.2.1.tif"
        url = f"{BASE_URL}/{source_name}"
        name = f"chelsa_{var}"
        source = args.download_dir / source_name
        cog = args.root / name / f"{name}_1981_2010_v2-1.tif"
        download(url, source)
        to_cog(source, cog)
        cog_hashes = hashes(cog)
        blob = bucket.blob(f"datasets/{name}/{VERSION}/{cog.name}")
        with contextlib.suppress(PreconditionFailed):
            blob.upload_from_filename(cog, content_type="image/tiff", if_generation_match=0)
        blob.reload()
        if blob.crc32c != cog_hashes["crc32c"]:
            raise RuntimeError(f"gs://{BUCKET}/{blob.name} exists with a different CRC32C")
        with rasterio.open(cog) as src:
            raster = {
                "dtype": src.dtypes[0],
                "nodata": src.nodata,
                "scale": src.scales[0],
                "offset": src.offsets[0],
                "shape": [src.height, src.width],
                "transform": list(src.transform)[:6],
            }
        records.append(
            {
                "variable": var,
                "source_url": url,
                "source_bytes": source.stat().st_size,
                "source_sha256": hashes(source)["sha256"],
                "gcs_uri": f"gs://{BUCKET}/{blob.name}",
                "gcs_generation": str(blob.generation),
                "cog_bytes": cog.stat().st_size,
                **{f"cog_{k}": v for k, v in cog_hashes.items()},
                **raster,
            }
        )
        print(json.dumps(records[-1]), flush=True)

    RECORD.mkdir(parents=True, exist_ok=True)
    (RECORD / "manifest.json").write_text(
        json.dumps(
            {
                "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "licence": "CC0 1.0",
                "dois": ["10.16904/envidat.228", "10.16904/envidat.332"],
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "layers": records,
            },
            indent=1,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
