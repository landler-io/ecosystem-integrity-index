#!/usr/bin/env python3
"""Publish v2 raster products as GCS COGs and register each as an Earth Engine collection.

Products (``PRODUCTS``), all float32 tiles on the CLMS NPP grid with nodata -9999:

- ``potential_npp``: ``<root>/tiles/<tile>.tif`` of ``predict_potential_npp_v2.py``
  (5 bands, gC m-2 yr-1); one image per tile.
- ``eii``: ``<root>/<year>/<tile>.tif`` of ``compute_eii_v2.py`` (9 bands); one image
  ``eii_<year>_<tile>`` per tile and reference year, with property ``reference_year`` and the
  calendar year as its time range.

Tiles without any valid pixel are skipped.

Stages:

- ``upload`` (application default credentials): rewrite each tile as a ZSTD COG with average
  overviews, upload it without overwriting to ``gs://<bucket>/<prefix>/<name>.tif``, record the
  generation, CRC32C, SHA-256, grid and ``SAMPLES`` pixel values per tile in the publish
  manifest (``--manifest``), and delete the local COG once its upload is verified.
- ``register`` (Earth Engine credentials with write access to the project): create the
  collection ``projects/<project>/assets/<prefix>`` (public read only with ``--public``), one
  image per tile that pins its GCS generation, then check every image's grid and the recorded
  pixel values via Earth Engine. Reruns skip registered images and check them again.

Both stages write to the Landler bucket and Earth Engine project (``BUCKET``, ``PROJECT``);
they record how the published collections were made. To publish elsewhere, change both.

Example:
    python pipelines/production/publish_v2.py upload potential_npp potential_npp_v2
    python pipelines/production/publish_v2.py upload eii eii_v2
    python pipelines/production/publish_v2.py register potential_npp --public
    python pipelines/production/publish_v2.py register eii --public
"""

from __future__ import annotations

import argparse
import base64
import contextlib
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

BUCKET = "landler-opendata-gee-assets"
PROJECT = "landler-open-data"
NODATA = -9999.0
SAMPLES = 3
POTENTIAL_NPP_PROPERTIES = {
    "title": "EII v2 potential NPP with prediction intervals",
    "units": "gC m-2 yr-1",
    "bands": "potential_npp, lower_68, upper_68, lower_90, upper_90",
    "model": "LightGBM Tweedie x global bias correction c = 1.0203",
    "training": "pristine reference sites plus gap fill: relaxed sites (HMI < 0.2, weight 0.5) "
    "where a 4 deg block x aridity class has < 40 pristine sites",
    "intervals": "asymmetric QRF-scaled, per-biome split conformal (norm_asym); lower >= 0",
    "predictors": "15: CHELSA V2.1 bio1/4/12/15, aridity, CHELSA NPP, MERIT/Geomorpho90m "
    "terrain, coast distance, EOT20 tidal range, WOA23 salinity, "
    "MERIT Hydro upstream area",
    "response": "CLMS DMP V2.0.2 2019-2025 (OLCI) mean at reference sites",
    "grid": "CLMS NPP v1_1 grid, EPSG:4326, 1/336 deg, 10 deg tiles",
    "producer": "The Landbanking Group (EII v2)",
    "code": "landler-io/ecosystem-integrity-index pipelines/training/fit_potential_npp_v2.py, "
    "pipelines/training/predict_potential_npp_v2.py, pipelines/production/publish_v2.py",
}
EII_PROPERTIES = {
    "title": "Ecosystem Integrity Index v2, annual",
    "bands": "eii, functional_integrity, structural_integrity, compositional_integrity "
    "(0-1); limiting_component (1 functional, 2 structural, 3 compositional); wetland_flag "
    "(0 none, 1 wetland, 2 peatland, 3 mangrove); with the change term fi_state, fi_change "
    "(0-1) and fi_change_sensors (0-2 sensors, 0 where functional_integrity = fi_state)",
    "fi_change": "v2.0: functional_integrity = fi_state x fi_change; fi_change from MODIS Terra "
    "and VIIRS EVI against each pixel's drift-corrected baseline; serving series (TV, TNc, Nc "
    "or none) and constants SHA-256 in each tile's GeoTIFF tags",
    "aggregation": "min_fuzzy_logic, masked where any component is missing",
    "functional": "observed CLMS DMP V2.0.2 NPP (3 years ending in reference_year, from 2019) "
    "vs v2 potential NPP recentred per biome; two-sided deviation tolerance with natural CV, "
    "deficit k 1.0, surplus k 1.5, surplus scores 0 at three times the potential; "
    "settings and factors in each tile's GeoTIFF tags",
    "structural": "GHM 2022 landscape quality (Hill et al. 2026): sqrt(w) times the "
    "exponential-kernel mean (5 km, cut at 26.5 km) of sqrt(w) opened by 300 m, with w the "
    "logistic habitat quality of HMI (midpoint 0.4, steepness 10)",
    "compositional": "Biodiversity Intactness v1.1 of reference_year, clamped to [0, 1]",
    "grid": "CLMS NPP v1_1 grid, EPSG:4326, 1/336 deg, 10 deg tiles",
    "producer": "The Landbanking Group (EII v2)",
    "code": "landler-io/ecosystem-integrity-index pipelines/production/compute_eii_v2.py, "
    "pipelines/production/publish_v2.py",
}
PRODUCTS = {
    "potential_npp": {
        "prefix": "products/eii/v2_0/intermediate/functional/potential_npp_olci_gapfill",
        "properties": POTENTIAL_NPP_PROPERTIES,
    },
    "eii": {
        "prefix": "products/eii/v2_0/release/global",
        "properties": EII_PROPERTIES,
    },
}


def product_tiles(product: str, root: Path) -> list[tuple[Path, str, dict]]:
    """(tile file, object name below the prefix, image properties) of a product's tiles."""
    if product == "potential_npp":
        return [(t, t.stem, {}) for t in sorted((root / "tiles").glob("*.tif"))]
    return [
        (t, f"{year.name}/{t.stem}", {"reference_year": int(year.name)} | fi_change_tags(t))
        for year in sorted(p for p in root.iterdir() if p.is_dir() and p.name.isdigit())
        for t in sorted(year.glob("*.tif"))
    ]


def fi_change_tags(tile: Path) -> dict:
    """Serving series and constants SHA-256 of the FI change term from the tile tags."""
    import rasterio

    with rasterio.open(tile) as src:
        tags = src.tags()
    keys = ("fi_change_series", "fi_change_k", "fi_change_constants_sha256")
    return {k: tags[k] for k in keys if k in tags}


def image_name(record: dict) -> str:
    """EE image ID below the collection: the tile, or ``eii_<year>_<tile>``."""
    year = record.get("properties", {}).get("reference_year")
    return record["tile"] if year is None else f"eii_{year}_{record['tile'].split('/')[1]}"


def hashes(path: Path) -> dict[str, str]:
    import google_crc32c

    sha, crc = hashlib.sha256(), google_crc32c.Checksum()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            sha.update(chunk)
            crc.update(chunk)
    return {"sha256": sha.hexdigest(), "crc32c": base64.b64encode(crc.digest()).decode()}


def to_cog(tile: Path, cog: Path, name: str, properties: dict) -> dict | None:
    """COG of one tile plus its grid and sampled pixels; None if the tile has no valid pixel."""
    import rasterio
    import rasterio.shutil

    with rasterio.open(tile) as src:
        valid = src.read(1) != src.nodata
        if not valid.any():
            return None
        rows, cols = np.nonzero(valid)
        pick = np.random.default_rng(len(tile.stem)).choice(len(rows), SAMPLES, replace=False)
        samples = []
        for i in pick:
            r, c = int(rows[i]), int(cols[i])
            lon, lat = src.xy(r, c)
            values = src.read(window=((r, r + 1), (c, c + 1)))[:, 0, 0]
            samples.append({"lon": lon, "lat": lat, "values": [float(v) for v in values]})
        grid = {
            "crs": src.crs.to_string(),
            "transform": list(src.transform)[:6],
            "shape": [src.height, src.width],
            "bands": list(src.descriptions),
        }
    if not cog.exists():
        part = cog.with_name(cog.name + ".part")
        with rasterio.open(tile) as src:
            rasterio.shutil.copy(
                src,
                part,
                driver="COG",
                compress="ZSTD",
                predictor="FLOATING_POINT",
                blocksize=512,
                overview_resampling="average",
                num_threads="ALL_CPUS",
            )
        part.rename(cog)
    return {
        "tile": name,
        "properties": properties,
        "grid": grid,
        "samples": samples,
        **hashes(cog),
    }


def upload(args: argparse.Namespace) -> None:
    from google.api_core.exceptions import PreconditionFailed
    from google.cloud import storage

    product = PRODUCTS[args.product]
    prefix = product["prefix"]
    cog_dir = args.root / "cog"
    records = []
    with ProcessPoolExecutor(args.workers) as pool:
        futures = []
        for tile, name, properties in product_tiles(args.product, args.root):
            cog = cog_dir / f"{name}.tif"
            cog.parent.mkdir(parents=True, exist_ok=True)
            futures.append(pool.submit(to_cog, tile, cog, name, properties))
        for future in as_completed(futures):
            if (record := future.result()) is not None:
                records.append(record)
    bucket = storage.Client(project=PROJECT).bucket(BUCKET, user_project=PROJECT)
    total = 0
    for i, record in enumerate(sorted(records, key=lambda r: r["tile"]), start=1):
        cog = cog_dir / f"{record['tile']}.tif"
        blob = bucket.blob(f"{prefix}/{record['tile']}.tif")
        with contextlib.suppress(PreconditionFailed):
            blob.upload_from_filename(cog, content_type="image/tiff", if_generation_match=0)
        blob.reload()
        if blob.crc32c != record["crc32c"]:
            raise RuntimeError(f"gs://{BUCKET}/{blob.name} exists with a different CRC32C")
        record |= {"gcs_uri": f"gs://{BUCKET}/{blob.name}", "gcs_generation": str(blob.generation)}
        total += blob.size
        cog.unlink()
        if i % 25 == 0:
            print(f"{i}/{len(records)} uploaded", flush=True)
    records.sort(key=lambda r: r["tile"])
    manifest = args.manifest
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(
        json.dumps(
            {
                "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "collection": f"projects/{PROJECT}/assets/{prefix}",
                "properties": product["properties"],
                "tiles": records,
            },
            indent=1,
        )
        + "\n"
    )
    print(f"{len(records)} tiles, {total / 1e9:.1f} GB -> gs://{BUCKET}/{prefix}/; {manifest}")


def register(public: bool, manifest_path: Path) -> None:
    import ee

    ee.Initialize(project=PROJECT)
    manifest = json.loads(manifest_path.read_text())
    collection = manifest["collection"]

    def exists(asset_id: str) -> dict | None:
        try:
            return ee.data.getAsset(asset_id)
        except ee.EEException:
            return None

    parts = collection.split("/")
    for length in range(4, len(parts) + 1):
        name = "/".join(parts[:length])
        kind = "IMAGE_COLLECTION" if length == len(parts) else "FOLDER"
        found = exists(name)
        if found is None:
            ee.data.createAsset({"type": kind}, name)
        elif found["type"] != kind:
            raise ValueError(f"{name} exists as {found['type']}, expected {kind}")
    ee.data.updateAsset(collection, {"properties": manifest["properties"]}, ["properties"])

    errors = []
    for n, record in enumerate(manifest["tiles"], start=1):
        image_id = f"{collection}/{image_name(record)}"
        year = record.get("properties", {}).get("reference_year")
        times = (
            {}
            if year is None
            else {
                "startTime": f"{year}-01-01T00:00:00Z",
                "endTime": f"{year + 1}-01-01T00:00:00Z",
            }
        )
        uri = f"{record['gcs_uri']}#{record['gcs_generation']}"
        if exists(image_id) is None:
            ee.data.startExternalImageIngestion(
                {
                    "name": image_id,
                    "tilesets": [{"id": "data", "sources": [{"uris": [uri]}]}],
                    "bands": [
                        {
                            "id": band,
                            "tilesetId": "data",
                            "tilesetBandIndex": b,
                            "missingData": {"values": [NODATA]},
                        }
                        for b, band in enumerate(record["grid"]["bands"])
                    ],
                    **times,
                    "properties": {
                        **record.get("properties", {}),
                        "tile": record["tile"].split("/")[-1],
                        "gcs_uri": record["gcs_uri"],
                        "gcs_generation": record["gcs_generation"],
                        "gcs_crc32c": record["crc32c"],
                        "gcs_sha256": record["sha256"],
                    },
                },
                allow_overwrite=False,
            )
        asset = ee.data.getAsset(image_id)
        uris = [u for t in asset.get("tilesets", []) for s in t["sources"] for u in s["uris"]]
        band = asset["bands"][0]["grid"]
        affine = band["affineTransform"]
        transform = [affine.get(k, 0) for k in ("scaleX", "shearX", "translateX")]
        transform += [affine.get(k, 0) for k in ("shearY", "scaleY", "translateY")]
        if uris != [uri]:
            errors.append(f"{record['tile']}: uris {uris}")
        if not np.allclose(transform, record["grid"]["transform"], rtol=0, atol=1e-12):
            errors.append(f"{record['tile']}: transform {transform}")
        if [band["dimensions"]["height"], band["dimensions"]["width"]] != record["grid"]["shape"]:
            errors.append(f"{record['tile']}: dimensions differ")
        points = ee.FeatureCollection(
            [
                ee.Feature(ee.Geometry.Point([s["lon"], s["lat"]]), {"i": i})
                for i, s in enumerate(record["samples"])
            ]
        )
        projection = ee.Projection(record["grid"]["crs"], record["grid"]["transform"])
        got = {
            f["properties"]["i"]: [f["properties"].get(b) for b in record["grid"]["bands"]]
            # Unmasked, so points with a nodata band are returned with NODATA, as recorded.
            for f in ee.Image(image_id)
            .unmask(NODATA, sameFootprint=False)
            .sampleRegions(collection=points, projection=projection, geometries=False)
            .getInfo()["features"]
        }
        for i, sample in enumerate(record["samples"]):
            values = got.get(i)
            if values is None or not np.allclose(values, sample["values"], rtol=1e-6):
                errors.append(f"{record['tile']}: sample {sample['values']} != {values}")
        if n % 25 == 0:
            print(f"{n}/{len(manifest['tiles'])} registered and checked", flush=True)
    acl = ee.data.getAssetAcl(collection)
    if public and not acl.get("all_users_can_read"):
        ee.data.setAssetAcl(collection, {**acl, "all_users_can_read": True})
    print(json.dumps({"images": len(manifest["tiles"]), "errors": errors}, indent=1))
    if errors:
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    u = sub.add_parser("upload")
    u.add_argument("product", choices=list(PRODUCTS))
    u.add_argument("root", type=Path)
    u.add_argument("--workers", type=int, default=8)
    r = sub.add_parser("register")
    r.add_argument("product", choices=list(PRODUCTS))
    r.add_argument("--public", action="store_true", help="Give all users read access")
    for p in (u, r):
        p.add_argument(
            "--manifest",
            type=Path,
            help="Publish manifest (default publish_manifest_<product>.json)",
        )
    args = parser.parse_args()
    args.manifest = args.manifest or Path(f"publish_manifest_{args.product}.json")
    if args.command == "upload":
        upload(args)
    else:
        register(args.public, args.manifest)


if __name__ == "__main__":
    main()
