#!/usr/bin/env python3
"""Build the v2 coastal-setting predictor rasters and register them as COG-backed EE images.

The potential-NPP model has four coastal predictors. Three come from local sources and are
built here, so training and production read the same values; the fourth (``log_upa_10km``)
comes from MERIT Hydro (``export_predictor_stack.py``).

- ``osm_coast_distance``: great-circle distance (m, UInt16) from the pixel centre to the
  nearest OSM coastline vertex, capped at ``COAST_KM``; land and sea alike.
- ``eot20_tidal_range``: 2 x (M2 + S2 + K1 + O1) EOT20 amplitude (m, float32), each at its
  nearest ocean cell within 60 km. Only pixels within ``COAST_KM`` of the coast hold a value.
- ``woa23_sss``: WOA23 1991-2020 annual surface salinity (psu, float32) at the nearest ocean
  cell within 60 km, same rule.

Pixels without a value are nodata (-9999). ``build_predictor_stack.py`` applies the model
rules: tidal range and salinity 0 where missing and tidal range capped at 16 m.

Grid: the 300 m EPSG:4326 grid that ``reproject(crs="EPSG:4326", scale=300)`` uses in Earth
Engine (pixel 300 / 111319.49 deg, aligned to (0, 0)), 512-pixel COG tiles with overviews.

Stages, from the repository root, with the sources of ``EOT20_URL``, ``WOA_URL`` and ``OSM_URL``
unpacked under ``data/external``::

    python -I pipelines/preprocessing/ingest_coastal_v2.py build
    python -I pipelines/preprocessing/ingest_coastal_v2.py register   # upload, register, verify

``register`` uploads with ``if_generation_match=0`` (no clobber), registers each COG as
``projects/landler-open-data/assets/datasets/<name>/v1_0``, makes it public and compares
sampled pixels of the local COG with Earth Engine. Reruns check existing objects and assets
and do not rewrite them. The inventory goes to ``data/ingest_records/coastal``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import ee
import google_crc32c
import numpy as np
import rasterio
import rasterio.shutil
from rasterio.windows import Window
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

import ingest_bii_v1_1 as bii  # noqa: E402

# Sources and nearest-ocean helpers.
EXTERNAL = REPO / "data" / "external"
EOT20_URL = "https://www.seanoe.org/data/00683/79489/data/85762.zip"
EOT20_CONSTITUENTS = ("M2", "S2", "K1", "O1")
WOA_URL = "https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/{0}/netcdf/decav91C0/0.25/woa23_decav91C0_{1}00_04.nc"
WOA_VARIABLES = {"sst": ("temperature", "t"), "sss": ("salinity", "s")}
OSM_URL = "https://osmdata.openstreetmap.de/download/coastlines-split-4326.zip"
COAST_KM = 50.0  # ocean values (tides, WOA) only for points this close to the coast
OCEAN_SEARCH_KM = 60.0  # nearest valid ocean cell must be this close

VERSION = "v1_0"
WORK = REPO / "data" / "coastal_ingest"
RECORD = REPO / "data" / "ingest_records" / "coastal"
ASSET_ROOT = f"projects/{bii.EE_PROJECT}/assets/datasets"
NODATA = -9999.0

PIXEL = 300 / 111319.49079327357  # EE's degrees for scale=300 in EPSG:4326
COLS, ROWS = 2 * 66792, 2 * 33395  # covers -180..180, stays inside +-90
TRANSFORM = rasterio.Affine(PIXEL, 0, -COLS / 2 * PIXEL, 0, -PIXEL, ROWS / 2 * PIXEL)
STRIP = 512
CELL = 64  # prefilter cell (pixels): skip cells whose centre has no coast within reach
CELL_REACH_KM = COAST_KM + CELL * 0.3 * math.sqrt(2) / 2  # + half-diagonal at the equator

SAMPLES = 24
OSM_CREDIT = "OpenStreetMap coastlines (c) OpenStreetMap contributors, ODbL 1.0"
COMMON = {
    "producer": "The Landbanking Group (EII v2), derived product",
    "author": "Benjamin Leutner",
    "code": "landler-io/ecosystem-integrity-index pipelines/preprocessing/ingest_coastal_v2.py",
    "grid": "EPSG:4326, 300 m (EE scale=300 grid, pixel 0.0026949459 deg, origin 0,0)",
    "dataset_version": VERSION,
}
DATASETS = {
    "osm_coast_distance": {
        "band": "dist_coast_m",
        "dtype": "uint16",
        "properties": {
            "title": "Distance to the OpenStreetMap coastline, capped at 50 km (EII v2)",
            "values": "great-circle metres from the pixel centre to the nearest OSM coastline "
            "vertex (lines segmentized to 0.005 deg); 50000 = 50 km or more; land and sea",
            "units": "m",
            "source": f"{OSM_URL} (server build 2026-10-27)",
            "attribution": OSM_CREDIT,
            "license": "ODbL-1.0",
            "license_url": "https://opendatacommons.org/licenses/odbl/1-0/",
            "share_alike": "Derived database under ODbL 1.0: redistribute under ODbL with "
            "attribution to OpenStreetMap contributors (https://www.openstreetmap.org/copyright)",
        },
    },
    "eot20_tidal_range": {
        "band": "tidal_range_m",
        "dtype": "float32",
        "properties": {
            "title": "EOT20 approximate spring tidal range near the coast (EII v2)",
            "values": "2 x (M2 + S2 + K1 + O1 amplitude) in m, each constituent at its "
            "nearest EOT20 ocean cell within 60 km; only pixels within 50 km of the OSM "
            "coastline, nodata elsewhere; not capped (EII caps at 16 m)",
            "units": "m",
            "source": EOT20_URL,
            "citation": "Hart-Davis, M. G., Piccioni, G., Dettmering, D., Schwatke, C., "
            "Passaro, M., Seitz, F. (2021). EOT20: a global ocean tide model from "
            "multi-mission satellite altimetry. ESSD 13, 3869-3884. "
            "doi:10.5194/essd-13-3869-2021; data doi:10.17882/79489",
            "license": "CC-BY-4.0",
            "license_url": "https://creativecommons.org/licenses/by/4.0/",
            "attribution": f"EOT20 (Hart-Davis et al. 2021), CC BY 4.0; coastal mask from "
            f"{OSM_CREDIT}",
        },
    },
    "woa23_sss": {
        "band": "sss_psu",
        "dtype": "float32",
        "properties": {
            "title": "WOA23 1991-2020 annual sea surface salinity near the coast (EII v2)",
            "values": "WOA23 decav91C0 0.25 deg annual objectively analysed surface s_an at "
            "the nearest ocean cell within 60 km; only pixels within 50 km of the OSM "
            "coastline, nodata elsewhere",
            "units": "psu",
            "source": WOA_URL.format("salinity", "s"),
            "citation": "Reagan, J. R. et al. (2024). World Ocean Atlas 2023, Volume 2: "
            "Salinity. NOAA Atlas NESDIS 90",
            "license": "public domain (NOAA NCEI)",
            "attribution": f"NOAA NCEI World Ocean Atlas 2023; coastal mask from {OSM_CREDIT}",
        },
    },
}


def folder(name: str) -> Path:
    return EXTERNAL / name


def unit_xyz(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lon, lat = np.radians(lon), np.radians(lat)
    return np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def chord_to_km(chord: np.ndarray) -> np.ndarray:
    return 2 * 6371.0088 * np.arcsin(np.clip(chord / 2, 0, 1))


def km_to_chord(km: float) -> float:
    return 2 * math.sin(km / 6371.0088 / 2)


def coast_tree() -> cKDTree:
    """KD-tree of OSM coastline vertices (segmentized to 0.005 deg) on the unit sphere."""
    import pyogrio
    import shapely

    path = folder("osm_coastlines") / "coastlines-split-4326.zip"
    inner = next(n for n in zipfile.ZipFile(path).namelist() if n.endswith(".shp"))
    geometry = pyogrio.read_dataframe(f"/vsizip/{path}/{inner}", columns=[]).geometry
    xy = shapely.get_coordinates(shapely.segmentize(geometry.to_numpy(), 0.005))
    return cKDTree(unit_xyz(xy[:, 0], xy[:, 1]))


def ocean_tree(grid: np.ndarray, transform: rasterio.Affine) -> tuple[cKDTree, np.ndarray]:
    """KD-tree of the finite cells of an ocean grid and their values."""
    rows, cols = np.nonzero(np.isfinite(grid))
    x, y = rasterio.transform.xy(transform, rows, cols)
    x = (np.asarray(x) + 180) % 360 - 180
    return cKDTree(unit_xyz(x, np.asarray(y))), grid[rows, cols]


def nearest_ocean_value(ocean: tuple[cKDTree, np.ndarray], xyz: np.ndarray) -> np.ndarray:
    """Value of the nearest ocean cell within ``OCEAN_SEARCH_KM`` (NaN otherwise)."""
    tree, values = ocean
    chord, index = tree.query(xyz, distance_upper_bound=km_to_chord(OCEAN_SEARCH_KM), workers=-1)
    found = np.isfinite(chord)
    out = np.full(len(xyz), np.nan)
    out[found] = values[index[found]]
    return out


def read_grid(path: str) -> tuple[np.ndarray, rasterio.Affine]:
    with rasterio.open(path) as ds:
        grid = ds.read(1, masked=True).astype("float64").filled(np.nan)
        grid = grid * (ds.scales[0] or 1) + (ds.offsets[0] or 0)
        return grid, ds.transform


def eot20_grids() -> list[tuple[np.ndarray, rasterio.Affine]]:
    """EOT20 amplitude grids (cm) of ``EOT20_CONSTITUENTS``."""
    outer = zipfile.ZipFile(folder("eot20") / "EOT20_85762.zip")
    archive = zipfile.ZipFile(outer.open("ocean_tides.zip"))
    grids = []
    for constituent in EOT20_CONSTITUENTS:
        member = f"ocean_tides/{constituent}_ocean_eot20.nc"
        target = folder("eot20") / member
        if not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(archive.read(member))
        grids.append(read_grid(f"netcdf:{target}:amplitude"))
    return grids


def woa_path(name: str) -> str:
    kind, letter = WOA_VARIABLES[name]
    path = folder("woa23") / WOA_URL.format(kind, letter).rsplit("/", 1)[1]
    return f"netcdf:{path}:{letter}_an"


def path(name: str) -> Path:
    return WORK / f"{name}_{VERSION}.tif"


def pixel_centres(rows: np.ndarray, cols: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lon = TRANSFORM.c + (cols + 0.5) * PIXEL
    return (lon + 180) % 360 - 180, TRANSFORM.f - (rows + 0.5) * PIXEL


def coastal_values(
    xyz: np.ndarray, coast: cKDTree, tides: list, salinity: tuple
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Distance (m, capped), tidal range (m) and salinity (psu) at unit vectors ``xyz``."""
    chord, _ = coast.query(xyz, distance_upper_bound=km_to_chord(COAST_KM), workers=-1)
    near = np.isfinite(chord)
    dist = np.where(near, chord_to_km(np.where(near, chord, 0)), COAST_KM)
    tidal = np.full(len(xyz), NODATA, dtype="float32")
    sss = np.full(len(xyz), NODATA, dtype="float32")
    if near.any():
        total = sum(nearest_ocean_value(ocean, xyz[near]) for ocean in tides)
        tidal[near] = np.nan_to_num(2 * total / 100, nan=NODATA)  # cm -> m
        sss[near] = np.nan_to_num(nearest_ocean_value(salinity, xyz[near]), nan=NODATA)
    return np.rint(dist * 1000).astype("uint16"), tidal, sss


def strip_values(row0: int, coast: cKDTree, tides: list, salinity: tuple) -> list[np.ndarray]:
    """The three rasters for rows ``row0`` to ``row0 + STRIP`` of the global grid."""
    height = min(STRIP, ROWS - row0)
    cell_rows = np.arange(math.ceil(height / CELL))
    cell_cols = np.arange(math.ceil(COLS / CELL))
    lon, lat = pixel_centres(
        row0 + (cell_rows[:, None] + 0.5) * CELL - 0.5, (cell_cols + 0.5) * CELL - 0.5
    )
    lon, lat = np.broadcast_arrays(lon, lat)
    chord, _ = coast.query(
        unit_xyz(lon.ravel(), lat.ravel()),
        distance_upper_bound=km_to_chord(CELL_REACH_KM),
        workers=-1,
    )
    cells = np.isfinite(chord).reshape(lon.shape)
    todo = cells[np.arange(height)[:, None] // CELL, np.arange(COLS)[None, :] // CELL]
    out = [
        np.full((height, COLS), round(COAST_KM * 1000), dtype="uint16"),
        np.full((height, COLS), NODATA, dtype="float32"),
        np.full((height, COLS), NODATA, dtype="float32"),
    ]
    rows, cols = np.nonzero(todo)
    if len(rows):
        values = coastal_values(unit_xyz(*pixel_centres(row0 + rows, cols)), coast, tides, salinity)
        for grid, value in zip(out, values, strict=True):
            grid[rows, cols] = value
    return out


def run_build() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    bii.log("building trees")
    coast = coast_tree()
    tides = [ocean_tree(grid, transform) for grid, transform in eot20_grids()]
    salinity = ocean_tree(*read_grid(woa_path("sss")))
    bii.log(
        f"coast vertices {coast.n:,}; EOT20 cells {tides[0][0].n:,}; WOA cells {salinity[0].n:,}"
    )
    profile = {
        "driver": "GTiff",
        "width": COLS,
        "height": ROWS,
        "count": 1,
        "crs": "EPSG:4326",
        "transform": TRANSFORM,
        "tiled": True,
        "blockxsize": 512,
        "blockysize": 512,
        "compress": "deflate",
        "bigtiff": "yes",
    }
    staging = {name: WORK / f"{name}.staging.tif" for name in DATASETS}
    sinks = {
        name: rasterio.open(
            staging[name],
            "w",
            **profile,
            dtype=spec["dtype"],
            nodata=None if spec["dtype"] == "uint16" else NODATA,
        )
        for name, spec in DATASETS.items()
    }
    for row0 in range(0, ROWS, STRIP):
        values = strip_values(row0, coast, tides, salinity)
        for sink, value in zip(sinks.values(), values, strict=True):
            sink.write(value, 1, window=Window(0, row0, COLS, len(value)))
        if row0 // STRIP % 10 == 0:
            bii.log(f"rows {row0:,}/{ROWS:,}")
    for name, sink in sinks.items():
        sink.close()
        rasterio.shutil.copy(
            staging[name],
            path(name),
            driver="COG",
            compress="DEFLATE",
            blocksize=512,
            predictor="2" if DATASETS[name]["dtype"] == "uint16" else "3",
            overview_resampling="average",
            bigtiff="YES",
            num_threads="ALL_CPUS",
        )
        staging[name].unlink()
        bii.log(f"{path(name)}: {path(name).stat().st_size:,} bytes")


def file_hashes(file: Path) -> dict:
    sha, crc = hashlib.sha256(), google_crc32c.Checksum()
    with file.open("rb") as stream:
        while chunk := stream.read(bii.READ_CHUNK):
            sha.update(chunk)
            crc.update(chunk)
    return {"sha256": sha.hexdigest(), "crc32c": bii.base64.b64encode(crc.digest()).decode()}


def upload(bucket: bii.storage.Bucket, name: str, hashes: dict) -> bii.storage.Blob:
    blob = bucket.blob(f"datasets/{name}/{VERSION}/{path(name).name}")
    try:
        blob.reload()
    except bii.NotFound:
        bii.log(f"uploading {path(name)} to gs://{bucket.name}/{blob.name}")
        blob.metadata = hashes
        blob.upload_from_filename(path(name), content_type="image/tiff", if_generation_match=0)
        blob.reload()
    if blob.crc32c != hashes["crc32c"]:
        raise ValueError(f"GCS object differs from the local COG: {blob.name}")
    return blob


def manifest(name: str, blob: bii.storage.Blob, hashes: dict) -> dict:
    spec = DATASETS[name]
    band = {"id": spec["band"], "tilesetId": "data", "tilesetBandIndex": 0}
    if spec["dtype"] != "uint16":
        band["missingData"] = {"values": [NODATA]}
    return {
        "name": f"{ASSET_ROOT}/{name}/{VERSION}",
        "tilesets": [
            {"id": "data", "sources": [{"uris": [f"gs://{blob.bucket.name}/{blob.name}"]}]}
        ],
        "bands": [band],
        "properties": {
            **COMMON,
            **spec["properties"],
            "gcs_uri": f"gs://{blob.bucket.name}/{blob.name}",
            "gcs_generation": str(blob.generation),
            "gcs_crc32c": blob.crc32c,
            "gcs_sha256": hashes["sha256"],
        },
    }


def sample_check(name: str, asset_id: str) -> list[dict]:
    """Compare local COG pixels (half near the coast) with EE ``sampleRegions``."""
    band = DATASETS[name]["band"]
    rng = np.random.default_rng(0)
    with rasterio.open(path("osm_coast_distance")) as dist:
        overview = dist.read(1, out_shape=(ROWS // 256, COLS // 256))
    coastal = np.argwhere(overview < COAST_KM * 1000)
    anywhere = np.argwhere(overview >= 0)
    cells = [*coastal[rng.choice(len(coastal), SAMPLES // 2, replace=False)]]
    cells += [*anywhere[rng.choice(len(anywhere), SAMPLES // 2, replace=False)]]
    pixels = [(int((r + rng.random()) * 256), int((c + rng.random()) * 256)) for r, c in cells]
    with rasterio.open(path(name)) as dataset:
        local = [float(dataset.read(1, window=Window(c, r, 1, 1))[0, 0]) for r, c in pixels]
    points = ee.FeatureCollection(
        [
            ee.Feature(ee.Geometry.Point(list(TRANSFORM * (c + 0.5, r + 0.5))), {"i": i})
            for i, (r, c) in enumerate(pixels)
        ]
    )
    sampled = (
        ee.Image(asset_id)
        .unmask(NODATA, sameFootprint=False)
        .sampleRegions(
            collection=points, projection=ee.Projection("EPSG:4326", list(TRANSFORM)[:6])
        )
        .getInfo()["features"]
    )
    remote = {f["properties"]["i"]: f["properties"][band] for f in sampled}
    rows = [
        {"row": r, "col": c, "local": local[i], "ee": remote.get(i)}
        for i, (r, c) in enumerate(pixels)
    ]
    bad = [
        row
        for row in rows
        if row["ee"] is None or not math.isclose(row["local"], row["ee"], rel_tol=1e-6)
    ]
    if bad:
        raise ValueError(f"{name}: EE/local sample mismatch {bad}")
    return rows


def ensure_folder(name: str) -> None:
    folder = f"{ASSET_ROOT}/{name}"
    existing = bii.asset_or_none(folder)
    if existing is None:
        ee.data.createAsset({"type": "FOLDER"}, folder)
    elif existing["type"] != "FOLDER":
        raise ValueError(f"Not a folder: {folder}")


def register(bucket: bii.storage.Bucket, name: str) -> dict:
    hashes = file_hashes(path(name))
    blob = upload(bucket, name, hashes)
    wanted = manifest(name, blob, hashes)
    asset = bii.asset_or_none(wanted["name"])
    if asset is None:
        ensure_folder(name)
        ee.data.startExternalImageIngestion(wanted, allow_overwrite=False)
        asset = ee.data.getAsset(wanted["name"])
    uris = [u for t in asset["tilesets"] for s in t["sources"] for u in s["uris"]]
    errors = [] if uris == [f"{wanted['properties']['gcs_uri']}#{blob.generation}"] else ["uris"]
    errors += [k for k, v in wanted["properties"].items() if asset["properties"].get(k) != v]
    if errors:
        raise ValueError(f"{wanted['name']} differs from its manifest: {errors}")
    acl = ee.data.getAssetAcl(wanted["name"])
    if not acl.get("all_users_can_read"):
        ee.data.setAssetAcl(wanted["name"], {**acl, "all_users_can_read": True})
    with rasterio.open(path(name)) as dataset:
        grid = {
            "dtype": dataset.dtypes[0],
            "nodata": dataset.nodata,
            "shape": [dataset.height, dataset.width],
            "transform": list(dataset.transform)[:6],
            "overviews": dataset.overviews(1),
            "layout": dataset.tags(ns="IMAGE_STRUCTURE").get("LAYOUT"),
        }
    record = {
        "ee_id": wanted["name"],
        "gcs_uri": wanted["properties"]["gcs_uri"],
        "gcs_generation": str(blob.generation),
        "size": blob.size,
        **hashes,
        "cog": grid,
        "properties": wanted["properties"],
        "samples": sample_check(name, wanted["name"]),
    }
    bii.log(f"{name} verified: {wanted['name']} gen {blob.generation}, {SAMPLES} samples match")
    return record


def run_register() -> None:
    bucket = bii.connect()
    records = [register(bucket, name) for name in DATASETS]
    RECORD.mkdir(parents=True, exist_ok=True)
    inventory = {
        "verified_at": datetime.now(timezone.utc).isoformat(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "software": {
            "earthengine": ee.__version__,
            "rasterio": rasterio.__version__,
            "gdal": rasterio.__gdal_version__,
        },
        "images": records,
    }
    (RECORD / "inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "register"))
    args = parser.parse_args()
    run_build() if args.command == "build" else run_register()


if __name__ == "__main__":
    main()
