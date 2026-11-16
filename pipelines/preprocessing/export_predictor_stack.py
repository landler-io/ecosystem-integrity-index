#!/usr/bin/env python3
"""Export the v2 potential-NPP terrain and upstream-area predictors from Earth Engine as COG tiles.

Each layer is exported on its own grid and resampled locally by ``build_predictor_stack.py``:

- ``terrain``: MERIT DEM elevation and Geomorpho90m slope, cti, tpi, northness, eastness,
  area-mean from native 3 arcsec straight onto the CLMS NPP v1_1 grid (1/336 deg).
- ``upa``: MERIT Hydro upstream area, max of the 10 x 10 native pixels (1/120 deg). The
  focal max and log of ``log_upa_10km`` are applied later.

Climate, CHELSA NPP and the coastal bands come from the bucket (``datasets/*/v1_0``).

``run`` writes to the Landler bucket from the Landler Earth Engine project
(``eii.config.GEE_PROJECT``) and records how the published export was made; its output is at
``gs://landler-opendata-gee-assets/products/eii/v2_0/predictor_stack/ee_export/``, and a
reproduction starts from ``sync_predictor_sources.py``. It is the global export: terrain per
land tile, upa as one task over 80 N to 60 S. It keeps ``--concurrency`` tasks running, cancels
a task past its layer cap, stops submitting past ``--budget-eecu-hours`` and is safe to
restart: its state is the project's EE operation list (description prefix
``predictor_stack_v2_``). ``check`` summarises a downloaded tile. ``land-tiles`` lists the
10 deg tiles between 80 N and 60 S that contain land (``export_predictor_stack_land_tiles.csv``,
also the tile list of ``build_predictor_stack.py``).

Example (Earth Engine credentials of the project):
    python pipelines/preprocessing/export_predictor_stack.py run
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path

import ee

from eii.config import GEE_PROJECT

BUCKET = "landler-opendata-gee-assets"
RUN_PREFIX = "products/eii/v2_0/predictor_stack/ee_export"
RUN_DESCRIPTION = "predictor_stack_v2_"
# Per-task EECU-h caps, about 3x the measured use (terrain <= 0.28 per tile; upa ~4 globally).
LAYER_CAP_EECU_H = {"terrain": 1.0, "upa": 12.0}
GLOBAL_LAYERS = ("upa",)
ACTIVE = {"PENDING", "RUNNING", "CANCELLING"}
MAX_ATTEMPTS = 3
LAND_TILES = Path(__file__).with_name("export_predictor_stack_land_tiles.csv")
NODATA = -9999
GEOMORPHO = "projects/sat-io/open-datasets/Geomorpho90m"
# (pixel size, left edge, top edge) of the target grids.
# CLMS NPP v1_1 (projects/landler-open-data/assets/datasets/clms_npp/v1_1/annual): pixel
# centres on multiples of 1/336 deg.
CLMS_GRID = (1 / 336, -180 - 1 / 672, 80 + 1 / 672)
# 1/120 deg blocks of exactly 10 x 10 MERIT Hydro pixels (MERIT centres on multiples of 1/1200).
MERIT_1KM_GRID = (1 / 120, -180 - 1 / 2400, 85 - 1 / 2400)


def native_mean(image: ee.Image, projection: ee.Projection | None = None) -> ee.Image:
    """Area-weighted mean of native pixels into whatever grid the export requests."""
    if projection is not None:
        image = image.setDefaultProjection(projection)
    return image.reduceResolution(ee.Reducer.mean(), maxPixels=64)


def geomorpho(name: str) -> ee.Image:
    collection = ee.ImageCollection(f"{GEOMORPHO}/{name}")
    return native_mean(collection.mosaic(), collection.first().projection()).rename(name)


def terrain() -> ee.Image:
    elevation = native_mean(ee.Image("MERIT/DEM/v1_0_3").select("dem")).rename("elevation")
    bands = [geomorpho(n) for n in ("slope", "cti", "tpi", "northness", "eastness")]
    return ee.Image.cat([elevation, *bands]).unmask(NODATA).toFloat()


def upa() -> ee.Image:
    image = ee.Image("MERIT/Hydro/v1_0_1").select("upa")
    return image.reduceResolution(ee.Reducer.max(), maxPixels=256).unmask(NODATA).toFloat()


LAYERS = {
    "terrain": (terrain, CLMS_GRID),
    "upa": (upa, MERIT_1KM_GRID),
}


def tile_grid(
    grid: tuple[float, float, float],
    lon: float,
    lat: float,
    width_deg: float,
    height_deg: float | None = None,
) -> tuple[list[float], str]:
    """crsTransform and dimensions from the grid edge nearest the tile corner (ties: west/north)."""
    pixel, x0, y0 = grid
    height_deg = width_deg if height_deg is None else height_deg
    col = math.floor((lon - x0) / pixel + 0.5 - 1e-6)
    row = math.floor((y0 - lat - height_deg) / pixel + 0.5 - 1e-6)
    width, height = round(width_deg / pixel), round(height_deg / pixel)
    return [pixel, 0, x0 + col * pixel, 0, -pixel, y0 - row * pixel], f"{width}x{height}"


def tile_name(lon: float, lat: float) -> str:
    return f"{'e' if lon >= 0 else 'w'}{abs(lon):03.0f}{'n' if lat >= 0 else 's'}{abs(lat):02.0f}"


def start(
    layer: str,
    name: str,
    lon: float,
    lat: float,
    width_deg: float,
    height_deg: float,
    prefix: str,
    description: str,
) -> tuple[str, list[float], str]:
    """Start one export task; returns its operation name, crsTransform and dimensions."""
    build, grid = LAYERS[layer]
    transform, dimensions = tile_grid(grid, lon, lat, width_deg, height_deg)
    task = ee.batch.Export.image.toCloudStorage(
        image=build(),
        description=description,
        bucket=BUCKET,
        fileNamePrefix=f"{prefix}/{name}",
        crs="EPSG:4326",
        crsTransform=transform,
        # With crs/crsTransform and no region, EE exports exactly this grid.
        dimensions=dimensions,
        maxPixels=1e13,
        fileFormat="GeoTIFF",
        formatOptions={"cloudOptimized": True, "noData": NODATA},
    )
    task.start()
    operation = task.operation_name or f"projects/{GEE_PROJECT}/operations/{task.id}"
    return operation, transform, dimensions


def run_jobs() -> list[tuple[str, str, float, float, float, float]]:
    """(layer, name, lon, lat, width, height): global coarse layers first, then land tiles."""
    jobs = [(layer, f"{layer}_global", -180.0, -60.0, 360.0, 140.0) for layer in GLOBAL_LAYERS]
    with LAND_TILES.open() as f:
        for row in csv.DictReader(f):
            lon, lat = float(row["lon"]), float(row["lat"])
            jobs.append(("terrain", f"terrain_{row['tile']}", lon, lat, 10.0, 10.0))
    return jobs


def run_operations() -> dict[str, list[dict]]:
    """This run's EE operations per description, oldest first."""
    by_description: dict[str, list[dict]] = {}
    for op in ee.data.listOperations():
        description = op.get("metadata", {}).get("description", "")
        if description.startswith(RUN_DESCRIPTION):
            by_description.setdefault(description, []).append(op)
    for ops in by_description.values():
        ops.sort(key=lambda op: op["metadata"].get("createTime", ""))
    return by_description


def run(args: argparse.Namespace) -> None:
    """Global export loop; restart-safe because all state is in the EE operation list."""
    jobs = run_jobs()
    while True:
        ops = run_operations()
        spent = sum(
            op["metadata"].get("batchEecuUsageSeconds", 0) / 3600
            for hist in ops.values()
            for op in hist
        )
        states: dict[str, str] = {}
        active = 0
        for layer, name, *_ in jobs:
            history = ops.get(RUN_DESCRIPTION + name, [])
            latest = history[-1] if history else None
            state = latest["metadata"].get("state", "PENDING") if latest else "TODO"
            if latest and state in ACTIVE:
                active += 1
                eecu_h = latest["metadata"].get("batchEecuUsageSeconds", 0) / 3600
                if state != "CANCELLING" and eecu_h > LAYER_CAP_EECU_H[layer]:
                    ee.data.cancelOperation(latest["name"])
                    print(json.dumps({"cancel": name, "eecu_h": round(eecu_h, 3)}), flush=True)
            elif state == "FAILED" and len(history) < MAX_ATTEMPTS:
                state = "RETRY"
            states[name] = state
        counts: dict[str, int] = {}
        for state in states.values():
            counts[state] = counts.get(state, 0) + 1
        print(
            json.dumps(
                {"t": time.strftime("%Y-%m-%d %H:%M:%S"), "eecu_h": round(spent, 2), **counts}
            ),
            flush=True,
        )
        over_budget = spent >= args.budget_eecu_hours
        pending = [job for job in jobs if states[job[1]] in ("TODO", "RETRY")]
        if not pending or over_budget:
            if active == 0:
                failed = [n for n, s in states.items() if s not in ("SUCCEEDED", "TODO", "RETRY")]
                print(
                    json.dumps({"done": not pending, "over_budget": over_budget, "not_ok": failed}),
                    flush=True,
                )
                return
        else:
            for layer, name, lon, lat, width, height in pending[
                : max(0, args.concurrency - active)
            ]:
                operation, _, dims = start(
                    layer, name, lon, lat, width, height, args.prefix, RUN_DESCRIPTION + name
                )
                print(json.dumps({"start": name, "operation": operation, "dims": dims}), flush=True)
        time.sleep(args.poll_seconds)


def check(args: argparse.Namespace) -> None:
    """Per-band ranges, nodata fraction, and pixels valid in all bands."""
    import numpy as np
    import rasterio

    with rasterio.open(args.path) as src:
        print(
            json.dumps(
                {
                    "crs": str(src.crs),
                    "transform": list(src.transform)[:6],
                    "shape": [src.count, src.height, src.width],
                    "dtype": src.dtypes[0],
                    "bands": src.descriptions,
                    "nodata": src.nodata,
                }
            )
        )
        all_valid = np.ones((src.height, src.width), bool)
        for i, band in enumerate(src.descriptions, start=1):
            data = src.read(i)
            valid = data != src.nodata
            all_valid &= valid
            v = data[valid] if valid.any() else np.array([np.nan])
            print(
                f"{band:20s} nodata={1 - valid.mean():.3f} "
                f"min={v.min():.4g} max={v.max():.4g} mean={v.mean():.4g}"
            )
        print(json.dumps({"pixels": int(all_valid.size), "all_valid": int(all_valid.sum())}))


def land_tiles(args: argparse.Namespace) -> None:
    """Tiles (80 N to 60 S) touching land, with land pixels on the CLMS grid estimated."""
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    from shapely.geometry import box

    land = pd.concat([gpd.read_file(p).to_crs("EPSG:4326") for p in args.land], ignore_index=True)
    # Land pixels: rasterise at 1/48 deg (7 x 7 CLMS pixels) by pixel centre.
    step = 7 * CLMS_GRID[0]
    width, height = 360 * 48, 140 * 48
    mask = rasterize(
        ((g, 1) for g in land.geometry),
        out_shape=(height, width),
        transform=from_origin(-180, 80, step, step),
        dtype="uint8",
    )
    per_tile = mask.reshape(14, 480, 36, 480).sum(axis=(1, 3)) * 49
    rows = []
    for i, lat in enumerate(range(70, -70, -10)):
        for j, lon in enumerate(range(-180, 180, 10)):
            if len(land.sindex.query(box(lon, lat, lon + 10, lat + 10), "intersects")):
                rows.append((tile_name(lon, lat), lon, lat, int(per_tile[i, j])))
    table = pd.DataFrame(rows, columns=["tile", "lon", "lat", "land_pixels_clms_est"])
    table.to_csv(LAND_TILES, index=False)
    print(json.dumps({"tiles": len(table), "land_pixels": int(table.iloc[:, 3].sum())}))
    print(f"wrote {LAND_TILES}; total grid cells with land: {int(np.count_nonzero(mask))}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run")
    r.add_argument("--prefix", default=RUN_PREFIX)
    r.add_argument("--concurrency", type=int, default=10)
    r.add_argument("--budget-eecu-hours", type=float, default=100)
    r.add_argument("--poll-seconds", type=int, default=60)
    c = sub.add_parser("check")
    c.add_argument("path")
    t = sub.add_parser("land-tiles")
    t.add_argument(
        "land",
        nargs="+",
        help="Land polygons, e.g. zip:// paths to ne_10m_land.zip and ne_10m_minor_islands.zip "
        "from https://naciscdn.org/naturalearth/10m/physical/",
    )
    args = parser.parse_args()
    if args.command == "run":
        ee.Initialize(project=GEE_PROJECT)
        run(args)
    elif args.command == "check":
        check(args)
    else:
        land_tiles(args)


if __name__ == "__main__":
    main()
