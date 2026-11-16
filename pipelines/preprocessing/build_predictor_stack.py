#!/usr/bin/env python3
"""Build the v2 potential-NPP predictor stack on the CLMS NPP grid, one 10 deg land tile per file.

Inputs come from ``sync_predictor_sources.py`` (``<root>/ee_export``, ``chelsa_npp`` and the
coastal COGs). Each output tile is a float32 GeoTIFF with one band per model-card predictor on
the CLMS NPP v1_1 grid (1/336 deg), nodata -9999, valid only where the terrain bands are (land):

- terrain (elevation, slope, cti, tpi, northness, eastness): the EE export is already on the
  grid (area mean of the native 90 m pixels);
- CHELSA V2.1 1981-2010 bio1/4/12/15, aridity = bio12 / (12 x pet_penman_mean) and CHELSA NPP
  (30 arcsec, scale/offset applied, e.g. bio1 in deg C) and coast distance, tidal range,
  salinity (EE 300 m grid): bilinear;
- the coastal rules of the v2 model card: dist_coast in km, tidal_range and sss 0 where
  missing, tidal_range capped at 16 m, and log_upa_10km =
  log10(1 + 10 km focal max of the 1 km max upstream area), nearest onto the grid, 0 where
  missing or 50 km or more from the coast.

CHELSA bio12 saturates at 6553.5 mm (uint16, scale 0.1) in the wettest places; aridity there
is low by the same factor. Gaps in the 30 arcsec layers (e.g. coast pixels where PET has no
value) are filled from the nearest valid source pixel within ``FILL_PX``. Sources:
``ingest_chelsa_v2_1.py`` (CHELSA bioclim and PET, CC0 1.0), CHELSA V2.1 NPP 1981-2010
(``chelsa_npp``, no ingest script; see ``sync_predictor_sources.py``), ``ingest_coastal_v2.py``
and ``export_predictor_stack.py``. The RESOLVE Ecoregions 2017 raster synced alongside is not a
predictor; it gives the biome of the interval calibration in ``predict_potential_npp_v2.py``.
Existing tiles are skipped, so the script can be re-run as more terrain tiles arrive. ``qa.csv`` gets one row per built tile.

Example:
    python pipelines/preprocessing/build_predictor_stack.py predictor_stack
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine
from rasterio.warp import reproject
from rasterio.windows import Window, from_bounds
from scipy.ndimage import distance_transform_edt, maximum_filter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from export_predictor_stack import CLMS_GRID, LAND_TILES, tile_grid  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from eii.compute.settings import COASTAL_SETTING_KM, MAX_TIDAL_RANGE_M  # noqa: E402
from eii.training.settings import PREDICTOR_VARIABLES  # noqa: E402

NODATA = -9999.0
CRS = "EPSG:4326"
TERRAIN = ["elevation", "slope", "cti", "tpi", "northness", "eastness"]
CHELSA_BIOCLIM = {
    "mean_annual_temp": "bio1",
    "temp_seasonality": "bio4",
    "annual_precip": "bio12",
    "precip_seasonality": "bio15",
}
CHELSA_VARIABLES = [*CHELSA_BIOCLIM.values(), "pet_penman_mean"]
# EE focalMax(10000, "circle", "meters") on a scale=1000 EPSG:4326 grid: 10 km at the equator.
UPA_RADIUS_DEG = 10000 / 111319.49079327357
MARGIN_PX = 16
FILL_PX = 8  # 30 arcsec pixels (~8 km at the equator)
PROFILE = {
    "driver": "GTiff",
    "dtype": "float32",
    "nodata": NODATA,
    "crs": CRS,
    "tiled": True,
    "blockxsize": 512,
    "blockysize": 512,
    "compress": "zstd",
    "predictor": 3,
    "bigtiff": "if_safer",
}


def mosaic_vrt(files: list[Path], out: Path) -> Path:
    """GDAL VRT mosaic of same-grid, same-band GeoTIFF pieces (no copying)."""
    metas = []
    for path in files:
        with rasterio.open(path) as src:
            metas.append((path, src.transform, src.width, src.height, src.count, src))
            dtype, nodata, descriptions = src.dtypes[0], src.nodata, src.descriptions
    pixel = metas[0][1].a
    left = min(t.c for _, t, *_ in metas)
    top = max(t.f for _, t, *_ in metas)
    right = max(t.c + w * pixel for _, t, w, *_ in metas)
    bottom = min(t.f - h * pixel for _, t, _, h, *_ in metas)
    width, height = round((right - left) / pixel), round((top - bottom) / pixel)
    gdal_type = {"float32": "Float32", "uint16": "UInt16"}[dtype]
    bands = []
    for b in range(1, metas[0][4] + 1):
        sources = []
        for path, t, w, h, *_ in metas:
            x, y = round((t.c - left) / pixel), round((top - t.f) / pixel)
            sources.append(
                f'<SimpleSource><SourceFilename relativeToVRT="0">{escape(str(path))}'
                f"</SourceFilename><SourceBand>{b}</SourceBand>"
                f'<SrcRect xOff="0" yOff="0" xSize="{w}" ySize="{h}"/>'
                f'<DstRect xOff="{x}" yOff="{y}" xSize="{w}" ySize="{h}"/></SimpleSource>'
            )
        nd = f"<NoDataValue>{nodata}</NoDataValue>" if nodata is not None else ""
        desc = f"<Description>{escape(descriptions[b - 1] or '')}</Description>"
        bands.append(
            f'<VRTRasterBand dataType="{gdal_type}" band="{b}">{desc}{nd}{"".join(sources)}'
            "</VRTRasterBand>"
        )
    out.write_text(
        f'<VRTDataset rasterXSize="{width}" rasterYSize="{height}">'
        f"<SRS>{CRS}</SRS><GeoTransform>{left}, {pixel}, 0, {top}, 0, {-pixel}</GeoTransform>"
        f"{''.join(bands)}</VRTDataset>\n"
    )
    return out


def sources(root: Path) -> dict[str, Path]:
    vrt_dir = root / "vrt"
    vrt_dir.mkdir(exist_ok=True)
    ee = root / "ee_export"

    def one(pattern: str, base: Path = ee) -> list[Path]:
        files = sorted(base.glob(pattern))
        if not files:
            raise FileNotFoundError(f"{base}/{pattern}")
        return files

    return {
        **{var: one("*.tif", root / f"chelsa_{var}")[0] for var in CHELSA_VARIABLES},
        "upa": mosaic_vrt(one("upa_global*.tif"), vrt_dir / "upa.vrt"),
        "chelsa_npp": mosaic_vrt(one("*/*.tif", root / "chelsa_npp"), vrt_dir / "chelsa.vrt"),
        "dist_coast": one("*.tif", root / "osm_coast_distance")[0],
        "tidal_range": one("*.tif", root / "eot20_tidal_range")[0],
        "sss": one("*.tif", root / "woa23_sss")[0],
    }


def read_window(
    path: Path, bounds: tuple[float, float, float, float], window: Window | None = None
) -> tuple[np.ndarray, Affine, Window]:
    """Bands over ``bounds`` plus ``MARGIN_PX``, float32, nodata -> NaN, scale/offset applied."""
    with rasterio.open(path) as src:
        if window is None:
            window = from_bounds(*bounds, transform=src.transform)
            window = Window(
                int(np.floor(window.col_off)) - MARGIN_PX,
                int(np.floor(window.row_off)) - MARGIN_PX,
                int(np.ceil(window.width)) + 2 * MARGIN_PX,
                int(np.ceil(window.height)) + 2 * MARGIN_PX,
            )
        fill = src.nodata if src.nodata is not None else 0
        data = src.read(window=window, boundless=True, fill_value=fill).astype("float32")
        if src.nodata is not None:
            data[data == src.nodata] = np.nan
        scales = np.asarray(src.scales, dtype="float32")[:, None, None]
        offsets = np.asarray(src.offsets, dtype="float32")[:, None, None]
        return data * scales + offsets, src.window_transform(window), window


def warp(
    data: np.ndarray,
    transform: Affine,
    dst_transform: Affine,
    size: int,
    resampling: Resampling = Resampling.bilinear,
) -> np.ndarray:
    out = np.full((data.shape[0], size, size), np.nan, dtype="float32")
    reproject(
        data,
        out,
        src_transform=transform,
        src_crs=CRS,
        src_nodata=np.nan,
        dst_transform=dst_transform,
        dst_crs=CRS,
        dst_nodata=np.nan,
        resampling=resampling,
    )
    return out


def fill_nearest(data: np.ndarray, max_px: float = FILL_PX) -> np.ndarray:
    """Per band, NaN pixels take the nearest valid pixel within ``max_px``."""
    out = data.copy()
    for band in out:
        missing = ~np.isfinite(band)
        if missing.any() and not missing.all():
            dist, (iy, ix) = distance_transform_edt(missing, return_indices=True)
            near = missing & (dist <= max_px)
            band[near] = band[iy[near], ix[near]]
    return out


def disk(radius: float) -> np.ndarray:
    r = int(np.floor(radius))
    y, x = np.mgrid[-r : r + 1, -r : r + 1]
    return x**2 + y**2 <= radius**2


def build_tile(tile: str, lon: float, lat: float, root: Path, src: dict[str, Path]) -> dict:
    out = root / "stack" / f"{tile}.tif"
    with rasterio.open(root / "ee_export" / f"terrain_{tile}.tif") as t:
        terrain = t.read().astype("float32")
        terrain[terrain == t.nodata] = np.nan
        dst_transform, size = t.transform, t.width
    expected, _ = tile_grid(CLMS_GRID, lon, lat, 10)
    assert np.allclose(list(dst_transform)[:6], expected, atol=1e-9), (tile, dst_transform)
    bounds = rasterio.transform.array_bounds(size, size, dst_transform)  # west, south, east, north
    land = np.isfinite(terrain).all(axis=0)
    bands: dict[str, np.ndarray] = dict(zip(TERRAIN, terrain, strict=True))

    chelsa = {}
    for var in CHELSA_VARIABLES:  # one 30 arcsec grid, so one window and transform
        data, ctransform, _ = read_window(src[var], bounds)
        chelsa[var] = data[0]
    for name, var in CHELSA_BIOCLIM.items():
        bands[name] = warp(fill_nearest(chelsa[var][None]), ctransform, dst_transform, size)[0]
    aridity = chelsa["bio12"] / (12 * chelsa["pet_penman_mean"])  # mm yr-1 / mean mm month-1
    bands["aridity"] = warp(fill_nearest(aridity[None]), ctransform, dst_transform, size)[0]
    chelsa, chtransform, _ = read_window(src["chelsa_npp"], bounds)
    bands["chelsa_npp"] = warp(fill_nearest(chelsa), chtransform, dst_transform, size)[0]

    dist, dtransform, _ = read_window(src["dist_coast"], bounds)
    bands["dist_coast"] = warp(dist, dtransform, dst_transform, size)[0] / 1000
    tidal, ttransform, _ = read_window(src["tidal_range"], bounds)
    tidal = np.nan_to_num(warp(tidal, ttransform, dst_transform, size)[0], nan=0.0)
    bands["tidal_range"] = np.minimum(tidal, MAX_TIDAL_RANGE_M)
    sss, stransform, _ = read_window(src["sss"], bounds)
    bands["sss"] = np.nan_to_num(warp(sss, stransform, dst_transform, size)[0], nan=0.0)

    upa, utransform, _ = read_window(src["upa"], bounds)
    focal = maximum_filter(
        np.where(np.isfinite(upa[0]), upa[0], -np.inf),
        footprint=disk(UPA_RADIUS_DEG / utransform.a),
        mode="constant",
        cval=-np.inf,
    )
    focal = np.where(np.isfinite(focal), focal, np.nan)[None]
    upa_10km = warp(focal, utransform, dst_transform, size, Resampling.nearest)[0]
    log_upa = np.nan_to_num(np.log10(1 + upa_10km), nan=0.0)
    log_upa[bands["dist_coast"] >= COASTAL_SETTING_KM] = 0.0
    bands["log_upa_10km"] = log_upa

    stack = np.stack([bands[name] for name in PREDICTOR_VARIABLES])
    stack[:, ~land] = np.nan
    qa = {"tile": tile, "land_px": int(land.sum())}
    for name, band in zip(PREDICTOR_VARIABLES, stack, strict=True):
        qa[f"missing_{name}"] = int((land & ~np.isfinite(band)).sum())
    qa["missing_any"] = int((land & ~np.isfinite(stack).all(axis=0)).sum())

    part = out.with_name(out.name + ".part")
    with rasterio.open(
        part, "w", width=size, height=size, count=len(stack), transform=dst_transform, **PROFILE
    ) as dst:
        dst.write(np.where(np.isfinite(stack), stack, NODATA).astype("float32"))
        dst.descriptions = tuple(PREDICTOR_VARIABLES)
    part.rename(out)
    return qa


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--tiles", nargs="*", help="Only these tile names (e.g. e000n40)")
    args = parser.parse_args()

    (args.root / "stack").mkdir(exist_ok=True)
    src = sources(args.root)
    with LAND_TILES.open() as f:
        tiles = [(r["tile"], float(r["lon"]), float(r["lat"])) for r in csv.DictReader(f)]
    if args.tiles:
        tiles = [t for t in tiles if t[0] in args.tiles]
    ready = [
        t
        for t in tiles
        if (args.root / "ee_export" / f"terrain_{t[0]}.tif").exists()
        and not (args.root / "stack" / f"{t[0]}.tif").exists()
    ]
    built = sum((args.root / "stack" / f"{t[0]}.tif").exists() for t in tiles)
    print(json.dumps({"tiles": len(tiles), "built": built, "to_build": len(ready)}), flush=True)

    qa_path = args.root / "stack" / "qa.csv"
    with ProcessPoolExecutor(args.workers) as pool:
        futures = {pool.submit(build_tile, *t, args.root, src): t[0] for t in ready}
        for i, future in enumerate(as_completed(futures), start=1):
            qa = future.result()
            new = not qa_path.exists()
            with qa_path.open("a", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(qa))
                if new:
                    writer.writeheader()
                writer.writerow(qa)
            print(
                f"{i}/{len(ready)} {qa['tile']} land={qa['land_px']:,} "
                f"missing_any={qa['missing_any']:,}",
                flush=True,
            )

    stack_tiles = sorted((args.root / "stack").glob("*.tif"))
    if stack_tiles:
        mosaic_vrt(stack_tiles, args.root / "stack.vrt")
        print(f"stack.vrt: {len(stack_tiles)} tiles", flush=True)


if __name__ == "__main__":
    main()
