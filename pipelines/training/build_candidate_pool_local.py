"""
Candidate pool of the potential-NPP training sample, built locally tile by tile.

Stages (``candidates`` resumable per 10x10-degree tile):

1. ``candidates``: on the 300 m HMI grid (tile + halo, wrapped at the antimeridian)
   compute the gap-filled max HMI 1990-2020, the natural-areas mask (``natural_areas_arrays``),
   distance to HMI > 0.4, WDPA high protection (``--wdpa``) and modal WorldCover 2021. Eligible:
   WorldCover not water/snow-ice/built/cropland and (HMI max <= 0.3 or WDPA). One random eligible 300 m
   pixel per 30-arc-second cell becomes a candidate; points sit a quarter pixel off the HMI
   pixel centre so they fall strictly inside one HMI pixel and one 300 m predictor pixel.
2. ``select`` (all tiles): Poisson sampling with ``inclusion_probabilities``: the biome pool
   (``pool_points``), an extra Europe quota and a minimum per HMI band.
3. ``points``: the selected candidates with their annotations and design columns as
   ``points.csv.gz`` (index ``point_id``).

The candidate pool read by ``sample_training_v2.py select`` is ``points.csv.gz`` plus, per point,
the predictors and annual CLMS NPP 2014-2025 (``sum_YEAR``, ``sd_YEAR``; -9999 where masked)
sampled at the exact point. The published pool
(``gs://landler-opendata-gee-assets/products/eii/v2_0/training/candidate_pool/
npp_candidate_pool_v2_local_seed42_20261027.csv.gz``, requester pays) holds values sampled in
Earth Engine from CLMS NPP v1.1 and the v1 predictor stack. They screen eligibility (valid years)
and define the climate clusters of the draw; ``extract_stack_at_points.py``,
``training_table_localstack.py`` and ``training_table_dmp_v2.py`` replace both in the training
table.

WDPA may not be redistributed, so ``candidates`` reads a WDPA raster you build from your own
WDPA download (``--wdpa``): a GeoTIFF on the GHM 300 m grid (``local_sources.GHM_PIXEL``,
origin ``GHM_X0``, ``GHM_Y0``), 1 where protected. The published pool used terrestrial WDPA
polygons (``MARINE`` 0) with IUCN category Ia, Ib or II or designation "National Park" and a
reported area above 5 km2, buffered 1000 m inward and rasterised in Earth Engine at
``scale=300`` in EPSG:4326. The RESOLVE ecoregions are read from the requester-pays bucket
(``EII_BILLING_PROJECT`` or ``--billing-project``).

Run from the repository root::

    python pipelines/training/build_candidate_pool_local.py candidates --workers 8 --download-ghm \
        --wdpa wdpa_high_protection.tif
    python pipelines/training/build_candidate_pool_local.py select
    python pipelines/training/build_candidate_pool_local.py points
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent / "preprocessing")]
import local_sources as src  # noqa: E402
from sample_training_v2 import POOL_MIN, POOL_TOTAL, land_shares, pool_points  # noqa: E402

from eii.training.sampling import (  # noqa: E402
    BUFFER_PIXELS,
    CANDIDATE_CELL_DEG,
    CANDIDATE_HMI_MAX,
    EXCLUDED_WORLDCOVER,
    HIGH_MOD_THRESHOLD,
    distance_to_mask_km,
    hmi_band,
    in_europe,
    inclusion_probabilities,
    natural_areas_arrays,
    one_per_cell,
)

TILE_DEG = 10
REGION = (-180, -60, 180, 80)  # CLMS NPP coverage
DIST_CAP_KM = 10.0
EUROPE_POOL = 20_000
BAND_POOL = 20_000
REFERENCE = src.REPO_ROOT / "pipelines/preprocessing/reference"
TILE_DIR = src.DATA_DIR / "tiles"


def tile_name(west: int, south: int) -> str:
    return f"{'e' if west >= 0 else 'w'}{abs(west):03d}{'n' if south >= 0 else 's'}{abs(south):02d}"


def tile_bounds(name: str) -> tuple[int, int, int, int]:
    west = int(name[1:4]) * (1 if name[0] == "e" else -1)
    south = int(name[5:7]) * (1 if name[4] == "n" else -1)
    return west, south, west + TILE_DEG, south + TILE_DEG


def all_tiles() -> list[str]:
    west, south, east, north = REGION
    return [
        tile_name(x, y) for y in range(south, north, TILE_DEG) for x in range(west, east, TILE_DEG)
    ]


def _lookup(name: str) -> dict[int, int]:
    return {int(k): int(v) for k, v in json.loads((REFERENCE / name).read_text()).items()}


# --------------------------------------------------------------------------------------
# Stage 1: candidates
# --------------------------------------------------------------------------------------


def _worldcover(vrt: Path, transform: rasterio.Affine, shape: tuple[int, int]) -> np.ndarray:
    """Modal class in parallel row strips (GDAL releases the GIL)."""
    strips = np.array_split(np.arange(shape[0]), 16)

    def run(rows: np.ndarray) -> np.ndarray:
        strip_transform = transform * rasterio.Affine.translation(0, int(rows[0]))
        return src.worldcover_mode(vrt, strip_transform, (len(rows), shape[1]))

    with ThreadPoolExecutor(len(strips)) as pool:
        return np.vstack(list(pool.map(run, strips)))


def core_pixels(bounds: tuple[float, ...]) -> tuple[np.ndarray, np.ndarray]:
    """
    Global HMI-grid rows and columns whose candidate point lies in ``bounds``.

    Points sit at pixel centre + (s/4, -s/4), so a tile owns the points in
    [west, east) x [south, north) and adjacent tiles never share a pixel.
    """
    west, south, east, north = bounds
    s, x0, y0 = src.GHM_PIXEL, src.GHM_X0, src.GHM_Y0
    cols = np.arange(math.floor((west - x0) / s) - 1, math.ceil((east - x0) / s) + 1)
    rows = np.arange(math.floor((y0 - north) / s) - 1, math.ceil((y0 - south) / s) + 1)
    lon, lat = x0 + (cols + 0.75) * s, y0 - (rows + 0.75) * s
    return rows[(lat >= south) & (lat < north)], cols[(lon >= west) & (lon < east)]


def build_candidates(tile: str, seed: int, wdpa_path: Path, metrics: dict) -> dict[str, np.ndarray]:
    west, south, east, north = tile_bounds(tile)
    s, x0, y0 = src.GHM_PIXEL, src.GHM_X0, src.GHM_Y0
    rows, cols = core_pixels((west, south, east, north))
    max_lat = min(max(abs(south), abs(north)), 85)
    halo = math.ceil(DIST_CAP_KM / (s * 111.32 * math.cos(math.radians(max_lat))))
    halo = max(halo, BUFFER_PIXELS + 3) + 1
    r0, c0 = int(rows[0]) - halo, int(cols[0]) - halo
    h, w = len(rows) + 2 * halo, len(cols) + 2 * halo

    files = [*src.GHM_SERIES.values(), src.GHM_2022]
    with src.gdal_env(gcs=False), ThreadPoolExecutor(len(files)) as pool:
        layers = list(pool.map(lambda f: src.read_ghm_window(f[0], f[1], r0, c0, h, w), files))
    hmi_max, natural = natural_areas_arrays(np.stack(layers[:-1]))
    lat_window = y0 - (np.arange(r0, r0 + h) + 0.5) * s
    dist = distance_to_mask_km(hmi_max > HIGH_MOD_THRESHOLD, lat_window, s, DIST_CAP_KM)
    core = (slice(halo, halo + len(rows)), slice(halo, halo + len(cols)))
    hmi_max, natural, dist, hmi_2022 = (a[core] for a in (hmi_max, natural, dist, layers[-1]))

    wdpa = src.read_wdpa_window(wdpa_path, int(rows[0]), int(cols[0]), len(rows), len(cols))

    with src.gdal_env(gcs=False):
        vrt = src.worldcover_vrt((west, south, east, north), TILE_DIR / tile / "wc.vrt")
        transform = rasterio.Affine(s, 0, x0 + cols[0] * s, 0, -s, y0 - rows[0] * s)
        cover = _worldcover(vrt, transform, (len(rows), len(cols)))

    excluded = np.isin(cover, list(EXCLUDED_WORLDCOVER)) | (cover == 0)
    with np.errstate(invalid="ignore"):
        eligible = ~excluded & ((hmi_max <= CANDIDATE_HMI_MAX) | wdpa)
    rr, cc = np.nonzero(eligible)
    lon = x0 + (cols[cc] + 0.75) * s
    lat = y0 - (rows[rr] + 0.75) * s
    cell = np.floor((90 - lat) / CANDIDATE_CELL_DEG).astype("int64") * 43200
    cell += np.floor((lon + 180) / CANDIDATE_CELL_DEG).astype("int64")
    rng = np.random.default_rng([seed, zlib.crc32(tile.encode())])
    pick = one_per_cell(cell, rng.random(len(cell)))
    rr, cc, lon, lat, cell = rr[pick], cc[pick], lon[pick], lat[pick], cell[pick]

    with src.gdal_env():
        eco = src.Shards("ecoregion").sample(lon, lat)[:, 0]
    to_biome, to_realm = _lookup("eco_id_to_biome.json"), _lookup("eco_id_to_realm.json")
    eco = np.nan_to_num(eco, nan=-1).astype("int64")
    hmi = hmi_max[rr, cc]
    metrics["n_eligible_pixels"] = int(eligible.sum())
    metrics["n_land_pixels"] = int((cover > 0).sum())
    return {
        "lon": lon,
        "lat": lat,
        "cell_id": cell,
        "worldcover": cover[rr, cc],
        "hmi_max": hmi.astype("float32"),
        "hmi_2022": np.nan_to_num(hmi_2022[rr, cc], nan=-1).astype("float32"),
        "dist_hmi04_km": dist[rr, cc],
        "natural": natural[rr, cc],
        "wdpa": wdpa[rr, cc],
        "biome": np.array([to_biome.get(e, -1) for e in eco], dtype="int8"),
        "realm": np.array([to_realm.get(e, 0) for e in eco], dtype="int8"),
        "europe": in_europe(lon, lat),
        "hmi_band": hmi_band(hmi).astype("int8"),
        "u": rng.random(len(lon)),
    }


def run_candidates(tile: str, seed: int, wdpa_path: Path) -> dict:
    out = TILE_DIR / tile / "candidates.npz"
    if out.exists():
        return {}
    out.parent.mkdir(parents=True, exist_ok=True)
    metrics: dict = {}
    arrays = build_candidates(tile, seed, wdpa_path, metrics)
    metrics["n_candidates"] = len(arrays["lon"])
    np.savez_compressed(out.with_suffix(".tmp.npz"), **arrays)
    out.with_suffix(".tmp.npz").rename(out)
    (out.parent / "candidates.json").write_text(json.dumps(metrics, indent=1))
    return metrics


# --------------------------------------------------------------------------------------
# Stage 2: selection
# --------------------------------------------------------------------------------------


def run_select(tiles: list[str], args: argparse.Namespace) -> dict:
    keys = ("biome", "europe", "hmi_band", "u")
    parts = {t: np.load(TILE_DIR / t / "candidates.npz") for t in tiles}
    parts = {t: {k: d[k] for k in keys} for t, d in parts.items() if len(d["u"])}
    cat = {k: np.concatenate([p[k] for p in parts.values()]) for k in keys}
    targets = pool_points(land_shares(), args.pool_total, args.pool_min)
    p = inclusion_probabilities(
        cat["biome"].astype("int64"),
        cat["europe"],
        cat["hmi_band"].astype("int64"),
        targets,
        europe_target=args.europe_pool,
        band_target=args.band_pool,
    )
    start = 0
    for tile, part in parts.items():
        n = len(part["u"])
        prob = p[start : start + n]
        chosen = np.flatnonzero(part["u"] < prob)
        np.savez(TILE_DIR / tile / "selected.npz", index=chosen, inclusion_prob=prob[chosen])
        start += n
    chosen = cat["u"] < p
    summary = {
        "n_candidates": int(len(p)),
        "n_selected": int(chosen.sum()),
        "expected_selected": float(p.sum()),
        "per_biome": pd.Series(cat["biome"][chosen]).value_counts().sort_index().to_dict(),
        "per_hmi_band": pd.Series(cat["hmi_band"][chosen]).value_counts().sort_index().to_dict(),
        "europe": int(cat["europe"][chosen].sum()),
    }
    (src.DATA_DIR / "selection.json").write_text(json.dumps(summary, indent=1, default=int))
    return summary


# --------------------------------------------------------------------------------------
# Stage 3: selected points
# --------------------------------------------------------------------------------------


def selected_points(tiles: list[str]) -> pd.DataFrame:
    """Selected candidates of ``tiles`` with their local annotations, indexed by ``point_id``."""
    frames = []
    for tile in tiles:
        if not (TILE_DIR / tile / "selected.npz").exists():
            continue
        cand = np.load(TILE_DIR / tile / "candidates.npz")
        sel = np.load(TILE_DIR / tile / "selected.npz")
        idx = sel["index"]
        natural, wdpa = cand["natural"][idx], cand["wdpa"][idx]
        frames.append(
            pd.DataFrame(
                {
                    "longitude_point": cand["lon"][idx],
                    "latitude_point": cand["lat"][idx],
                    "pool_biome": cand["biome"][idx],
                    "biome": cand["biome"][idx],
                    "realm": cand["realm"][idx],
                    "worldcover": cand["worldcover"][idx],
                    "hmi_2022": cand["hmi_2022"][idx],
                    "pristine_mask": (natural | wdpa).astype(int),
                    "hmi_natural_mask": natural.astype(int),
                    "wdpa_high_protection": wdpa.astype(int),
                    "hmi_max": cand["hmi_max"][idx],
                    "hmi_band": cand["hmi_band"][idx],
                    "dist_hmi04_km": cand["dist_hmi04_km"][idx],
                    "europe": cand["europe"][idx].astype(int),
                    "inclusion_prob": sel["inclusion_prob"],
                    "design_weight": 1 / sel["inclusion_prob"],
                    "cell_id": cand["cell_id"][idx],
                    "tile": tile,
                }
            )
        )
    return pd.concat(frames, ignore_index=True).rename_axis("point_id")


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def _parallel(fn, tiles: list[str], workers: int, *args) -> None:  # noqa: ANN001
    with ProcessPoolExecutor(workers) as pool:
        for tile, result in zip(
            tiles, pool.map(fn, tiles, *[[a] * len(tiles) for a in args]), strict=True
        ):
            print(tile, json.dumps(result), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["candidates", "select", "points"])
    parser.add_argument("--tiles", nargs="*", help="default: all tiles in REGION")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--pool-total", type=int, default=POOL_TOTAL)
    parser.add_argument("--pool-min", type=int, default=POOL_MIN)
    parser.add_argument("--europe-pool", type=int, default=EUROPE_POOL)
    parser.add_argument("--band-pool", type=int, default=BAND_POOL)
    parser.add_argument(
        "--download-ghm", action="store_true", help="candidates: whole GHM files first"
    )
    parser.add_argument(
        "--wdpa", type=Path, help="candidates: WDPA high-protection GeoTIFF on the GHM grid"
    )
    parser.add_argument(
        "--billing-project",
        help="Google Cloud project billed for the bucket reads (default: EII_BILLING_PROJECT)",
    )
    parser.add_argument("--output", type=Path, default=src.DATA_DIR / "points.csv.gz")
    args = parser.parse_args()
    if args.billing_project:
        os.environ[src.BILLING_ENV] = args.billing_project
    tiles = args.tiles or all_tiles()
    if args.stage == "candidates":
        if args.wdpa is None:
            parser.error("candidates needs --wdpa")
        if args.download_ghm:
            src.download_ghm()
        _parallel(run_candidates, tiles, args.workers, args.seed, args.wdpa)
    elif args.stage == "select":
        print(json.dumps(run_select(tiles, args), indent=1, default=int))
    else:
        points = selected_points(tiles)
        points.to_csv(args.output, compression="gzip")
        print(f"wrote {len(points):,} rows to {args.output}")


if __name__ == "__main__":
    main()
