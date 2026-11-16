#!/usr/bin/env python3
"""Calibration of the structural-integrity habitat-quality curve against WorldCover 2021.

Structural integrity is ``S_i = sqrt(q_i) * kernel mean of sqrt(q'_j)`` with habitat quality
``q = habitat_quality(HMI)``: a logistic curve in HMI, rescaled to 1 at HMI 0 and 0 at HMI 1
(``eii.compute.structural``), the exponential kernel (5 km scale, cut at 26.5 km) and the 300 m
grey opening ``q'`` of the landscape term.

There is no ground truth for structural integrity. The soft target is WorldCover 2021: the share
``n`` of natural classes (tree, shrub, grass, bare, snow and ice, wetland, mangrove, moss and
lichen; water excluded) in each 300 m HMI pixel, scored like S,
``T = sqrt(n_i) * kernel mean of sqrt(n_j)``. WorldCover counts pasture and plantations as
natural, so T is an upper bound for a habitat-based score: S above T credits structure the
land cover does not have.

The steepness is fixed (``QUALITY_STEEPNESS``, 10 per unit HMI); the midpoint is calibrated
over a set of candidates (default 0.3, 0.4, 0.5). Acceptance criteria, at sample points with
their RESOLVE biome and anthrome class (Ellis et al. 2010, anthromes 2000):

- C1 upper bound: per biome (n >= 100) and anthrome group, median S <= median T + 0.05.
- C2 wilderness kept: wild anthromes (61, 62) whose kernel natural share is >= 0.95 keep
  median S >= 0.95 and p10 S >= 0.85.

The pick is the gentlest midpoint (highest mean S) passing C1 and C2. On the 29,585-point
equal-area global land sample, midpoint 0.5 fails C1 in village anthromes (median S 0.13 above
T), 0.4 and 0.3 pass both criteria, and 0.4 is the pick (``QUALITY_MIDPOINT``). The sample
(``lon``, ``lat``, ``biome``, ``anthrome``) is published as ``si_calibration_sample.csv.gz``
(see ``pipelines/README.md``); the code that drew it is not in this repository.

Stages:

- ``natural``: WorldCover 2021 natural share on the HMI grid, one GeoTIFF per WorldCover tile
  from the 80 m overviews of the ESA COGs on AWS open data, mosaicked as ``natural.vrt``
  (80 m overview: bias -0.007, MAE 0.02 at 480 m against the native 10 m).
- ``calibrate``: S per candidate midpoint, T and the kernel natural share at the sample points
  (``points.csv.gz``), C1/C2 per midpoint (``criteria.csv``) and the pick.

Example:
    python pipelines/training/calibrate_structural.py natural si_calibration/natural \\
        --hmi HMv20240801_2022s_AA_300.tif
    python pipelines/training/calibrate_structural.py calibrate si_calibration_sample.csv.gz \\
        si_calibration --hmi HMv20240801_2022s_AA_300.tif --natural si_calibration/natural/natural.vrt
"""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd

from eii.compute.settings import (
    BIOME_NAMES,
    LANDSCAPE_DECAY_M,
    LANDSCAPE_OPENING_M,
    LANDSCAPE_RADIUS_M,
    METERS_PER_DEGREE,
    QUALITY_MIDPOINT,
    QUALITY_STEEPNESS,
)
from eii.compute.structural import exponential_kernel, ground_disk, habitat_quality

WORLDCOVER = "https://esa-worldcover.s3.eu-central-1.amazonaws.com"
WORLDCOVER_PREFIX = "v200/2021/map/"
WORLDCOVER_OVERVIEW = 2
NATURAL_CLASSES = (10, 20, 30, 60, 70, 90, 95, 100)
MIDPOINTS = (0.3, 0.4, 0.5)
ANTHROME_GROUPS = {
    "settlements": (11, 12),
    "villages": (21, 22, 23, 24),
    "croplands": (31, 32, 33, 34),
    "rangelands": (41, 42, 43),
    "seminatural": (51, 52, 53, 54),
    "wild": (61, 62),
}
C1_MAX_EXCESS = 0.05
C1_MIN_BIOME_ROWS = 100
C2_MIN_NATURAL = 0.95
C2_MIN_MEDIAN = 0.95
C2_MIN_P10 = 0.85


def worldcover_tiles() -> list[str]:
    """Object keys of the WorldCover 2021 map COGs in the public S3 bucket."""
    keys, token = [], None
    ns = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
    while True:
        url = f"{WORLDCOVER}/?list-type=2&prefix={WORLDCOVER_PREFIX}"
        if token:
            url += "&continuation-token=" + urllib.parse.quote(token)
        root = ET.fromstring(urllib.request.urlopen(url, timeout=60).read())
        keys += [k.text for k in root.iterfind("s3:Contents/s3:Key", ns)]
        nxt = root.find("s3:NextContinuationToken", ns)
        if nxt is None:
            return [k for k in keys if k.endswith("_Map.tif")]
        token = nxt.text


def natural_tile(key: str, out_dir: Path, hmi_path: Path) -> str:
    """Natural share of one WorldCover tile, area-averaged onto the HMI grid."""
    import rasterio
    from rasterio.warp import Resampling, reproject

    out = out_dir / (re.search(r"_([NS]\d\d[EW]\d{3})_Map", key).group(1) + ".tif")
    if out.exists():
        return out.name
    with rasterio.open(hmi_path) as hmi:
        grid, res = hmi.transform, hmi.transform.a
    with rasterio.open(f"/vsicurl/{WORLDCOVER}/{key}", OVERVIEW_LEVEL=WORLDCOVER_OVERVIEW) as src:
        lc = src.read(1)
        src_transform, crs, bounds = src.transform, src.crs, src.bounds
    land = (lc != 80) & (lc != 0)
    share = np.where(land, np.isin(lc, NATURAL_CLASSES), np.nan).astype("float32")
    c0 = math.floor((bounds.left - grid.c) / res)
    c1 = math.ceil((bounds.right - grid.c) / res)
    r0 = math.floor((grid.f - bounds.top) / res)
    r1 = math.ceil((grid.f - bounds.bottom) / res)
    dst_transform = rasterio.Affine(res, 0, grid.c + c0 * res, 0, -res, grid.f - r0 * res)
    dst = np.full((r1 - r0, c1 - c0), np.nan, dtype="float32")
    # Pixels straddling WorldCover tile edges (~0.1 %) keep one side's partial mean.
    reproject(share, dst, src_transform=src_transform, src_crs=crs, src_nodata=np.nan,
              dst_transform=dst_transform, dst_crs=crs, dst_nodata=np.nan,
              resampling=Resampling.average)  # fmt: skip
    profile = {"driver": "GTiff", "width": dst.shape[1], "height": dst.shape[0], "count": 1,
               "dtype": "float32", "crs": crs, "transform": dst_transform, "nodata": np.nan,
               "tiled": True, "compress": "deflate", "predictor": 3}  # fmt: skip
    tmp = out.with_suffix(".part.tif")
    with rasterio.open(tmp, "w", **profile) as dst_ds:
        dst_ds.write(dst, 1)
    tmp.rename(out)
    return out.name


def natural(args: argparse.Namespace) -> None:
    os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    args.out.mkdir(parents=True, exist_ok=True)
    keys = worldcover_tiles()
    print(f"{len(keys)} WorldCover tiles", flush=True)
    work = partial(natural_tile, out_dir=args.out, hmi_path=args.hmi)
    with ProcessPoolExecutor(args.workers) as pool:
        for i, name in enumerate(pool.map(work, keys), 1):
            if i % 100 == 0:
                print(i, name, flush=True)
    tiles = sorted(str(p) for p in args.out.glob("[NS]*.tif") if ".part" not in p.name)
    subprocess.run(["gdalbuildvrt", str(args.out / "natural.vrt"), *tiles], check=True)


_DATASETS: dict[str, object] = {}


def read_window(path: Path, lon: float, lat: float, half_r: int, half_c: int) -> np.ndarray:
    """Band 1 around the pixel holding ``lon, lat``, ``+- half_r`` rows, ``+- half_c`` columns.

    Pixels outside the raster and nodata are NaN.
    """
    import rasterio
    from rasterio.windows import Window

    if str(path) not in _DATASETS:
        _DATASETS[str(path)] = rasterio.open(path)
    src = _DATASETS[str(path)]
    t = src.transform
    row, col = int((t.f - lat) // -t.e), int((lon - t.c) // t.a)
    out = np.full((2 * half_r + 1, 2 * half_c + 1), np.nan, dtype="float32")
    r0, c0 = row - half_r, col - half_c
    top, left = max(r0, 0), max(c0, 0)
    bottom, right = min(row + half_r + 1, src.height), min(col + half_c + 1, src.width)
    if bottom > top and right > left:
        data = src.read(1, window=Window(left, top, right - left, bottom - top)).astype("float32")
        if src.nodata is not None and not np.isnan(src.nodata):
            data[data == src.nodata] = np.nan
        out[top - r0 : bottom - r0, left - c0 : right - c0] = data
    return out


def centre_mean(values: np.ndarray, valid: np.ndarray, kernel: np.ndarray) -> float:
    """Kernel-weighted mean of ``values`` over ``valid`` pixels, kernel centred on the window."""
    hr, hc = kernel.shape[0] // 2, kernel.shape[1] // 2
    cr, cc = values.shape[0] // 2, values.shape[1] // 2
    sl = np.s_[cr - hr : cr + hr + 1, cc - hc : cc + hc + 1]
    k = kernel * valid[sl]
    total = k.sum()
    return float((k * np.where(valid[sl], values[sl], 0)).sum() / total) if total > 0 else np.nan


def point_scores(
    hmi: np.ndarray,
    nat: np.ndarray | None,
    lat: float,
    res: float,
    midpoints: tuple[float, ...] = MIDPOINTS,
) -> dict[str, float]:
    """S per midpoint, T and the kernel natural share at the centre pixel of the windows.

    Args:
        hmi: HMI window (NaN outside land), centred on the point.
        nat: Natural-share window on the same grid, or ``None``.
        lat: Latitude of the point, for the ground kernels.
        res: Pixel size in degrees.
        midpoints: Logistic midpoints to score.

    Returns:
        ``hmi``, ``s_<midpoint>`` per midpoint and, with ``nat``, ``n``, ``nmean`` and ``t``.
    """
    from scipy.ndimage import grey_opening

    valid = np.isfinite(hmi)
    cr, cc = hmi.shape[0] // 2, hmi.shape[1] // 2
    row = {"hmi": float(hmi[cr, cc])}
    if not valid[cr, cc]:
        return row
    clipped = np.where(valid, np.clip(np.nan_to_num(hmi, nan=0.0), 0, 1), 0).astype("float64")
    kernel = exponential_kernel(LANDSCAPE_RADIUS_M, LANDSCAPE_DECAY_M, lat, res).astype("float64")
    w = 1 - clipped
    opened = 1 - grey_opening(w, footprint=ground_disk(LANDSCAPE_OPENING_M, lat, res),
                              mode="nearest")  # fmt: skip
    for mid in midpoints:
        site = math.sqrt(float(habitat_quality(clipped[cr, cc], mid, QUALITY_STEEPNESS)))
        landscape = np.sqrt(habitat_quality(opened, mid, QUALITY_STEEPNESS))
        row[f"s_{mid}"] = site * centre_mean(landscape, valid, kernel)
    if nat is not None:
        nvalid = np.isfinite(nat)
        row["n"] = float(nat[cr, cc])
        row["nmean"] = centre_mean(nat, nvalid, kernel)
        land = centre_mean(np.sqrt(nat), nvalid, kernel)
        row["t"] = math.sqrt(nat[cr, cc]) * land if nvalid[cr, cc] else np.nan
    return row


def points_chunk(
    chunk: pd.DataFrame, hmi_path: Path, natural_path: Path, midpoints: tuple[float, ...]
) -> list[dict]:
    import rasterio

    with rasterio.open(hmi_path) as src:
        res = src.transform.a
    reach_m = LANDSCAPE_RADIUS_M + 2 * LANDSCAPE_OPENING_M
    half_r = math.ceil(reach_m / (res * METERS_PER_DEGREE)) + 1
    rows = []
    for idx, lon, lat in zip(chunk.index, chunk["lon"], chunk["lat"], strict=True):
        half_c = math.ceil(half_r / max(math.cos(math.radians(lat)), 0.01)) + 1
        hmi = read_window(hmi_path, lon, lat, half_r, half_c)
        hmi[(hmi < 0) | (hmi > 1)] = np.nan
        nat = read_window(natural_path, lon, lat, half_r, half_c)
        rows.append({"idx": idx, **point_scores(hmi, nat, lat, res, midpoints)})
    return rows


def criteria(frame: pd.DataFrame, midpoints: tuple[float, ...]) -> pd.DataFrame:
    """C1 and C2 per midpoint.

    Args:
        frame: One row per point with ``biome``, ``anthrome``, ``t``, ``nmean`` and
            ``s_<midpoint>`` columns.
        midpoints: Midpoints with an ``s_<midpoint>`` column.

    Returns:
        One row per midpoint: mean and median S, the largest median excess of S over T and its
        stratum (C1), the wilderness median and p10 (C2) and the pass flags.
    """
    group = {code: g for g, codes in ANTHROME_GROUPS.items() for code in codes}
    frame = frame.assign(
        biome_name=frame["biome"].map(BIOME_NAMES), anthrome_group=frame["anthrome"].map(group)
    )
    biomes = frame.groupby("biome_name").filter(lambda g: len(g) >= C1_MIN_BIOME_ROWS)
    wild = frame[(frame["anthrome_group"] == "wild") & (frame["nmean"] >= C2_MIN_NATURAL)]
    rows = []
    for mid in midpoints:
        s = f"s_{mid}"
        excess = pd.concat([
            biomes.groupby("biome_name")[s].median() - biomes.groupby("biome_name")["t"].median(),
            frame.groupby("anthrome_group")[s].median()
            - frame.groupby("anthrome_group")["t"].median(),
        ])  # fmt: skip
        wild_median, wild_p10 = wild[s].median(), wild[s].quantile(0.1)
        rows.append({
            "midpoint": mid, "steepness": QUALITY_STEEPNESS, "n": len(frame),
            "mean_s": frame[s].mean(), "median_s": frame[s].median(), "mean_t": frame["t"].mean(),
            "max_excess": excess.max(), "worst_stratum": excess.idxmax(),
            "c1": bool(excess.max() <= C1_MAX_EXCESS),
            "wild_median": wild_median, "wild_p10": wild_p10,
            "c2": bool(wild_median >= C2_MIN_MEDIAN and wild_p10 >= C2_MIN_P10),
            "spearman_t": frame[s].corr(frame["t"], method="spearman"),
        })  # fmt: skip
    return pd.DataFrame(rows)


def pick(table: pd.DataFrame) -> float | None:
    """Gentlest midpoint (highest mean S) passing C1 and C2, or ``None``."""
    passing = table[table["c1"] & table["c2"]]
    return None if passing.empty else float(passing.loc[passing["mean_s"].idxmax(), "midpoint"])


def calibrate(args: argparse.Namespace) -> None:
    midpoints = tuple(args.midpoints)
    sample = pd.read_csv(args.sample, usecols=["lon", "lat", "biome", "anthrome"])
    order = np.lexsort((sample["lon"].to_numpy(), np.floor(90.0 - sample["lat"].to_numpy())))
    ordered = sample.iloc[order]
    chunks = [ordered.iloc[i] for i in np.array_split(np.arange(len(ordered)), args.workers * 8)]
    work = partial(points_chunk, hmi_path=args.hmi, natural_path=args.natural, midpoints=midpoints)
    rows = []
    with ProcessPoolExecutor(args.workers) as pool:
        for i, part in enumerate(pool.map(work, chunks), 1):
            rows += part
            print(f"chunk {i}/{len(chunks)}", flush=True)
    scores = pd.DataFrame(rows).set_index("idx").sort_index()
    args.out.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.out / "points.csv.gz", float_format="%.5g")

    frame = sample.join(scores)
    frame = frame[frame[[f"s_{midpoints[0]}", "t"]].notna().all(axis=1)]
    table = criteria(frame, midpoints)
    table.to_csv(args.out / "criteria.csv", index=False, float_format="%.4f")
    chosen = pick(table)
    print(table.to_string(index=False, float_format="%.4f"))
    print(f"pick: {chosen}; configured QUALITY_MIDPOINT: {QUALITY_MIDPOINT}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("natural", help="WorldCover 2021 natural share on the HMI grid")
    p.add_argument("out", type=Path, help="Output folder for the tiles and natural.vrt")
    p.add_argument("--hmi", type=Path, required=True, help="GHM 2022 300 m GeoTIFF (the grid)")
    p.add_argument("--workers", type=int, default=32)
    p.set_defaults(func=natural)
    p = sub.add_parser("calibrate", help="S, T and C1/C2 per midpoint at the sample points")
    p.add_argument("sample", type=Path, help="CSV with lon, lat, biome (RESOLVE), anthrome")
    p.add_argument("out", type=Path, help="Output folder for points.csv.gz and criteria.csv")
    p.add_argument("--hmi", type=Path, required=True, help="GHM 2022 300 m GeoTIFF")
    p.add_argument("--natural", type=Path, required=True, help="natural.vrt of the natural stage")
    p.add_argument("--midpoints", type=float, nargs="+", default=list(MIDPOINTS))
    p.add_argument("--workers", type=int, default=48)
    p.set_defaults(func=calibrate)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
