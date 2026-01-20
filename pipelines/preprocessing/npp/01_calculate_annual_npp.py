#!/usr/bin/env python3
"""
This script calculates total annual Net Primary Production (NPP) from 10-day
(dekadal) time-series rasters. It is designed to be robust to data gaps by
performing linear interpolation and applying a quality mask from QFLAG files.

This version uses Dask and xarray for scalable, parallel, out-of-core processing.

Workflow:
1.  For each year, find all 10-day DMP and corresponding QFLAG raster files.
2.  Open all rasters for a year as a single Dask-backed xarray DataArray.
3.  Apply a quality mask using the QFLAG layer to remove low-quality pixels.
4.  For each pixel:
    a. Create a complete daily time axis for the year (1-365 or 1-366).
    b. Linearly interpolate the 10-day data across this daily axis to fill gaps.
    c. Sum the resulting daily values to get the total annual production.
    d. Compute the standard deviation of daily values (intra-annual variability).
    e. Normalize the sum to a 365-day equivalent to ensure all years are comparable.
5.  Convert results from native units (kg DM/ha/year) to (g C/m²/year).
6.  Scale by 100 and save as multi-band INT32 GeoTIFF (band 1: sum, band 2: std).
"""

import re
import sys
import traceback
from datetime import datetime
from pathlib import Path

# Add preprocessing root to path for utils
sys.path.insert(0, str(Path(__file__).parent.parent))
import dask  # noqa: E402
import dask.array as da_dask  # noqa: E402
import numpy as np  # noqa: E402
import xarray as xr  # noqa: E402
from dask.distributed import Client, LocalCluster  # noqa: E402
from numba import njit, prange  # noqa: E402
from utils.config_utils import load_config  # noqa: E402

dask.config.set({"optimization.fuse.active": False})

# --- CONFIGURATION ---
SCRIPT_DIR = Path(__file__).parent
config = load_config(SCRIPT_DIR / "config.cfg")

CLMS_RAW_DATA_DIR = Path(config["CLMS_RAW_DATA_DIR"])
NPP_OUTPUT_DIR = Path(config["NPP_OUTPUT_DIR"])

START_YEAR = config["START_YEAR"]
END_YEAR = config["END_YEAR"]
OVERWRITE = config["OVERWRITE"]

FILENAME_DATE_PATTERN = r"(\d{8})"
FILENAME_DMP_KEYWORD = "-DMP-"
FILENAME_QFLAG_KEYWORD = "-QFLAG-"

DASK_N_WORKERS = config["DASK_N_WORKERS"]
DASK_MEMORY_LIMIT = config["DASK_MEMORY_LIMIT"]

NODATA_VALUE = config["NODATA_VALUE"]
CHUNK_SIZE = config["CHUNK_SIZE"]

ANNUAL_SUM_TO_GC_M2_YEAR_FACTOR = 0.05
SCALE_FACTOR = 100
# --- END CONFIGURATION ---


def get_day_of_year_from_filename(filepath: Path) -> int:
    match = re.search(FILENAME_DATE_PATTERN, filepath.name)
    if match:
        date_str = match.group(1)
        dt = datetime.strptime(date_str, "%Y%m%d")
        return dt.timetuple().tm_yday
    raise ValueError(f"Could not parse date from filename: {filepath.name}")


@njit(parallel=True, cache=True)
def _interp_sum_std_numba(block_flat, observed_sorted, target_days, n_pixels):
    """Single-pass computation of sum and std using Welford's algorithm."""
    result_sum = np.empty(n_pixels, dtype=np.float32)
    result_std = np.empty(n_pixels, dtype=np.float32)
    n_target = len(target_days)

    for pix_idx in prange(n_pixels):
        pixel_ts = block_flat[:, pix_idx]

        valid_count = 0
        for i in range(len(pixel_ts)):
            if not np.isnan(pixel_ts[i]):
                valid_count += 1

        if valid_count < 1:
            result_sum[pix_idx] = np.nan
            result_std[pix_idx] = np.nan
            continue

        valid_days = np.empty(valid_count, dtype=np.float64)
        valid_vals = np.empty(valid_count, dtype=np.float64)
        j = 0
        for i in range(len(pixel_ts)):
            if not np.isnan(pixel_ts[i]):
                valid_days[j] = observed_sorted[i]
                valid_vals[j] = pixel_ts[i]
                j += 1

        total = 0.0
        mean = 0.0
        m2 = 0.0

        for t in range(n_target):
            day = target_days[t]

            if valid_count == 1 or day <= valid_days[0]:
                val = valid_vals[0]
            elif day >= valid_days[valid_count - 1]:
                val = valid_vals[valid_count - 1]
            else:
                lo, hi = 0, valid_count - 1
                while hi - lo > 1:
                    mid = (lo + hi) // 2
                    if valid_days[mid] <= day:
                        lo = mid
                    else:
                        hi = mid
                t_frac = (day - valid_days[lo]) / (valid_days[hi] - valid_days[lo])
                val = valid_vals[lo] + t_frac * (valid_vals[hi] - valid_vals[lo])

            total += val
            delta = val - mean
            mean += delta / (t + 1)
            delta2 = val - mean
            m2 += delta * delta2

        result_sum[pix_idx] = total
        result_std[pix_idx] = np.sqrt(m2 / n_target) if n_target > 0 else 0.0

    return result_sum, result_std


def _interpolate_and_stats(block, observed_days, target_days):
    """
    Interpolate and compute sum + std in single pass.
    Returns 3D array (2, y, x) with [sum, std].
    """
    n_time, n_y, n_x = block.shape

    sort_idx = np.argsort(observed_days)
    observed_sorted = observed_days[sort_idx].astype(np.float64)
    target_days_f = target_days.astype(np.float64)

    block_flat = block[sort_idx].reshape(n_time, -1).astype(np.float64)
    n_pixels = block_flat.shape[1]

    sum_flat, std_flat = _interp_sum_std_numba(block_flat, observed_sorted, target_days_f, n_pixels)

    return np.stack([sum_flat.reshape(n_y, n_x), std_flat.reshape(n_y, n_x)], axis=0)


def process_year(year: int, file_pairs: list[tuple[Path, Path]]) -> Path | None:
    if not file_pairs:
        raise ValueError(f"No file pairs provided for year {year}")

    dmp_files, qflag_files = zip(*file_pairs, strict=True)
    days_of_year = np.array([get_day_of_year_from_filename(f) for f in dmp_files])

    with xr.open_dataarray(dmp_files[0], engine="rasterio") as ref:
        src_crs = ref.rio.crs
        src_transform = ref.rio.transform()

    dmp_da = (
        xr.open_mfdataset(
            dmp_files,
            engine="rasterio",
            chunks={"x": CHUNK_SIZE, "y": CHUNK_SIZE, "band": 1},
            concat_dim="time",
            combine="nested",
            mask_and_scale=False,
            parallel=True,
        )
        .to_array()
        .squeeze("variable", drop=True)
        .assign_coords(time=days_of_year)
    )
    if "band" in dmp_da.dims:
        dmp_da = dmp_da.squeeze("band", drop=True)
    dmp_da = dmp_da.chunk({"time": -1, "x": CHUNK_SIZE, "y": CHUNK_SIZE})

    qflag_da = (
        xr.open_mfdataset(
            qflag_files,
            engine="rasterio",
            chunks={"x": CHUNK_SIZE, "y": CHUNK_SIZE, "band": 1},
            concat_dim="time",
            combine="nested",
            mask_and_scale=False,
            parallel=True,
        )
        .to_array()
        .squeeze("variable", drop=True)
        .assign_coords(time=days_of_year)
    )
    if "band" in qflag_da.dims:
        qflag_da = qflag_da.squeeze("band", drop=True)
    qflag_da = qflag_da.chunk({"time": -1, "x": CHUNK_SIZE, "y": CHUNK_SIZE})

    da = dmp_da.where(dmp_da != NODATA_VALUE)
    quality_mask = (qflag_da.astype(np.int32) & 3) == 0
    da = da.where(quality_mask)

    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    num_days_in_year = 366 if is_leap else 365
    target_days = np.arange(1, num_days_in_year + 1)
    observed_days = days_of_year

    da_transposed = da.transpose("time", "y", "x")
    dask_data = da_transposed.data

    stats_dask = da_dask.map_blocks(
        _interpolate_and_stats,
        dask_data,
        observed_days=observed_days,
        target_days=target_days,
        dtype=np.float32,
        drop_axis=0,
        new_axis=0,
        chunks=(2, dask_data.chunks[1], dask_data.chunks[2]),
    )

    stats_da = xr.DataArray(
        stats_dask,
        dims=["band", "y", "x"],
        coords={"y": da.coords["y"], "x": da.coords["x"], "band": ["sum", "std"]},
    )
    stats_da = stats_da.rio.write_crs(src_crs)
    stats_da = stats_da.rio.write_transform(src_transform)

    annual_sum = stats_da.sel(band="sum")
    annual_std = stats_da.sel(band="std")

    annual_sum_normalized = annual_sum * (365.0 / num_days_in_year)
    annual_sum_final = (
        annual_sum_normalized * ANNUAL_SUM_TO_GC_M2_YEAR_FACTOR * SCALE_FACTOR
    ).round()
    annual_sum_final = annual_sum_final.fillna(NODATA_VALUE)

    annual_std_final = (annual_std * ANNUAL_SUM_TO_GC_M2_YEAR_FACTOR * SCALE_FACTOR).round()
    annual_std_final = annual_std_final.fillna(NODATA_VALUE)

    output_da = xr.concat([annual_sum_final, annual_std_final], dim="band")
    output_da = output_da.rio.write_crs(src_crs)
    output_da = output_da.rio.write_transform(src_transform)

    output_path = NPP_OUTPUT_DIR / f"annual_npp_{year}.tif"
    output_da.rio.to_raster(
        output_path,
        compress="deflate",
        dtype=np.int32,
        tiled=True,
        blockxsize=CHUNK_SIZE,
        blockysize=CHUNK_SIZE,
        nodata=NODATA_VALUE,
        BIGTIFF="YES",
        tags={
            "SCALE_FACTOR": str(SCALE_FACTOR),
            "UNITS": "g C / m² / year (scaled)",
            "BAND_1": "annual_sum",
            "BAND_2": "annual_std",
        },
    )

    return output_path


def main():
    if not CLMS_RAW_DATA_DIR or str(CLMS_RAW_DATA_DIR) == "":
        raise ValueError("CLMS_RAW_DATA_DIR not set in config.yaml")

    if not CLMS_RAW_DATA_DIR.exists():
        raise FileNotFoundError(f"Input directory not found: {CLMS_RAW_DATA_DIR}")

    if not NPP_OUTPUT_DIR.exists():
        NPP_OUTPUT_DIR.mkdir(parents=True)

    cluster = LocalCluster(
        n_workers=DASK_N_WORKERS,
        memory_limit=DASK_MEMORY_LIMIT,
    )
    client = Client(cluster)

    for year in range(START_YEAR, END_YEAR + 1):
        output_path = NPP_OUTPUT_DIR / f"annual_npp_{year}.tif"

        if output_path.exists() and not OVERWRITE:
            continue

        year_dmp_files = sorted(CLMS_RAW_DATA_DIR.glob(f"{year}/**/*{FILENAME_DMP_KEYWORD}*.tif*"))

        if not year_dmp_files:
            continue

        year_file_pairs = []
        for dmp_path in year_dmp_files:
            qflag_name = dmp_path.name.replace(FILENAME_DMP_KEYWORD, FILENAME_QFLAG_KEYWORD)
            qflag_path = dmp_path.with_name(qflag_name)

            if qflag_path.exists():
                year_file_pairs.append((dmp_path, qflag_path))

        if not year_file_pairs:
            continue

        try:
            process_year(year, year_file_pairs)
        except Exception:
            traceback.print_exc()

    client.close()
    cluster.close()


if __name__ == "__main__":
    main()
