#!/usr/bin/env python3
"""
This script calculates total annual Net Primary Production (NPP) from 10-day
(dekadal) time-series rasters. It is designed to be robust to data gaps by
performing linear interpolation and applying a quality mask from QFLAG files.

This version streams spatial windows through one process with bounded memory.

Workflow:
1.  For each year, find all 10-day DMP and corresponding QFLAG raster files.
2.  Read one spatial window across all observations at a time.
3.  Mask missing (-1) and water (-2) DMP values, and pixels whose QFLAG marks
    invalid FAPAR or meteorology (bits 1 and 2; 255 is missing).
4.  For each pixel:
    a. Create a complete daily time axis for the year (1-365 or 1-366).
    b. Linearly interpolate the 10-day data across this daily axis to fill gaps.
       Each dekad is a mean daily rate, anchored at the middle of its period.
    c. Sum the resulting daily values to get the total annual production.
    d. Compute the standard deviation of daily values (intra-annual variability).
    e. Normalize the sum to a 365-day equivalent to ensure all years are comparable.
5.  Convert results from native units (kg DM/ha/year) to (g C/m²/year).
6.  Save as multi-band INT32 GeoTIFF (band 1: sum, band 2: std). The DMP digital
    number scale (0.01) is not applied and the result is multiplied by 100, so
    stored values are physical values times 10,000 (``scale_factor=0.0001``).

Each year must contain every dekad (``REQUIRED_DEKADS``) of a single product
version. Use ``--start-year``/``--end-year`` to run years in separate processes.
"""

import argparse
import re
import sys
from contextlib import ExitStack
from datetime import datetime, timedelta
from pathlib import Path

# Add preprocessing root to path for utils
sys.path.insert(0, str(Path(__file__).parent.parent))
import numpy as np  # noqa: E402
import rasterio  # noqa: E402
from affine import Affine  # noqa: E402
from numba import njit  # noqa: E402
from rasterio.windows import Window  # noqa: E402
from utils.config_utils import load_config  # noqa: E402

# --- CONFIGURATION ---
SCRIPT_DIR = Path(__file__).parent
config = load_config(SCRIPT_DIR / "config.cfg")

CLMS_RAW_DATA_DIR = Path(config["CLMS_RAW_DATA_DIR"])
NPP_OUTPUT_DIR = Path(config["NPP_OUTPUT_DIR"])

START_YEAR = config["START_YEAR"]
END_YEAR = config["END_YEAR"]
OVERWRITE = config["OVERWRITE"]

FILENAME_DATE_PATTERN = r"(\d{8})"
FILENAME_VERSION_PATTERN = r"_GLOBE_([A-Z0-9]+)_(V\d+\.\d+\.\d+)\.tiff?$"
FILENAME_DMP_KEYWORD = "-DMP-"
FILENAME_QFLAG_KEYWORD = "-QFLAG-"

NODATA_VALUE = config["NODATA_VALUE"]
CHUNK_SIZE = config["CHUNK_SIZE"]
REQUIRED_DEKADS = 36

# QFLAG bits 1 (FAPAR invalid) and 2 (meteo invalid); identical in V1.x and V2.0.2.
QFLAG_INVALID_MASK = 3

# Grid of the existing annual NPP (datasets/clms_npp/v1_1). CLMS V2.0.2 files differ
# from it only by float rounding of the pixel height (< 1e-9 degrees over the full
# raster), so outputs are written on this exact transform to keep grids identical.
CLMS_300M_TRANSFORM = Affine(
    0.0029761904762040103, 0.0, -180.0014880952381, 0.0, -0.0029761904762040103, 80.00148809523809
)
GRID_TOLERANCE_DEGREES = 1e-9

ANNUAL_SUM_TO_GC_M2_YEAR_FACTOR = 0.05
SCALE_FACTOR = 100
# --- END CONFIGURATION ---


def get_day_of_year_from_filename(filepath: Path) -> float:
    """Day of year at the middle of the dekad that the file covers.

    CLMS names a dekad by its last day (10th, 20th or month end); the period
    starts on the 1st, 11th or 21st. A dekad value is the mean daily rate over that
    period, so it is anchored at the period's midpoint (e.g. 5.5 for Jan 1-10).
    """
    match = re.search(FILENAME_DATE_PATTERN, filepath.name)
    if not match:
        raise ValueError(f"Could not parse date from filename: {filepath.name}")
    end = datetime.strptime(match.group(1), "%Y%m%d")
    if end.day not in (10, 20) and (end + timedelta(days=1)).day != 1:
        raise ValueError(f"Not a dekad end date: {filepath.name}")
    start_day = 21 if end.day > 20 else end.day - 9
    start = end.replace(day=start_day)
    return (start.timetuple().tm_yday + end.timetuple().tm_yday) / 2


def get_product_version(filepath: Path) -> tuple[str, str]:
    """Return (sensor, version), e.g. ("OLCI", "V2.0.2"), from a CLMS filename."""
    match = re.search(FILENAME_VERSION_PATTERN, filepath.name)
    if not match:
        raise ValueError(f"Could not parse sensor and version from filename: {filepath.name}")
    return match.group(1), match.group(2)


@njit
def _interp_sum_std_numba(block_flat, observed_sorted, target_days, n_pixels):
    """Single-pass computation of sum and std using Welford's algorithm."""
    result_sum = np.empty(n_pixels, dtype=np.float32)
    result_std = np.empty(n_pixels, dtype=np.float32)
    n_target = len(target_days)

    for pix_idx in range(n_pixels):
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

    block_flat = block[sort_idx].reshape(n_time, -1)
    n_pixels = block_flat.shape[1]

    sum_flat, std_flat = _interp_sum_std_numba(block_flat, observed_sorted, target_days_f, n_pixels)

    return np.stack([sum_flat.reshape(n_y, n_x), std_flat.reshape(n_y, n_x)], axis=0)


def process_year(year: int, file_pairs: list[tuple[Path, Path]]) -> Path | None:
    if not file_pairs:
        raise ValueError(f"No file pairs provided for year {year}")

    file_pairs = sorted(file_pairs, key=lambda pair: get_day_of_year_from_filename(pair[0]))
    observed_days = np.array([get_day_of_year_from_filename(pair[0]) for pair in file_pairs])
    if len(np.unique(observed_days)) != len(observed_days):
        raise ValueError(f"Duplicate observation dates for {year}")
    if len(observed_days) != REQUIRED_DEKADS:
        raise ValueError(
            f"{year} has {len(observed_days)} dekads, expected {REQUIRED_DEKADS}; "
            "an incomplete year biases the annual sum"
        )
    products = sorted({get_product_version(dmp) for dmp, _ in file_pairs})
    if len({version for _, version in products}) != 1:
        raise ValueError(f"Mixed product versions for {year}: {products}")
    is_leap = year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
    num_days = 366 if is_leap else 365
    target_days = np.arange(1, num_days + 1)
    output_path = NPP_OUTPUT_DIR / f"annual_npp_{year}.tif"
    temporary_path = output_path.with_suffix(".partial.tif")

    # GDAL's default cache is proportional to host RAM; explicitly bound it.
    with rasterio.Env(GDAL_CACHEMAX=64 * 1024 * 1024), ExitStack() as stack:
        sources = [
            (stack.enter_context(rasterio.open(dmp)), stack.enter_context(rasterio.open(qflag)))
            for dmp, qflag in file_pairs
        ]
        ref = sources[0][0]
        for pair in sources:
            for src in pair:
                if (
                    src.count != 1
                    or src.shape != ref.shape
                    or src.crs != ref.crs
                    or src.transform != ref.transform
                ):
                    raise ValueError(f"Raster grid mismatch: {src.name}")
        profile = ref.profile.copy()
        if ref.shape == (47040, 120960):
            corner_offsets = np.subtract(
                (ref.transform * (ref.width, ref.height)),
                (CLMS_300M_TRANSFORM * (ref.width, ref.height)),
            )
            if not ref.transform.almost_equals(CLMS_300M_TRANSFORM, precision=1e-12) or np.any(
                np.abs(corner_offsets) > GRID_TOLERANCE_DEGREES
            ):
                raise ValueError(f"{year} grid differs from the CLMS 300 m grid: {ref.transform}")
            profile["transform"] = CLMS_300M_TRANSFORM
        profile.update(
            driver="GTiff",
            count=2,
            dtype="int32",
            nodata=NODATA_VALUE,
            compress="deflate",
            tiled=True,
            blockxsize=CHUNK_SIZE,
            blockysize=CHUNK_SIZE,
            BIGTIFF="YES",
            NUM_THREADS="1",
        )
        with rasterio.open(temporary_path, "w", **profile) as dst:
            sensors = "+".join(sensor for sensor, _ in products)
            dst.update_tags(
                scale_factor="0.0001",
                source=f"CLMS DMP 300m 10-daily {products[0][1]} RT6 ({sensors})",
                processing="annual_sum_std_linear_interpolation_dekad_midpoint",
                n_dekads=str(len(file_pairs)),
                BAND_1="annual_sum",
                BAND_2="annual_std",
            )
            dst.set_band_description(1, "annual_sum")
            dst.set_band_description(2, "annual_std")
            dst.scales = (1e-4, 1e-4)
            dst.units = ("gC m-2 yr-1", "gC m-2 day-1")
            for row in range(0, ref.height, CHUNK_SIZE):
                for col in range(0, ref.width, CHUNK_SIZE):
                    height = min(CHUNK_SIZE, ref.height - row)
                    width = min(CHUNK_SIZE, ref.width - col)
                    window = Window(col, row, width, height)
                    block = np.empty((len(sources), height, width), dtype=np.float32)
                    for i, (dmp, qflag) in enumerate(sources):
                        values = dmp.read(1, window=window)
                        quality = qflag.read(1, window=window)
                        block[i] = values
                        # Negative DMP is missing (-1) or water (-2).
                        block[i][(values < 0) | ((quality & QFLAG_INVALID_MASK) != 0)] = np.nan
                    stats = _interpolate_and_stats(block, observed_days, target_days)
                    stats[0] *= 365.0 / num_days
                    stats *= ANNUAL_SUM_TO_GC_M2_YEAR_FACTOR * SCALE_FACTOR
                    np.rint(stats, out=stats)
                    stats[np.isnan(stats)] = NODATA_VALUE
                    dst.write(stats.astype(np.int32), window=window)
                print(
                    f"{year}: processed {min(row + CHUNK_SIZE, ref.height)}/{ref.height} rows",
                    flush=True,
                )
    # A failed or interrupted run leaves only a .partial.tif, never a completed name.
    temporary_path.replace(output_path)
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--start-year", type=int, default=START_YEAR)
    parser.add_argument("--end-year", type=int, default=END_YEAR)
    return parser.parse_args()


def main():
    args = parse_args()
    if not CLMS_RAW_DATA_DIR or str(CLMS_RAW_DATA_DIR) == "":
        raise ValueError("CLMS_RAW_DATA_DIR not set in config.cfg")

    if not CLMS_RAW_DATA_DIR.exists():
        raise FileNotFoundError(f"Input directory not found: {CLMS_RAW_DATA_DIR}")

    if not NPP_OUTPUT_DIR.exists():
        NPP_OUTPUT_DIR.mkdir(parents=True)

    if not isinstance(CHUNK_SIZE, int) or not 16 <= CHUNK_SIZE <= 512 or CHUNK_SIZE % 16:
        raise ValueError("CHUNK_SIZE must be a multiple of 16 between 16 and 512")

    for year in range(args.start_year, args.end_year + 1):
        output_path = NPP_OUTPUT_DIR / f"annual_npp_{year}.tif"

        if output_path.exists() and not OVERWRITE:
            continue

        # Exact suffixes skip temporary files of interrupted downloads (NAME.tiff.1a2B3c4D).
        year_dmp_files = sorted(
            path
            for path in CLMS_RAW_DATA_DIR.glob(f"{year}/**/*{FILENAME_DMP_KEYWORD}*")
            if path.suffix in (".tif", ".tiff")
        )

        if not year_dmp_files:
            continue

        year_file_pairs = []
        for dmp_path in year_dmp_files:
            qflag_name = dmp_path.name.replace(FILENAME_DMP_KEYWORD, FILENAME_QFLAG_KEYWORD)
            qflag_path = dmp_path.with_name(qflag_name)

            if not qflag_path.exists():
                raise FileNotFoundError(f"Missing QFLAG for {dmp_path.name}")
            year_file_pairs.append((dmp_path, qflag_path))

        if not year_file_pairs:
            continue

        print(f"Processing {year}: {len(year_file_pairs)} observations", flush=True)
        process_year(year, year_file_pairs)


if __name__ == "__main__":
    main()
