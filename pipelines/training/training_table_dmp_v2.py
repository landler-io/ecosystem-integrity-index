#!/usr/bin/env python3
"""Training table with the NPP response sampled from the CLMS DMP V2.0.2 annual sums.

Takes the local-stack table (``training_table_localstack.py``: design, flag and predictor
columns) and samples every annual NPP value anew from ``pipelines/preprocessing/npp``:
``sum_YEAR`` and ``sd_YEAR`` for the observed years and the long-term response columns over the
OLCI years ``POTENTIAL_NPP_RESPONSE_YEARS``. Points, predictors and flags stay. Rows with fewer
than ``POTENTIAL_NPP_MIN_VALID_YEARS`` valid years are dropped; the run stops if more than
``MAX_DROPPED_SHARE`` of the rows go or the new response correlates with the table's previous
one below ``MIN_RESPONSE_R``.

Values are the stored integers of ``pipelines/preprocessing/npp`` (physical times 10,000), read
from the pixel that contains each point; -1 marks an invalid year.

Example:
    python pipelines/training/training_table_dmp_v2.py \\
        npp_training_samples_v2_wp11_localstack_20261028.csv.gz npp_annual_sums_v2 \\
        npp_training_samples_v2_wp11_localstack_dmpv2_olci_20261104.csv.gz
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sample_training_v2 import BANDS, long_term_response  # noqa: E402

from eii.compute.settings import (  # noqa: E402
    NPP_STORED_SCALE,
    OBSERVED_NPP_YEARS,
    POTENTIAL_NPP_MIN_VALID_YEARS,
    POTENTIAL_NPP_RESPONSE_YEARS,
)

NODATA = -1
BLOCK = 512
RESPONSE_SOURCE = "CLMS DMP 300m V2.0.2 RT6 annual sums, OLCI 2019-2025 (preprocessing/npp)"
MAX_DROPPED_SHARE = 0.02
MIN_RESPONSE_R = 0.9


def sample_pixels(path: Path | str, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Band values at the pixels containing the points, shape (bands, points), int64.

    Points are read block by block, so scattered points do not reread a block. Points outside
    the raster and nodata pixels get ``NODATA``.
    """
    with rasterio.open(path) as src:
        inverse = ~src.transform
        col = np.floor(inverse.a * lon + inverse.b * lat + inverse.c).astype("int64")
        row = np.floor(inverse.d * lon + inverse.e * lat + inverse.f).astype("int64")
        inside = (col >= 0) & (col < src.width) & (row >= 0) & (row < src.height)
        out = np.full((src.count, lon.size), NODATA, dtype="int64")
        indices = np.flatnonzero(inside)
        blocks = (row[indices] // BLOCK) * (src.width // BLOCK + 1) + col[indices] // BLOCK
        order = np.argsort(blocks, kind="stable")
        starts = np.flatnonzero(np.diff(blocks[order], prepend=-1))
        for points in np.split(indices[order], starts[1:]):
            if not points.size:
                continue
            r0, c0 = (row[points[0]] // BLOCK) * BLOCK, (col[points[0]] // BLOCK) * BLOCK
            window = Window(c0, r0, min(BLOCK, src.width - c0), min(BLOCK, src.height - r0))
            raw = src.read(window=window)[:, row[points] - r0, col[points] - c0]
            values = raw.astype("int64")
            if src.nodata is not None:
                values[raw == src.nodata] = NODATA
            out[:, points] = values
    return out


def annual_columns(
    annual_dir: Path, lon: np.ndarray, lat: np.ndarray, workers: int
) -> pd.DataFrame:
    """``sum_YEAR`` and ``sd_YEAR`` for the observed years from ``annual_npp_YEAR.tif``."""
    years = range(OBSERVED_NPP_YEARS[0], OBSERVED_NPP_YEARS[1] + 1)
    paths = [annual_dir / f"annual_npp_{year}.tif" for year in years]
    with ProcessPoolExecutor(workers) as pool:
        sampled = pool.map(sample_pixels, paths, [lon] * len(paths), [lat] * len(paths))
        columns: dict[str, np.ndarray] = {}
        for year, values in zip(years, sampled, strict=True):
            columns[f"sum_{year}"], columns[f"sd_{year}"] = values
    return pd.DataFrame(columns)


def rebuild(table: pd.DataFrame, annual: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Table with the new response and a report comparing it with the old one."""
    years = range(POTENTIAL_NPP_RESPONSE_YEARS[0], POTENTIAL_NPP_RESPONSE_YEARS[1] + 1)
    response = long_term_response(annual, years)
    old_response = table["longterm_avg_npp_sum"].to_numpy("float64")
    out = table.copy()
    out[annual.columns] = annual.to_numpy()
    out[response.columns] = response.to_numpy()
    out["n_valid_years"] = response[[f"npp_valid_years_{b}" for b in BANDS]].min(axis=1)
    out["response_source"] = RESPONSE_SOURCE
    keep = (out["n_valid_years"] >= POTENTIAL_NPP_MIN_VALID_YEARS).to_numpy()

    new = response["longterm_avg_npp_sum"].to_numpy("float64")[keep] * NPP_STORED_SCALE
    old = old_response[keep] * NPP_STORED_SCALE
    positive = old > 50
    report = {
        "rows_in": len(table),
        "rows_out": int(keep.sum()),
        "dropped": int((~keep).sum()),
        "dropped_share": float((~keep).mean()),
        "dropped_selected": int((~keep & table["selected"].astype(bool).to_numpy()).sum()),
        "dropped_reference": int((~keep & table["reference"].astype(bool).to_numpy()).sum()),
        "response_pearson_r_old_new": float(np.corrcoef(new, old)[0, 1]),
        "response_median_ratio_new_over_old": float(np.median(new[positive] / old[positive])),
        "response_mean_old_gc": float(old.mean()),
        "response_mean_new_gc": float(new.mean()),
    }
    return out[keep].reset_index(drop=True), report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("table", type=Path)
    parser.add_argument("annual_dir", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()

    table = pd.read_csv(args.table)
    lon, lat = table["longitude"].to_numpy(), table["latitude"].to_numpy()
    annual = annual_columns(args.annual_dir, lon, lat, args.workers)
    rebuilt, report = rebuild(table, annual)

    print(json.dumps(report, indent=1))
    if report["dropped_share"] > MAX_DROPPED_SHARE:
        raise SystemExit(f"Too many rows dropped: {report['dropped_share']:.3%}")
    if report["response_pearson_r_old_new"] < MIN_RESPONSE_R:
        raise SystemExit(f"Response r too low: {report['response_pearson_r_old_new']:.3f}")
    rebuilt.to_csv(args.out, index=False)
    args.out.with_suffix(".report.json").write_text(json.dumps(report, indent=1) + "\n")
    print(f"sha256 {hashlib.sha256(args.out.read_bytes()).hexdigest()}  {args.out}")


if __name__ == "__main__":
    main()
