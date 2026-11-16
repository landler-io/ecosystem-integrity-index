#!/usr/bin/env python3
"""Sample the local predictor stack (``build_predictor_stack.py``) at table points.

Reads ``point_id``, ``longitude`` and ``latitude`` from a CSV, samples every band of
``<root>/stack.vrt`` at the pixel that contains each point, and writes ``point_id`` plus one
column per stack band (NaN where the stack has nodata). Points are sorted by position first, so
reads stay within few tiles at a time.

Example:
    python pipelines/training/extract_stack_at_points.py predictor_stack \\
        npp_training_samples_v2_wp11_blockcap_seed42_20261027.csv.gz stack_values_wp11.csv.gz
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("points", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()

    points = pd.read_csv(args.points, usecols=["point_id", "longitude", "latitude"])
    points = points.sort_values(["latitude", "longitude"], ascending=[False, True])
    with rasterio.open(args.root / "stack.vrt") as src:
        values = np.array(
            list(src.sample(zip(points["longitude"], points["latitude"], strict=True)))
        )
        names = list(src.descriptions)
        nodata = src.nodata
    values = np.where(values == nodata, np.nan, values)
    out = pd.concat(
        [points[["point_id"]].reset_index(drop=True), pd.DataFrame(values, columns=names)],
        axis=1,
    ).sort_values("point_id")
    out.to_csv(args.out, index=False)
    complete = out[names].notna().all(axis=1)
    print(f"{len(out):,} points, {complete.sum():,} complete -> {args.out}")


if __name__ == "__main__":
    main()
