#!/usr/bin/env python3
"""Training table with the potential-NPP predictors sampled from the local predictor stack.

Takes the sample of ``sample_training_v2.py`` (design, response and flag columns) and replaces
its predictors with the values ``extract_stack_at_points.py`` sampled from the local stack
(CHELSA climate, terrain and coastal bands on the CLMS NPP grid). Pool predictors that the stack
does not hold (ph, sand, clay, tri, tpi_regional) are dropped. Rows keep their ``selected`` and
``reference`` flags; rows with a missing predictor keep NaN.

Example:
    python pipelines/training/training_table_localstack.py \\
        npp_training_samples_v2_wp11_blockcap_seed42_20261027.csv.gz \\
        stack_values_wp11.csv.gz \\
        npp_training_samples_v2_wp11_localstack_20261028.csv.gz
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from eii.training.settings import PREDICTOR_VARIABLES  # noqa: E402

DROPPED = ["ph", "sand", "clay", "tri", "tpi_regional"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path)
    parser.add_argument("stack_values", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()

    table = pd.read_csv(args.table)
    values = pd.read_csv(args.stack_values, usecols=["point_id", *PREDICTOR_VARIABLES])
    assert values["point_id"].is_unique and set(values["point_id"]) == set(table["point_id"])
    table = table.drop(columns=[c for c in [*PREDICTOR_VARIABLES, *DROPPED] if c in table])
    table = table.merge(values, on="point_id", how="left", validate="one_to_one")
    table.to_csv(args.out, index=False)

    complete = table[PREDICTOR_VARIABLES].notna().all(axis=1)
    reference = table["reference"].astype(bool)
    print(
        f"{len(table):,} rows, {complete.sum():,} complete; reference {reference.sum():,} "
        f"({(reference & complete).sum():,} complete); selected {table['selected'].sum():,} "
        f"({(table['selected'].astype(bool) & complete).sum():,} complete)"
    )
    print(f"sha256 {hashlib.sha256(args.out.read_bytes()).hexdigest()}  {args.out}")


if __name__ == "__main__":
    main()
