"""Offline test of the training-sample selection on a synthetic candidate pool."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipelines" / "training"))
import sample_training_v2 as sample  # noqa: E402


def synthetic_pool(n: int = 600, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    pool = pd.DataFrame(
        {
            "longitude_point": rng.uniform(-60, 60, n),
            "latitude_point": rng.uniform(-40, 60, n),
            "biome": rng.choice([1, 6, 13], n),
            "worldcover": rng.choice([10, 20, 30, 60, 40], n, p=[0.3, 0.2, 0.2, 0.2, 0.1]),
            "pristine_mask": (rng.random(n) < 0.7).astype(int),
            "hmi_natural_mask": 1,
            "wdpa_high_protection": 0,
            "design_weight": rng.uniform(1, 3, n),
        }
    )
    for column in sample.ENV_CLUSTER_COLUMNS:
        pool[column] = rng.normal(size=n)
    for year in range(2014, 2026):
        pool[f"sum_{year}"] = rng.integers(1_000_000, 9_000_000, n)
        pool[f"sd_{year}"] = rng.integers(1_000, 9_000, n)
    pool.loc[:9, [f"sum_{y}" for y in range(2014, 2020)]] = -9999
    return pool


def test_select_keeps_relaxed_rows_and_draws_only_reference_rows():
    pool = sample.prepare_pool(synthetic_pool(), sample.SCREEN_YEARS, 8)
    assert pool.loc[:9, "exclusion_reason"].str.contains("few_valid_years").all()
    assert (pool.loc[pool["worldcover"] == 40, "exclusion_reason"] != "").all()
    shares = {1: 0.5, 6: 0.3, 13: 0.2}
    table, accounting = sample.select(
        pool, shares, total=120, min_quota=20, max_quota=None, min_distance_km=2.0,
        max_per_block=40, seed=42, keep_relaxed=True, min_cluster_share=0.01,
    )  # fmt: skip
    assert set(table["exclusion_reason"]) <= {"", "not_pristine"}
    assert table["reference"].eq(table["exclusion_reason"] == "").all()
    assert not (table["selected"] & ~table["reference"]).any()
    assert table["selected"].sum() == accounting["n_selected"].sum() <= 120
    assert table["block_4deg"].eq(sample.block_id(table["latitude"], table["longitude"], 4.0)).all()
    out = sample.add_provenance(table, accounting, sample.SCREEN_YEARS, 8, 2.0)
    assert out["response_window"].eq("2014-2024").all()
