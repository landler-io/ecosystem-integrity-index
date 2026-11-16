"""Tests for pipelines/training/build_fi_change_controls.py on synthetic data."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipelines" / "training"))
import build_fi_change_controls as bfc  # noqa: E402

from eii.compute.settings import REFERENCE_YEARS  # noqa: E402
from eii.training.settings import PREDICTOR_VARIABLES  # noqa: E402


def test_draw_cells_snaps_to_unique_grid_cells() -> None:
    cells = bfc.draw_cells(50_000, np.random.default_rng(0))
    assert cells["cell_id"].is_unique
    assert cells["lat"].between(-60, 80).all()
    row, col = np.divmod(cells["cell_id"].to_numpy(), bfc.NCOLS)
    assert np.allclose(cells["lon"], -180 + col / 336) and np.allclose(cells["lat"], 80 - row / 336)


def test_within_km_uses_great_circle_distance() -> None:
    ref = pd.DataFrame({"lon": [10.0], "lat": [60.0]})
    deg_lon = 1 / (111.195 * np.cos(np.radians(60.0)))  # one km east at 60 N
    hit = bfc.within_km(10 + np.array([4.9, 5.1]) * deg_lon, np.array([60.0, 60.0]), ref, 5.0)
    assert hit.tolist() == [True, False]


def test_pick_controls_purity_quota_and_order() -> None:
    rng = np.random.default_rng(1)
    n = 400
    cands = pd.DataFrame({
        "cell_id": np.arange(n),
        "biome": np.r_[np.full(10, 1), np.full(n - 10, 2)],
        "draw_key": rng.random(n),
        "worldcover": rng.choice([10, 20, 30, 40], n),
        "wc_share": rng.uniform(0.6, 1.0, n),
    })  # fmt: skip
    out = bfc.pick_controls(cands, n=60)
    pure = cands["worldcover"].isin(bfc.NATURAL_WORLDCOVER) & (cands["wc_share"] >= bfc.PURE)
    n_small = int((pure & (cands["biome"] == 1)).sum())
    assert out["worldcover"].isin(bfc.NATURAL_WORLDCOVER).all() and (out["wc_share"] >= 0.75).all()
    assert len(out) == 60 and (out["biome"] == 1).sum() == n_small  # biome 1 short: 2 takes rest
    assert out["draw_key"].is_monotonic_increasing
    kept2 = cands[pure & (cands["biome"] == 2)].nsmallest(60 - n_small, "draw_key")["cell_id"]
    assert set(out.loc[out["biome"] == 2, "cell_id"]) == set(kept2)


def test_annual_series_counts_dates_of_other_tiles_as_invalid() -> None:
    days = [f"2020_{m:02d}_01" for m in range(1, 13)] + ["2020_12_15", "2020_12_31"]
    codes = pd.DataFrame(1000, index=[1, 2], columns=days)
    codes.loc[1, "2020_06_01"] = 5000  # point 1's tile has the June composite
    codes.loc[2, "2020_06_01"] = np.nan  # point 2's tile lacks it: interpolated between May/July
    codes[[f"2020_{m:02d}_02" for m in (2, 3)]] = 3000
    out = bfc.annual_series(codes, "T")
    assert list(out.columns) == [f"evi_T_{y}" for y in bfc.YEARS] and out["evi_T_2019"].isna().all()
    assert out.loc[1, "evi_T_2020"] == round((13 * 0.1 + 0.5 + 2 * 0.3) / 16, 5)
    assert out.loc[2, "evi_T_2020"] == round((14 * 0.1 + 2 * 0.3) / 16, 5)


def test_control_table_layout() -> None:
    controls = pd.DataFrame({"point_id": [7, 3], "biome": [1, 2], "worldcover": [10, 30]})
    pred = pd.DataFrame(np.ones((2, len(PREDICTOR_VARIABLES))), columns=list(PREDICTOR_VARIABLES))
    state = {y: np.array([0.5, 1.0]) for y in REFERENCE_YEARS}
    evi = pd.DataFrame({"evi_T_2001": [0.3], "evi_N_2025": [0.4]}, index=[3])
    table = bfc.control_table(controls, pred, state, evi)
    expected = ["point_id", "biome", "worldcover", *PREDICTOR_VARIABLES,
                *[f"fi_state_{y}" for y in REFERENCE_YEARS],
                *[f"evi_{s}_{y}" for s in "TVN" for y in bfc.YEARS]]  # fmt: skip
    assert list(table.columns) == expected and len(expected) == 98
    assert np.isnan(table.loc[0, "evi_T_2001"]) and table.loc[1, "evi_T_2001"] == 0.3
    assert table.loc[1, "evi_N_2025"] == 0.4 and table["evi_V_2010"].isna().all()


def test_pixel_of_finds_the_pixel_centre() -> None:
    h, v, n = 18, 4, 2400
    size = bfc.es.TILE_M / n
    x = bfc.es.X0 + h * bfc.es.TILE_M + (np.array([0, 1234]) + 0.5) * size
    y = bfc.es.Y0 - v * bfc.es.TILE_M - (np.array([5, 2399]) + 0.5) * size
    lat = np.degrees(y / bfc.es.EARTH_R)
    lon = np.degrees(x / (bfc.es.EARTH_R * np.cos(np.radians(lat))))
    row, col = bfc.pixel_of(lon, lat, h, v, n)
    assert row.tolist() == [5, 2399] and col.tolist() == [0, 1234]
