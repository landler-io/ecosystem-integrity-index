"""Offline test of the potential-NPP fit on a small synthetic training table."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("lightgbm")
pytest.importorskip("quantile_forest")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipelines" / "training"))
import fit_potential_npp_v2 as fit  # noqa: E402

from eii.training.settings import PREDICTOR_VARIABLES  # noqa: E402


def test_conformal_factors_cover_the_requested_share():
    rng = np.random.default_rng(0)
    scores = np.abs(rng.normal(size=2000))
    biomes = np.where(np.arange(2000) < 1980, "Boreal", "Mangrove")
    k = fit.biome_factors(scores, biomes, alpha=0.32)
    assert np.mean(scores[biomes == "Boreal"] <= k["Boreal"]) >= 0.68
    assert k["Mangrove"] == k["__pooled__"]  # 20 rows: below MIN_CALIBRATION_ROWS


def test_interval_scores_are_asymmetric_and_floored():
    quantiles = np.array([[80.0, 100.0, 140.0], [100.0, 100.0, 100.0]])
    potential = np.array([100.0, 100.0])
    scores = fit.interval_scores(np.array([60.0, 103.0]), potential, quantiles)
    np.testing.assert_allclose(scores, [2.0, 3.0])  # (100-60)/20; (103-100)/EPS


def test_natural_cv_needs_valid_years_and_a_minimum_mean():
    years = range(fit.POTENTIAL_NPP_RESPONSE_YEARS[0], fit.POTENTIAL_NPP_RESPONSE_YEARS[1] + 1)
    table = pd.DataFrame({f"sum_{y}": [5_000_000, 100_000, -1] for y in years})
    table.loc[0, "sum_2019"] = 7_000_000
    table.loc[2, ["sum_2019", "sum_2020"]] = 5_000_000
    cv, valid = fit.natural_cv(table)
    assert valid.tolist() == [True, False, False]  # 500 gC; 10 gC < floor; 2 valid years
    assert cv[0] > 0


def synthetic_table(n: int = 1500, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    table = pd.DataFrame(
        rng.normal(size=(n, len(PREDICTOR_VARIABLES))), columns=PREDICTOR_VARIABLES
    )
    table["aridity"] = rng.uniform(0.02, 1.0, n)
    npp = 400 + 150 * table["mean_annual_temp"] + 200 * table["aridity"] + rng.normal(0, 20, n)
    years = range(2014, 2026)
    for y in years:
        table[f"sum_{y}"] = (np.maximum(npp, 1) * rng.uniform(0.9, 1.1, n) * 1e4).astype(int)
    table["longterm_avg_npp_sum"] = table[[f"sum_{y}" for y in range(2019, 2026)]].mean(axis=1)
    table["latitude"] = rng.uniform(-40, 60, n)
    table["longitude"] = rng.uniform(-100, 100, n)
    table["biome"] = rng.choice([1, 6, 13], n)
    table["point_id"] = np.arange(n)
    table["hmi_2022"] = rng.uniform(0, 0.4, n)
    table["reference"] = rng.random(n) < 0.8
    table["selected"] = table["reference"] & (rng.random(n) < 0.5)
    table["exclusion_reason"] = np.where(table["reference"], "", "not_pristine")
    return table


def test_fit_writes_the_model_files_predict_reads(tmp_path, monkeypatch):
    path = tmp_path / "table.csv"
    synthetic_table().to_csv(path, index=False)
    monkeypatch.setitem(fit.MEAN_PARAMS, "n_estimators", 50)
    monkeypatch.setitem(fit.CV_PARAMS, "n_estimators", 50)
    monkeypatch.setitem(fit.QRF_PARAMS, "n_estimators", 20)
    monkeypatch.setattr(sys, "argv", ["fit", str(path), str(tmp_path / "model")])
    fit.main()
    card = json.loads((tmp_path / "model" / "model_card.json").read_text())
    for spec in (card["mean_model"], card["interval_model"], card["interannual_cv_model"]):
        assert (tmp_path / "model" / spec["file"]).exists()
    assert card["predictors"] == list(PREDICTOR_VARIABLES)
    assert 0.8 < card["calibration"]["c"] < 1.2
    k = card["calibration"]["designs"][fit.DESIGN]
    assert set(k) == {"68", "90"} and k["90"]["__pooled__"] > k["68"]["__pooled__"]
