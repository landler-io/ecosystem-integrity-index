"""Tests for pipelines/training/derive_fi_change_constants.py on synthetic controls."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipelines" / "training"))
import derive_fi_change_constants as dfc  # noqa: E402

from eii.compute.settings import REFERENCE_YEARS  # noqa: E402
from eii.training.settings import PREDICTOR_VARIABLES  # noqa: E402


def controls(n: int = 600, seed: int = 0) -> pd.DataFrame:
    """Controls in two biomes with a known yearly drift of +2 % per year in biome 6."""
    rng = np.random.default_rng(seed)
    years = np.array(dfc.YEARS)
    biome = np.where(np.arange(n) < n // 2, 4, 6)
    drift = np.where(biome[:, None] == 6, 0.02 * (years - years.mean()), 0.0)
    level = rng.uniform(0.2, 0.6, n)[:, None]
    frame = {"point_id": np.arange(n), "biome": biome, "worldcover": rng.choice([10, 20, 30], n)}
    frame |= {p: rng.normal(size=n) for p in PREDICTOR_VARIABLES}
    frame |= {f"fi_state_{y}": np.ones(n) for y in REFERENCE_YEARS}
    for sensor, first in (("T", 2001), ("V", 2012), ("N", 2018)):
        noise = rng.normal(0, 0.01, (n, len(years)))
        evi = np.where(years >= first, level * np.exp(drift + noise), np.nan)
        frame |= {f"evi_{sensor}_{y}": evi[:, j] for j, y in enumerate(years)}
    return pd.DataFrame(frame)


def test_drift_table_recovers_biome_drift() -> None:
    c = controls()
    x = c[[f"evi_T_{y}" for y in dfc.YEARS]].to_numpy()
    table = dfc.drift_table(x, c["biome"].to_numpy(), np.ones(len(c), int))
    l1 = table[table["level"] == "L1"].pivot(index="group", columns="year", values="log_factor")
    expected = 0.02 * (np.array(dfc.YEARS) - np.mean(dfc.YEARS))
    np.testing.assert_allclose(l1.loc[6].to_numpy(), expected, atol=0.005)
    np.testing.assert_allclose(l1.loc[4].to_numpy(), 0, atol=0.005)
    corrected = dfc.corrected(x, table, c["biome"].to_numpy(), np.ones(len(c), int))
    trend = np.polyfit(dfc.YEARS, np.log(corrected[c["biome"] == 6]).mean(0), 1)[0]
    assert abs(trend) < 0.002


@pytest.mark.slow
def test_derive_writes_constants(tmp_path: Path) -> None:
    pytest.importorskip("lightgbm")
    dfc.derive(controls(), tmp_path)
    kc = json.loads((tmp_path / "k_c.json").read_text())
    assert kc["n_controls"] == 600
    for series in ("TV", "TNc", "Nc"):
        assert (tmp_path / f"h_model_{series}.txt").exists()
        entry = kc["series"][series]
        assert entry["k_c"] == 1.0 and entry["met_all_controls"]
        assert min(entry["p10_all_controls_production"].values()) >= dfc.P10_BAR
    drift = pd.read_csv(tmp_path / "drift_factors.csv")
    assert set(drift["sensor"]) == {"T", "V", "Nc"}
    v = drift[(drift["sensor"] == "V") & drift["log_factor"].notna()]
    assert v["year"].min() == 2012
