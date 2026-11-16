#!/usr/bin/env python3
"""Derive the stored constants of the FI change term from the near-pristine control points.

Writes the files of ``src/eii/compute/fi_change_constants/``:

- ``drift_factors.csv``: log drift factor per sensor (T MODIS Terra, V VIIRS Suomi NPP, Nc
  NOAA-20 VIIRS chained to Suomi NPP), level, group and year. The anomaly of a point is its log
  annual EVI minus its mean log EVI over its years; the factor is the median anomaly of the
  controls in the group, with at least ``MIN_GROUP_POINTS`` values (empty otherwise; the
  lookup then falls back to the next coarser level). Levels: ``all`` controls
  (group 0), ``L1`` biome, ``L2`` biome x natural class (``biome * 10 + class``);
- ``h_model_<series>.txt``: LightGBM model of the tolerance h per series (``TV``, ``TNc``,
  ``Nc``), the 68 % quantile of |r - 1| on the controls with the reference-year windows pooled;
  combined series add the number of contributing sensors as 16th predictor;
- ``k_c.json``: k_c per series, the smallest ``K_GRID`` value at which the 10th percentile of
  FI = FI_state x FI_change on the controls is at least ``P10_BAR`` in every reference year.

Input is one table of the control points (fresh near-pristine points: HMI < 0.05, WorldCover
tree, shrubland or grassland) with the columns ``point_id``, ``biome`` (RESOLVE number),
``worldcover`` (ESA WorldCover code), the 15 potential-NPP predictors
(``eii.training.settings.PREDICTOR_VARIABLES``), ``fi_state_<year>`` (production FI_state per
reference year) and ``evi_<sensor>_<year>`` (annual mean EVI for sensors T, V and N, years
``YEARS``; empty where missing).

The table behind the stored constants (19,984 controls) is published as
``gs://landler-opendata-gee-assets/products/eii/v2_0/training/fi_change_controls.csv.gz``
(requester pays). The code that selected the controls and extracted their EVI is not in this
repository, and ``fi_state_<year>`` comes from a ``compute_eii_v2.py`` run, which itself needs the
constants: production uses the stored copies in ``src/eii/compute/fi_change_constants/``, and
this script re-derives and checks them.

Example:
    python pipelines/training/derive_fi_change_constants.py fi_change_controls.csv.gz out/
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from eii.compute import baseline_median, chain_factor, chained, change_ratio, natural_class
from eii.compute import fi_change as fi_change_score
from eii.compute import functional_with_change as product
from eii.compute.settings import (
    FI_CHANGE_CONSTANTS_SHA256,
    FI_CHANGE_H_FLOOR,
    FI_CHANGE_SERIES,
    FI_CHANGE_ZERO_R,
    REFERENCE_YEARS,
)
from eii.training.settings import PREDICTOR_VARIABLES

YEARS = list(range(2001, 2026))
"""Columns of the annual EVI series."""
DRIFT_SENSORS = ("T", "V", "Nc")
MIN_GROUP_POINTS = 10
K_GRID = (1.0, 1.5, 2.0, 2.5, 3.0)
P10_BAR = 0.90
MIN_WINDOW_POINTS = 30
"""Controls with r and h a reference year needs to count in the k_c selection."""
H_QUANTILE = 0.68
H_MAX_ROWS = 100_000
H_SEED = 35
H_PARAMS = {
    "objective": "quantile",
    "alpha": H_QUANTILE,
    "n_estimators": 200,
    "learning_rate": 0.05,
    "num_leaves": 31,
    "min_child_samples": 200,
    "subsample": 0.8,
    "subsample_freq": 1,
    "colsample_bytree": 0.8,
    "num_threads": 12,
    "verbose": -1,
}
"""LightGBM settings of h; ``num_threads`` is fixed because it is written into the model file."""


def drift_table(x: np.ndarray, biome: np.ndarray, cls: np.ndarray) -> pd.DataFrame:
    """Log drift factors of one sensor per level, group and year.

    Args:
        x: Annual EVI of the controls (points, ``YEARS``), NaN where missing.
        biome: RESOLVE biome number per point.
        cls: Natural class per point (1 tree, 2 shrubland, 3 grassland).

    Returns:
        Columns ``level``, ``group``, ``year``, ``log_factor`` for groups with enough values.
    """
    if (cls == 0).any():
        raise ValueError("controls must be tree, shrubland or grassland")
    fit = (np.nan_to_num(x) > 0).any(1)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        logx = np.log(np.where(x > 0, x, np.nan))
        anom = logx - np.nanmean(logx, 1, keepdims=True)
    keys = {"all": np.zeros(len(x), int), "L1": biome, "L2": biome * 10 + cls}
    rows = []
    for level, key in keys.items():
        g = pd.DataFrame(anom[fit], columns=YEARS).groupby(key[fit])
        med = g.median().where(g.count() >= MIN_GROUP_POINTS).rename_axis("group").reset_index()
        long = med.melt("group", var_name="year", value_name="log_factor")
        rows.append(long.sort_values(["group", "year"], kind="stable").assign(level=level))
    return pd.concat(rows)[["level", "group", "year", "log_factor"]].reset_index(drop=True)


def corrected(x: np.ndarray, table: pd.DataFrame, biome: np.ndarray, cls: np.ndarray) -> np.ndarray:
    """Series divided by its drift factor: biome x class, else biome, else all; 1 where none."""
    wide = table.pivot_table("log_factor", ["level", "group"], "year").reindex(columns=YEARS)
    every = wide.loc["all"].to_numpy()[0]
    l1 = wide.loc["L1"].reindex(biome).to_numpy()
    l1 = np.where(np.isnan(l1), every, l1)
    l2 = wide.loc["L2"].reindex(biome * 10 + cls).to_numpy()
    factor = np.nan_to_num(np.where(np.isnan(l2), l1, l2))
    return (x / np.exp(factor)).astype("float32")


def combine(rs: list[dict[int, np.ndarray]]) -> tuple[dict[int, np.ndarray], dict[int, np.ndarray]]:
    """Mean of the sensors' r with a value, and the number of such sensors, per reference year."""
    out, count = {}, {}
    for y in REFERENCE_YEARS:
        stack = np.stack([r[y] for r in rs])
        count[y] = np.isfinite(stack).sum(0).astype("float64")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            out[y] = np.nanmean(stack, 0)
    return out, count


def fit_h(dev: dict[int, np.ndarray], count: dict[int, np.ndarray] | None, x: np.ndarray):  # noqa: ANN201
    """LightGBM booster of h on |r - 1| of all controls, reference years pooled.

    Args:
        dev: r - 1 per reference year.
        count: Sensors with a value per reference year (combined series), else None.
        x: The 15 predictors (points, 15).

    Returns:
        The fitted ``lightgbm.Booster``.
    """
    import lightgbm as lgb

    years = list(dev)
    if count is None:
        parts, xs = [dev[y] for y in years], x
        idx = np.concatenate([np.flatnonzero(np.isfinite(d)) for d in parts])
        y = np.concatenate([np.abs(d[np.isfinite(d)]) for d in parts])
    else:
        d = np.concatenate([dev[y] for y in years])
        xs = np.vstack([np.column_stack([x, count[y]]) for y in years])
        idx, y = np.flatnonzero(np.isfinite(d)), np.abs(d[np.isfinite(d)])
    use = np.arange(len(idx))
    if len(use) > H_MAX_ROWS:
        use = np.random.default_rng(H_SEED).choice(use, H_MAX_ROWS, replace=False)
    return lgb.LGBMRegressor(**H_PARAMS).fit(xs[idx[use]], y[use]).booster_


def predict_h(booster, x: np.ndarray, count: np.ndarray | None) -> np.ndarray:  # noqa: ANN001
    """h, floored, where the predictors (and a positive sensor count, if used) are present."""
    feats = (x if count is None else np.column_stack([x, count])).astype("float64")
    ok = np.isfinite(feats).all(1) & (True if count is None else count > 0)
    h = np.full(len(x), np.nan)
    if ok.any():
        h[ok] = np.maximum(booster.predict(feats[ok]), FI_CHANGE_H_FLOOR)
    return h


def p10(state: np.ndarray, r: np.ndarray, h: np.ndarray, k: float) -> float:
    """10th percentile of FI_state x FI_change over the controls with a value."""
    fi = product(state, fi_change_score(r, h, k, FI_CHANGE_ZERO_R))
    return float(np.quantile(fi[np.isfinite(fi)], 0.1))


def select_k(
    state: dict[int, np.ndarray], r: dict[int, np.ndarray], h: dict[int, np.ndarray]
) -> tuple[float, bool]:
    """Smallest ``K_GRID`` value with p10 >= ``P10_BAR`` in every served year, and if one met it.

    Years count where at least ``MIN_WINDOW_POINTS`` controls have r and h. Without such a value
    the largest is returned with False.
    """
    years = [y for y in r if (np.isfinite(r[y]) & np.isfinite(h[y])).sum() >= MIN_WINDOW_POINTS]
    for k in K_GRID:
        if years and min(p10(state[y], r[y], h[y], k) for y in years) >= P10_BAR:
            return k, True
    return K_GRID[-1], False


def derive(controls: pd.DataFrame, out: Path) -> None:
    """Write the drift table, the h models and ``k_c.json`` of ``controls`` to ``out``."""
    biome = controls["biome"].to_numpy(int)
    cls = natural_class(controls["worldcover"].to_numpy())
    evi = {s: controls[[f"evi_{s}_{y}" for y in YEARS]].to_numpy("float64") for s in "TVN"}
    raw = {"T": evi["T"], "V": evi["V"]}
    raw["Nc"] = chained(evi["V"], evi["N"], YEARS, chain_factor(evi["V"], evi["N"], YEARS))
    tables, rs = [], {}
    for s in DRIFT_SENSORS:
        table = drift_table(raw[s], biome, cls)
        xc = corrected(raw[s], table, biome, cls)
        median = baseline_median(xc, YEARS, s)
        rs[s] = {y: change_ratio(xc, YEARS, y, median) for y in REFERENCE_YEARS}
        tables.append(table.assign(sensor=s))
    x = controls[list(PREDICTOR_VARIABLES)].to_numpy("float64")
    state = {y: controls[f"fi_state_{y}"].to_numpy("float64") for y in REFERENCE_YEARS}
    out.mkdir(parents=True, exist_ok=True)
    kc: dict = {"n_controls": len(controls), "series": {}}
    for name, members in FI_CHANGE_SERIES.items():
        r, count = combine([rs[s] for s in members])
        n = None if len(members) == 1 else count
        booster = fit_h({y: r[y] - 1 for y in REFERENCE_YEARS}, n, x)
        h = {y: predict_h(booster, x, None if n is None else n[y]) for y in REFERENCE_YEARS}
        k, met = select_k(state, r, h)
        booster.save_model(str(out / f"h_model_{name}.txt"))
        kc["series"][name] = {
            "sensors": list(members),
            "k_c": k,
            "k_c_all_controls": k,
            "met_all_controls": met,
            "p10_all_controls_production": {
                str(y): round(p10(state[y], r[y], h[y], k), 4) for y in REFERENCE_YEARS
            },
        }
    drift = pd.concat(tables)[["sensor", "level", "group", "year", "log_factor"]]
    drift.to_csv(out / "drift_factors.csv", index=False, float_format="%.6f")
    (out / "k_c.json").write_text(json.dumps(kc, indent=1) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("controls", type=Path, help="Control-point table (CSV)")
    parser.add_argument("out", type=Path, help="Output folder")
    args = parser.parse_args()
    derive(pd.read_csv(args.controls), args.out)
    stale = 0
    for name, digest in FI_CHANGE_CONSTANTS_SHA256.items():
        found = hashlib.sha256((args.out / name).read_bytes()).hexdigest()
        same = found == digest
        stale += not same
        print(f"{found}  {name}  {'matches' if same else 'differs from'} the stored constants")
    sys.exit(1 if stale else 0)


if __name__ == "__main__":
    main()
