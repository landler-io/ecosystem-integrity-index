#!/usr/bin/env python3
"""Fit the v2 potential-NPP model, its prediction intervals and the natural inter-annual CV model.

Input is the training table of ``training_table_dmp_v2.py`` (published as
``gs://landler-opendata-gee-assets/products/eii/v2_0/training/
npp_training_samples_v2_wp11_localstack_dmpv2_olci_20261104.csv.gz``, requester pays). Rows with
a negative response or a missing predictor, coordinate or biome are dropped. Response: the mean
annual CLMS DMP V2.0.2 NPP of the OLCI years, in gC m-2 yr-1.

Training rows are the reference rows (near-pristine sites) plus gap-fill rows
(``eii.training.sampling.gap_fill``: low-HMI rows topping up 4-degree x aridity strata short of
reference rows), weighted 0.5. The ``selected`` rows are the scored rows of the spatial
cross-validation: 5 folds of whole 4-degree blocks, training rows within 1 degree of a test
block left out.

- Mean model: LightGBM Tweedie (``MEAN_PARAMS``) on all training rows, ``lgbm_mean.txt``.
- Bias correction: potential NPP = c x LightGBM, with c = sum(observed) / sum(predicted) over
  the out-of-fold predictions of the scored rows.
- Intervals: a quantile regression forest (``QRF_PARAMS``) on the same rows, ``qrf.joblib``.
  The interval of level L is ``[P_c - k_b (q50 - q16), P_c + k_b (q84 - q50)]``; k_b is the
  per-biome split-conformal factor of the out-of-fold scores (pooled where a biome has fewer
  than ``MIN_CALIBRATION_ROWS`` scored rows), scales are floored at ``EPS`` gC.
- Natural inter-annual CV: LightGBM Tweedie (``CV_PARAMS``) on the CV (sample std / mean) of the
  annual sums at reference rows with enough valid years and a mean of at least ``P_FLOOR``,
  times its own out-of-fold bias correction (scored rows with a CV response),
  ``lgbm_interannual_cv.txt``.

``model_card.json`` records the predictors, parameters and calibration that
``predict_potential_npp_v2.py`` reads. Fitting takes a few minutes on 64 cores.

Example:
    python pipelines/training/fit_potential_npp_v2.py \\
        data/npp_training_samples_v2_wp11_localstack_dmpv2_olci_20261104.csv.gz \\
        models/potential_npp_v2
"""

from __future__ import annotations

import argparse
import json
import time
from collections.abc import Callable
from pathlib import Path

import joblib
import lightgbm as lgb
import numpy as np
import pandas as pd
from quantile_forest import RandomForestQuantileRegressor

from eii.compute.settings import (
    BIOME_NAMES,
    NPP_STORED_SCALE,
    POTENTIAL_NPP_MIN_VALID_YEARS,
    POTENTIAL_NPP_RESPONSE_YEARS,
)
from eii.training import validation as protocol
from eii.training.sampling import GAP_HMI_MAX, GAP_MIN_REFERENCE, GAP_WEIGHT, gap_fill
from eii.training.settings import (
    PREDICTOR_VARIABLES,
    RF_BAG_FRACTION,
    RF_MIN_LEAF_POPULATION,
    RF_NUM_TREES,
    RF_VARIABLES_PER_SPLIT,
)

RESPONSE = "longterm_avg_npp_sum"
SEED = 42
BLOCK_DEG = 4.0
MEAN_PARAMS = {
    "objective": "tweedie",
    "n_estimators": 3000,
    "learning_rate": 0.03,
    "tweedie_variance_power": 1.35,
    "num_leaves": 31,
    "min_child_samples": 20,
    "colsample_bytree": 0.8,
    "subsample": 1.0,
    "reg_lambda": 0.0,
}
CV_PARAMS = {
    "objective": "tweedie",
    "n_estimators": 1000,
    "learning_rate": 0.03,
    "tweedie_variance_power": 1.35,
    "num_leaves": 15,
    "min_child_samples": 100,
    "colsample_bytree": 0.8,
    "subsample": 0.8,
    "reg_lambda": 1.0,
}
QRF_PARAMS = {
    "n_estimators": RF_NUM_TREES,
    "min_samples_leaf": RF_MIN_LEAF_POPULATION,
    "max_features": RF_VARIABLES_PER_SPLIT,
    "max_samples": RF_BAG_FRACTION,
}
QUANTILES = (0.05, 0.16, 0.5, 0.84, 0.95)
"""Quantiles listed in the model card; the intervals use q16, q50 and q84."""
LEVELS = {"68": 0.32, "90": 0.10}
"""Interval level -> miscoverage alpha."""
DESIGN = "norm_asym"
P_FLOOR = 50.0
"""gC m-2 yr-1; minimum mean NPP of a row with a CV response."""
EPS = 1.0
"""gC m-2 yr-1; floor of the interval scales."""
MIN_CALIBRATION_ROWS = 50
MEAN_FILE = "lgbm_mean.txt"
QRF_FILE = "qrf.joblib"
CV_FILE = "lgbm_interannual_cv.txt"


def load_table(path: Path) -> pd.DataFrame:
    """Usable rows of the training table with ``npp_gc``, ``block_4deg`` and ``biome_name``."""
    table = pd.read_csv(path)
    required = [*PREDICTOR_VARIABLES, "latitude", "longitude", "biome", RESPONSE]
    keep = (table[RESPONSE] >= 0) & table[required].notna().all(axis=1)
    table = table[keep].reset_index(drop=True)
    table["npp_gc"] = table[RESPONSE] * NPP_STORED_SCALE
    table["block_4deg"] = protocol.spatial_block_ids(table["latitude"], table["longitude"], 4)
    table["biome_name"] = table["biome"].map(BIOME_NAMES).fillna("unknown")
    return table


def block_folds(table: pd.DataFrame, scored: np.ndarray) -> np.ndarray:
    """Fold of each row: whole 4-degree blocks of the scored rows; -1 (always training) elsewhere."""
    blocks = table.loc[scored, "block_4deg"]
    folds = protocol.group_kfold_folds(blocks, protocol.N_FOLDS)
    mapping = dict(zip(blocks, folds, strict=True))
    return table["block_4deg"].map(mapping).fillna(-1).astype(int).to_numpy()


def lgbm(params: dict) -> lgb.LGBMRegressor:
    return lgb.LGBMRegressor(subsample_freq=1, random_state=SEED, n_jobs=-1, verbose=-1, **params)


def qrf() -> RandomForestQuantileRegressor:
    return RandomForestQuantileRegressor(n_jobs=-1, random_state=SEED, **QRF_PARAMS)


def out_of_fold(
    table: pd.DataFrame,
    folds: np.ndarray,
    test_rows: np.ndarray,
    train_rows: np.ndarray,
    fit_predict: Callable[[np.ndarray, np.ndarray], np.ndarray],
) -> np.ndarray:
    """Buffered block-CV predictions ``(rows, k)`` for ``test_rows`` (NaN elsewhere).

    ``fit_predict(train, test)`` fits on the boolean ``train`` rows and returns the predictions
    of the ``test`` rows.
    """
    lat, lon = table["latitude"].to_numpy(), table["longitude"].to_numpy()
    out = None
    for fold in np.unique(folds[test_rows & (folds >= 0)]):
        test = test_rows & (folds == fold)
        train = train_rows & (folds != fold)
        train &= ~protocol.spatial_buffer_mask(lat, lon, test, BLOCK_DEG)
        values = np.asarray(fit_predict(train, test), dtype="float64").reshape(int(test.sum()), -1)
        if out is None:
            out = np.full((len(table), values.shape[1]), np.nan)
        out[test] = values
        print(f"  fold {fold}: {int(train.sum()):,} training, {int(test.sum()):,} test", flush=True)
    return out


def conformal_quantile(scores: np.ndarray, alpha: float) -> float:
    """Split-conformal quantile of ``scores`` for miscoverage ``alpha``."""
    n = len(scores)
    level = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n)
    return float(np.quantile(scores, level, method="higher"))


def biome_factors(scores: np.ndarray, biomes: np.ndarray, alpha: float) -> dict[str, float]:
    """Conformal factor per biome name, plus ``__pooled__`` for small or unknown biomes."""
    pooled = conformal_quantile(scores, alpha)
    factors = {"__pooled__": pooled}
    for biome in np.unique(biomes):
        rows = biomes == biome
        factors[biome] = (
            conformal_quantile(scores[rows], alpha)
            if rows.sum() >= MIN_CALIBRATION_ROWS
            else pooled
        )
    return factors


def interval_scores(
    observed: np.ndarray, potential: np.ndarray, quantiles: np.ndarray
) -> np.ndarray:
    """Normalised exceedance of the asymmetric interval: the k needed to cover each row."""
    lower = np.maximum(quantiles[:, 1] - quantiles[:, 0], EPS)
    upper = np.maximum(quantiles[:, 2] - quantiles[:, 1], EPS)
    return np.maximum((potential - observed) / lower, (observed - potential) / upper)


def natural_cv(table: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """CV (sample std / mean) of the annual sums of the response years, and its valid rows."""
    years = range(POTENTIAL_NPP_RESPONSE_YEARS[0], POTENTIAL_NPP_RESPONSE_YEARS[1] + 1)
    annual = table[[f"sum_{year}" for year in years]].to_numpy("float64")
    annual = np.where(annual >= 0, annual * NPP_STORED_SCALE, np.nan)
    n_valid = np.isfinite(annual).sum(axis=1)
    with np.errstate(all="ignore"):
        mean = np.nanmean(annual, axis=1)
        cv = np.nanstd(annual, axis=1, ddof=1) / mean
    valid = (n_valid >= POTENTIAL_NPP_MIN_VALID_YEARS) & (mean >= P_FLOOR) & np.isfinite(cv)
    return cv, valid


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("table", type=Path, help="Training table (csv or csv.gz)")
    parser.add_argument("model_dir", type=Path, help="Output directory of the model files")
    args = parser.parse_args()
    args.model_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    table = load_table(args.table)
    scored = table["selected"].astype(bool).to_numpy()
    reference = table["reference"].astype(bool).to_numpy() | scored
    fill, _ = gap_fill(table)
    fill &= ~reference
    train_rows = reference | fill
    weight = np.where(fill, GAP_WEIGHT, 1.0)
    x = table[PREDICTOR_VARIABLES].to_numpy(dtype="float64")
    y = table["npp_gc"].to_numpy()
    folds = block_folds(table, scored)
    print(
        f"{len(table):,} rows: {int(scored.sum()):,} scored, {int(train_rows.sum()):,} training "
        f"({int(fill.sum()):,} gap fill)",
        flush=True,
    )

    def fit_mean(train: np.ndarray, test: np.ndarray) -> np.ndarray:
        return (
            lgbm(MEAN_PARAMS).fit(x[train], y[train], sample_weight=weight[train]).predict(x[test])
        )

    def fit_qrf(train: np.ndarray, test: np.ndarray) -> np.ndarray:
        model = qrf().fit(x[train], y[train], sample_weight=weight[train])
        return model.predict(x[test], quantiles=[0.16, 0.5, 0.84])

    print("mean model, block CV", flush=True)
    p_oof = out_of_fold(table, folds, scored, train_rows, fit_mean)[scored, 0]
    print("interval model, block CV", flush=True)
    q_oof = out_of_fold(table, folds, scored, train_rows, fit_qrf)[scored]
    observed = y[scored]
    c = float(observed.sum() / p_oof.sum())
    scores = interval_scores(observed, c * p_oof, q_oof)
    ok = np.isfinite(scores)
    biomes = table.loc[scored, "biome_name"].to_numpy()
    k = {level: biome_factors(scores[ok], biomes[ok], a) for level, a in LEVELS.items()}

    print("final mean and interval models", flush=True)
    mean_model = lgbm(MEAN_PARAMS).fit(
        x[train_rows], y[train_rows], sample_weight=weight[train_rows]
    )
    mean_model.booster_.save_model(str(args.model_dir / MEAN_FILE))
    interval_model = qrf().fit(x[train_rows], y[train_rows], sample_weight=weight[train_rows])
    joblib.dump(interval_model, args.model_dir / QRF_FILE, compress=3)

    print("natural CV model, block CV", flush=True)
    cv, cv_ok = natural_cv(table)
    cv_train = reference & cv_ok

    def fit_cv(train: np.ndarray, test: np.ndarray) -> np.ndarray:
        return lgbm(CV_PARAMS).fit(x[train], cv[train]).predict(x[test])

    # Folds hold out all scored rows, so the training-side buffer matches the mean model's.
    cv_oof = out_of_fold(table, folds, scored, cv_train, fit_cv)[:, 0]
    cv_scored = scored & cv_ok
    c_cv = float(cv[cv_scored].sum() / cv_oof[cv_scored].sum())
    lgbm(CV_PARAMS).fit(x[cv_train], cv[cv_train]).booster_.save_model(
        str(args.model_dir / CV_FILE)
    )

    first, last = POTENTIAL_NPP_RESPONSE_YEARS
    card = {
        "created": pd.Timestamp.now(tz="UTC").isoformat(timespec="seconds"),
        "table": args.table.name,
        "n_train": int(train_rows.sum()),
        "gap_fill": {
            "min_reference": GAP_MIN_REFERENCE,
            "hmi_max": GAP_HMI_MAX,
            "weight": GAP_WEIGHT,
            "n_rows": int(fill.sum()),
        },
        "predictors": list(PREDICTOR_VARIABLES),
        "mean_model": {"file": MEAN_FILE, **MEAN_PARAMS},
        "interval_model": {"file": QRF_FILE, "quantiles": list(QUANTILES), **QRF_PARAMS},
        "calibration": {"c": c, "p_floor": P_FLOOR, "eps": EPS, "designs": {DESIGN: k}},
        "interannual_cv_model": {
            "file": CV_FILE,
            "response": f"CV (sample std / mean) of annual CLMS DMP V2.0.2 NPP sums {first}-{last}, "
            f">= {POTENTIAL_NPP_MIN_VALID_YEARS} valid years, mean >= {P_FLOOR:g} gC m-2 yr-1",
            "predictors": list(PREDICTOR_VARIABLES),
            **CV_PARAMS,
            "c": c_cv,
            "n_train": int(cv_train.sum()),
        },
    }
    (args.model_dir / "model_card.json").write_text(json.dumps(card, indent=1) + "\n")
    print(
        json.dumps(
            {
                "c": round(c, 4),
                "c_cv": round(c_cv, 4),
                "k68_pooled": k["68"]["__pooled__"],
                "minutes": round((time.perf_counter() - t0) / 60, 1),
            }
        ),  # fmt: skip
        flush=True,
    )


if __name__ == "__main__":
    main()
