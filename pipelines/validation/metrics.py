"""Correlation and class metrics of the external validation (pure numpy/scipy, no I/O)."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import stats


def pearson_spearman(x: Any, y: Any) -> dict[str, float]:
    """Pearson r and Spearman rho with p-values over rows where both are finite.

    Args:
        x: First variable, array-like.
        y: Second variable, array-like of the same length.

    Returns:
        ``n``, ``pearson_r``, ``pearson_p``, ``spearman_rho`` and ``spearman_p``; NaN
        statistics when fewer than two rows are valid or either variable is constant.
    """
    x = np.asarray(x, dtype="float64")
    y = np.asarray(y, dtype="float64")
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    out = {
        "n": int(len(x)),
        "pearson_r": np.nan,
        "pearson_p": np.nan,
        "spearman_rho": np.nan,
        "spearman_p": np.nan,
    }
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return out
    pearson = stats.pearsonr(x, y)
    spearman = stats.spearmanr(x, y)
    out.update(
        pearson_r=float(pearson.statistic),
        pearson_p=float(pearson.pvalue),
        spearman_rho=float(spearman.statistic),
        spearman_p=float(spearman.pvalue),
    )
    return out


def locke_condition(
    hf: Any,
    intensive_frac: Any,
    hf_max: float = 4.0,
    wild_intensive_max: float = 0.005,
    cities_min: float = 0.5,
) -> np.ndarray:
    """Locke et al. (2019) global condition per cell.

    ``cities_farms`` if intensive land use exceeds 50 %; else ``large_wild_areas`` if the
    human footprint is below 4 and intensive land use below 0.5 %; else ``shared_lands``.
    Cities take precedence over wild areas, so the three classes are mutually exclusive.

    Args:
        hf: Human footprint values.
        intensive_frac: Fraction (0-1) of cropland plus built-up cover.
        hf_max: Footprint below which a cell can be a large wild area.
        wild_intensive_max: Intensive fraction below which a cell can be a large wild area.
        cities_min: Intensive fraction above which a cell is cities and farms.

    Returns:
        Object array of condition names with the broadcast shape of the inputs.
    """
    hf = np.asarray(hf, dtype="float64")
    intensive = np.asarray(intensive_frac, dtype="float64")
    hf, intensive = np.broadcast_arrays(hf, intensive)
    out = np.full(hf.shape, "shared_lands", dtype=object)
    out[(hf < hf_max) & (intensive < wild_intensive_max)] = "large_wild_areas"
    out[intensive > cities_min] = "cities_farms"
    return out
