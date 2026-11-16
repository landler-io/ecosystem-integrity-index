"""
FI change term (v2.0): vegetation loss against the pixel's own EVI baseline.

FI = FI_state x FI_change. FI_change compares the mean EVI of the 3-year window ending in the
reference year with the median of a fixed, drift-corrected baseline of the same sensor
(r = window / baseline); it is 1 while the deficit 1 - r stays within k_c h and falls linearly
to 0 at r = 0. Sensors are scored on their own grid and combined after resampling to the EII
grid (``combine_ratios``). The constants (drift factors, h models, k_c) are derived by
``pipelines/training/derive_fi_change_constants.py``.

Arrays are (pixels, years) unless noted, NaN for nodata. Annual series carry a ``years``
sequence that names their columns.
"""

from __future__ import annotations

import hashlib
import json
import warnings
from collections.abc import Mapping, Sequence
from functools import lru_cache
from importlib import resources
from typing import Any

import numpy as np

from .settings import (
    EVI_COMPOSITES_PER_YEAR,
    EVI_INVALID_CODE,
    EVI_MAX_GAP,
    EVI_MIN_COMPOSITES,
    EVI_SNOW_CODE,
    EVI_SNOW_PERCENTILE,
    FI_CHANGE_BASELINES,
    FI_CHANGE_CHAIN_MIN_YEARS,
    FI_CHANGE_CHAIN_SWITCH,
    FI_CHANGE_CHAIN_YEARS,
    FI_CHANGE_CONSTANTS_SHA256,
    FI_CHANGE_H_FLOOR,
    FI_CHANGE_K,
    FI_CHANGE_MIN_BASELINE_YEARS,
    FI_CHANGE_MIN_WINDOW_YEARS,
    FI_CHANGE_SERIES,
    FI_CHANGE_WINDOW_YEARS,
    FI_CHANGE_ZERO_R,
)

MODIS_EVI_FILL = -3000
WORLDCOVER_NATURAL = {10: 1, 20: 2, 30: 3}
"""WorldCover tree, shrubland and grassland -> natural class of the drift groups."""

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------


def constants_sha256() -> dict[str, str]:
    """SHA-256 of the stored constant files, after checking them against the settings.

    Raises:
        ValueError: A stored file differs from ``FI_CHANGE_CONSTANTS_SHA256``, or ``k_c.json``
            disagrees with ``FI_CHANGE_K``.
    """
    return dict(_checked_constants())


@lru_cache(maxsize=1)
def _checked_constants() -> tuple[tuple[str, str], ...]:
    folder = resources.files("eii.compute") / "fi_change_constants"
    for name, digest in FI_CHANGE_CONSTANTS_SHA256.items():
        found = hashlib.sha256((folder / name).read_bytes()).hexdigest()
        if found != digest:
            raise ValueError(f"FI change constant {name} has SHA-256 {found}, expected {digest}")
    stored = json.loads((folder / "k_c.json").read_text())["series"]
    if {s: v["k_c"] for s, v in stored.items()} != FI_CHANGE_K:
        raise ValueError("k_c.json does not match FI_CHANGE_K")
    return tuple(FI_CHANGE_CONSTANTS_SHA256.items())


def _constant_path(name: str) -> Any:
    _checked_constants()
    return resources.files("eii.compute") / "fi_change_constants" / name


@lru_cache(maxsize=1)
def drift_table() -> dict[str, dict[tuple[str, int], dict[int, float]]]:
    """Stored log drift factors: sensor -> (level, group) -> year -> factor (NaN if none)."""
    import csv

    out: dict[str, dict[tuple[str, int], dict[int, float]]] = {}
    with _constant_path("drift_factors.csv").open() as f:
        for row in csv.DictReader(f):
            key = (row["level"], int(row["group"]))
            value = float(row["log_factor"]) if row["log_factor"] else np.nan
            out.setdefault(row["sensor"], {}).setdefault(key, {})[int(row["year"])] = value
    return out


@lru_cache(maxsize=3)
def h_model(series: str) -> Any:
    """Stored LightGBM booster of the tolerance h of ``series`` (needs ``lightgbm``)."""
    import lightgbm as lgb

    return lgb.Booster(model_str=_constant_path(f"h_model_{series}.txt").read_text())


# --------------------------------------------------------------------------------------
# Annual EVI per sensor
# --------------------------------------------------------------------------------------


def encode_modis_evi(evi: np.ndarray, reliability: np.ndarray) -> np.ndarray:
    """Encode MOD13Q1 composites: EVI x 1e4 where ``pixel_reliability`` is 0-1, snow where 2.

    Args:
        evi: ``250m 16 days EVI`` (x 1e4, fill -3000).
        reliability: ``250m 16 days pixel reliability``.

    Returns:
        int16 codes: EVI, ``EVI_SNOW_CODE`` or ``EVI_INVALID_CODE``.
    """
    evi, rel = np.asarray(evi, "int64"), np.asarray(reliability, "int64")
    out = np.where((evi != MODIS_EVI_FILL) & (rel >= 0) & (rel <= 1), evi, EVI_INVALID_CODE)
    return np.where(rel == 2, EVI_SNOW_CODE, out).astype("int16")


def encode_viirs_evi(evi: np.ndarray, reliability: np.ndarray, qa: np.ndarray) -> np.ndarray:
    """Encode VNP13A1 / VJ113A1 composites: EVI where MODLAND_QA <= 1 and reliability <= 7.

    Args:
        evi: ``16 days EVI`` (x 1e4).
        reliability: ``16 days pixel reliability``; 8 is snow.
        qa: ``16 days VI Quality``; bits 0-1 are MODLAND_QA.

    Returns:
        int16 codes: EVI, ``EVI_SNOW_CODE`` or ``EVI_INVALID_CODE``.
    """
    evi, rel, qa = (np.asarray(a, "int64") for a in (evi, reliability, qa))
    valid = ((qa & 3) <= 1) & (rel >= 0) & (rel <= 7) & (np.abs(evi) <= 10000)
    out = np.where(rel == 8, EVI_SNOW_CODE, np.where(valid, evi, EVI_INVALID_CODE))
    return out.astype("int16")


def _row_percentile(x: np.ndarray, q: float) -> np.ndarray:
    """``np.nanpercentile(x, q, axis=1)`` (linear), by one sort; NaN where a row is empty."""
    s = np.sort(x, axis=1)
    k = np.isfinite(x).sum(1)
    pos = q / 100 * np.maximum(k - 1, 0)
    lo = np.floor(pos).astype("int64")
    hi = np.minimum(lo + 1, np.maximum(k - 1, 0))
    a = np.take_along_axis(s, lo[:, None], 1)[:, 0]
    b = np.take_along_axis(s, hi[:, None], 1)[:, 0]
    with np.errstate(invalid="ignore"):
        return np.where(k > 0, a + (b - a) * (pos - lo), np.nan)


def _fill_short_gaps(x: np.ndarray, max_gap: int) -> np.ndarray:
    """Linear interpolation along axis 1 of NaN runs of at most ``max_gap`` inside the record."""
    m = x.shape[1]
    valid = np.isfinite(x)
    idx = np.arange(m)
    prev = np.maximum.accumulate(np.where(valid, idx, -1), axis=1)
    nxt = np.minimum.accumulate(np.where(valid, idx, m)[:, ::-1], axis=1)[:, ::-1]
    rows, cols = np.nonzero(~valid & (prev >= 0) & (nxt < m) & (nxt - prev - 1 <= max_gap))
    p, q = prev[rows, cols], nxt[rows, cols]
    out: np.ndarray = np.array(x, copy=True)
    out[rows, cols] = x[rows, p] + (x[rows, q] - x[rows, p]) * (cols - p) / (q - p)
    return out


def annual_evi(
    codes: np.ndarray, composite_years: Sequence[int], years: Sequence[int]
) -> tuple[np.ndarray, np.ndarray]:
    """Annual mean EVI and the snow value of each pixel from its encoded 16-day composites.

    Snow composites take the pixel's ``EVI_SNOW_PERCENTILE`` percentile of valid EVI over all
    composites given (the sensor's record); runs of up to ``EVI_MAX_GAP`` missing composites
    are interpolated; a year needs ``EVI_MIN_COMPOSITES`` valid or snow composites.

    Args:
        codes: (pixels, composites) codes of ``encode_modis_evi`` / ``encode_viirs_evi``,
            composites in date order.
        composite_years: Calendar year of each composite.
        years: Years to return.

    Returns:
        float32 annual means (pixels, years) and the float32 snow value (pixels,).
    """
    codes = np.asarray(codes)
    evi = np.where(codes > -10001, codes / 1e4, np.nan)
    snow = codes == EVI_SNOW_CODE
    p5 = _row_percentile(evi, EVI_SNOW_PERCENTILE)
    filled = _fill_short_gaps(np.where(snow, p5[:, None], evi), EVI_MAX_GAP)
    counted = np.isfinite(evi) | snow
    year = np.asarray(composite_years)
    mean = np.full((len(codes), len(years)), np.nan)
    with np.errstate(all="ignore"):
        for j, y in enumerate(years):
            cols = year == y
            if cols.any():
                ok = counted[:, cols].sum(1) >= EVI_MIN_COMPOSITES
                total = np.nansum(filled[:, cols], 1)
                n = np.isfinite(filled[:, cols]).sum(1)
                mean[:, j] = np.where(ok & (n > 0), total / n, np.nan)
    return mean.astype("float32"), p5.astype("float32")


def complete_years(composite_dates: Sequence[str]) -> set[int]:
    """Calendar years with all ``EVI_COMPOSITES_PER_YEAR`` composites (dates ``YYYY_MM_DD``)."""
    years = [int(d[:4]) for d in set(composite_dates)]
    return {y for y in set(years) if years.count(y) >= EVI_COMPOSITES_PER_YEAR}


# --------------------------------------------------------------------------------------
# Drift correction, chain and ratio
# --------------------------------------------------------------------------------------


def natural_class(worldcover: np.ndarray) -> np.ndarray:
    """1 tree, 2 shrubland, 3 grassland, 0 otherwise, from ESA WorldCover codes."""
    wc = np.asarray(worldcover)
    return np.select([wc == c for c in WORLDCOVER_NATURAL], list(WORLDCOVER_NATURAL.values()), 0)


def _lookup(table: Mapping[tuple[str, int], Mapping[int, float]], level: str, year: int,
            key: np.ndarray) -> np.ndarray:  # fmt: skip
    """Factor of ``level`` and ``year`` for each group in ``key``; NaN where there is none."""
    lut = {g: v.get(year, np.nan) for (lev, g), v in table.items() if lev == level}
    arr = np.full(max(lut) + 1, np.nan)
    arr[list(lut)] = list(lut.values())
    ok = (key >= 0) & (key < len(arr))
    return np.where(ok, arr[np.where(ok, key, 0)], np.nan)


def drift_log_factor(
    sensor: str, years: Sequence[int], biome: np.ndarray, cls: np.ndarray
) -> np.ndarray:
    """Stored log drift factor per pixel and year: biome x class, else biome, else all.

    Years without any stored factor of the sensor (before its record) get 0.

    Args:
        sensor: ``T``, ``V`` or ``Nc``.
        years: Years of the series, none after the last stored year.
        biome: RESOLVE biome number (1-14; anything else has no biome group).
        cls: ``natural_class`` of the pixel.

    Returns:
        float64 (pixels, years) log factors.
    """
    table = drift_table()[sensor]
    last = max(table[("all", 0)])
    if max(years) > last:
        raise ValueError(f"no stored drift factors for {sensor} after {last}")
    biome, cls = np.asarray(biome, "int64"), np.asarray(cls, "int64")
    l2_key = np.where(cls > 0, biome * 10 + cls, -1)
    out = np.empty((len(biome), len(years)))
    for j, y in enumerate(years):
        f = np.full(len(biome), table[("all", 0)].get(y, np.nan))
        for level, key in (("L1", biome), ("L2", l2_key)):
            g = _lookup(table, level, y, key)
            f = np.where(np.isnan(g), f, g)
        out[:, j] = np.nan_to_num(f)
    return out


def correct_drift(
    annual: np.ndarray, sensor: str, years: Sequence[int], biome: np.ndarray, cls: np.ndarray
) -> np.ndarray:
    """Annual series divided by the stored drift factor of its sensor, group and year."""
    factor = drift_log_factor(sensor, years, biome, cls)
    return (np.asarray(annual, "float64") / np.exp(factor)).astype("float32")


def chain_factor(v: np.ndarray, n: np.ndarray, years: Sequence[int]) -> np.ndarray:
    """NOAA-20 chain factor f: median of annual V / N over the chain years (>= 2 valid)."""
    cols = [list(years).index(y) for y in range(FI_CHANGE_CHAIN_YEARS[0],
                                                 FI_CHANGE_CHAIN_YEARS[1] + 1)]  # fmt: skip
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN rows give NaN
        ratio = v[:, cols] / np.where(n[:, cols] > 0, n[:, cols], np.nan)
        ok = np.isfinite(ratio).sum(1) >= FI_CHANGE_CHAIN_MIN_YEARS
        return np.where(ok, np.nanmedian(ratio, 1), np.nan)


def chained(v: np.ndarray, n: np.ndarray, years: Sequence[int], f: np.ndarray) -> np.ndarray:
    """Nc annual series: V before ``FI_CHANGE_CHAIN_SWITCH``, N x f from it on."""
    out = np.array(v, "float64")
    j = list(years).index(FI_CHANGE_CHAIN_SWITCH)
    out[:, j:] = n[:, j:] * np.asarray(f)[:, None]
    return out


def baseline_median(x: np.ndarray, years: Sequence[int], sensor: str) -> np.ndarray:
    """Median of the sensor's baseline years (>= ``FI_CHANGE_MIN_BASELINE_YEARS`` valid)."""
    first, last = FI_CHANGE_BASELINES[sensor]
    ys = list(years)
    base = x[:, ys.index(first) : ys.index(last) + 1]
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ok = np.isfinite(base).sum(1) >= FI_CHANGE_MIN_BASELINE_YEARS
        return np.where(ok, np.nanmedian(base, 1), np.nan)


def window_mean(x: np.ndarray, years: Sequence[int], year: int) -> np.ndarray:
    """Mean of the 3-year window ending in ``year`` (>= ``FI_CHANGE_MIN_WINDOW_YEARS`` valid)."""
    j = list(years).index(year)
    win = x[:, j - FI_CHANGE_WINDOW_YEARS + 1 : j + 1]
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ok = np.isfinite(win).sum(1) >= FI_CHANGE_MIN_WINDOW_YEARS
        return np.where(ok, np.nanmean(win, 1), np.nan)


def change_ratio(x: np.ndarray, years: Sequence[int], year: int, median: np.ndarray) -> np.ndarray:
    """r of the window ending in ``year``: window mean / baseline median (NaN if median <= 0).

    Args:
        x: Drift-corrected annual series of one sensor.
        years: Years of the columns of ``x``.
        year: Reference year.
        median: ``baseline_median`` of the same corrected series.
    """
    with np.errstate(all="ignore"):
        return window_mean(x, years, year) / np.where(median > 0, median, np.nan)


def combine_ratios(ratios: Sequence[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Unweighted mean of the sensors' r that have a value, and the number of sensors."""
    stack = np.stack([np.asarray(r, "float64") for r in ratios])
    count = np.isfinite(stack).sum(0)
    with np.errstate(all="ignore"):
        mean = np.where(count > 0, np.nansum(stack, 0) / np.maximum(count, 1), np.nan)
    return mean, count.astype("uint8")


# --------------------------------------------------------------------------------------
# Tolerance, score and product
# --------------------------------------------------------------------------------------


def tolerance_h(
    series: str, predictors: np.ndarray, count: np.ndarray | None = None, threads: int = 0
) -> np.ndarray:
    """Tolerance h of ``series`` from its stored LightGBM model, floored at 0.02.

    Args:
        series: ``TV``, ``TNc`` or ``Nc``.
        predictors: (pixels, 15) potential-NPP predictors in model-card ``v2`` order
            (``eii.training.settings.PREDICTOR_VARIABLES``).
        count: Sensors contributing to r (combined series only; the 16th predictor).
        threads: LightGBM threads (0: its default).

    Returns:
        float64 h, NaN where a predictor is missing or no sensor contributes.
    """
    single = len(FI_CHANGE_SERIES[series]) == 1
    if single != (count is None):
        raise ValueError(f"{series}: pass the sensor count exactly for combined series")
    feats = np.asarray(predictors, "float64")
    if count is not None:
        feats = np.column_stack([feats, np.asarray(count, "float64")])
    ok = np.isfinite(feats).all(1)
    if count is not None:
        ok &= feats[:, -1] > 0
    h = np.full(len(feats), np.nan)
    if ok.any():
        pred = h_model(series).predict(feats[ok], num_threads=threads)
        h[ok] = np.maximum(pred, FI_CHANGE_H_FLOOR)
    return h


def fi_change(
    r: np.ndarray, h: np.ndarray, k: float, zero_r: float = FI_CHANGE_ZERO_R
) -> np.ndarray:
    """FI_change: 1 where 1 - r <= k h, then linear to 0 at r = ``zero_r``.

    Example:
        >>> fi_change(np.array([1.1, 0.95, 0.5, 0.0]), np.full(4, 0.1), 1.0).round(3).tolist()
        [1.0, 1.0, 0.556, 0.0]
    """
    r, h = np.asarray(r, "float64"), np.asarray(h, "float64")
    z = 1 - zero_r
    deficit = np.nan_to_num(np.maximum(1 - r, 0))
    t = np.clip(np.nan_to_num(k * h), 0, 0.999 * z)
    score = np.clip(1 - (deficit - t) / (z - t), 0, 1)
    return np.where(np.isfinite(r) & np.isfinite(h), score, np.nan)


def functional_with_change(fi_state: np.ndarray, change: np.ndarray) -> np.ndarray:
    """FI = FI_state x FI_change; FI_state where FI_change is missing."""
    return fi_state * np.where(np.isfinite(change), change, 1.0)


def series_for(year: int, complete: Mapping[str, set[int]]) -> str | None:
    """Series that serves reference ``year``, or None (then FI = FI_state everywhere).

    Args:
        year: Reference year; the window is ``year - 2`` to ``year``.
        complete: Per sensor (``T`` Terra, ``V`` Suomi NPP, ``N`` NOAA-20) the calendar
            years with all composites (``complete_years``).

    A sensor serves only if every year of the window is complete. The VIIRS slot is Suomi NPP
    if it serves, else NOAA-20 chained (Nc); the two are never mixed. A combination without
    stored constants (Terra alone) serves nothing.

    Example:
        >>> full = set(range(2009, 2026))
        >>> series_for(2025, {"T": full, "V": full, "N": full})
        'TV'
    """
    window = set(range(year - FI_CHANGE_WINDOW_YEARS + 1, year + 1))
    terra = "T" if window <= complete.get("T", set()) else ""
    viirs = next((s for s, k in (("V", "V"), ("Nc", "N")) if window <= complete.get(k, set())), "")
    name = terra + viirs
    return name if name in FI_CHANGE_SERIES else None
