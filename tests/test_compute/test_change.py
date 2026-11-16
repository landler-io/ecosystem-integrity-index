"""Offline tests of the FI change term rules and stored constants in ``eii.compute``."""

from __future__ import annotations

import numpy as np
import pytest

from eii.compute import (
    baseline_median,
    change_ratio,
    combine_ratios,
    complete_years,
    constants_sha256,
    fi_change,
    functional_with_change,
    series_for,
    tolerance_h,
)
from eii.compute.settings import FI_CHANGE_CONSTANTS_SHA256, FI_CHANGE_K

YEARS = list(range(2001, 2026))


def test_stored_constants_match_their_checksums():
    assert constants_sha256() == FI_CHANGE_CONSTANTS_SHA256
    assert FI_CHANGE_K == {"TV": 1.0, "TNc": 1.0, "Nc": 1.0}


def test_series_schedule():
    full = set(range(2009, 2026))
    complete = {"T": full | {2026}, "V": full, "N": set(range(2019, 2028))}
    assert [series_for(y, complete) for y in (2021, 2025, 2026, 2027)] == [
        "TV", "TV", "TNc", "Nc"]  # fmt: skip
    assert series_for(2026, complete | {"T": full}) == "Nc"
    assert series_for(2026, complete | {"N": {2025, 2026}}) is None  # V and N never mixed
    assert series_for(2028, complete | {"N": set()}) is None  # T alone has no constants


def test_complete_years_need_all_23_composites():
    days = [f"2020_{k:02d}" for k in range(23)] + [f"2021_{k:02d}" for k in range(22)]
    assert complete_years(days + days[:3]) == {2020}


def test_fi_change_band_and_zero_point():
    r = np.array([1.2, 0.95, 0.9, 0.5, 0.0, -0.1, np.nan, 0.8])
    h = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, np.nan])
    expected = [1.0, 1.0, 1.0, 1 - 0.4 / 0.9, 0.0, 0.0, np.nan, np.nan]
    np.testing.assert_allclose(fi_change(r, h, 1.0), expected)


def test_product_keeps_state_without_change():
    state = np.array([0.8, 0.8, np.nan])
    np.testing.assert_allclose(
        functional_with_change(state, np.array([0.5, np.nan, 0.5])), [0.4, 0.8, np.nan]
    )


def test_combine_mean_of_available_sensors():
    r, count = combine_ratios([np.array([0.8, np.nan, np.nan]), np.array([1.0, 0.9, np.nan])])
    np.testing.assert_allclose(r, [0.9, 0.9, np.nan])
    assert count.tolist() == [2, 1, 0]


def test_tolerance_h_floor_and_missing_inputs():
    x = np.tile(np.linspace(1, 2, 15), (3, 1))
    x[2, 4] = np.nan
    h = tolerance_h("TV", x, np.array([2, 0, 2]))
    assert h[0] >= 0.02 and np.isnan(h[1]) and np.isnan(h[2])
    with pytest.raises(ValueError):
        tolerance_h("Nc", x, np.array([1, 1, 1]))


def test_ratio_window_and_baseline_rules():
    x = np.full((3, len(YEARS)), 0.5)
    x[1, YEARS.index(2009) : YEARS.index(2012)] = np.nan  # 7 baseline years: ok
    x[2, YEARS.index(2009) : YEARS.index(2013)] = np.nan  # 6 valid of 10 ... still ok
    x[2, YEARS.index(2013)] = np.nan  # 5 valid: no baseline
    x[0, YEARS.index(2024) :] = np.nan  # window 2025: 1 valid year -> NaN
    med = baseline_median(x, YEARS, "T")
    assert np.isfinite(med[:2]).all() and np.isnan(med[2])
    r = change_ratio(x, YEARS, 2025, med)
    assert np.isnan(r[0]) and r[1] == pytest.approx(1.0)
