"""Offline tests of the numpy EII v2 calculation."""

from __future__ import annotations

import numpy as np
import pytest

from eii.compute import (
    combine_components,
    compositional_integrity,
    deviation_tolerance_integrity,
    limiting_component,
    observed_npp,
    observed_npp_window,
    recentre_potential,
    recentring_factors,
    recentring_sums,
    structural_integrity,
    wetland_flag,
)
from eii.compute.structural import (
    exponential_kernel,
    ground_disk,
    habitat_quality,
    landscape_quality,
)

GHM_PIXEL = 0.0026949458523585646  # Earth Engine "300 m" on EPSG:4326


def test_observed_npp_masks_sentinel_and_counts_valid_years():
    annual = np.array([[5e6, 5e6, 0.0], [5e6, 5e6, -1], [5e6, -1, -1]]).T
    assert observed_npp(annual).tolist() == pytest.approx([1000 / 3, 500, np.nan], nan_ok=True)
    assert np.isnan(observed_npp(annual, min_valid_years=3)[1])


def test_observed_npp_window():
    assert observed_npp_window(2024) == [2022, 2023, 2024]
    assert observed_npp_window(2021) == [2019, 2020, 2021]
    assert observed_npp_window(2020) == [2019, 2020]
    for bad in (2015, 2019, 2026, 2024.0, True):
        with pytest.raises(ValueError):
            observed_npp_window(bad)


@pytest.mark.parametrize(
    ("actual", "expected"),
    # P = 100, interval 80-110: tolerance 0.2 below, 0.1 above; 0 at |A - P| = P.
    [(100, 1.0), (85, 1.0), (80, 1.0), (60, 0.75), (0, 0.0), (110, 1.0), (155, 0.5), (250, 0.0)],
)
def test_deviation_tolerance_penalises_deficit_and_surplus(actual, expected):
    score = deviation_tolerance_integrity(
        np.array(actual),
        np.array(100.0),
        np.array(80.0),
        np.array(110.0),
        k_surplus=None,
        surplus_zero=1.0,
    )
    assert float(score) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("actual", "expected"),
    # Defaults, P = 100, interval 80-110: deficit tolerance 1.0 x 0.2, surplus 1.5 x 0.1,
    # surplus scoring 0 at 2 x D (three times the potential).
    [(80, 1.0), (60, 0.75), (115, 1.0), (155, 1 - 0.4 / 1.85), (250, 1 - 1.35 / 1.85), (300, 0.0)],
)
def test_deviation_tolerance_defaults_widen_only_the_surplus_side(actual, expected):
    args = (np.array(actual), np.array(100.0), np.array(80.0), np.array(110.0))
    assert float(deviation_tolerance_integrity(*args)) == pytest.approx(expected, rel=1e-6)
    assert float(deviation_tolerance_integrity(*args, k=1.0, k_surplus=1.5)) == pytest.approx(
        expected, rel=1e-6
    )


def test_deviation_tolerance_surplus_k_leaves_deficits_alone():
    actual = np.array([40.0, 60.0, 140.0, 170.0])
    args = (actual, np.array(100.0), np.array(70.0), np.array(130.0))
    cv = np.array(0.4)
    kw = {"natural_cv": cv, "surplus_zero": 1.0}
    same = deviation_tolerance_integrity(*args, k=1.0, k_surplus=None, **kw)
    assert np.array_equal(same, deviation_tolerance_integrity(*args, k_surplus=1.0, **kw))
    wide = deviation_tolerance_integrity(*args, k=1.0, k_surplus=1.5, **kw)
    assert np.array_equal(wide[:2], same[:2])
    # surplus: u = 0.3, v = 0.4 -> t = 0.75; A = 170 is d = 0.7 -> inside the tolerance
    assert wide[3] == 1 and same[3] == pytest.approx(1 - 0.2 / 0.5)


def test_deviation_tolerance_scores_deserts_on_the_floor_and_masks_missing():
    # P = 10 is below the 50 gC floor: deviations are relative to 50; a surplus of 2 x 50
    # scores 0.
    lower, upper = np.array(0.0), np.array(15.0)
    desert = deviation_tolerance_integrity(np.array([0.0, 110.0]), np.array(10.0), lower, upper)
    assert desert == pytest.approx([1.0, 0.0])
    missing = deviation_tolerance_integrity(
        np.array([np.nan, 100.0]), np.array([100.0, np.nan]), lower, upper
    )
    assert np.isnan(missing).all()
    # A zero lower bound gives the widest tolerance; total loss still scores 0.
    assert (
        deviation_tolerance_integrity(
            np.array(0.0), np.array(100.0), np.array(0.0), np.array(200.0)
        )
        == 0
    )


def test_deviation_tolerance_combines_prediction_and_interannual_uncertainty():
    # u = 0.3, v = 0.4 -> t = 0.5; A = 40 is d = 0.6 -> 1 - 0.1 / 0.5 = 0.8
    args = (np.array([60.0, 40.0]), np.array(100.0), np.array(70.0), np.array(130.0))
    assert deviation_tolerance_integrity(*args, natural_cv=np.array(0.4)) == pytest.approx([1, 0.8])
    assert deviation_tolerance_integrity(*args, natural_cv=np.array(np.nan)) == pytest.approx(
        deviation_tolerance_integrity(*args)
    )


def test_recentring_factor_is_a_ratio_of_means_per_biome():
    n = 600
    biome = np.repeat([6, 13], n // 2)
    potential = np.full(n, 100.0)
    actual = np.where(biome == 6, 120.0, 0.0)
    actual[n // 2 : n // 2 + 120] = 200.0  # 40 % of the desert observes 200, the rest 0
    factors = recentring_factors(recentring_sums(actual, potential, biome, np.ones(n, bool)))
    assert factors.shape == (15,)
    assert factors[6] == pytest.approx(1.2)
    assert factors[13] == pytest.approx(0.8)  # mean A / mean P, not the median ratio (0)
    assert np.all(np.delete(factors, [6, 13]) == 1)


def test_recentring_sums_add_up_over_tiles_and_skip_non_reference_pixels():
    rng = np.random.default_rng(42)
    n = 2000
    actual, potential = rng.uniform(0, 900, n), rng.uniform(0, 600, n)
    potential[:5] = 0  # A / P undefined: counted in the overflow bin
    actual[5:10] = np.nan
    biome = rng.integers(0, 16, n).astype("float64")
    biome[10:15] = np.nan
    reference = rng.random(n) < 0.8
    whole = recentring_sums(actual, potential, biome, reference)
    parts = sum(
        recentring_sums(actual[s], potential[s], biome[s], reference[s])
        for s in (slice(0, 700), slice(700, None))
    )
    assert np.allclose(whole, parts)
    counted = reference & np.isfinite(actual) & (biome >= 1) & (biome <= 14)
    assert whole[..., 0].sum() == counted.sum()
    assert whole[..., 1].sum() == pytest.approx(actual[counted].sum())
    assert whole[0].sum() == 0


def test_recentring_keeps_factor_one_for_small_biomes():
    sums = recentring_sums(np.full(200, 150.0), np.full(200, 100.0), np.full(200, 4), True)
    assert recentring_factors(sums)[4] == 1
    assert recentring_factors(sums, min_pixels=199)[4] == pytest.approx(1.5)


def test_recentring_is_a_ratio_of_means():
    actual = np.concatenate([np.full(960, 100.0), np.zeros(20), np.full(20, 1000.0)])
    sums = recentring_sums(actual, np.full(1000, 100.0), np.full(1000, 1), True)
    assert recentring_factors(sums)[1] == pytest.approx((96000 + 20000) / 100000)


def test_recentre_potential_scales_interval_and_keeps_unknown_biomes():
    factors = np.ones(15)
    factors[6], factors[0] = 1.2, 5.0  # biome 0 (unknown) is never scaled
    biome = np.array([6, 0, np.nan, 99, 13])
    p, lo, hi = recentre_potential(
        np.full(5, 100.0), np.full(5, 80.0), np.full(5, 130.0), biome, factors
    )
    assert p.dtype == np.float32
    assert p.tolist() == pytest.approx([120, 100, 100, 100, 100])
    assert lo.tolist() == pytest.approx([96, 80, 80, 80, 80])
    assert hi.tolist() == pytest.approx([156, 130, 130, 130, 130])


def test_recentring_removes_a_biome_offset_but_keeps_pixel_loss():
    n = 1000
    potential, lower, upper = np.full(n, 400.0), np.full(n, 360.0), np.full(n, 440.0)
    actual = np.full(n, 480.0)  # the whole biome sits 20 % above potential in this window
    actual[:10] = 240.0  # 50 % below the window's level
    biome, reference = np.full(n, 6), np.ones(n, bool)
    plain = deviation_tolerance_integrity(actual, potential, lower, upper, k_surplus=1.0)
    factors = recentring_factors(recentring_sums(actual, potential, biome, reference))
    recentred = deviation_tolerance_integrity(
        actual, *recentre_potential(potential, lower, upper, biome, factors)
    )
    assert plain[10:].max() < 1 and np.all(recentred[10:] == 1)
    assert np.all(recentred[:10] < plain[:10]) and recentred[0] == pytest.approx(0.558, abs=1e-3)


def test_wetland_flag_classes_and_cover_thresholds():
    classes = np.array([28, 28, 22, 22, 10, 5, 33, np.nan])
    cover = np.array([60, 40, 33, 32, 51, 100, 100, 100])
    assert wetland_flag(classes, cover).tolist() == [3, 0, 2, 0, 1, 0, 0, 0]


def test_ground_disk_matches_earth_engine_pixels_and_widens_with_latitude():
    assert ground_disk(300, 0, GHM_PIXEL).astype(int).tolist() == [
        [0, 1, 0],
        [1, 1, 1],
        [0, 1, 0],
    ]
    disk = ground_disk(5000, 60, GHM_PIXEL)
    assert disk.shape == (33, 67)
    assert disk.sum() == pytest.approx(ground_disk(5000, 0, GHM_PIXEL).sum() * 2, rel=0.02)


def test_exponential_kernel_decays_and_widens_with_latitude():
    k0 = exponential_kernel(26500, 5000, 0, GHM_PIXEL)
    assert k0.shape == ground_disk(26500, 0, GHM_PIXEL).shape
    centre = k0.shape[0] // 2
    assert k0[centre, centre] == 1
    step = GHM_PIXEL * 111319.49079327357
    assert k0[centre, centre + 10] == pytest.approx(np.exp(-10 * step / 5000), rel=1e-5)
    k60 = exponential_kernel(26500, 5000, 60, GHM_PIXEL)
    assert k60.sum() == pytest.approx(k0.sum() * 2, rel=0.02)


@pytest.mark.parametrize("hmi", [0.0, 0.15, 0.35, 0.5, 0.9])
def test_structural_uniform_landscape_scores_its_habitat_quality(hmi):
    out = structural_integrity(np.full((60, 60), hmi), top_lat=0.1, pixel_deg=GHM_PIXEL)
    assert np.allclose(out, habitat_quality(hmi), atol=1e-5)


def test_habitat_quality_is_a_rescaled_logistic_through_one_half_at_the_midpoint():
    h = np.linspace(0, 1, 101)
    q = habitat_quality(h)
    assert q[0] == pytest.approx(1) and q[-1] == pytest.approx(0)
    assert np.all(np.diff(q) < 0)
    assert habitat_quality(0.4) == pytest.approx(0.5, abs=0.01)
    assert habitat_quality(0.1) > 1 - 0.1 and habitat_quality(0.6) < 1 - 0.6


def test_structural_ramp_has_no_plateaus_and_is_monotone():
    hmi = np.tile(np.linspace(0, 1, 400), (40, 1))
    out = landscape_quality(hmi, top_lat=0.1, pixel_deg=GHM_PIXEL)[20]
    assert np.all(np.diff(out) < 0)


def test_structural_islet_adds_no_landscape_quality_but_keeps_its_site_term():
    hmi = np.ones((300, 300))
    hmi[100:140, 100:140] = 0.0  # a 12 km block of natural land
    without = landscape_quality(hmi, top_lat=0.1, pixel_deg=GHM_PIXEL)
    hmi[150, 150] = 0.0  # plus a single natural pixel in fully modified land
    out = landscape_quality(hmi, top_lat=0.1, pixel_deg=GHM_PIXEL)
    no_open = landscape_quality(hmi, top_lat=0.1, pixel_deg=GHM_PIXEL, opening_m=0)
    others = np.ones(hmi.shape, bool)
    others[150, 150] = False
    assert np.allclose(out[others], without[others], atol=1e-6)  # no landscape contribution
    assert not np.allclose(no_open[others], without[others], atol=1e-6)
    assert 0 < out[150, 150] < no_open[150, 150]  # its site term sees the block, not itself
    assert out[200, 200] == 0  # modified land scores 0 whatever its surroundings


def test_structural_sea_neither_counts_nor_opens_coastal_land():
    hmi = np.zeros((120, 120))
    hmi[:, :60] = np.nan  # sea
    hmi[:, 60] = 0.0
    out = structural_integrity(hmi, top_lat=0.1, pixel_deg=GHM_PIXEL)
    assert np.isnan(out[:, :60]).all()
    assert np.allclose(out[:, 60:], 1, atol=1e-5)


def test_structural_matches_a_direct_sum():
    rng = np.random.default_rng(3)
    hmi = rng.uniform(0, 1, (90, 90))
    hmi[:10, :10] = np.nan
    out = landscape_quality(hmi, top_lat=0.1, pixel_deg=GHM_PIXEL, radius_m=3000, opening_m=0)
    kernel = exponential_kernel(3000, 5000, 0.1 - 45 * GHM_PIXEL, GHM_PIXEL)
    r = kernel.shape[0] // 2
    i, j = 45, 50
    patch = hmi[i - r : i + r + 1, j - r : j + r + 1]
    valid = np.isfinite(patch)
    quality = np.where(valid, habitat_quality(np.nan_to_num(patch)), 0)
    landscape = (kernel * np.sqrt(quality)).sum() / kernel[valid].sum()
    assert out[i, j] == pytest.approx(np.sqrt(habitat_quality(hmi[i, j])) * landscape, rel=1e-4)


def test_combine_components_methods_and_missing_component():
    f, s, c = np.array([0.2, 0.5]), np.array([0.5, np.nan]), np.array([0.8, 1.0])
    expected = {
        "minimum": 0.2,
        "product": 0.08,
        "geometric_mean": 0.08 ** (1 / 3),
        "min_fuzzy_logic": 0.2 * (0.5 + 0.8 - 0.4),
    }
    for method, value in expected.items():
        eii = combine_components(f, s, c, method=method)
        assert eii[0] == pytest.approx(value, rel=1e-6) and np.isnan(eii[1])
    with pytest.raises(ValueError):
        combine_components(f, s, c, method="mean")  # type: ignore[arg-type]


def test_limiting_component_codes_and_ties():
    f = np.array([0.2, 0.5, 0.9, 0.5, np.nan])
    s = np.array([0.5, 0.2, 0.9, 0.5, 0.1])
    c = np.array([0.8, 0.9, 0.1, 0.5, 0.1])
    assert limiting_component(f, s, c).tolist() == [1, 2, 3, 1, 0]


def test_compositional_clamps_and_masks_off_land():
    out = compositional_integrity(np.array([-0.1, 0.5, 1.2, 0.0]), np.array([1, 1, 1, 0], bool))
    assert out[:3].tolist() == [0, 0.5, 1] and np.isnan(out[3])


@pytest.mark.parametrize(
    ("actual", "expected"),
    # P = 100, interval 80-110, k_surplus 1.5: surplus tolerance 0.15, zero at 3 x D.
    [(115, 1.0), (200, 1 - 0.85 / 2.85), (300, 1 - 1.85 / 2.85), (400, 0.0), (60, 0.75)],
)
def test_deviation_tolerance_surplus_zero_softens_only_the_surplus_ramp(actual, expected):
    args = (np.array(actual), np.array(100.0), np.array(80.0), np.array(110.0))
    score = deviation_tolerance_integrity(*args, surplus_zero=3.0)
    assert float(score) == pytest.approx(expected, rel=1e-6)


def test_deviation_tolerance_surplus_zero_must_be_positive():
    with pytest.raises(ValueError):
        deviation_tolerance_integrity(
            np.array(1.0), np.array(1.0), np.array(1.0), np.array(1.0), surplus_zero=0
        )
