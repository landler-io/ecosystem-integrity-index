"""
Functional integrity: observed NPP relative to the modelled natural potential.

All functions take numpy arrays in gC m-2 yr-1 with NaN as nodata and return float32
arrays with NaN where integrity is not defined.
"""

from __future__ import annotations

import numpy as np

from .settings import (
    N_BIOMES,
    NPP_STORED_SCALE,
    OBSERVED_NPP_FIRST_YEAR,
    OBSERVED_NPP_MIN_VALID_YEARS,
    OBSERVED_NPP_WINDOW_YEARS,
    OBSERVED_NPP_YEARS,
    POTENTIAL_NPP_FLOOR,
    RECENTRE_MIN_REFERENCE_PIXELS,
    RECENTRE_RATIO_BINS,
    RECENTRE_RATIO_MAX,
    SURPLUS_ZERO,
    TOLERANCE_K,
    TOLERANCE_K_SURPLUS,
    WETLAND_FLAG_CLASSES,
    WETLAND_FLAG_MIN_COVER_PCT,
)


def observed_npp_window(reference_year: int) -> list[int]:
    """Return the annual NPP years averaged for one reference year.

    Args:
        reference_year: Last calendar year of the window.

    Returns:
        The ``OBSERVED_NPP_WINDOW_YEARS`` years ending in ``reference_year``, without years
        before ``OBSERVED_NPP_FIRST_YEAR``.

    Example:
        >>> observed_npp_window(2024)
        [2022, 2023, 2024]
        >>> observed_npp_window(2020)
        [2019, 2020]
    """
    first = OBSERVED_NPP_FIRST_YEAR + OBSERVED_NPP_MIN_VALID_YEARS - 1
    if (
        isinstance(reference_year, bool)
        or not isinstance(reference_year, int)
        or not first <= reference_year <= OBSERVED_NPP_YEARS[1]
    ):
        raise ValueError(
            f"reference_year must be an integer from {first} to {OBSERVED_NPP_YEARS[1]}"
        )
    start = max(reference_year - OBSERVED_NPP_WINDOW_YEARS + 1, OBSERVED_NPP_FIRST_YEAR)
    return list(range(start, reference_year + 1))


def observed_npp(
    annual: np.ndarray, min_valid_years: int = OBSERVED_NPP_MIN_VALID_YEARS
) -> np.ndarray:
    """Mean annual NPP in gC m-2 yr-1 from stored CLMS annual sums.

    Args:
        annual: Stored annual sums, shape ``(years, ...)``. Negative values (the -1
            sentinel) and NaN are nodata.
        min_valid_years: Minimum valid years per pixel; fewer gives NaN.

    Returns:
        float32 mean of the valid years times ``NPP_STORED_SCALE``.
    """
    values = np.asarray(annual, dtype="float64")
    valid = values >= 0
    count = valid.sum(axis=0)
    total = np.where(valid, values, 0).sum(axis=0)
    mean = total / np.maximum(count, 1) * NPP_STORED_SCALE
    return np.where(count >= min_valid_years, mean, np.nan).astype("float32")


def deviation_tolerance_integrity(
    actual: np.ndarray,
    potential: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    k: float = TOLERANCE_K,
    floor: float = POTENTIAL_NPP_FLOOR,
    natural_cv: np.ndarray | None = None,
    k_surplus: float | None = TOLERANCE_K_SURPLUS,
    surplus_zero: float = SURPLUS_ZERO,
) -> np.ndarray:
    """Functional integrity as a deficit or surplus beyond a tolerance.

    Deviations are relative to ``D = max(P, floor)``: ``d = |A - P| / D``. The relative
    prediction uncertainty is ``u = (P - lower) / D`` for a deficit and ``(upper - P) / D``
    for a surplus; with the optional modelled natural inter-annual CV ``v`` the tolerance
    is ``t = k sqrt(u^2 + v^2)`` (``t = k u`` without ``natural_cv``), with ``k`` for a
    deficit and ``k_surplus`` for a surplus. Integrity is 1 for ``d <= t`` and falls
    linearly to 0 at ``d = z``: ``z = 1`` for a deficit (total loss) and
    ``z = surplus_zero`` for a surplus (with 1, a surplus of ``D``). The floor keeps deserts
    scored: where ``P`` is near 0, a surplus of ``surplus_zero * floor`` (e.g. irrigation)
    scores 0.

    Args:
        actual: Observed NPP A.
        potential: Potential NPP P.
        lower: Lower bound of the potential-NPP prediction interval (68 % by default).
        upper: Upper bound of the same interval.
        k: Deficit tolerance multiplier.
        floor: Smallest denominator of the relative deviation, in gC m-2 yr-1.
        natural_cv: Natural inter-annual coefficient of variation of NPP; NaN or None
            leaves only the prediction uncertainty.
        k_surplus: Surplus tolerance multiplier; None uses ``k`` on both sides.
        surplus_zero: Relative surplus scoring 0; must be positive.

    Returns:
        float32 integrity in [0, 1], NaN where A or P is missing.

    Example:
        >>> a, p = np.array([70.0, 140.0]), np.array(100.0)
        >>> deviation_tolerance_integrity(a, p, np.array(80.0), np.array(120.0)).round(3)
        array([0.875, 0.941], dtype=float32)
    """
    scale = np.maximum(potential, floor)
    surplus = actual > potential
    deviation = np.abs(actual - potential) / scale
    relative = np.where(surplus, upper - potential, potential - lower) / scale
    if natural_cv is not None:
        relative = np.hypot(relative, np.nan_to_num(natural_cv, nan=0.0))
    if surplus_zero <= 0:
        raise ValueError("surplus_zero must be positive")
    multiplier = np.where(surplus, k if k_surplus is None else k_surplus, k)
    zero = np.where(surplus, surplus_zero, 1.0)
    # Capped below the zero point so that a deviation of z still scores 0 for very wide
    # intervals.
    tolerance = np.clip(multiplier * np.nan_to_num(relative, nan=0.0), 0, 0.999 * zero)
    score: np.ndarray = np.clip(1 - (deviation - tolerance) / (zero - tolerance), 0, 1)
    return score.astype("float32")


def recentring_sums(
    actual: np.ndarray, potential: np.ndarray, biome: np.ndarray, reference: np.ndarray
) -> np.ndarray:
    """Per-biome sums of observed and potential NPP on reference pixels, binned by A / P.

    Sums from separate tiles add up, so a global run accumulates them over all tiles of
    one observed window and passes the total to ``recentring_factors``. Sums are binned by
    A / P (width ``RECENTRE_RATIO_MAX / RECENTRE_RATIO_BINS``, plus one bin for ratios at or
    above ``RECENTRE_RATIO_MAX`` and for ``P = 0``).

    Args:
        actual: Observed NPP A for the window.
        potential: Potential NPP P, same shape.
        biome: RESOLVE biome number per pixel (1-14); 0, NaN or other values are unknown
            and never counted.
        reference: True on reference pixels (natural land chosen independently of the
            run's structural and compositional scores).

    Returns:
        float64 array ``(N_BIOMES, RECENTRE_RATIO_BINS + 1, 3)``: pixel count, sum of A and
        sum of P per biome and ratio bin.
    """
    index = np.nan_to_num(np.asarray(biome, dtype="float64"), nan=0).astype("int64")
    ok = np.asarray(reference, dtype=bool) & np.isfinite(actual) & np.isfinite(potential)
    ok &= (index > 0) & (index < N_BIOMES)
    a = np.asarray(actual, dtype="float64")[ok]
    p = np.asarray(potential, dtype="float64")[ok]
    ratio = np.divide(a, p, out=np.where(a > 0, np.inf, 0.0), where=p > 0)
    scaled = np.clip(ratio, 0, RECENTRE_RATIO_MAX) / RECENTRE_RATIO_MAX * RECENTRE_RATIO_BINS
    bins = np.minimum(scaled.astype("int64"), RECENTRE_RATIO_BINS)
    flat = index[ok] * (RECENTRE_RATIO_BINS + 1) + bins
    size = N_BIOMES * (RECENTRE_RATIO_BINS + 1)
    sums = [np.bincount(flat, weights=w, minlength=size) for w in (None, a, p)]
    return np.stack(sums, axis=-1).astype("float64").reshape(N_BIOMES, -1, 3)


def recentring_factors(
    sums: np.ndarray, min_pixels: int = RECENTRE_MIN_REFERENCE_PIXELS
) -> np.ndarray:
    """Recentring factor per biome: sum of observed over sum of potential NPP.

    A ratio of means, not a median of ratios: about half of the near-pristine desert
    pixels observe NPP 0.

    Args:
        sums: Output of ``recentring_sums``, summed over tiles.
        min_pixels: Biomes with this many reference pixels or fewer keep factor 1.

    Returns:
        float64 factors indexed by biome number (length ``N_BIOMES``; index 0 is 1).

    Example:
        >>> a, p = np.full(300, 120.0), np.full(300, 100.0)
        >>> sums = recentring_sums(a, p, np.full(300, 6), np.ones(300, bool))
        >>> float(recentring_factors(sums)[6].round(6))
        1.2
    """
    total, observed, potential = sums.sum(axis=1).T
    use = (total > min_pixels) & (potential > 0)
    factors = np.ones(len(sums))
    factors[use] = observed[use] / potential[use]
    return factors


def recentre_potential(
    potential: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    biome: np.ndarray,
    factors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Scale potential NPP and its interval bounds by the factor of each pixel's biome.

    Args:
        potential: Potential NPP P.
        lower: Lower bound of its prediction interval.
        upper: Upper bound of the same interval.
        biome: RESOLVE biome number per pixel; unknown biomes (0, NaN, out of range) keep
            factor 1.
        factors: Factors indexed by biome number, from ``recentring_factors``.

    Returns:
        float32 recentred potential, lower and upper bound.
    """
    index = np.nan_to_num(np.asarray(biome, dtype="float64"), nan=0).astype("int64")
    lookup = np.append(np.asarray(factors, dtype="float64"), 1.0)
    lookup[0] = 1.0
    index = np.where((index > 0) & (index < len(factors)), index, len(factors))
    factor = lookup[index]
    return tuple(  # type: ignore[return-value]
        (np.asarray(x) * factor).astype("float32") for x in (potential, lower, upper)
    )


def wetland_flag(glwd_class: np.ndarray, glwd_cover: np.ndarray) -> np.ndarray:
    """Flag GLWD v2 wetlands where potential NPP is known to be overestimated.

    Codes follow ``WETLAND_FLAG_CLASSES``: 0 none, 1 wetland, 2 peatland, 3 mangrove.
    A pixel is flagged where its dominant GLWD class belongs to a family and the GLWD
    cover reaches that family's ``WETLAND_FLAG_MIN_COVER_PCT``. The flag annotates
    functional integrity and never changes its values.

    Args:
        glwd_class: GLWD v2 dominant class (0-33; NaN or other values are no wetland).
        glwd_cover: GLWD v2 combined cover in percent of the cell.

    Returns:
        uint8 flag codes.
    """
    code = np.zeros(256, dtype="uint8")
    minimum = np.full(256, np.inf)
    for flag, classes in WETLAND_FLAG_CLASSES.items():
        code[list(classes)] = flag
        minimum[list(classes)] = WETLAND_FLAG_MIN_COVER_PCT[flag]
    index = np.nan_to_num(glwd_class, nan=0).astype("int64")
    index = np.where((index >= 0) & (index < 256), index, 0)
    with np.errstate(invalid="ignore"):
        flagged = np.nan_to_num(glwd_cover, nan=-1) >= minimum[index]
    return np.where(flagged, code[index], 0).astype("uint8")
