"""
Core modulation logic for Natural Capital adjustment of EII.

This module provides functions to:
1. Calculate the Natural Capital (NC) score from KPI components
2. Apply the NC modulation to base EII scores
"""

from __future__ import annotations

import ee

from .settings import DEFAULT_KPI_WEIGHTS, MODULATION_RANGE


def calculate_nc_score(
    biodiversity_kpi: float,
    soil_kpi: float,
    water_kpi: float,
    weights: dict[str, float] | None = None,
) -> float:
    """
    Calculate the Natural Capital score from component KPIs.

    The NC score is a weighted average of the three KPI dimensions.
    Default weighting is equal (1/3 each).

    Args:
        biodiversity_kpi: Normalized biodiversity KPI (0-1)
        soil_kpi: Normalized soil KPI (0-1)
        water_kpi: Normalized water KPI (0-1)
        weights: Optional custom weights. Keys must be "biodiversity", "soil", "water".
            Values must sum to 1.0.

    Returns:
        Natural Capital score (0-1)

    Raises:
        ValueError: If weights don't sum to 1.0 (within tolerance)

    Example:
        >>> nc = calculate_nc_score(0.45, 0.68, 0.52)
        >>> print(f"NC Score: {nc:.3f}")
        NC Score: 0.550
    """
    kpis = {
        "biodiversity": biodiversity_kpi,
        "soil": soil_kpi,
        "water": water_kpi,
    }

    if weights is None:
        weights = DEFAULT_KPI_WEIGHTS

    weight_sum = sum(weights.get(k, 0) for k in kpis)
    if abs(weight_sum - 1.0) > 0.001:
        raise ValueError(f"Weights must sum to 1.0, got {weight_sum:.4f}")

    nc_score = sum(kpis[k] * weights.get(k, 0) for k in kpis)

    return nc_score


def calculate_nc_score_image(
    biodiversity_kpi: ee.Image,
    soil_kpi: ee.Image,
    water_kpi: ee.Image,
    weights: dict[str, float] | None = None,
) -> ee.Image:
    """
    Calculate the Natural Capital score from KPI raster layers.

    Args:
        biodiversity_kpi: Normalized biodiversity KPI image (0-1)
        soil_kpi: Normalized soil KPI image (0-1)
        water_kpi: Normalized water KPI image (0-1)
        weights: Optional custom weights (must sum to 1.0)

    Returns:
        ee.Image with NC score (0-1)
    """
    if weights is None:
        weights = DEFAULT_KPI_WEIGHTS

    weight_sum = sum(weights.values())
    if abs(weight_sum - 1.0) > 0.001:
        raise ValueError(f"Weights must sum to 1.0, got {weight_sum:.4f}")

    nc_score = (
        biodiversity_kpi.multiply(weights["biodiversity"])
        .add(soil_kpi.multiply(weights["soil"]))
        .add(water_kpi.multiply(weights["water"]))
    )

    return nc_score.rename("nc_score")


def apply_modulation(
    eii: float,
    nc_score: float,
    modulation_range: float = MODULATION_RANGE,
) -> float:
    """
    Apply Natural Capital modulation to a base EII score.

    The modulation formula shifts EII by up to ±(modulation_range/2)
    based on the NC score:
    - NC = 0 -> EII decreases by modulation_range/2 (default: -0.05)
    - NC = 0.5 -> EII unchanged
    - NC = 1 -> EII increases by modulation_range/2 (default: +0.05)

    Formula: EII_final = EII + (NC - 0.5) * modulation_range

    Args:
        eii: Base EII score (0-1)
        nc_score: Natural Capital score (0-1)
        modulation_range: Total range of modulation (default 0.1 = ±0.05)

    Returns:
        Modulated EII score, clamped to [0, 1]

    Example:
        >>> eii_base = 0.532
        >>> nc = 0.62
        >>> eii_mod = apply_modulation(eii_base, nc)
        >>> print(f"Base: {eii_base:.3f}, Modulated: {eii_mod:.3f}")
        Base: 0.532, Modulated: 0.544
    """
    modulation = (nc_score - 0.5) * modulation_range
    eii_final = eii + modulation

    return max(0.0, min(1.0, eii_final))


def apply_modulation_image(
    eii: ee.Image,
    nc_score: ee.Image,
    modulation_range: float = MODULATION_RANGE,
) -> ee.Image:
    """
    Apply Natural Capital modulation to an EII raster.

    Args:
        eii: Base EII image (0-1)
        nc_score: Natural Capital score image (0-1)
        modulation_range: Total range of modulation (default 0.1)

    Returns:
        ee.Image with modulated EII score (0-1), clamped
    """
    modulation = nc_score.subtract(0.5).multiply(modulation_range)
    eii_final = eii.add(modulation).clamp(0, 1)

    return eii_final.rename("eii_modulated")
