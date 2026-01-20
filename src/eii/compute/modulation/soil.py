"""
Soil KPI calculation based on soil organic carbon relative to climate potential.

Uses SoilGrids 250m for actual SOC and Köppen climate zones for
reference SOC values representing natural/undisturbed conditions.
"""

from __future__ import annotations

import ee

from .settings import (
    DEFAULT_SOC_REFERENCE,
    KOPPEN_ASSET,
    KOPPEN_CODE_TO_NAME,
    SOC_REFERENCE_BY_CLIMATE,
    SOIL_MIN_THRESHOLD,
    SOILGRIDS_SOC_ASSET,
)


def calculate_soil_kpi(
    aoi: ee.Geometry,
    depth: str = "0-30cm",
) -> ee.Image:
    """
    Calculate soil KPI as SOC relative to climate-typical potential.

    This function:
    1. Loads SoilGrids SOC at specified depth
    2. Loads Köppen climate zones
    3. Creates a reference SOC image from climate zone lookup
    4. Calculates actual/reference ratio with a minimum threshold and clamps to [0, 1]

    Args:
        aoi: Area of interest geometry
        depth: Soil depth for SOC extraction. Options: "0-5cm", "5-15cm",
            "15-30cm", "30-60cm", "60-100cm", "100-200cm", or "0-30cm" (default)

    Returns:
        ee.Image with soil KPI values (0-1) at 250m resolution

    Note:
        - SoilGrids SOC is stored as g/kg * 10, converted to g/kg here
        - Reference values represent typical SOC in undisturbed soils per climate zone
        - Ratio > 1.0 (well-managed or pristine soils) is clamped to 1.0
        - Ratio <= SOIL_MIN_THRESHOLD is mapped to 0

    Example:
        >>> aoi = ee.Geometry.Rectangle([10.5, 47.5, 11.0, 48.0])
        >>> soil_kpi = calculate_soil_kpi(aoi)
        >>> # Result is 0-1 normalized at 250m resolution
    """
    soc = _get_soc(depth)

    reference_soc = _get_reference_soc()

    ratio = soc.divide(reference_soc)
    relative_soc = ratio.subtract(SOIL_MIN_THRESHOLD).divide(1 - SOIL_MIN_THRESHOLD).clamp(0, 1)

    return relative_soc.rename("soil_kpi").clip(aoi)


def _get_soc(depth: str = "0-30cm") -> ee.Image:
    """
    Get soil organic carbon from SoilGrids.

    Args:
        depth: Soil depth layer to extract

    Returns:
        ee.Image with SOC in g/kg (global, unclipped)
    """
    soilgrids = ee.Image(SOILGRIDS_SOC_ASSET)

    depth_map = {
        "0-5cm": "soc_0-5cm_mean",
        "5-15cm": "soc_5-15cm_mean",
        "15-30cm": "soc_15-30cm_mean",
        "30-60cm": "soc_30-60cm_mean",
        "60-100cm": "soc_60-100cm_mean",
        "100-200cm": "soc_100-200cm_mean",
        "0-30cm": "soc_0-5cm_mean",
    }

    if depth == "0-30cm":
        soc_0_5 = soilgrids.select("soc_0-5cm_mean")
        soc_5_15 = soilgrids.select("soc_5-15cm_mean")
        soc_15_30 = soilgrids.select("soc_15-30cm_mean")

        soc = soc_0_5.multiply(5).add(soc_5_15.multiply(10)).add(soc_15_30.multiply(15)).divide(30)
    elif depth in depth_map:
        soc = soilgrids.select(depth_map[depth])
    else:
        raise ValueError(f"Unknown depth: {depth}. Use one of {list(depth_map.keys())}")

    soc = soc.divide(10).rename("soc")

    return soc


def _get_reference_soc() -> ee.Image:
    """
    Create reference SOC image from Köppen climate zones.

    Returns:
        ee.Image with reference SOC values (g/kg) based on climate zone (global, unclipped)
    """
    koppen = ee.Image(KOPPEN_ASSET).select("b1")

    codes = list(KOPPEN_CODE_TO_NAME.keys())
    reference_values = [
        SOC_REFERENCE_BY_CLIMATE.get(KOPPEN_CODE_TO_NAME[code], DEFAULT_SOC_REFERENCE)
        for code in codes
    ]

    reference_soc = koppen.remap(codes, reference_values, DEFAULT_SOC_REFERENCE)

    return reference_soc.rename("reference_soc")


def get_soc(
    aoi: ee.Geometry,
    depth: str = "0-30cm",
) -> ee.Image:
    """
    Get raw soil organic carbon without normalization.

    Useful for inspection or custom normalization.

    Args:
        aoi: Area of interest geometry
        depth: Soil depth for SOC extraction

    Returns:
        ee.Image with SOC in g/kg at 250m resolution
    """
    return _get_soc(depth).clip(aoi)


def get_reference_soc(aoi: ee.Geometry) -> ee.Image:
    """
    Get reference SOC image based on climate zones.

    Useful for inspection of reference values.

    Args:
        aoi: Area of interest geometry

    Returns:
        ee.Image with reference SOC values (g/kg)
    """
    return _get_reference_soc().clip(aoi)
