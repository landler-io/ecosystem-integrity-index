"""
Water KPI calculation based on available water capacity (AWC).

Uses SoilGrids sand/clay and SOC with the Saxton & Rawls (2006)
pedotransfer function to estimate field capacity and wilting point,
then calculates AWC, normalized relative to texture-typical maximum.
"""

from __future__ import annotations

import ee

from .settings import (
    AWC_BASE_MAX,
    AWC_BASE_MIN,
    AWC_CLAY_PENALTY,
    AWC_MAX_BOUND,
    AWC_MIN_BOUND,
    AWC_MIN_MAX_BOUND,
    AWC_MIN_MIN_BOUND,
    AWC_OPTIMAL_CLAY,
    AWC_OPTIMAL_SAND,
    AWC_SAND_PENALTY,
    SOILGRIDS_CLAY_ASSET,
    SOILGRIDS_SAND_ASSET,
    SOILGRIDS_SOC_ASSET,
)


def calculate_water_kpi(
    aoi: ee.Geometry,
    depth: str = "0-30cm",
) -> ee.Image:
    """
    Calculate water KPI as AWC relative to texture-typical maximum.

    This function:
    1. Calculates AWC using the Saxton & Rawls (2006) pedotransfer function
       (sand, clay, SOC -> SOM -> field capacity and wilting point)
    2. Estimates maximum AWC based on soil texture (sand/clay content)
    3. Calculates actual/maximum ratio, clamped to [0, 1]

    The texture-based normalization accounts for the fact that:
    - Loamy soils have highest AWC potential (~25 vol%)
    - Sandy soils have lower AWC potential (~10-15 vol%)
    - Clayey soils have moderate AWC potential (~15-20 vol%)

    Args:
        aoi: Area of interest geometry
        depth: Soil depth for AWC calculation. Options: "0-5cm", "5-15cm",
            "15-30cm", "30-60cm", "60-100cm", "100-200cm", or "0-30cm" (default)

    Returns:
        ee.Image with water KPI values (0-1) at 250m resolution

    Note:
        - AWC = Field Capacity (33 kPa) - Permanent Wilting Point (1500 kPa)
        - SOC is converted to SOM using som = soc * 1.9 (Saxton & Rawls)
        - AWC is returned in vol% (m3/m3 * 100)
        - High values indicate soil is performing well for its texture class

    Example:
        >>> aoi = ee.Geometry.Rectangle([10.5, 47.5, 11.0, 48.0])
        >>> water_kpi = calculate_water_kpi(aoi)
        >>> # Result is 0-1 normalized at 250m resolution
    """
    awc = _get_awc(depth)

    max_awc = _get_max_awc_by_texture(depth)
    min_awc = _get_min_awc_by_texture(depth)

    normalized = awc.subtract(min_awc).divide(max_awc.subtract(min_awc)).clamp(0, 1)

    return normalized.rename("water_kpi").clip(aoi)


def _get_awc(depth: str = "0-30cm") -> ee.Image:
    """
    Calculate Available Water Capacity from SoilGrids using Saxton & Rawls (2006).

    Args:
        depth: Soil depth layer

    Returns:
        ee.Image with AWC in vol% (global, unclipped)
    """
    sand_source = ee.Image(SOILGRIDS_SAND_ASSET)
    clay_source = ee.Image(SOILGRIDS_CLAY_ASSET)
    soc_source = ee.Image(SOILGRIDS_SOC_ASSET)

    if depth == "0-30cm":
        sand_0_5 = sand_source.select("sand_0-5cm_mean")
        sand_5_15 = sand_source.select("sand_5-15cm_mean")
        sand_15_30 = sand_source.select("sand_15-30cm_mean")
        sand = (
            sand_0_5.multiply(5).add(sand_5_15.multiply(10)).add(sand_15_30.multiply(15)).divide(30)
        )

        clay_0_5 = clay_source.select("clay_0-5cm_mean")
        clay_5_15 = clay_source.select("clay_5-15cm_mean")
        clay_15_30 = clay_source.select("clay_15-30cm_mean")
        clay = (
            clay_0_5.multiply(5).add(clay_5_15.multiply(10)).add(clay_15_30.multiply(15)).divide(30)
        )

        soc_0_5 = soc_source.select("soc_0-5cm_mean")
        soc_5_15 = soc_source.select("soc_5-15cm_mean")
        soc_15_30 = soc_source.select("soc_15-30cm_mean")
        soc = soc_0_5.multiply(5).add(soc_5_15.multiply(10)).add(soc_15_30.multiply(15)).divide(30)
    else:
        depth_suffix = depth.replace("cm", "cm_mean")
        sand = sand_source.select(f"sand_{depth_suffix}")
        clay = clay_source.select(f"clay_{depth_suffix}")
        soc = soc_source.select(f"soc_{depth_suffix}")

    sand = sand.divide(10)
    clay = clay.divide(10)
    soc = soc.divide(10)

    sand_frac = sand.divide(100)
    clay_frac = clay.divide(100)
    soc_percent = soc.divide(10)
    som_percent = soc_percent.multiply(1.9)

    theta_33t = (
        sand_frac.multiply(-0.251)
        .add(clay_frac.multiply(0.195))
        .add(som_percent.multiply(0.011))
        .add(sand_frac.multiply(som_percent).multiply(0.006))
        .subtract(clay_frac.multiply(som_percent).multiply(0.027))
        .add(sand_frac.multiply(clay_frac).multiply(0.452))
        .add(0.299)
    )
    theta_33 = theta_33t.add(
        theta_33t.pow(2).multiply(1.283).subtract(theta_33t.multiply(0.374)).subtract(0.015)
    )

    theta_1500t = (
        sand_frac.multiply(-0.024)
        .add(clay_frac.multiply(0.487))
        .add(som_percent.multiply(0.006))
        .add(sand_frac.multiply(som_percent).multiply(0.005))
        .subtract(clay_frac.multiply(som_percent).multiply(0.013))
        .add(sand_frac.multiply(clay_frac).multiply(0.068))
        .add(0.031)
    )
    theta_1500 = theta_1500t.add(theta_1500t.multiply(0.14).subtract(0.02))

    awc = theta_33.subtract(theta_1500).multiply(100).rename("awc")

    return awc


def _get_min_awc_by_texture(depth: str = "0-30cm") -> ee.Image:
    """
    Estimate minimum AWC based on soil texture.

    This provides a texture-specific lower bound for normalization.
    """
    sand_source = ee.Image(SOILGRIDS_SAND_ASSET)
    clay_source = ee.Image(SOILGRIDS_CLAY_ASSET)

    if depth == "0-30cm":
        sand_0_5 = sand_source.select("sand_0-5cm_mean")
        sand_5_15 = sand_source.select("sand_5-15cm_mean")
        sand_15_30 = sand_source.select("sand_15-30cm_mean")
        sand = (
            sand_0_5.multiply(5).add(sand_5_15.multiply(10)).add(sand_15_30.multiply(15)).divide(30)
        )

        clay_0_5 = clay_source.select("clay_0-5cm_mean")
        clay_5_15 = clay_source.select("clay_5-15cm_mean")
        clay_15_30 = clay_source.select("clay_15-30cm_mean")
        clay = (
            clay_0_5.multiply(5).add(clay_5_15.multiply(10)).add(clay_15_30.multiply(15)).divide(30)
        )
    else:
        depth_suffix = depth.replace("cm", "cm_mean")
        sand = sand_source.select(f"sand_{depth_suffix}")
        clay = clay_source.select(f"clay_{depth_suffix}")

    sand = sand.divide(10)
    clay = clay.divide(10)

    sand_penalty = sand.subtract(AWC_OPTIMAL_SAND).pow(2).multiply(AWC_SAND_PENALTY)
    clay_penalty = clay.subtract(AWC_OPTIMAL_CLAY).pow(2).multiply(AWC_CLAY_PENALTY)

    min_awc = (
        ee.Image(AWC_BASE_MIN)
        .subtract(sand_penalty)
        .subtract(clay_penalty)
        .clamp(AWC_MIN_MIN_BOUND, AWC_MIN_MAX_BOUND)
    )

    return min_awc.rename("min_awc")


def _get_max_awc_by_texture(depth: str = "0-30cm") -> ee.Image:
    """
    Estimate maximum AWC based on soil texture.

    Uses a parabolic function that peaks around loam texture
    (optimal sand ~35%, optimal clay ~25%) and decreases for
    sandy or clayey extremes.

    Args:
        depth: Soil depth layer

    Returns:
        ee.Image with maximum AWC (vol%) for each pixel's texture (global, unclipped)
    """
    sand_source = ee.Image(SOILGRIDS_SAND_ASSET)
    clay_source = ee.Image(SOILGRIDS_CLAY_ASSET)

    if depth == "0-30cm":
        sand_0_5 = sand_source.select("sand_0-5cm_mean")
        sand_5_15 = sand_source.select("sand_5-15cm_mean")
        sand_15_30 = sand_source.select("sand_15-30cm_mean")
        sand = (
            sand_0_5.multiply(5).add(sand_5_15.multiply(10)).add(sand_15_30.multiply(15)).divide(30)
        )

        clay_0_5 = clay_source.select("clay_0-5cm_mean")
        clay_5_15 = clay_source.select("clay_5-15cm_mean")
        clay_15_30 = clay_source.select("clay_15-30cm_mean")
        clay = (
            clay_0_5.multiply(5).add(clay_5_15.multiply(10)).add(clay_15_30.multiply(15)).divide(30)
        )
    else:
        depth_suffix = depth.replace("cm", "cm_mean")
        sand = sand_source.select(f"sand_{depth_suffix}")
        clay = clay_source.select(f"clay_{depth_suffix}")

    sand = sand.divide(10)
    clay = clay.divide(10)

    sand_penalty = sand.subtract(AWC_OPTIMAL_SAND).pow(2).multiply(AWC_SAND_PENALTY)
    clay_penalty = clay.subtract(AWC_OPTIMAL_CLAY).pow(2).multiply(AWC_CLAY_PENALTY)

    max_awc = (
        ee.Image(AWC_BASE_MAX)
        .subtract(sand_penalty)
        .subtract(clay_penalty)
        .clamp(AWC_MIN_BOUND, AWC_MAX_BOUND)
    )

    return max_awc.rename("max_awc")


def get_awc(
    aoi: ee.Geometry,
    depth: str = "0-30cm",
) -> ee.Image:
    """
    Get raw Available Water Capacity without normalization.

    Useful for inspection or custom normalization.

    Args:
        aoi: Area of interest geometry
        depth: Soil depth for AWC calculation

    Returns:
        ee.Image with AWC in vol% at 250m resolution
    """
    return _get_awc(depth).clip(aoi)


def get_max_awc_by_texture(
    aoi: ee.Geometry,
    depth: str = "0-30cm",
) -> ee.Image:
    """
    Get texture-based maximum AWC.

    Useful for inspection of reference values.

    Args:
        aoi: Area of interest geometry
        depth: Soil depth for texture extraction

    Returns:
        ee.Image with maximum AWC (vol%) based on texture
    """
    return _get_max_awc_by_texture(depth).clip(aoi)
