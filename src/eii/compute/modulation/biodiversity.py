"""
Biodiversity KPI calculation based on natural habitat fraction.

Uses ESA WorldCover 10m land cover classification to calculate the
fraction of natural and semi-natural habitats within an area of interest.
"""

from __future__ import annotations

import ee

from .settings import (
    BIODIVERSITY_AGGREGATION_SCALE,
    BIODIVERSITY_MAX_THRESHOLD,
    BIODIVERSITY_MIN_THRESHOLD,
    NATURAL_LANDCOVER_CLASSES,
    WORLDCOVER_ASSET,
)


def calculate_biodiversity_kpi(
    aoi: ee.Geometry,
    natural_classes: list[int] | None = None,
    aggregation_scale: int = BIODIVERSITY_AGGREGATION_SCALE,
    min_threshold: float = BIODIVERSITY_MIN_THRESHOLD,
    max_threshold: float = BIODIVERSITY_MAX_THRESHOLD,
) -> ee.Image:
    """
    Calculate biodiversity KPI as natural habitat fraction.

    This function:
    1. Loads ESA WorldCover 10m land cover classification
    2. Creates a binary mask of natural/semi-natural classes
    3. Aggregates to specified resolution (default 100m = 1ha) to get fraction
    4. Normalizes using min-max thresholds

    Args:
        aoi: Area of interest geometry
        natural_classes: List of WorldCover class codes considered "natural".
            Default includes forests, shrubland, grassland, wetlands, etc.
        aggregation_scale: Resolution for fraction calculation in meters.
            Default 100m (1 hectare grid cells).
        min_threshold: Lower bound for normalization (default 0.0 = 0%)
        max_threshold: Upper bound for normalization (default 0.5 = 50%)

    Returns:
        ee.Image with biodiversity KPI values (0-1) at aggregation_scale resolution

    Note:
        WorldCover classes:
        - 10: Tree cover (includes plantations)
        - 20: Shrubland
        - 30: Grassland
        - 40: Cropland (excluded)
        - 50: Built-up (excluded)
        - 60: Bare/sparse vegetation
        - 70: Snow and ice
        - 80: Permanent water bodies
        - 90: Herbaceous wetland
        - 95: Mangroves
        - 100: Moss and lichen

    Example:
        >>> aoi = ee.Geometry.Rectangle([10.5, 47.5, 11.0, 48.0])
        >>> biodiv_kpi = calculate_biodiversity_kpi(aoi)
        >>> # Result is 0-1 normalized at 100m resolution
    """
    if natural_classes is None:
        natural_classes = NATURAL_LANDCOVER_CLASSES

    worldcover = ee.ImageCollection(WORLDCOVER_ASSET).first().select("Map")

    natural_mask = worldcover.remap(
        natural_classes,
        [1] * len(natural_classes),
        defaultValue=0,
    )

    natural_fraction = (
        natural_mask.reduceResolution(
            reducer=ee.Reducer.mean(),
            maxPixels=1024,
        )
        .reproject(crs="EPSG:4326", scale=aggregation_scale)
        .rename("natural_fraction")
    )

    normalized = (
        natural_fraction.subtract(min_threshold).divide(max_threshold - min_threshold).clamp(0, 1)
    )

    return normalized.rename("biodiversity_kpi").clip(aoi)


def get_natural_habitat_fraction(
    aoi: ee.Geometry,
    natural_classes: list[int] | None = None,
    aggregation_scale: int = BIODIVERSITY_AGGREGATION_SCALE,
) -> ee.Image:
    """
    Get raw natural habitat fraction without normalization.

    Useful for inspection or custom normalization.

    Args:
        aoi: Area of interest geometry
        natural_classes: List of WorldCover class codes considered "natural"
        aggregation_scale: Resolution for fraction calculation in meters

    Returns:
        ee.Image with natural habitat fraction (0-1) at aggregation_scale resolution
    """
    if natural_classes is None:
        natural_classes = NATURAL_LANDCOVER_CLASSES

    worldcover = ee.ImageCollection(WORLDCOVER_ASSET).first().select("Map")

    natural_mask = worldcover.remap(
        natural_classes,
        [1] * len(natural_classes),
        defaultValue=0,
    )

    natural_fraction = (
        natural_mask.reduceResolution(
            reducer=ee.Reducer.mean(),
            maxPixels=1024,
        )
        .reproject(crs="EPSG:4326", scale=aggregation_scale)
        .rename("natural_fraction")
    )

    return natural_fraction.clip(aoi)
