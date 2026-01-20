"""
Combined Ecosystem Integrity Index calculation.
"""

from __future__ import annotations

from typing import Literal

import ee

from .settings import (
    DEFAULT_AGGREGATION_METHOD,
    OBSERVED_NPP_YEAR_RANGE,
)


def calculate_eii(
    aoi: ee.Geometry,
    method: Literal[
        "minimum", "product", "min_fuzzy_logic", "geometric_mean"
    ] = DEFAULT_AGGREGATION_METHOD,
    year_range: list[str] | None = None,
    include_seasonality: bool = True,
) -> dict[str, ee.Image]:
    """
    Calculate the full Ecosystem Integrity Index for an area.

    Args:
        aoi: Earth Engine geometry defining the area of interest.
        method: Method to combine components. Options:
            - "minimum": Simple minimum of components (most conservative)
            - "product": Product of components
            - "min_fuzzy_logic": Minimum with fuzzy compensation (default)
            - "geometric_mean": Geometric mean of components
        year_range: Date range for actual NPP [start_date, end_date].
            Defaults to OBSERVED_NPP_YEAR_RANGE.
        include_seasonality: Include seasonality in functional integrity.

    Returns:
        Dictionary containing functional_integrity, structural_integrity,
        compositional_integrity, and combined eii images.
    """
    from .compositional import calculate_compositional_integrity
    from .npp import calculate_functional_integrity
    from .structural import calculate_structural_integrity

    if year_range is None:
        year_range = OBSERVED_NPP_YEAR_RANGE

    # Buffer AOI for edge calculations
    buffered_aoi = aoi.buffer(20000).bounds()

    # Calculate components
    functional_result = calculate_functional_integrity(
        aoi=buffered_aoi,
        year_range=year_range,
        include_seasonality=include_seasonality,
    )
    functional = functional_result["functional_integrity"]

    structural = calculate_structural_integrity(aoi=buffered_aoi)

    compositional = calculate_compositional_integrity(aoi=buffered_aoi)

    # Combine components
    eii = combine_components(
        functional=functional,
        structural=structural,
        compositional=compositional,
        method=method,
    )

    # Clip to original AOI
    return {
        "functional_integrity": functional.clip(aoi),
        "structural_integrity": structural.clip(aoi),
        "compositional_integrity": compositional.clip(aoi),
        "eii": eii.clip(aoi),
    }


def combine_components(
    functional: ee.Image,
    structural: ee.Image,
    compositional: ee.Image,
    method: Literal[
        "minimum", "product", "min_fuzzy_logic", "geometric_mean"
    ] = DEFAULT_AGGREGATION_METHOD,
) -> ee.Image:
    """
    Combine integrity components into a single EII score.

    Args:
        functional: Functional integrity image (0-1).
        structural: Structural integrity image (0-1).
        compositional: Compositional integrity image (0-1).
        method: Aggregation method.

    Returns:
        Combined EII image (0-1).
    """
    components = ee.Image.cat([functional, structural, compositional])

    if method == "minimum":
        eii = components.reduce(ee.Reducer.min()).rename("eii")

    elif method == "product":
        eii = components.reduce(ee.Reducer.product()).rename("eii")

    elif method == "min_fuzzy_logic":
        # Fuzzy logic combination:
        # EII = M * F, where M = min, F = fuzzy sum of other two
        M = components.reduce(ee.Reducer.min())
        Med = components.reduce(ee.Reducer.median())
        Max = components.reduce(ee.Reducer.max())

        # Fuzzy sum: F = y + z - y*z
        F = Med.add(Max).subtract(Med.multiply(Max))

        eii = M.multiply(F).rename("eii")

    elif method == "geometric_mean":
        # Geometric mean = (a * b * c)^(1/n) where n=3 for three components
        product = components.reduce(ee.Reducer.product())
        eii = product.pow(1.0 / 3.0).rename("eii")

    else:
        raise ValueError(f"Unknown method '{method}'")

    return eii
