"""
Full EII computation module.

This module provides functions to calculate the Ecosystem Integrity Index
from scratch, including all component calculations.

Requires additional dependencies: pip install eii[compute]

Example:
    >>> import ee
    >>> from eii.compute import calculate_eii
    >>>
    >>> ee.Initialize()
    >>> polygon = ee.Geometry.Rectangle([-60, -10, -55, -5])
    >>> results = calculate_eii(polygon)

For Natural Capital modulation KPIs:
    >>> from eii.compute import modulation
    >>> biodiv = modulation.calculate_biodiversity_kpi(polygon)
"""

from . import modulation
from .compositional import calculate_compositional_integrity
from .integrity import (
    calculate_eii,
    combine_components,
)
from .npp import (
    calculate_functional_integrity,
    load_natural_npp,
    load_natural_npp_tiles,
    load_npp_diff_percentiles,
    setup_predictor_stack,
    setup_response,
)
from .structural import calculate_structural_integrity

__all__ = [
    # Core EII computation
    "calculate_eii",
    "combine_components",
    "calculate_functional_integrity",
    "calculate_structural_integrity",
    "calculate_compositional_integrity",
    # NPP utilities
    "load_npp_diff_percentiles",
    "load_natural_npp",
    "load_natural_npp_tiles",
    "setup_predictor_stack",
    "setup_response",
    # Natural Capital modulation
    "modulation",
]
