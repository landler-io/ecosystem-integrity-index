"""
Lightweight client for retrieving pre-computed EII data from GEE assets.

This module provides simple functions to access the Ecosystem Integrity Index
without needing to run the full computation pipeline.

Example:
    >>> import ee
    >>> from eii.client import get_stats
    >>>
    >>> ee.Initialize(project='your-project')
    >>> polygon = ee.Geometry.Rectangle([-60, -10, -55, -5])
    >>> stats = get_stats(polygon)

For Natural Capital modulated EII:
    >>> from eii.client import get_modulated_eii, get_default_kpis
    >>> result = get_modulated_eii(polygon, compute_default_kpis=True)
"""

from .analysis import compare_methods, get_zonal_stats
from .assets import ASSETS, get_asset_path
from .modulation import (
    get_default_kpis,
    get_kpi_layers,
    get_modulated_eii,
    get_nc_score,
)
from .retrieve import get_layers, get_raster, get_stats
from .utils import normalize_client_input, quicklook

__all__ = [
    # Core retrieval
    "get_stats",
    "get_raster",
    "get_layers",
    "get_zonal_stats",
    "compare_methods",
    # Natural Capital modulation
    "get_modulated_eii",
    "get_default_kpis",
    "get_kpi_layers",
    "get_nc_score",
    # Utilities
    "quicklook",
    "normalize_client_input",
    "ASSETS",
    "get_asset_path",
]
