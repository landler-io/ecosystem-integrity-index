"""
Natural Capital modulation for EII.

This module provides functions to calculate the Natural Capital (NC) score
from three KPI dimensions (biodiversity, soil, water) and apply it as a
modulation to the base EII score.

The NC modulation shifts EII by up to ±0.05 based on plot-level performance:
- NC = 0 -> EII decreases by 0.05
- NC = 0.5 -> EII unchanged
- NC = 1 -> EII increases by 0.05

Example:
    >>> from eii.compute.modulation import (
    ...     calculate_biodiversity_kpi,
    ...     calculate_soil_kpi,
    ...     calculate_water_kpi,
    ...     calculate_nc_score_image,
    ...     apply_modulation_image,
    ... )
    >>> import ee
    >>> ee.Initialize()
    >>> aoi = ee.Geometry.Rectangle([10.5, 47.5, 11.0, 48.0])
    >>> biodiv = calculate_biodiversity_kpi(aoi)
    >>> soil = calculate_soil_kpi(aoi)
    >>> water = calculate_water_kpi(aoi)
    >>> nc = calculate_nc_score_image(biodiv, soil, water)
    >>> eii_modulated = apply_modulation_image(eii_base, nc)
"""

from .biodiversity import (
    calculate_biodiversity_kpi,
    get_natural_habitat_fraction,
)
from .core import (
    apply_modulation,
    apply_modulation_image,
    calculate_nc_score,
    calculate_nc_score_image,
)
from .settings import (
    BIODIVERSITY_MAX_THRESHOLD,
    DEFAULT_KPI_WEIGHTS,
    MODULATION_RANGE,
)
from .soil import (
    calculate_soil_kpi,
    get_reference_soc,
    get_soc,
)
from .water import (
    calculate_water_kpi,
    get_awc,
    get_max_awc_by_texture,
)

__all__ = [
    # Core modulation functions
    "calculate_nc_score",
    "calculate_nc_score_image",
    "apply_modulation",
    "apply_modulation_image",
    # Biodiversity KPI
    "calculate_biodiversity_kpi",
    "get_natural_habitat_fraction",
    # Soil KPI
    "calculate_soil_kpi",
    "get_soc",
    "get_reference_soc",
    # Water KPI
    "calculate_water_kpi",
    "get_awc",
    "get_max_awc_by_texture",
    # Settings
    "DEFAULT_KPI_WEIGHTS",
    "MODULATION_RANGE",
    "BIODIVERSITY_MAX_THRESHOLD",
]
