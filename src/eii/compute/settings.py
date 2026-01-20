"""
Settings for EII computation.

Contains configuration for computing EII from scratch, including
model paths, resolution settings, and pre-computed decile breaks.
"""

from typing import Literal

# Spatial resolution for computation (meters)
SPATIAL_RESOLUTION = 300

# Predictor settings
INCLUDE_LAT_LON_PREDICTORS = True  # Whether to include lat/lon as predictors

# Default NPP model asset paths
DEFAULT_NPP_MODEL_PATH = (
    "projects/landler-open-data/assets/eii/models/potential_npp_mean_classifier"
)
DEFAULT_NPP_STD_MODEL_PATH = (
    "projects/landler-open-data/assets/eii/models/potential_npp_std_classifier"
)

# Pre-computed natural NPP
NATURAL_NPP_ASSET_PATH = (
    "projects/landler-open-data/assets/eii/intermediate/functional/predictions/natural_npp_v1"
)

# Observed NPP from CLMS (band 0 = sum, band 1 = std)
OBSERVED_NPP_ASSET_PATH = "projects/landler-open-data/assets/datasets/clms/npp/annual"

# NPP difference percentiles asset path (for functional integrity scoring)
# Computed by: pipelines/modeling/eii_calculation/compute_npp_decile_breaks.ipynb
NPP_DIFF_PERCENTILES_ASSET_PATH = (
    "projects/landler-open-data/assets/eii/intermediate/functional/npp_diff_decile_breaks"
)

# Default EII aggregation method
DEFAULT_AGGREGATION_METHOD: Literal[
    "minimum",
    "product",
    "min_fuzzy_logic",
    "geometric_mean",
] = "min_fuzzy_logic"

# Reference dates
NPP_YEAR_RANGE = ["2014-01-01", "2025-01-01"]  # For model training (historical)
OBSERVED_NPP_YEAR_RANGE = [
    "2022-01-01",
    "2025-01-01",
]  # 3-year rolling window for current state
EII_REFERENCE_DATE = "2024-12-31"
LANDCOVER_REFERENCE_DATE = "2023-01-01"

# Functional integrity weights
MAGNITUDE_WEIGHT = 2 / 3  # ~0.67
SEASONALITY_WEIGHT = 1 / 3  # ~0.33
