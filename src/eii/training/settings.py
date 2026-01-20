"""
Settings for NPP model training.

These settings control the training data generation and model training process.
"""

from ..client.assets import DATASETS_BASE, EII_ASSET_BASE

# Training asset paths
TRAINING_ASSET_PATH = f"{EII_ASSET_BASE}/training/samples"
FINAL_TRAINING_ASSET_PATH = f"{EII_ASSET_BASE}/training"
MODEL_ASSETS_PATH = f"{EII_ASSET_BASE}/models"
VALIDATION_ASSET_PATH = f"{EII_ASSET_BASE}/validation"
GLOBAL_GRID_ASSET_PATH = f"{EII_ASSET_BASE}/processing/grids"

# HMI natural areas mask collection (pre-computed raster tiles)
HMI_MASKS_COLLECTION = f"{DATASETS_BASE}/natural_lands/hmi_masks/v1"

# WDPA raster asset (pre-computed raster mask of high protection areas)
WDPA_RASTER_ASSET = f"{DATASETS_BASE}/protected_areas/high_protection_status_areas"

# Sampling settings
TRAINING_GRID_SIZE_DEG = 20  # Grid cell size for spatial sampling (tile generation)
CV_GRID_SIZE_DEG = 2  # Grid cell size for spatial cross-validation
CV_BUFFER_DEG = 0.5  # Buffer to exclude around spatial blocks (negative buffer)
SAMPLES_PER_ECOREGION = 500  # Target samples per ecoregion
TOTAL_POINTS_PER_GRID_CELL = 10000  # Total points per grid cell (area-proportional)

# Ecoregions
ECOREGIONS_ASSET = "RESOLVE/ECOREGIONS/2017"
ECOREGIONS_RASTER_ASSET = f"{DATASETS_BASE}/ecoregions/resolve_2017_raster"

# Model hyperparameters
RF_NUM_TREES = 200
RF_MIN_LEAF_POPULATION = 5
RF_VARIABLES_PER_SPLIT = 5
RF_BAG_FRACTION = 0.7
RF_SEED = 42

# Predictor settings
INCLUDE_LAT_LON_PREDICTORS = True  # Whether to include lat/lon as predictors

# Base predictors used for model training
PREDICTOR_VARIABLES = [
    "mean_annual_temp",
    "temp_seasonality",
    "annual_precip",
    "precip_seasonality",
    "aridity",
    "elevation",
    "slope",
    "tpi",
    "tri",
    "cti",
    "northness",
    "eastness",
    "sand",
    "clay",
    "ph",
    "chelsa_npp",
    "tpi_regional",
]

if INCLUDE_LAT_LON_PREDICTORS:
    PREDICTOR_VARIABLES.extend(["latitude", "longitude"])

# GEE performance tuning
TILE_SCALE = 4  # tileScale for sampling operations
SAMPLE_SCALE = 300  # Scale for sampling in meters

# Train/test split
TRAIN_TEST_SPLIT_RATIO = 0.9
