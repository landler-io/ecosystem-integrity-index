import logging
from pathlib import Path

import ee
import yaml
from IPython.core.getipython import get_ipython

ipython = get_ipython()
if ipython:
    ipython.run_line_magic("load_ext", "autoreload")
    ipython.run_line_magic("autoreload", "2")

from eii.utils import create_assets_folder  # noqa: E402

# Load project-level config (asset paths, GEE project)
try:
    project_config_path = Path(__file__).parent / "../../../config.yaml"
    with open(project_config_path) as f:
        project_config = yaml.safe_load(f)
except FileNotFoundError:
    print("Warning: project config.yaml not found.")
    raise

# Load pipeline-specific config (structural integrity parameters)
try:
    pipeline_config_path = Path(__file__).parent / "config.yaml"
    with open(pipeline_config_path) as f:
        pipeline_config = yaml.safe_load(f)
except FileNotFoundError:
    print("Warning: pipeline config.yaml not found, using defaults.")
    pipeline_config = {}

# Project settings
GEE_PROJECT = project_config["gee"]["project"]
GEE_ASSET_ROOT = project_config["gee"]["asset_root"]
GEE_MODEL_VERSION = project_config["gee"]["version"]

# Asset paths
INTERMEDIATE_ROOT = f"{GEE_ASSET_ROOT}/{project_config['structure']['intermediate']}"
PRODUCTS_ROOT = f"{GEE_ASSET_ROOT}/{project_config['structure']['products']}/{GEE_MODEL_VERSION}"

STRUCTURAL_ASSET_PATH = f"{PRODUCTS_ROOT}/structural_integrity"
GLOBAL_GRID_ASSET_PATH = f"{INTERMEDIATE_ROOT}/{project_config['modules']['functional']}/{project_config['folders']['grid']}"

FINAL_PRODUCT_ASSET_PATH = f"{PRODUCTS_ROOT}/structural_integrity/structural_integrity_core_area"

# Pipeline parameters (from pipeline config)
EDGE_DEPTH_M = pipeline_config.get("edge_depth_m", 300)
NEIGHBORHOOD_M = pipeline_config.get("neighborhood_m", 5000)
HMI_THRESHOLD = pipeline_config.get("hmi_threshold", 0.4)
SCALE_M = pipeline_config.get("scale_m", 300)
TILE_SIZE_DEG = pipeline_config.get("tile_size_deg", 10.0)
MIN_LAT = pipeline_config.get("min_lat", -60)
MAX_LAT = pipeline_config.get("max_lat", 80)
MAX_TILES_PER_BATCH = pipeline_config.get("max_tiles_per_batch", 3000)

print(f"Asset Root: {GEE_ASSET_ROOT}")
print(f"Model Version: {GEE_MODEL_VERSION}")

logging.getLogger("ee").setLevel(logging.INFO)

try:
    ee.Initialize(project=GEE_PROJECT)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=GEE_PROJECT)

folders_to_create = [
    STRUCTURAL_ASSET_PATH,
    FINAL_PRODUCT_ASSET_PATH,
]

for folder in folders_to_create:
    create_assets_folder(folder)
