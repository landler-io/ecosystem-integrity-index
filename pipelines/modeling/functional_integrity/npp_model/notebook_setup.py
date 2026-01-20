import logging
from pathlib import Path

import ee
import yaml
from IPython import get_ipython

from eii.utils import create_assets_folder

ipython = get_ipython()
if ipython:
    ipython.run_line_magic("load_ext", "autoreload")
    ipython.run_line_magic("autoreload", "2")


try:
    config_path = Path(__file__).parent / "../../../../config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    print("Configuration loaded from config.yaml")
except FileNotFoundError:
    print("Warning: config.yaml not found. Please ensure it exists in the project root.")
    raise

GEE_PROJECT = config["gee"]["project"]
GEE_ASSET_ROOT = config["gee"]["asset_root"]
GEE_MODEL_VERSION = config["gee"]["version"]

INTERMEDIATE_ROOT = f"{GEE_ASSET_ROOT}/{config['structure']['intermediate']}"
PRODUCTS_ROOT = f"{GEE_ASSET_ROOT}/{config['structure']['products']}/{GEE_MODEL_VERSION}"

FUNCTIONAL_PATH = f"{INTERMEDIATE_ROOT}/{config['modules']['functional']}"
TRAINING_BASE_PATH = f"{FUNCTIONAL_PATH}/{config['folders']['training']}"
TRAINING_TILES_PATH = f"{TRAINING_BASE_PATH}/tiles"
TRAINING_ASSET_PATH = TRAINING_BASE_PATH
FINAL_TRAINING_ASSET_PATH = TRAINING_BASE_PATH

MODEL_ASSETS_PATH = f"{FUNCTIONAL_PATH}/{config['folders']['models']}"
NPP_PREDICTION_ASSET_PATH = f"{FUNCTIONAL_PATH}/{config['folders']['predictions']}"
GLOBAL_GRID_ASSET_PATH = f"{FUNCTIONAL_PATH}/{config['folders']['grid']}"
VALIDATION_ASSET_PATH = f"{FUNCTIONAL_PATH}/{config['folders']['validation']}"

FINAL_PRODUCT_ASSET_PATH = f"{PRODUCTS_ROOT}/functional_integrity"

print(f"Asset Root: {GEE_ASSET_ROOT}")
print(f"Model Version: {GEE_MODEL_VERSION}")
print(f"Training Tiles: {TRAINING_TILES_PATH}")
print(f"Model Assets: {MODEL_ASSETS_PATH}")

logging.getLogger("ee").setLevel(logging.INFO)

try:
    ee.Initialize(project=GEE_PROJECT)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=GEE_PROJECT)

folders_to_create = [
    GEE_ASSET_ROOT,
    INTERMEDIATE_ROOT,
    FUNCTIONAL_PATH,
    TRAINING_BASE_PATH,
    TRAINING_TILES_PATH,
    MODEL_ASSETS_PATH,
    NPP_PREDICTION_ASSET_PATH,
    GLOBAL_GRID_ASSET_PATH,
    VALIDATION_ASSET_PATH,
    PRODUCTS_ROOT,
    FINAL_PRODUCT_ASSET_PATH,
]

for folder in folders_to_create:
    create_assets_folder(folder)
