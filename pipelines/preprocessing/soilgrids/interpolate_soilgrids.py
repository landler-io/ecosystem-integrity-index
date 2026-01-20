"""
SoilGrids gap-filling with Gaussian convolution.
Fills NA pixels using weighted average from valid neighbors
"""

import sys
from pathlib import Path

import ee

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.config_utils import load_config  # noqa: E402
from utils.utils_gee_assets import make_asset_dir  # noqa: E402

# --- CONFIGURATION ---
SCRIPT_DIR = Path(__file__).parent
config = load_config(SCRIPT_DIR / "config.cfg")

GEE_PROJECT = config["GEE_PROJECT"]
OUTPUT_ASSET_PREFIX = config["OUTPUT_ASSET_PREFIX"]
SEARCH_DISTANCE_M = config["SEARCH_DISTANCE_M"]
SOILGRIDS_RESOLUTION = config["SOILGRIDS_RESOLUTION"]

variables_raw = config["VARIABLES"]
VARIABLES = []
for item in variables_raw.split(";"):
    item = item.strip()
    if not item:
        continue
    parts = item.split(",")
    if len(parts) == 2:
        VARIABLES.append((parts[0].strip(), parts[1].strip()))
# --- END CONFIGURATION ---


def fill_gaps_gaussian(image, kernel_radius, sigma):
    valid_mask = image.mask()

    kernel = ee.Kernel.gaussian(radius=kernel_radius, sigma=sigma, units="pixels", normalize=False)

    image_zeros = image.unmask(0)
    weighted_sum = image_zeros.multiply(valid_mask.unmask(0)).convolve(kernel)
    weight_sum = valid_mask.unmask(0).convolve(kernel)
    interpolated = weighted_sum.divide(weight_sum.max(1e-10))

    filled = image.unmask(interpolated)
    has_valid_neighbors = weight_sum.gt(0)
    return filled.updateMask(valid_mask.Or(has_valid_neighbors))


def process_variable(var_name, depth, kernel_radius, sigma):
    band_name = f"{var_name}_{depth}_mean"
    make_asset_dir(f"{OUTPUT_ASSET_PREFIX}/{var_name}")

    image = ee.Image(f"projects/soilgrids-isric/{var_name}_mean").select(band_name)

    filled = fill_gaps_gaussian(image, kernel_radius, sigma).rename(band_name)

    print(f"{var_name}_{depth}_gapfill")
    task = ee.batch.Export.image.toAsset(
        image=filled.toInt16(),
        description=f"{var_name}_{depth}_gapfill",
        assetId=f"{OUTPUT_ASSET_PREFIX}/{var_name}/{var_name}_{depth}_mean_gapfilled",
        crs="EPSG:4326",
        scale=SOILGRIDS_RESOLUTION,
        maxPixels=1e13,
        pyramidingPolicy={band_name: "mean"},
    )
    return task


def main():
    ee.Initialize(project=GEE_PROJECT)

    sigma = (SEARCH_DISTANCE_M / SOILGRIDS_RESOLUTION) / 3
    radius = int(SEARCH_DISTANCE_M / SOILGRIDS_RESOLUTION)

    for var_name, depth in VARIABLES:
        task = process_variable(var_name, depth, radius, sigma)
        task.start()


if __name__ == "__main__":
    main()
