#!/usr/bin/env python3
"""
Rasterize RESOLVE Ecoregions 2017 to a GEE Image asset.

This pre-rasterization improves performance when ecoregions are used
for stratified sampling, avoiding repeated reduceToImage() calls.

Usage:
    python rasterize_ecoregions.py
"""

import sys
from pathlib import Path

import ee

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.config_utils import load_config  # noqa: E402
from utils.utils_gee_assets import make_asset_dir  # noqa: E402


def rasterize_ecoregions(config: dict) -> ee.batch.Task:
    """Rasterize RESOLVE Ecoregions 2017 to a GEE Image asset."""
    asset_id = config["asset_id"]
    make_asset_dir("/".join(asset_id.split("/")[:-1]))

    raw_ecoregions = ee.FeatureCollection(config["source_asset"])

    # --- FIX 1: Sanitize Data Type ---
    # Ensure ECO_ID is cast to a Number.
    # 'reduceToImage' will silently mask pixels if the property is a String.
    def sanitize(feature):
        return feature.set("numeric_eco_id", ee.Number(feature.get("ECO_ID")))

    ecoregions = raw_ecoregions.map(sanitize)

    # --- FIX 2: Rasterize ---
    # Use the sanitized numeric property
    ecoregion_raster = (
        ecoregions.reduceToImage(properties=["numeric_eco_id"], reducer=ee.Reducer.first())
        .rename("ecoregion")
        .toInt()
        .unmask(0)  # masked (oceans) -> 0
    )

    # --- FIX 3: Correct Geometry Definition ---
    # Use 'geodesic=False' to prevent -180/180 from collapsing into a zero-width line.
    # Alternatively, use [-179.999, -60, 179.999, 90]
    global_bounds = ee.Geometry.Rectangle([-180, -60, 180, 90], "EPSG:4326", False)

    # Do not .clip() here. Let the Export task handle the region cropping.
    # Clipping an infinite image to a complex geometry can sometimes cause artifacts.

    task = ee.batch.Export.image.toAsset(
        image=ecoregion_raster,
        description="Rasterize_Ecoregions",
        assetId=config["asset_id"],
        region=global_bounds,
        scale=config["scale"],  # 5000 meters
        crs="EPSG:4326",
        maxPixels=1e13,
    )
    task.start()
    return task


def main():
    ee.Initialize()
    config = load_config(Path(__file__).parent / "config.cfg")
    task = rasterize_ecoregions(config)
    print(f"Export task started: {task.id}")


if __name__ == "__main__":
    main()
