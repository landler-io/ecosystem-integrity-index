#!/usr/bin/env python3
"""
Generate HMI-based natural areas mask and export to Earth Engine assets.

Usage:
    python calculate_hmi_natural_areas.py
"""

import sys
from pathlib import Path

import ee

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.config_utils import load_config  # noqa: E402
from utils.utils_gee_assets import make_asset_dir  # noqa: E402


def get_natural_areas_mask(
    hmi_threshold: float,
    high_mod_threshold: float,
    buffer_pixels: int,
    aoi: ee.Geometry | None = None,
) -> ee.Image:
    """Generate binary natural areas mask from HMI time-series."""
    hm_collection = ee.ImageCollection(
        "projects/sat-io/open-datasets/GHM/HM_1990_2020_OVERALL_300M"
    )
    hmi_series = hm_collection.filter(ee.Filter.eq("threat_code", "AA"))

    if aoi:
        buffer_m = buffer_pixels * 300
        processing_geom = aoi.buffer(buffer_m)
        hmi_series = hmi_series.filterBounds(processing_geom)

    hmi_series = hmi_series.map(lambda img: img.select([0]).rename("hmi"))
    hmi_max = hmi_series.max()

    if aoi:
        hmi_max = hmi_max.clip(processing_geom)

    # Gap-fill with 5x5 focal max
    kernel_5x5 = ee.Kernel.square(radius=2, units="pixels")
    hmi_filled = hmi_max.focalMax(kernel=kernel_5x5)

    # Buffer high-modification areas
    high_mod_mask = hmi_filled.gt(high_mod_threshold)
    buffer_kernel = ee.Kernel.circle(radius=buffer_pixels, units="pixels")
    high_mod_buffered = high_mod_mask.focalMax(kernel=buffer_kernel)

    # Create mask: low HMI AND outside buffer
    low_hmi_mask = hmi_filled.lt(hmi_threshold)
    outside_buffer = high_mod_buffered.eq(0)
    natural_mask = low_hmi_mask.And(outside_buffer).rename("natural_mask")

    if aoi:
        natural_mask = natural_mask.clip(aoi)

    return natural_mask.set(
        {
            "hmi_threshold": hmi_threshold,
            "high_mod_threshold": high_mod_threshold,
            "buffer_pixels": buffer_pixels,
        }
    )


def create_grid_tiles(tile_size_deg: float) -> list[dict]:
    """Create global grid of tiles."""
    tiles = []
    lat = -90
    while lat < 90:
        lon = -180
        while lon < 180:
            geom = ee.Geometry.Rectangle([lon, lat, lon + tile_size_deg, lat + tile_size_deg])
            tile_id = f"lat_{lat:+07.2f}_lon_{lon:+08.2f}".replace(".", "_").replace("+", "_")
            tiles.append({"id": tile_id, "geometry": geom})
            lon += tile_size_deg
        lat += tile_size_deg
    return tiles


def export_tile(
    tile: dict,
    config: dict,
) -> ee.batch.Task:
    """Export a single tile mask to EE asset."""
    mask = get_natural_areas_mask(
        hmi_threshold=config["hmi_threshold"],
        high_mod_threshold=config["high_mod_threshold"],
        buffer_pixels=config["buffer_pixels"],
        aoi=tile["geometry"],
    )

    # Add all metadata properties (matching upload.sh)
    mask = mask.set(
        {
            "tile_id": tile["id"],
            "min_area_km2": config.get("min_area_km2", 0),
        }
    )

    asset_id = f"{config['asset_collection']}/{tile['id']}"

    task = ee.batch.Export.image.toAsset(
        image=mask.toInt8(),
        description=tile["id"],
        assetId=asset_id,
        region=tile["geometry"],
        scale=config["export_scale"],
        maxPixels=1e12,
    )
    task.start()
    return task


def main():
    ee.Initialize()
    config = load_config(Path(__file__).parent / "config.cfg")

    # Ensure parent folder exists and create collection
    collection_path = config["asset_collection"]
    make_asset_dir("/".join(collection_path.split("/")[:-1]))
    try:
        ee.data.getAsset(collection_path)
    except ee.EEException:
        ee.data.createAsset({"type": "ImageCollection"}, collection_path)

    tiles = create_grid_tiles(config["grid_tile_size_deg"])

    start_idx = config.get("start_tile_index") or 0
    end_idx = config.get("end_tile_index") or len(tiles)
    tiles = tiles[start_idx:end_idx]

    for tile in tiles:
        export_tile(tile, config)


if __name__ == "__main__":
    main()
