"""
Inference utilities for NPP model prediction.
"""

from __future__ import annotations

import datetime

import ee

from .._utils.gee import wait_for_tasks
from .grid import load_precomputed_grid
from .npp import setup_predictor_stack, setup_response


def predict_npp_tiled_batch_optimized(
    aoi: ee.Geometry | None = None,
    tile_size_deg: float = 5.0,
    model_version: str = "v10",
    resolution: int = 300,
    max_tiles_per_batch: int = 20,
    min_lat: float | None = None,
    max_lat: float | None = None,
    model_asset_path_mean: str | None = None,
    model_asset_path_std: str | None = None,
    tiled_folder: str | None = None,
    grid_asset_path: str | None = None,
) -> dict | None:
    """
    OPTIMIZED VERSION: Run NPP predictions using pre-computed grid with land/ocean classification.

    Predicts both Potential NPP Mean and Potential NPP Seasonality (Std Dev) to calculate
    Functional Integrity (Z-score).

    Args:
        aoi: Area of Interest (optional).
        tile_size_deg: Tile size in degrees.
        model_version: Version string.
        resolution: Spatial resolution in meters.
        max_tiles_per_batch: Batch size for exports.
        min_lat: Minimum latitude filter.
        max_lat: Maximum latitude filter.
        model_asset_path_mean: Full asset path for the Mean NPP model.
        model_asset_path_std: Full asset path for the Seasonality NPP model.
        tiled_folder: Folder to export tiles to.
        grid_asset_path: Path to pre-computed grid asset.
    """

    print("=== OPTIMIZED NPP Tiled Batch Prediction Pipeline ===")
    print(f"Tile size: {tile_size_deg}°")

    if (
        model_asset_path_mean is None
        or model_asset_path_std is None
        or tiled_folder is None
        or grid_asset_path is None
    ):
        raise ValueError(
            "model_asset_path_mean, model_asset_path_std, tiled_folder, and grid_asset_path must be provided"
        )

    # Check/Create tiled predictions folder
    try:
        ee.data.getAsset(tiled_folder)
        print(f"Using existing tiled folder: {tiled_folder}")
    except ee.EEException:
        print(f"Creating tiled folder: {tiled_folder}")
        ee.data.createAsset({"type": "Folder"}, tiled_folder)

    print("Loading trained models...")
    try:
        model_mean = ee.Classifier.load(model_asset_path_mean)
        model_std = ee.Classifier.load(model_asset_path_std)
    except Exception as e:
        print(f"Error loading models: {e}")
        return None

    print("Loading tiles from pre-computed grid...")
    tiles = load_precomputed_grid(
        tile_size_deg=tile_size_deg,
        land_tiles_only=True,
        min_lat=min_lat,
        max_lat=max_lat,
        grid_asset_path=grid_asset_path,
    )

    if not tiles:
        print(
            "Pre-computed grid not found. Need to run create_global_grid_with_land_classification first."
        )
        return None

    # Filter tiles by AOI if provided
    if aoi is not None:
        print("Filtering tiles by AOI...")
        aoi_bounds = aoi.bounds().getInfo()["coordinates"][0]
        min_lon = min([coord[0] for coord in aoi_bounds])
        max_lon = max([coord[0] for coord in aoi_bounds])
        min_lat = min([coord[1] for coord in aoi_bounds])
        max_lat = max([coord[1] for coord in aoi_bounds])

        # Filter tiles that intersect with AOI. Tile 'bounds' are just [minlon, minlat, maxlon, maxlat] usually
        # But our load_precomputed_grid might return different format.
        # Let's rely on parsed 'bounds' format [ [x,y], [x,y]... ] from getInfo
        # Actually load_precomputed_grid returns the result of bounds().coordinates()[0] which is a list of points.

        filtered_tiles = []
        for tile in tiles:
            # Simple bbox check using tile properties if available, or geometry
            # Let's use the geometry we loaded
            tile_geom = tile["geometry"]
            # Client side intersection check is hard without extra libs.
            # Simple bounds check:
            # We assume tile["bounds"] is the polygon ring.
            t_coords = tile["bounds"]
            t_min_lon = min(p[0] for p in t_coords)
            t_max_lon = max(p[0] for p in t_coords)
            t_min_lat = min(p[1] for p in t_coords)
            t_max_lat = max(p[1] for p in t_coords)

            if (
                t_min_lon < max_lon
                and t_max_lon > min_lon
                and t_min_lat < max_lat
                and t_max_lat > min_lat
            ):
                filtered_tiles.append(tile)

        original_count = len(tiles)
        tiles = filtered_tiles
        print(f"Filtered from {original_count} to {len(tiles)} tiles within AOI")

    print(f"Processing {len(tiles)} land tiles")
    # Process tiles in batches
    export_tasks = []
    timestamp = datetime.datetime.now().strftime("%Y_%m_%d")

    for i, tile in enumerate(tiles):
        tile_name = tile["name"]
        tile_geom = tile["geometry"]

        # Check if tile asset already exists
        tile_asset_path = f"{tiled_folder}/npp_{tile_name}"

        try:
            ee.data.getAsset(tile_asset_path)
            print(f"Tile {tile_name} already exists, skipping...")
            continue
        except ee.EEException:
            pass  # Asset doesn't exist, proceed with export

        print(
            f"Processing tile {i + 1}/{len(tiles)}: {tile_name} (Land: {tile.get('land_percentage', 0):.1f}%)"
        )

        try:
            # Setup predictors for this tile
            predictors_tile = setup_predictor_stack(resolution=resolution).clip(tile_geom)

            # Generate predictions for tile
            predicted_npp_mean = (
                predictors_tile.classify(model_mean).rename("potential_npp_mean").toFloat()
            )

            predicted_npp_std = (
                predictors_tile.classify(model_std).rename("potential_npp_std").toFloat()
            )

            # Get actual NPP for tile
            actual_npp_tile = (
                setup_response(
                    product="clms",
                    year_range=["2020-01-01", "2025-01-01"],
                )
                .clip(tile_geom)
                .rename("actual_npp")
                .toFloat()
            )

            # Calculate derived metrics for tile
            # Functional Integrity (Z-score): (Actual - Potential Mean) / Potential Std

            # Avoid division by zero by clamping std dev
            safe_std = predicted_npp_std.max(0.1)

            functional_integrity = (
                actual_npp_tile.subtract(predicted_npp_mean)
                .divide(safe_std)
                .rename("functional_integrity")
                .toFloat()
            )

            # Relative NPP (Ratio)
            relative_npp_tile = (
                actual_npp_tile.divide(predicted_npp_mean)
                .rename("relative_npp")
                .clamp(0, 5)  # Reasonable upper bound for ratio
                .toFloat()
            )

            # Raw Difference
            npp_difference_tile = (
                actual_npp_tile.subtract(predicted_npp_mean).rename("npp_difference").toFloat()
            )

            # Create composite for tile
            tile_composite = ee.Image.cat(
                [
                    predicted_npp_mean,
                    predicted_npp_std,
                    actual_npp_tile,
                    functional_integrity,
                    relative_npp_tile,
                    npp_difference_tile,
                ]
            ).set(
                {
                    "tile_name": tile_name,
                    "tile_id": tile.get("tile_id", i),
                    "land_percentage": tile.get("land_percentage", 0),
                    "model_version": model_version,
                    "resolution": resolution,
                    "prediction_date": datetime.datetime.now().strftime("%Y-%m-%d"),
                }
            )

            # Export tile
            task = ee.batch.Export.image.toAsset(
                image=tile_composite,
                description=f"NPP_Tile_{tile_name}",
                assetId=tile_asset_path,
                region=tile_geom,
                scale=resolution,
                maxPixels=1e13,
                crs="EPSG:4326",
            )
            task.start()

            export_tasks.append(task.id)

            # Limit concurrent tasks
            if len(export_tasks) >= max_tiles_per_batch:
                print(
                    f"\\nReached batch limit ({max_tiles_per_batch}). Waiting for current batch to complete..."
                )
                wait_for_tasks(task_ids=export_tasks[-max_tiles_per_batch:])

        except Exception as e:
            print(f"Error processing tile {tile_name}: {e}")

    print(f"\nStarted {len(export_tasks)} tile export tasks")

    return {
        "export_tasks": export_tasks,
        "tiled_folder": tiled_folder,
        "total_tiles": len(tiles),
        "processed_tiles": len(export_tasks),
        "tile_size_deg": tile_size_deg,
        "timestamp": timestamp,
    }
