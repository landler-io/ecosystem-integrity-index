"""
Grid management utilities for global tiling.
"""

from __future__ import annotations

import ee

from .._utils.gee import wait_for_tasks


def create_global_grid_with_land_classification(
    tile_size_deg: float = 10.0,
    min_land_percentage: float = 1.0,
    overwrite: bool = False,
    grid_asset_path: str | None = None,
) -> str:
    """
    Create a pre-computed global grid with land/ocean classification and store as EE asset.
    """

    print("=== Creating Pre-computed Global Grid ===")
    print(f"Tile size: {tile_size_deg}°")
    print(f"Minimum land percentage: {min_land_percentage}%")

    if grid_asset_path is None:
        raise ValueError("grid_asset_path must be provided")

    # Setup paths
    full_grid_asset_path = f"{grid_asset_path}/global_grid_{int(tile_size_deg)}deg"

    # Check if grid already exists
    if not overwrite:
        try:
            ee.data.getAsset(full_grid_asset_path)
            print(f"Grid already exists at: {full_grid_asset_path}")
            print("Use overwrite=True to recreate it.")
            return full_grid_asset_path
        except ee.EEException:
            pass  # Asset doesn't exist, proceed with creation

    # Generate global grid
    print("Generating global tile grid...")
    tiles = []

    # Global bounds
    min_lon, max_lon = -180, 180
    min_lat, max_lat = -80, 80

    tile_id = 0
    for lon in range(int(min_lon), int(max_lon), int(tile_size_deg)):
        for lat in range(int(min_lat), int(max_lat), int(tile_size_deg)):
            tile_bounds = [
                lon,
                lat,
                min(lon + tile_size_deg, max_lon),
                min(lat + tile_size_deg, max_lat),
            ]

            tile_geom = ee.Geometry.Rectangle(tile_bounds)
            tile_name = f"tile_{lon}_{lat}"

            tiles.append(
                ee.Feature(
                    tile_geom,
                    {
                        "tile_id": tile_id,
                        "tile_name": tile_name,
                        "min_lon": lon,
                        "min_lat": lat,
                        "max_lon": min(lon + tile_size_deg, max_lon),
                        "max_lat": min(lat + tile_size_deg, max_lat),
                        "tile_size_deg": tile_size_deg,
                    },
                )
            )

            tile_id += 1

    print(f"Generated {len(tiles)} tiles globally")

    # Create FeatureCollection
    tiles_fc = ee.FeatureCollection(tiles)

    # Add land classification using MODIS Land Cover
    print("Classifying tiles by land coverage using MODIS Land Cover...")

    # Load MODIS Land Cover
    modis_lc = ee.ImageCollection("MODIS/061/MCD12Q1").first()
    landcover = modis_lc.select("LC_Type1")

    # Create land mask (all classes except 17 = water)
    land_mask = landcover.neq(17).rename("land")

    def classify_tile_land_coverage(feature):
        """Add land coverage percentage and classification to each tile"""
        tile_geom = feature.geometry()

        # Calculate land percentage
        land_stats = land_mask.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=tile_geom,
            scale=500,  # MODIS native resolution
            maxPixels=1e7,
        )

        # Handle potential null values from reduceRegion
        land_value = land_stats.get("land", 0)
        land_percentage = ee.Algorithms.If(
            ee.Algorithms.IsEqual(land_value, None),
            0,
            ee.Number(land_value).multiply(100),
        )
        is_land_tile = ee.Number(land_percentage).gte(min_land_percentage)

        return feature.set(
            {
                "land_percentage": land_percentage,
                "is_land_tile": is_land_tile,
                "is_ocean_tile": ee.Number(is_land_tile).Not(),
                "classification_method": "MODIS_LC",
                "min_land_threshold": min_land_percentage,
                "created_date": "2025-01-24",
            }
        )

    # Apply classification to all tiles
    classified_tiles = tiles_fc.map(classify_tile_land_coverage)

    # Export grid to asset
    print(f"Exporting grid to asset: {full_grid_asset_path}")

    # Ensure parent folder exists (handled by notebook logic mostly, but good practice)
    # create_assets_folder(grid_asset_path) # Assuming caller handled this or will fail

    export_task = ee.batch.Export.table.toAsset(
        collection=classified_tiles,
        description=f"Global_Grid_{int(tile_size_deg)}deg",
        assetId=full_grid_asset_path,
    )
    export_task.start()

    print(f"Export task started: {export_task.id}")
    print("Waiting for grid export to complete...")
    wait_for_tasks(task_ids=[export_task.id])

    print(f"Global grid created and stored at: {full_grid_asset_path}")
    return full_grid_asset_path


def load_precomputed_grid(
    tile_size_deg: float = 10.0,
    land_tiles_only: bool = True,
    min_lat: float | None = None,
    max_lat: float | None = None,
    grid_asset_path: str | None = None,
) -> list[ee.Feature]:
    """
    Load the pre-computed global grid with land/ocean classification.
    """

    print("=== Loading Pre-computed Global Grid ===")

    if grid_asset_path is None:
        raise ValueError("grid_asset_path must be provided")

    full_grid_asset_path = f"{grid_asset_path}/global_grid_{int(tile_size_deg)}deg"

    try:
        # Load grid from asset
        print(f"Loading grid from: {full_grid_asset_path}")
        grid_fc = ee.FeatureCollection(full_grid_asset_path)

        # Filter for land tiles if requested
        if land_tiles_only:
            # Use land_percentage threshold (more reliable than boolean filter)
            grid_fc = grid_fc.filter(ee.Filter.gte("land_percentage", 1.0))
            print("Filtered for land tiles only")

        # Filter by latitude range if requested
        if min_lat is not None:
            grid_fc = grid_fc.filter(ee.Filter.gte("min_lat", min_lat))
            print(f"Filtered for tiles >= {min_lat}° latitude")

        if max_lat is not None:
            grid_fc = grid_fc.filter(ee.Filter.lte("max_lat", max_lat))
            print(f"Filtered for tiles <= {max_lat}° latitude")

        # Get all tiles
        tiles = grid_fc.getInfo()["features"]

        # Parse features into uniform dictionaries expected by logic
        parsed_tiles = []
        for feat in tiles:
            props = feat["properties"]
            geom = ee.Geometry(feat["geometry"])
            # Ensure required fields are present
            parsed_tiles.append(
                {
                    "name": props["tile_name"],
                    "geometry": geom,
                    "bounds": geom.bounds().getInfo()["coordinates"][0],  # Approx bounds from geom
                    "land_percentage": props.get("land_percentage", 0),
                    "tile_id": props.get("tile_id"),
                    "properties": props,
                }
            )

        print(f"Loaded {len(parsed_tiles)} tiles")
        return parsed_tiles

    except Exception as e:
        print(f"Error loading grid: {e}")
        return []
