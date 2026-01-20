"""
Training data sampling utilities.

Ecoregion-stratified sampling from pristine reference areas for NPP model training.
Uses HMI-based natural areas masks and high-category IUCN protected areas.
"""

from __future__ import annotations

import zlib

import ee

from .settings import (
    ECOREGIONS_RASTER_ASSET,
    HMI_MASKS_COLLECTION,
    SAMPLE_SCALE,
    SAMPLES_PER_ECOREGION,
    TILE_SCALE,
    TOTAL_POINTS_PER_GRID_CELL,
    TRAINING_GRID_SIZE_DEG,
    WDPA_RASTER_ASSET,
)


def get_pristine_mask(
    aoi: ee.Geometry | None = None,
) -> ee.Image:
    """
    Create a binary mask of pristine/natural areas.

    Combines:
    - HMI-based natural areas (pre-computed raster mask)
    - IUCN protected areas (pre-computed raster mask, already buffered)

    Args:
        aoi: Optional area of interest to constrain computation.

    Returns:
        Binary ee.Image (1 = pristine, 0 = not pristine).
    """

    # Load HMI natural areas mask (ImageCollection of tiles)
    hmi_masks = ee.ImageCollection(HMI_MASKS_COLLECTION)
    if aoi is not None:
        hmi_masks = hmi_masks.filterBounds(aoi)
    hmi_mask = hmi_masks.mosaic().unmask(0).gt(0)

    # Load Pre-computed Protected Areas Raster
    pa_mask = ee.Image(WDPA_RASTER_ASSET).unmask(0).gt(0)

    # Combine: pristine if HMI natural OR protected area
    pristine_mask = hmi_mask.Or(pa_mask).rename("pristine")

    if aoi is not None:
        pristine_mask = pristine_mask.clip(aoi)

    return pristine_mask


def get_ecoregion_image() -> ee.Image:
    """
    Load ecoregions as a raster image with ECO_ID as values.

    Attempts to load pre-rasterized version first, falls back to on-the-fly
    rasterization if not available.
    """
    return ee.Image(ECOREGIONS_RASTER_ASSET).rename("ecoregion")


def sample_pristine_areas_stratified(
    predictor_stack: ee.Image,
    samples_per_ecoregion: int = SAMPLES_PER_ECOREGION,
    aoi: ee.Geometry | None = None,
    seed: int = 42,
    tile_scale: int = TILE_SCALE,
) -> ee.FeatureCollection:
    """
    Sample predictor values from pristine areas, stratified by ecoregion.

    Args:
        predictor_stack: Image containing all predictor bands + response.
        samples_per_ecoregion: Target number of samples per ecoregion.
        aoi: Optional area of interest.
        seed: Random seed for reproducibility.
        tile_scale: GEE tileScale parameter for memory optimization.

    Returns:
        FeatureCollection of sampled points with predictor values.
    """

    pristine_mask = get_pristine_mask(aoi)
    ecoregion_img = ee.Image(ECOREGIONS_RASTER_ASSET).rename("ecoregion")

    sample_image = predictor_stack.addBands(ecoregion_img)
    sample_image = sample_image.updateMask(pristine_mask)

    region = aoi if aoi is not None else ee.Geometry.Rectangle([-180, -60, 180, 90])

    # Stratified sampling by ecoregion
    samples = sample_image.stratifiedSample(
        numPoints=samples_per_ecoregion,
        classBand="ecoregion",
        region=region,
        scale=SAMPLE_SCALE,
        seed=seed,
        geometries=True,
        tileScale=tile_scale,
    )

    return samples


def sample_grid_cell_stratified(
    predictor_stack: ee.Image,
    grid_cell: ee.Geometry,
    grid_cell_name: str,
    _total_points: int = TOTAL_POINTS_PER_GRID_CELL,
    seed: int = 42,
    tile_scale: int = TILE_SCALE,
) -> ee.FeatureCollection:
    """
    Sample from a grid cell with points stratified by ecoregion within pristine areas.

    This replaces area-proportional sampling to ensure better ecological representation
    within each spatial grid cell throughout the pristine mask.

    Args:
        predictor_stack: Image containing predictors + response.
        grid_cell: Geometry of the grid cell.
        grid_cell_name: Name/ID of the grid cell.
        total_points: Total target points for the cell.
        seed: Random seed.
        tile_scale: GEE tileScale for memory optimization.

    Returns:
        FeatureCollection of samples with grid_cell property.
    """

    pristine_mask = get_pristine_mask(aoi=grid_cell)
    ecoregion_img = get_ecoregion_image().clip(grid_cell)

    sample_image = predictor_stack.addBands(ecoregion_img)
    sample_image = sample_image.updateMask(pristine_mask)

    cell_seed = zlib.adler32(grid_cell_name.encode("utf-8")) % 10000 + seed
    samples = sample_image.stratifiedSample(
        numPoints=SAMPLES_PER_ECOREGION,
        classBand="ecoregion",
        region=grid_cell,
        scale=SAMPLE_SCALE,
        seed=cell_seed,
        geometries=True,
        tileScale=tile_scale,
        dropNulls=True,
    )

    samples = samples.map(lambda f: f.set("grid_cell", grid_cell_name))

    return samples


def sample_grid_cell_area_proportional(
    predictor_stack: ee.Image,
    grid_cell: ee.Geometry,
    grid_cell_name: str,
    total_points: int = TOTAL_POINTS_PER_GRID_CELL,
    seed: int = 42,
    tile_scale: int = TILE_SCALE,
) -> ee.FeatureCollection:
    """
    Sample from a grid cell with points distributed proportionally to available natural area.

    Uses efficient raster-based sampling (random sampling within valid mask).

    Args:
        predictor_stack: Image containing predictors + response.
        grid_cell: Geometry of the grid cell.
        grid_cell_name: Name/ID of the grid cell.
        total_points: Total target points for the cell.
        seed: Random seed.
        tile_scale: GEE tileScale for memory optimization.

    Returns:
        FeatureCollection of samples with grid_cell property.
    """

    pristine_mask = get_pristine_mask(aoi=grid_cell)
    sample_image = predictor_stack.updateMask(pristine_mask)
    cell_seed = zlib.adler32(grid_cell_name.encode("utf-8")) % 10000 + seed

    samples = sample_image.sample(
        region=grid_cell,
        scale=SAMPLE_SCALE,
        numPixels=total_points,
        seed=cell_seed,
        geometries=True,
        tileScale=tile_scale,
        dropNulls=True,
    )

    samples = samples.map(lambda f: f.set("grid_cell", grid_cell_name))

    return samples


def setup_training_grid(grid_size_deg: int = TRAINING_GRID_SIZE_DEG) -> dict:
    """
    Create a global grid for spatial cross-validation.

    Args:
        grid_size_deg: Size of grid cells in degrees.

    Returns:
        Dictionary mapping cell names to ee.Geometry objects.
    """

    lon_min, lon_max = -180, 180
    lat_min, lat_max = -60, 90  # Exclude Antarctica

    grid_cells = {}
    for lon in range(lon_min, lon_max, grid_size_deg):
        for lat in range(lat_min, lat_max, grid_size_deg):
            cell_name = f"grid_{lon}_{lat}"
            cell = ee.Geometry.Rectangle([lon, lat, lon + grid_size_deg, lat + grid_size_deg])
            grid_cells[cell_name] = cell

    return grid_cells


def sample_all_grid_cells(
    predictor_stack: ee.Image,
    grid_cells: dict | None = None,
    _total_points_per_cell: int = TOTAL_POINTS_PER_GRID_CELL,
    seed: int = 42,
) -> dict[str, ee.FeatureCollection]:
    """
    Sample from all grid cells using area-proportional sampling.

    Args:
        predictor_stack: Image containing predictors + response.
        grid_cells: Dictionary of grid cell geometries. If None, creates default.
        total_points_per_cell: Total target points per cell.
        seed: Random seed.

    Returns:
        Dictionary mapping cell names to sample FeatureCollections.
    """
    if grid_cells is None:
        grid_cells = setup_training_grid()

    samples_by_cell = {}
    for cell_name, cell_geom in grid_cells.items():
        samples = sample_grid_cell_stratified(
            predictor_stack=predictor_stack,
            grid_cell=cell_geom,
            grid_cell_name=cell_name,
            seed=seed,
        )
        samples_by_cell[cell_name] = samples

    return samples_by_cell
