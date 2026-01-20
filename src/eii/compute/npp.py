"""
NPP-based functional integrity calculation.

This module handles:
- Predictor stack setup for the NPP model
- NPP model inference
- Functional integrity score calculation
"""

from __future__ import annotations

import ee

from .settings import (
    DEFAULT_NPP_MODEL_PATH,
    INCLUDE_LAT_LON_PREDICTORS,
    MAGNITUDE_WEIGHT,
    NATURAL_NPP_ASSET_PATH,
    NPP_DIFF_PERCENTILES_ASSET_PATH,
    NPP_YEAR_RANGE,
    OBSERVED_NPP_ASSET_PATH,
    OBSERVED_NPP_YEAR_RANGE,
    SEASONALITY_WEIGHT,
    SPATIAL_RESOLUTION,
)

# Module-level cache for NPP diff percentiles (loaded once per session)
_NPP_DIFF_PERCENTILES_CACHE: dict | None = None


def load_npp_diff_percentiles(
    asset_path: str = NPP_DIFF_PERCENTILES_ASSET_PATH,
) -> dict[str, float]:
    """
    Load NPP difference percentiles from GEE asset.

    The asset is a FeatureCollection with a single feature containing
    percentile breaks as properties (p05, p10, ..., p95).

    Args:
        asset_path: Path to the percentiles asset.

    Returns:
        Dictionary with percentile values.

    Raises:
        RuntimeError: If the percentiles asset is not available.
            Run compute_npp_decile_breaks.ipynb to create it.
    """
    global _NPP_DIFF_PERCENTILES_CACHE

    if _NPP_DIFF_PERCENTILES_CACHE is not None:
        return _NPP_DIFF_PERCENTILES_CACHE

    try:
        fc = ee.FeatureCollection(asset_path)
        props = fc.first().getInfo()["properties"]

        percentiles = {f"p{p}": props[f"p{p}"] for p in range(5, 100, 5)}
        _NPP_DIFF_PERCENTILES_CACHE = percentiles
        return percentiles

    except Exception as e:
        raise RuntimeError(
            f"NPP diff percentiles asset not found at: {asset_path}\n"
            "Run pipelines/modeling/eii_calculation/compute_npp_decile_breaks.ipynb "
            "to create it."
        ) from e


def load_npp_decile_breaks(
    asset_path: str = NPP_DIFF_PERCENTILES_ASSET_PATH,
) -> dict[str, float]:
    """
    Backward-compatible alias for loading NPP difference percentiles.
    """
    return load_npp_diff_percentiles(asset_path=asset_path)


def setup_predictor_stack(
    resolution: int = SPATIAL_RESOLUTION,
    include_lat_lon: bool | None = None,
    include_regional_tpi: bool = True,
) -> ee.Image:
    """
    Set up the predictor stack for NPP model inference.

    Args:
        resolution: Output resolution in meters.
        include_lat_lon: Whether to include lat/lon as predictors.
            If None, uses INCLUDE_LAT_LON_PREDICTORS from settings.
        include_regional_tpi: Whether to calculate and include a regional TPI
            (Topographic Position Index) based on a coarser DEM.
            Defaults to True.

    Returns:
        Multi-band image containing all predictor variables.
    """
    if include_lat_lon is None:
        include_lat_lon = INCLUDE_LAT_LON_PREDICTORS

    # WorldClim bioclimatic variables
    worldclim_bio = ee.Image("WORLDCLIM/V1/BIO")
    bio01 = worldclim_bio.select("bio01").rename("mean_annual_temp")
    bio04 = worldclim_bio.select("bio04").rename("temp_seasonality")
    bio12 = worldclim_bio.select("bio12").rename("annual_precip")
    bio15 = worldclim_bio.select("bio15").rename("precip_seasonality")

    # Aridity index
    aridity = (
        ee.Image("projects/sat-io/open-datasets/global_ai/global_ai_yearly")
        .select("b1")
        .rename("aridity")
    )

    # Topography
    dem = ee.Image("MERIT/DEM/v1_0_3")
    elevation = dem.select("dem").rename("elevation")
    slope = (
        ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/slope")
        .mosaic()
        .rename("slope")
    )
    tpi = (
        ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/tpi").mosaic().rename("tpi")
    )
    tri = (
        ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/tri").mosaic().rename("tri")
    )
    cti = (
        ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/cti").mosaic().rename("cti")
    )
    northness = (
        ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/northness")
        .mosaic()
        .rename("northness")
    )
    eastness = (
        ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/eastness")
        .mosaic()
        .rename("eastness")
    )
    sand = ee.Image(
        "projects/landler-open-data/assets/datasets/soilgrids/sand/sand_15-30cm_mean_gapfilled"
    ).rename("sand")
    clay = ee.Image(
        "projects/landler-open-data/assets/datasets/soilgrids/clay/clay_15-30cm_mean_gapfilled"
    ).rename("clay")
    ph = ee.Image(
        "projects/landler-open-data/assets/datasets/soilgrids/phh2o/phh2o_15-30cm_mean_gapfilled"
    ).rename("ph")

    # CHELSA NPP (climatological reference - summarizes climatic NPP potential)
    chelsa_npp = ee.Image(
        "projects/landler-open-data/assets/datasets/chelsa/npp/chelsa_npp_1981_2010_v2-1"
    ).rename("chelsa_npp")

    # Base predictors
    predictor_bands = [
        bio01,
        bio04,
        bio12,
        bio15,
        aridity,
        elevation,
        slope,
        tpi,
        tri,
        cti,
        northness,
        eastness,
        chelsa_npp,
        sand,
        clay,
        ph,
    ]

    if include_regional_tpi:
        dem_coarse = dem.reproject(crs="EPSG:4326", scale=500)
        regional_mean = dem_coarse.focal_mean(2000, "circle", "meters")
        tpi_regional = elevation.subtract(regional_mean).rename("tpi_regional")
        predictor_bands.append(tpi_regional)

    if include_lat_lon:
        lat_lon = ee.Image.pixelLonLat().rename(["longitude", "latitude"])
        predictor_bands.append(lat_lon)

    continuous_predictors = ee.Image.cat(predictor_bands)
    continuous_predictors = continuous_predictors.setDefaultProjection(crs="EPSG:4326", scale=90)
    continuous_predictors = continuous_predictors.reduceResolution(
        reducer=ee.Reducer.mean(), maxPixels=1024
    ).reproject(crs="EPSG:4326", scale=resolution)

    return continuous_predictors


def setup_response(
    product: str = "clms",
    year_range: list[str] = NPP_YEAR_RANGE,
    include_std: bool = False,
) -> ee.Image:
    """
    Set up the response variable(s) (observed NPP).

    We model the "Average Annual NPP Sum".
    This is calculated by averaging the 'yearly NPP sum' assets over the provided range.

    Args:
        product: NPP data source ("clms" or "modis").
        year_range: Date range [start, end]. Defaults to ["2014-01-01", "2025-01-01"].
        include_std: If True, also compute inter-annual standard deviation.

    Returns:
        Image with "longterm_avg_npp_sum" band (representing the multi-year average of annual sums),
        and optionally "longterm_avg_npp_sd" band.
    """

    if product == "modis":
        # MOD17A3HGF is already Annual Sum
        npp_collection = (
            ee.ImageCollection("MODIS/061/MOD17A3HGF")
            .filter(ee.Filter.date(year_range[0], year_range[1]))
            .select("Npp")
        )
    else:
        from .._utils.gee import load_tiled_collection

        raw_collection = load_tiled_collection(
            "projects/landler-open-data/assets/datasets/clms/npp/annual"
        ).filter(ee.Filter.date(year_range[0], year_range[1]))

        npp_collection = raw_collection.select([0])

    npp_mean = npp_collection.mean().rename("longterm_avg_npp_sum")

    if include_std:
        if product == "modis":
            npp_std = npp_collection.reduce(ee.Reducer.stdDev()).rename("longterm_avg_npp_sd")
        else:
            npp_std = raw_collection.select([1]).mean().rename("longterm_avg_npp_sd")

        return npp_mean.addBands(npp_std)

    return npp_mean


def setup_training_image(
    product: str = "clms",
    year_range: list[str] = NPP_YEAR_RANGE,
    include_lat_lon: bool | None = None,
    include_regional_tpi: bool = True,
) -> ee.Image:
    """
    Set up combined predictor + response image for training.

    Includes both longterm_avg_npp_sum (mean) and longterm_avg_npp_sd as response variables.

    Args:
        product: NPP data source.
        year_range: Date range for NPP observations.
        include_lat_lon: Whether to include lat/lon as predictors.
        include_regional_tpi: Whether to include regional TPI logic.

    Returns:
        Image with predictor bands + longterm_avg_npp_sum + longterm_avg_npp_sd bands.
    """
    predictors = setup_predictor_stack(
        include_lat_lon=include_lat_lon,
        include_regional_tpi=include_regional_tpi,
    )
    response = setup_response(product=product, year_range=year_range, include_std=True)

    return predictors.addBands(response)


def load_natural_npp(
    aoi: ee.Geometry | None = None,
    asset_path: str = NATURAL_NPP_ASSET_PATH,
    use_tiled_collection: bool = False,
) -> dict[str, ee.Image]:
    """Load pre-computed natural NPP mean/std from an image or tiled collection."""
    if use_tiled_collection:
        from .._utils.gee import load_tiled_collection

        collection = load_tiled_collection(asset_path)
        natural_npp_image = collection.mosaic()
    else:
        natural_npp_image = ee.Image(asset_path)

    natural_npp_mean = natural_npp_image.select("natural_npp_mean")
    natural_npp_std = natural_npp_image.select("natural_npp_std")

    if aoi is not None:
        natural_npp_mean = natural_npp_mean.clip(aoi)
        natural_npp_std = natural_npp_std.clip(aoi)

    return {"natural_npp_mean": natural_npp_mean, "natural_npp_std": natural_npp_std}


def load_natural_npp_tiles(
    aoi: ee.Geometry | None = None,
    asset_path: str = NATURAL_NPP_ASSET_PATH,
    use_tiled_collection: bool = False,
) -> dict[str, ee.Image]:
    """Backward-compatible alias for loading natural NPP."""
    return load_natural_npp(
        aoi=aoi,
        asset_path=asset_path,
        use_tiled_collection=use_tiled_collection,
    )


def _load_observed_npp(
    aoi: ee.Geometry | None = None,
    year_range: list[str] = OBSERVED_NPP_YEAR_RANGE,
    asset_path: str = OBSERVED_NPP_ASSET_PATH,
) -> dict[str, ee.Image]:
    """Load observed NPP mean and std from CLMS annual assets."""
    from .._utils.gee import load_tiled_collection

    annual_collection = load_tiled_collection(asset_path).filterDate(year_range[0], year_range[1])
    observed_mean = annual_collection.select([0]).mean().rename("observed_npp_mean")
    observed_std = annual_collection.select([1]).mean().rename("observed_npp_std")

    if aoi is not None:
        observed_mean = observed_mean.clip(aoi)
        observed_std = observed_std.clip(aoi)

    return {"observed_npp_mean": observed_mean, "observed_npp_std": observed_std}


def calculate_functional_integrity(
    aoi: ee.Geometry | None = None,
    year_range: list[str] = OBSERVED_NPP_YEAR_RANGE,
    include_seasonality: bool = True,
    use_precomputed: bool = True,
    model_asset_path: str = DEFAULT_NPP_MODEL_PATH,
    natural_npp_asset_path: str = NATURAL_NPP_ASSET_PATH,
    observed_npp_asset_path: str = OBSERVED_NPP_ASSET_PATH,
    natural_npp_use_tiled_collection: bool = False,
    absolute_diff_percentile: str = "p95",
) -> dict[str, ee.Image]:
    """
    Calculate NPP-based functional integrity with magnitude and seasonality dimensions.

    Magnitude integrity: symmetric deviation score + absolute diff score (weight: 2/3)
    Seasonality integrity: comparison of observed vs natural intra-annual std (weight: 1/3)

    Args:
        aoi: Area of interest. If None, returns unclipped global image.
        year_range: Date range for observed NPP (3-year rolling window).
        include_seasonality: Include seasonality dimension in scoring.
        use_precomputed: Use pre-computed natural NPP tiles (recommended).
        model_asset_path: Path to NPP model (only used if use_precomputed=False).
        natural_npp_asset_path: Asset path for pre-computed natural NPP.
        observed_npp_asset_path: Asset path for observed NPP annual tiles.
        natural_npp_use_tiled_collection: Use tiled collection for natural NPP.
        absolute_diff_percentile: Percentile key for absolute NPP penalty (e.g. "p80").

    Returns:
        Dictionary with intermediate layers and final functional_integrity score.
    """
    if use_precomputed:
        natural_npp = load_natural_npp(
            aoi=aoi,
            asset_path=natural_npp_asset_path,
            use_tiled_collection=natural_npp_use_tiled_collection,
        )
        potential_npp = natural_npp["natural_npp_mean"].rename("potential_npp")
        natural_npp_std = natural_npp["natural_npp_std"]
    else:
        if include_seasonality:
            raise ValueError("include_seasonality=True requires use_precomputed=True")
        npp_model = ee.Classifier.load(model_asset_path)
        predictors = setup_predictor_stack()
        if aoi is not None:
            predictors = predictors.clip(aoi)
        potential_npp = predictors.classify(npp_model).rename("potential_npp")
        natural_npp_std = None

    observed = _load_observed_npp(
        aoi=aoi,
        year_range=year_range,
        asset_path=observed_npp_asset_path,
    )
    actual_npp = observed["observed_npp_mean"].rename("actual_npp")

    relative_npp = actual_npp.divide(potential_npp).rename("relative_npp")
    npp_difference = actual_npp.subtract(potential_npp).abs().rename("npp_difference")

    smooth_kernel = ee.Kernel.square(1, "pixels")
    relative_npp = relative_npp.focal_mean(kernel=smooth_kernel, iterations=1)
    npp_difference = npp_difference.focal_mean(kernel=smooth_kernel, iterations=1)
    proportional_score = _calculate_symmetric_deviation_score(relative_npp)
    absolute_score = _apply_npp_absolute_diff_scaling(
        npp_difference,
        percentile_key=absolute_diff_percentile,
    )
    magnitude_integrity = (
        proportional_score.add(absolute_score).divide(2.0).rename("magnitude_integrity")
    )

    if include_seasonality and natural_npp_std is not None:
        observed_std = observed["observed_npp_std"]
        std_ratio = observed_std.divide(natural_npp_std).rename("std_ratio")
        std_deviation = std_ratio.subtract(1).abs()
        seasonality_integrity = (
            ee.Image(1).divide(ee.Image(1).add(std_deviation)).rename("seasonality_integrity")
        )
        functional_integrity = (
            magnitude_integrity.multiply(MAGNITUDE_WEIGHT)
            .add(seasonality_integrity.multiply(SEASONALITY_WEIGHT))
            .rename("functional_integrity")
        )
    else:
        seasonality_integrity = None
        functional_integrity = magnitude_integrity.rename("functional_integrity")

    return {
        "potential_npp": potential_npp,
        "natural_npp_std": natural_npp_std,
        "actual_npp": actual_npp,
        "relative_npp": relative_npp,
        "npp_difference": npp_difference,
        "proportional_score": proportional_score,
        "absolute_score": absolute_score,
        "magnitude_integrity": magnitude_integrity,
        "seasonality_integrity": seasonality_integrity,
        "functional_integrity": functional_integrity,
    }


def _calculate_symmetric_deviation_score(relative_npp: ee.Image) -> ee.Image:
    """
    Calculate integrity score based on symmetric deviation from potential NPP.

    Both under-performance and over-performance are penalized equally.
    Score = 1 / (1 + |relative_npp - 1|)

    Returns:
        Score from 0-1, where 1 = at potential.
    """
    deviation = relative_npp.subtract(1).abs()
    score = ee.Image(1).divide(ee.Image(1).add(deviation))

    return score.rename("symmetric_npp_score")


def _apply_npp_absolute_diff_scaling(
    npp_difference: ee.Image,
    percentiles: dict | None = None,
    percentile_key: str = "p95",
) -> ee.Image:
    """
    Scale NPP difference into a truncated linear penalty score.

    Args:
        npp_difference: Absolute NPP difference image.
        percentiles: Dictionary with p05..p95 thresholds.
            If None, loads from GEE asset (uses p95).

    Returns:
        Score from 0-1 based on linear scaling to p95.
    """
    if percentiles is None:
        percentiles = load_npp_diff_percentiles()

    breaks = ee.Dictionary(percentiles)
    if percentile_key not in percentiles:
        available = ", ".join(sorted(percentiles.keys()))
        raise ValueError(f"Percentile '{percentile_key}' not found. Available: {available}")

    p95 = ee.Number(breaks.get(percentile_key))

    # Truncate at p95 and scale linearly
    diff_clamped = npp_difference.clamp(0, p95)
    score = ee.Image(1).subtract(diff_clamped.divide(p95))

    return score.rename("npp_difference_penalty")
