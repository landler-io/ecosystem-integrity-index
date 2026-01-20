"""
Analysis utilities for EII data.

Provides functions for zonal statistics and basic spatial analysis
using pre-computed EII assets.
"""

from __future__ import annotations

from typing import Literal

import ee

from eii.compute.settings import OBSERVED_NPP_YEAR_RANGE

from .retrieve import _build_reducer, _validate_stats_params, get_layers
from .utils import normalize_client_input

AggregationMethod = Literal["min_fuzzy_logic", "minimum", "product", "geometric_mean"]
ComputeMode = Literal["precomputed", "on_the_fly"]


def _features_to_geodataframe(features: ee.FeatureCollection):
    try:
        import geopandas as gpd
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("geopandas is required for geodataframe output") from exc

    info = features.getInfo()
    return gpd.GeoDataFrame.from_features(info["features"], crs="EPSG:4326")


def get_zonal_stats(
    features: ee.FeatureCollection,
    aggregation_method: AggregationMethod = "min_fuzzy_logic",
    include_components: bool = True,
    scale: int = 300,
    stats: list[str] | None = None,
    percentiles: list[int] | None = None,
    compute_mode: ComputeMode = "precomputed",
    npp_year_range: list[str] | None = None,
    include_seasonality: bool = True,
    bii_year: int | None = None,
    absolute_diff_percentile: str = "p95",
    keep_columns: list[str] | None = None,
):
    """
    Calculate EII zonal statistics for each feature in a collection.

    Args:
        features: FeatureCollection or GeoPandas/shapely inputs containing zones.
        aggregation_method: EII aggregation method.
        include_components: Include component statistics.
        scale: Resolution in meters for the reduction. Defaults to 300.
        stats: Statistics to compute. Supported: "mean", "median", "min",
            "max", "std". Defaults to ["mean"].
        percentiles: Optional list of percentiles (0-100) to include.
        compute_mode: "precomputed" (default) or "on_the_fly".
        npp_year_range: Date range for observed NPP (on-the-fly mode).
        include_seasonality: Include seasonality in functional integrity.
        bii_year: BII year to select for compositional integrity.
        absolute_diff_percentile: Percentile key for NPP absolute difference penalty.
        keep_columns: Optional list of input feature properties to preserve.
            Extracted stats are always included.

    Returns:
        GeoDataFrame with EII statistics for each feature. Columns use MultiIndex
        format: (metric, stat), e.g., (eii, mean), (functional_integrity, min).

    Example:
        >>> import ee
        >>> ee.Initialize()
        >>> countries = ee.FeatureCollection("FAO/GAUL/2015/level0")
        >>> gdf = get_zonal_stats(countries.limit(10), stats=["mean", "min", "max"])
    """
    layers = get_layers(
        layers="all" if include_components else "eii",
        aggregation_method=aggregation_method,
        compute_mode=compute_mode,
        geometry=None,
        npp_year_range=npp_year_range or OBSERVED_NPP_YEAR_RANGE,
        include_seasonality=include_seasonality,
        bii_year=bii_year or 2020,
        absolute_diff_percentile=absolute_diff_percentile,
    )
    image_stack = layers["eii"]
    if include_components:
        for name in ("functional", "structural", "compositional"):
            image_stack = image_stack.addBands(layers[name].rename(f"{name}_integrity"))

    features = normalize_client_input(features, target="features")
    stats = _validate_stats_params(stats, percentiles)
    reducer = _build_reducer(stats, percentiles)

    def add_stats(feature: ee.Feature) -> ee.Feature:
        stats_result = image_stack.reduceRegion(
            reducer=reducer,
            geometry=feature.geometry(),
            scale=scale,
            maxPixels=1e12,
        )
        return feature.set(stats_result)

    result = features.map(add_stats)

    try:
        import geopandas as gpd
        import pandas as pd
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("geopandas and pandas are required for get_zonal_stats") from exc

    gdf = _features_to_geodataframe(result)
    original_props = features.first().propertyNames().getInfo()
    original_set = set(original_props or [])

    if keep_columns:
        if isinstance(keep_columns, str):
            keep_columns = [keep_columns]
        keep_set = set(keep_columns)
        stat_columns = [col for col in gdf.columns if col not in original_set and col != "geometry"]
        preserved = [col for col in gdf.columns if col in keep_set or col in stat_columns]
        if "geometry" not in preserved:
            preserved.append("geometry")
        gdf = gdf[preserved]

    geometry = gdf.geometry
    data = gdf.drop(columns="geometry")
    tuples: list[tuple[str, str]] = []
    for col in data.columns:
        if col in original_set:
            tuples.append((col, "value"))
            continue
        parts = col.rsplit("_", 1)
        if len(parts) == 2:
            metric, stat = parts
            if stat == "stdDev":
                stat = "std"
            tuples.append((metric, stat))
        else:
            tuples.append((col, "value"))
    data.columns = pd.MultiIndex.from_tuples(tuples)

    metric_order = [
        "eii",
        "functional_integrity",
        "structural_integrity",
        "compositional_integrity",
    ]
    ordered_cols: list[tuple[str, str]] = []
    if keep_columns:
        ordered_cols.extend(
            [(col, "value") for col in keep_columns if (col, "value") in data.columns]
        )
    for metric in metric_order:
        metric_cols = [col for col in data.columns if col[0] == metric]
        ordered_cols.extend([col for col in metric_cols if col not in ordered_cols])
    remaining = [col for col in data.columns if col not in ordered_cols]
    ordered_cols.extend(remaining)
    if ordered_cols:
        data = data.loc[:, ordered_cols]
    geometry.name = ("geometry", "")
    return gpd.GeoDataFrame(data, geometry=geometry, crs=gdf.crs)


def compare_methods(
    geometry: ee.Geometry,
    scale: int = 300,
    compute_mode: ComputeMode = "precomputed",
    npp_year_range: list[str] | None = None,
    include_seasonality: bool = True,
    bii_year: int | None = None,
    absolute_diff_percentile: str = "p95",
) -> ee.Dictionary:
    """
    Compare EII values across different aggregation methods for a polygon.

    Args:
        geometry: AOI geometry (ee.Geometry/Feature/FeatureCollection, GeoPandas,
            shapely, or bbox tuple).
        scale: Resolution in meters.
        compute_mode: "precomputed" (default) or "on_the_fly".
        npp_year_range: Date range for observed NPP (on-the-fly mode).
        include_seasonality: Include seasonality in functional integrity.
        bii_year: BII year to select for compositional integrity.
        absolute_diff_percentile: Percentile key for NPP absolute difference penalty.

    Returns:
        Dictionary with mean EII for each aggregation method.

    Example:
        >>> comparison = compare_methods(my_polygon)
        >>> print(comparison.getInfo())
        {'minimum': 0.65, 'product': 0.58, 'min_fuzzy_logic': 0.72, 'geometric_mean': 0.68}
    """
    geometry = normalize_client_input(geometry, target="geometry")
    methods: list[AggregationMethod] = [
        "minimum",
        "product",
        "min_fuzzy_logic",
        "geometric_mean",
    ]
    result = ee.Dictionary({})

    for method in methods:
        if compute_mode == "precomputed":
            eii = get_layers(
                layers="eii",
                aggregation_method=method,
                compute_mode=compute_mode,
            )["eii"]
        elif compute_mode == "on_the_fly":
            eii = get_layers(
                layers="eii",
                aggregation_method=method,
                compute_mode=compute_mode,
                geometry=geometry,
                npp_year_range=npp_year_range or OBSERVED_NPP_YEAR_RANGE,
                include_seasonality=include_seasonality,
                bii_year=bii_year or 2020,
                absolute_diff_percentile=absolute_diff_percentile,
            )["eii"]
        else:
            raise ValueError("compute_mode must be one of: 'precomputed', 'on_the_fly'")
        mean_val = eii.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=geometry,
            scale=scale,
            maxPixels=1e12,
        ).get("eii")
        result = result.set(method, mean_val)

    return result
