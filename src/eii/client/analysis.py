"""
Analysis utilities for EII data.

Provides functions for zonal statistics and basic spatial analysis
using pre-computed EII assets.
"""

from __future__ import annotations

import ee

from eii.compute.settings import AggregationMethod

from .assets import resolve_version
from .retrieve import _build_reducer, _validate_stats_params, get_layers
from .utils import normalize_client_input


def _features_to_geodataframe(features: ee.FeatureCollection):
    try:
        import geopandas as gpd
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("geopandas is required for geodataframe output") from exc

    info = features.getInfo()
    return gpd.GeoDataFrame.from_features(info["features"], crs="EPSG:4326")


def get_zonal_stats(
    features: ee.FeatureCollection,
    include_components: bool = True,
    scale: int = 300,
    stats: list[str] | None = None,
    percentiles: list[int] | None = None,
    keep_columns: list[str] | None = None,
    version: str | None = None,
    reference_year: int | None = None,
    exclude_wetlands: bool = False,
):
    """
    Calculate EII zonal statistics for each feature in a collection.

    Args:
        features: FeatureCollection or GeoPandas/shapely inputs containing zones.
        include_components: Include component statistics.
        scale: Resolution in meters for the reduction. Defaults to 300.
        stats: Statistics to compute. Supported: "mean", "median", "min",
            "max", "std". Defaults to ["mean"].
        percentiles: Optional list of percentiles (0-100) to include.
        reference_year: v2 reference year, defaults to the latest; see ``get_layers``.
        keep_columns: Optional list of input feature properties to preserve.
            Extracted stats are always included.
        version: EII dataset version (e.g. "v1", "v2"). If None, uses the active
            session version (see eii.client.set_version), defaulting to config.yaml.
        exclude_wetlands: v2 only. Mask pixels with a nonzero ``wetland_flag``
            before the reduction. By default flagged pixels are only reported.

    Returns:
        GeoDataFrame with EII statistics for each feature. Columns use MultiIndex
        format: (metric, stat), e.g., (eii, mean), (functional_integrity, min).
        Version v2 adds (wetland_flag, share): the flagged share of valid EII
        pixels, computed before any exclusion.

    Example:
        >>> import ee
        >>> ee.Initialize()
        >>> countries = ee.FeatureCollection("FAO/GAUL/2015/level0")
        >>> gdf = get_zonal_stats(countries.limit(10), stats=["mean", "min", "max"])
    """
    version = resolve_version(version)
    if exclude_wetlands and version != "v2":
        raise ValueError("exclude_wetlands is available only for version='v2'")
    layers = get_layers(
        layers="all" if include_components else "eii",
        version=version,
        reference_year=reference_year,
    )
    image_stack = layers["eii"]
    if include_components:
        for name in ("functional", "structural", "compositional"):
            image_stack = image_stack.addBands(layers[name].rename(f"{name}_integrity"))
    flagged = None
    if version == "v2":
        flag = get_layers("wetland_flag", version, reference_year)["wetland_flag"]
        flagged = flag.gt(0).updateMask(layers["eii"].mask()).rename("wetland_flag_share")
        if exclude_wetlands:
            image_stack = image_stack.updateMask(flag.eq(0))

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
        if flagged is not None:
            share = flagged.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=feature.geometry(),
                scale=scale,
                maxPixels=10**12,
            )
            stats_result = stats_result.combine(share)
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
    version: str | None = None,
    reference_year: int | None = None,
) -> ee.Dictionary:
    """
    Compare EII values across different aggregation methods for a polygon.

    Each method recombines the pre-computed components; only ``min_fuzzy_logic`` is
    published as the EII band.

    Args:
        geometry: AOI geometry (ee.Geometry/Feature/FeatureCollection, GeoPandas,
            shapely, or bbox tuple).
        scale: Resolution in meters.
        reference_year: v2 reference year, defaults to the latest; see ``get_layers``.

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
    components = get_layers("components", version, reference_year)
    result = ee.Dictionary({})
    for method in methods:
        eii = combine_components(
            components["functional"], components["structural"], components["compositional"], method
        )
        mean_val = eii.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=geometry,
            scale=scale,
            maxPixels=1e12,
        ).get("eii")
        result = result.set(method, mean_val)

    return result


def combine_components(
    functional: ee.Image,
    structural: ee.Image,
    compositional: ee.Image,
    method: AggregationMethod = "min_fuzzy_logic",
) -> ee.Image:
    """Earth Engine counterpart of ``eii.compute.combine_components`` for component layers.

    Args:
        functional: Functional integrity image (0-1).
        structural: Structural integrity image (0-1).
        compositional: Compositional integrity image (0-1).
        method: "minimum", "product", "geometric_mean" or "min_fuzzy_logic".

    Returns:
        EII image named "eii", masked wherever any component is missing.
    """
    components = ee.Image.cat([functional, structural, compositional])
    if method == "minimum":
        eii = components.reduce(ee.Reducer.min())
    elif method == "product":
        eii = components.reduce(ee.Reducer.product())
    elif method == "geometric_mean":
        eii = components.reduce(ee.Reducer.product()).pow(1.0 / 3.0)
    elif method == "min_fuzzy_logic":
        low = components.reduce(ee.Reducer.min())
        mid = components.reduce(ee.Reducer.median())
        high = components.reduce(ee.Reducer.max())
        eii = low.multiply(mid.add(high).subtract(mid.multiply(high)))
    else:
        raise ValueError(f"Unknown method '{method}'")
    return eii.updateMask(components.mask().reduce(ee.Reducer.min())).rename("eii")
