"""
Client API for Natural Capital modulated EII.

This module provides user-facing functions to retrieve EII scores with
Natural Capital modulation applied. The modulation adjusts EII by up to
±0.05 based on three KPI dimensions: biodiversity, soil, and water.
"""

from __future__ import annotations

from typing import Any, Literal

import ee

from eii.compute.modulation import (
    BIODIVERSITY_MAX_THRESHOLD,
    MODULATION_RANGE,
    apply_modulation,
    apply_modulation_image,
    calculate_biodiversity_kpi,
    calculate_nc_score,
    calculate_nc_score_image,
    calculate_soil_kpi,
    calculate_water_kpi,
)

from .retrieve import get_layers, get_stats
from .utils import normalize_client_input

AggregationMethod = Literal["min_fuzzy_logic", "minimum", "product", "geometric_mean"]
ComputeMode = Literal["precomputed", "on_the_fly"]
OutputFormat = Literal["dict", "geodataframe"]


def get_modulated_eii(
    geometry: ee.Geometry | Any,
    kpis: dict[str, float] | None = None,
    kpi_layers: dict[str, ee.Image] | None = None,
    compute_default_kpis: bool = False,
    aggregation_method: AggregationMethod = "min_fuzzy_logic",
    compute_mode: ComputeMode = "precomputed",
    include_components: bool = True,
    kpi_weights: dict[str, float] | None = None,
    modulation_range: float = MODULATION_RANGE,
    biodiversity_max_threshold: float = BIODIVERSITY_MAX_THRESHOLD,
    scale: int = 300,
    stats: list[str] | None = None,
    percentiles: list[int] | None = None,
    output_format: OutputFormat = "dict",
) -> dict | Any:
    """
    Get EII with Natural Capital modulation applied.

    The Natural Capital (NC) score is calculated from three KPI dimensions
    (biodiversity, soil, water) and used to modulate the base EII by ±0.05.

    Three ways to provide KPIs:
    1. `kpis`: Pre-computed normalized values (0-1) - fastest
    2. `kpi_layers`: User-provided ee.Image layers (0-1) - flexible
    3. `compute_default_kpis=True`: Compute from GEE datasets - demo/fallback

    Args:
        geometry: AOI geometry (ee.Geometry, Feature, FeatureCollection,
            GeoPandas, shapely, or bbox tuple)
        kpis: Pre-computed KPI values as dict with keys "biodiversity",
            "soil", "water" and float values 0-1
        kpi_layers: User-provided ee.Image layers (must be normalized 0-1)
            with keys "biodiversity", "soil", "water"
        compute_default_kpis: If True, compute KPIs from built-in GEE datasets
        aggregation_method: EII aggregation method
        compute_mode: "precomputed" or "on_the_fly" for base EII
        include_components: Include EII component layers in output
        kpi_weights: Custom KPI weights (must sum to 1.0)
        modulation_range: Total modulation range (default 0.1 = ±0.05)
        biodiversity_max_threshold: Upper threshold for biodiversity KPI
            normalization (default 0.5 = 50% natural habitat)
        scale: Resolution in meters for statistics
        stats: Statistics to compute ("mean", "median", "min", "max", "std")
        percentiles: Optional list of percentiles (0-100)
        output_format: "dict" or "geodataframe"

    Returns:
        Dictionary or GeoDataFrame with:
        - eii: Base EII statistics
        - eii_modulated: Modulated EII statistics
        - nc_score: Natural Capital score
        - kpis: Individual KPI values (biodiversity, soil, water)
        - Component statistics if include_components=True

    Example:
        >>> import ee
        >>> from eii.client import get_modulated_eii
        >>> ee.Initialize()
        >>> aoi = ee.Geometry.Rectangle([10.5, 47.5, 11.0, 48.0])
        >>> result = get_modulated_eii(aoi, compute_default_kpis=True)
        >>> print(f"Base EII: {result['values']['eii']['mean']:.3f}")
        >>> print(f"Modulated: {result['values']['eii_modulated']['mean']:.3f}")
    """
    if stats is None:
        stats = ["mean"]

    kpi_source_count = sum([kpis is not None, kpi_layers is not None, compute_default_kpis])
    if kpi_source_count == 0:
        raise ValueError("Must provide one of: kpis, kpi_layers, or compute_default_kpis=True")
    if kpi_source_count > 1:
        raise ValueError("Provide only one of: kpis, kpi_layers, or compute_default_kpis=True")

    geometry = normalize_client_input(geometry, target="geometry")

    base_result = get_stats(
        geometry=geometry,
        aggregation_method=aggregation_method,
        include_components=include_components,
        scale=scale,
        stats=stats,
        percentiles=percentiles,
        compute_mode=compute_mode,
        output_format="dict",
    )

    if kpis is not None:
        kpi_values = _validate_kpi_dict(kpis)
        nc_score = calculate_nc_score(
            biodiversity_kpi=kpi_values["biodiversity"],
            soil_kpi=kpi_values["soil"],
            water_kpi=kpi_values["water"],
            weights=kpi_weights,
        )

        base_eii_mean = base_result["values"]["eii"].get("mean")
        if base_eii_mean is None:
            base_eii_mean = base_result["values"]["eii"].get("value")

        eii_modulated = apply_modulation(base_eii_mean, nc_score, modulation_range)

        result = _build_result_from_scalar_kpis(
            base_result=base_result,
            kpi_values=kpi_values,
            nc_score=nc_score,
            eii_modulated=eii_modulated,
            stats=stats,
        )

    elif kpi_layers is not None:
        kpi_images = _validate_kpi_layers(kpi_layers)
        result = _compute_modulation_from_layers(
            geometry=geometry,
            base_result=base_result,
            kpi_images=kpi_images,
            kpi_weights=kpi_weights,
            modulation_range=modulation_range,
            scale=scale,
            stats=stats,
            percentiles=percentiles,
        )

    else:
        kpi_images = get_kpi_layers(
            geometry=geometry,
            biodiversity_max_threshold=biodiversity_max_threshold,
        )
        result = _compute_modulation_from_layers(
            geometry=geometry,
            base_result=base_result,
            kpi_images=kpi_images,
            kpi_weights=kpi_weights,
            modulation_range=modulation_range,
            scale=scale,
            stats=stats,
            percentiles=percentiles,
        )

    if output_format == "geodataframe":
        return _to_geodataframe(geometry, result)

    return result


def get_default_kpis(
    geometry: ee.Geometry | Any,
    biodiversity_max_threshold: float = BIODIVERSITY_MAX_THRESHOLD,
    scale: int = 100,
) -> dict[str, float]:
    """
    Compute default KPIs using built-in GEE datasets.

    Uses:
    - ESA WorldCover 10m for biodiversity (natural habitat fraction)
    - SoilGrids 250m + Köppen zones for soil (SOC relative to reference)
    - SoilGrids 250m for water (AWC relative to texture maximum)

    Args:
        geometry: AOI geometry
        biodiversity_max_threshold: Upper threshold for biodiversity
            normalization (default 0.5 = 50% natural habitat)
        scale: Resolution for reduction (default 100m)

    Returns:
        Dictionary with keys "biodiversity", "soil", "water" and
        float values 0-1

    Example:
        >>> kpis = get_default_kpis(my_polygon)
        >>> print(kpis)
        {'biodiversity': 0.42, 'soil': 0.61, 'water': 0.55}
    """
    geometry = normalize_client_input(geometry, target="geometry")

    kpi_images = get_kpi_layers(geometry, biodiversity_max_threshold)

    kpi_values = {}
    for name, image in kpi_images.items():
        value = image.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=geometry,
            scale=scale,
            maxPixels=1e12,
        ).getInfo()

        kpi_values[name] = value.get(f"{name}_kpi") or value.get(name)

    return kpi_values


def get_kpi_layers(
    geometry: ee.Geometry | Any,
    biodiversity_max_threshold: float = BIODIVERSITY_MAX_THRESHOLD,
) -> dict[str, ee.Image]:
    """
    Get default KPI raster layers for an AOI.

    Returns normalized (0-1) ee.Image for each KPI:
    - biodiversity: Natural habitat fraction (100m resolution)
    - soil: SOC relative to climate reference (250m resolution)
    - water: AWC relative to texture maximum (250m resolution)

    Args:
        geometry: AOI geometry
        biodiversity_max_threshold: Upper threshold for biodiversity
            normalization (default 0.5 = 50% natural habitat)

    Returns:
        Dictionary with keys "biodiversity", "soil", "water" and
        ee.Image values (0-1)

    Example:
        >>> layers = get_kpi_layers(my_polygon)
        >>> biodiv_image = layers["biodiversity"]
        >>> # Use for visualization or custom aggregation
    """
    geometry = normalize_client_input(geometry, target="geometry")

    biodiversity = calculate_biodiversity_kpi(
        aoi=geometry,
        max_threshold=biodiversity_max_threshold,
    )

    soil = calculate_soil_kpi(aoi=geometry)

    water = calculate_water_kpi(aoi=geometry)

    return {
        "biodiversity": biodiversity,
        "soil": soil,
        "water": water,
    }


def get_nc_score(
    geometry: ee.Geometry | Any,
    kpis: dict[str, float] | None = None,
    kpi_layers: dict[str, ee.Image] | None = None,
    compute_default_kpis: bool = False,
    kpi_weights: dict[str, float] | None = None,
    biodiversity_max_threshold: float = BIODIVERSITY_MAX_THRESHOLD,
    scale: int = 100,
) -> float:
    """
    Calculate Natural Capital score from KPI inputs.

    Args:
        geometry: AOI geometry
        kpis: Pre-computed KPI values (0-1)
        kpi_layers: User-provided ee.Image layers (0-1)
        compute_default_kpis: Compute from GEE datasets
        kpi_weights: Custom weights (must sum to 1.0)
        biodiversity_max_threshold: Upper threshold for biodiversity
        scale: Resolution for reduction

    Returns:
        Natural Capital score (0-1)
    """
    if kpis is not None:
        kpi_values = _validate_kpi_dict(kpis)
    elif kpi_layers is not None:
        geometry = normalize_client_input(geometry, target="geometry")
        kpi_values = {}
        for name, image in kpi_layers.items():
            value = image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=geometry,
                scale=scale,
                maxPixels=1e12,
            ).getInfo()
            kpi_values[name] = list(value.values())[0]
    elif compute_default_kpis:
        kpi_values = get_default_kpis(
            geometry=geometry,
            biodiversity_max_threshold=biodiversity_max_threshold,
            scale=scale,
        )
    else:
        raise ValueError("Must provide one of: kpis, kpi_layers, or compute_default_kpis=True")

    return calculate_nc_score(
        biodiversity_kpi=kpi_values["biodiversity"],
        soil_kpi=kpi_values["soil"],
        water_kpi=kpi_values["water"],
        weights=kpi_weights,
    )


def _validate_kpi_dict(kpis: dict[str, float]) -> dict[str, float]:
    """Validate pre-computed KPI dictionary."""
    required_keys = {"biodiversity", "soil", "water"}
    missing = required_keys - set(kpis.keys())
    if missing:
        raise ValueError(f"Missing required KPI keys: {missing}")

    for key in required_keys:
        value = kpis[key]
        if not isinstance(value, int | float):
            raise ValueError(f"KPI '{key}' must be numeric, got {type(value)}")
        if not 0 <= value <= 1:
            raise ValueError(f"KPI '{key}' must be in [0, 1], got {value}")

    return {k: kpis[k] for k in required_keys}


def _validate_kpi_layers(kpi_layers: dict[str, ee.Image]) -> dict[str, ee.Image]:
    """Validate user-provided KPI layers."""
    required_keys = {"biodiversity", "soil", "water"}
    missing = required_keys - set(kpi_layers.keys())
    if missing:
        raise ValueError(f"Missing required KPI layer keys: {missing}")

    return {k: kpi_layers[k] for k in required_keys}


def _compute_modulation_from_layers(
    geometry: ee.Geometry,
    base_result: dict,
    kpi_images: dict[str, ee.Image],
    kpi_weights: dict[str, float] | None,
    modulation_range: float,
    scale: int,
    stats: list[str],
    percentiles: list[int] | None,
) -> dict:
    """Compute modulation from KPI raster layers."""
    from .retrieve import _build_reducer, _format_stats

    nc_image = calculate_nc_score_image(
        biodiversity_kpi=kpi_images["biodiversity"],
        soil_kpi=kpi_images["soil"],
        water_kpi=kpi_images["water"],
        weights=kpi_weights,
    )

    layers = get_layers(layers="eii", compute_mode="precomputed")
    eii_base = layers["eii"]

    eii_modulated = apply_modulation_image(
        eii=eii_base,
        nc_score=nc_image,
        modulation_range=modulation_range,
    )

    reducer = _build_reducer(stats, percentiles)

    kpi_stats = {}
    for name, image in kpi_images.items():
        raw = image.reduceRegion(
            reducer=reducer,
            geometry=geometry,
            scale=scale,
            maxPixels=1e12,
        ).getInfo()
        kpi_stats[name] = _format_stats(f"{name}_kpi", raw, False, stats, percentiles)

    nc_raw = nc_image.reduceRegion(
        reducer=reducer,
        geometry=geometry,
        scale=scale,
        maxPixels=1e12,
    ).getInfo()

    mod_raw = eii_modulated.reduceRegion(
        reducer=reducer,
        geometry=geometry,
        scale=scale,
        maxPixels=1e12,
    ).getInfo()

    values = dict(base_result.get("values", {}))
    values["eii_modulated"] = _format_stats(
        "eii_modulated", mod_raw, False, stats, percentiles
    ).get("eii_modulated", {})
    values["nc_score"] = _format_stats("nc_score", nc_raw, False, stats, percentiles).get(
        "nc_score", {}
    )

    for name, stat_dict in kpi_stats.items():
        values[f"{name}_kpi"] = stat_dict.get(f"{name}_kpi", {})

    return {
        "geometry_type": base_result.get("geometry_type"),
        "values": values,
    }


def _build_result_from_scalar_kpis(
    base_result: dict,
    kpi_values: dict[str, float],
    nc_score: float,
    eii_modulated: float,
    stats: list[str],
) -> dict:
    """Build result dict when KPIs are provided as scalars."""
    values = dict(base_result.get("values", {}))

    if "mean" in stats:
        values["eii_modulated"] = {"mean": eii_modulated}
        values["nc_score"] = {"mean": nc_score}
        for name, value in kpi_values.items():
            values[f"{name}_kpi"] = {"mean": value}
    else:
        stat_key = stats[0]
        values["eii_modulated"] = {stat_key: eii_modulated}
        values["nc_score"] = {stat_key: nc_score}
        for name, value in kpi_values.items():
            values[f"{name}_kpi"] = {stat_key: value}

    return {
        "geometry_type": base_result.get("geometry_type"),
        "values": values,
    }


def _to_geodataframe(geometry: ee.Geometry, result: dict):
    """Convert result to GeoDataFrame."""
    try:
        import geopandas as gpd
        import pandas as pd
        from shapely.geometry import shape
    except ImportError as exc:
        raise RuntimeError("geopandas and shapely required for geodataframe output") from exc

    geom_info = geometry.getInfo()
    geom_shape = shape(geom_info)

    values = result.get("values", {})
    data_rows = []
    columns = []

    for metric, metric_value in values.items():
        if isinstance(metric_value, dict):
            for stat, value in metric_value.items():
                columns.append((metric, stat))
                data_rows.append(value)
        else:
            columns.append((metric, "value"))
            data_rows.append(metric_value)

    data = pd.DataFrame([data_rows], columns=pd.MultiIndex.from_tuples(columns))

    return gpd.GeoDataFrame(data, geometry=[geom_shape], crs="EPSG:4326")
