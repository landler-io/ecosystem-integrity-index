"""
Core retrieval functions for EII data.

This module provides the main interface for extracting Ecosystem Integrity
Index data from pre-computed GEE assets.
"""

from __future__ import annotations

import os
import shutil
import tempfile
import urllib.request
import zipfile
from pathlib import Path
from typing import Literal

import ee

from eii.compute.compositional import calculate_compositional_integrity
from eii.compute.integrity import combine_components
from eii.compute.npp import calculate_functional_integrity
from eii.compute.settings import (
    NATURAL_NPP_ASSET_PATH,
    OBSERVED_NPP_ASSET_PATH,
    OBSERVED_NPP_YEAR_RANGE,
)

from .assets import get_asset_path
from .utils import normalize_client_input

DEFAULT_BII_YEAR = 2020


AggregationMethod = Literal["min_fuzzy_logic", "minimum", "product", "geometric_mean"]
ComputeMode = Literal["precomputed", "on_the_fly"]
OutputFormat = Literal["dict", "geodataframe"]
LayerSelection = Literal["eii", "components", "all"]


def get_layers(
    layers: LayerSelection = "all",
    aggregation_method: AggregationMethod = "min_fuzzy_logic",
    compute_mode: ComputeMode = "precomputed",
    geometry: ee.Geometry | None = None,
    npp_year_range: list[str] = OBSERVED_NPP_YEAR_RANGE,
    include_seasonality: bool = True,
    natural_npp_asset_path: str = NATURAL_NPP_ASSET_PATH,
    observed_npp_asset_path: str = OBSERVED_NPP_ASSET_PATH,
    structural_asset_path: str | None = None,
    bii_year: int = DEFAULT_BII_YEAR,
    absolute_diff_percentile: str = "p95",
) -> dict[str, ee.Image]:
    """
    Get EII and/or component layers for a given geometry and computation mode.

    Args:
        layers: Which layers to return: "eii", "components", or "all".
        aggregation_method: EII aggregation method.
        compute_mode: "precomputed" (default) or "on_the_fly".
        geometry: Optional geometry for on-the-fly component calculation.
        npp_year_range: Date range for observed NPP (on-the-fly mode).
        include_seasonality: Include seasonality in functional integrity.
        natural_npp_asset_path: Asset for pre-computed natural NPP.
        observed_npp_asset_path: Asset for observed NPP annual tiles.
        structural_asset_path: Asset for structural integrity core area.
        bii_year: BII year to select for compositional integrity.
        absolute_diff_percentile: Percentile key for NPP absolute difference penalty.

    Returns:
        Dictionary of Earth Engine images. Keys include "eii" and/or component
        names ("functional", "structural", "compositional").

    Example:
        >>> layers = get_layers(layers="all")
        >>> eii = layers["eii"].clip(my_polygon)
    """
    if geometry is not None:
        geometry = normalize_client_input(geometry, target="geometry")

    if compute_mode == "precomputed" and aggregation_method != "min_fuzzy_logic":
        raise ValueError(
            "Only the pre-computed default method is available. "
            "Use compute_mode='on_the_fly' for other aggregation methods."
        )

    include_eii = layers in ("eii", "all")
    include_components = layers in ("components", "all")

    result: dict[str, ee.Image] = {}
    components: dict[str, ee.Image] | None = None

    if include_components:
        if compute_mode == "precomputed":
            mosaic = ee.Image(get_asset_path("components"))
            components = {
                "functional": mosaic.select("functional_integrity"),
                "structural": mosaic.select("structural_integrity"),
                "compositional": mosaic.select("compositional_integrity"),
            }
        elif compute_mode == "on_the_fly":
            functional_result = calculate_functional_integrity(
                aoi=geometry,
                year_range=npp_year_range,
                include_seasonality=include_seasonality,
                natural_npp_asset_path=natural_npp_asset_path,
                observed_npp_asset_path=observed_npp_asset_path,
                absolute_diff_percentile=absolute_diff_percentile,
            )
            structural = _get_structural_integrity_precomputed(
                aoi=geometry,
                asset_path=structural_asset_path,
            )
            compositional = calculate_compositional_integrity(aoi=geometry, year=bii_year)
            components = {
                "functional": functional_result["functional_integrity"],
                "structural": structural,
                "compositional": compositional,
            }
        else:
            raise ValueError("compute_mode must be one of: 'precomputed', 'on_the_fly'")

        result.update(components)

    if include_eii:
        if compute_mode == "precomputed":
            result["eii"] = ee.Image(get_asset_path("eii")).select("eii")
        else:
            if components is None:
                components = get_layers(
                    layers="components",
                    aggregation_method=aggregation_method,
                    compute_mode=compute_mode,
                    geometry=geometry,
                    npp_year_range=npp_year_range,
                    include_seasonality=include_seasonality,
                    natural_npp_asset_path=natural_npp_asset_path,
                    observed_npp_asset_path=observed_npp_asset_path,
                    structural_asset_path=structural_asset_path,
                    bii_year=bii_year,
                    absolute_diff_percentile=absolute_diff_percentile,
                )
            result["eii"] = combine_components(
                functional=components["functional"],
                structural=components["structural"],
                compositional=components["compositional"],
                method=aggregation_method,
            ).rename("eii")

    return result


def _get_structural_integrity_precomputed(
    aoi: ee.Geometry | None = None,
    asset_path: str | None = None,
) -> ee.Image:
    """Load pre-computed structural integrity (core area) image."""
    if asset_path is None:
        asset_path = get_asset_path("structural_core_area")

    structural = ee.Image(asset_path).rename("structural_integrity")
    if aoi is not None:
        structural = structural.clip(aoi)

    return structural


def _is_point_geometry(geometry: ee.Geometry) -> ee.Boolean:
    """Return True if geometry is a point or multipoint."""
    return ee.List(["Point", "MultiPoint"]).contains(geometry.type())


def _reduce_point_values(
    image: ee.Image,
    geometry: ee.Geometry,
    scale: int,
) -> ee.Dictionary:
    """Extract single values for point geometry."""
    return image.reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=geometry,
        scale=scale,
        maxPixels=1e12,
    )


VALID_STATS = {"mean", "median", "min", "max", "std"}


def _build_reducer(
    stats: list[str],
    percentiles: list[int] | None = None,
) -> ee.Reducer:
    """Build a combined reducer from requested stats and percentiles."""
    stats_set = set(stats)
    reducers: list[ee.Reducer] = []

    if "mean" in stats_set:
        reducers.append(ee.Reducer.mean())
    if "median" in stats_set:
        reducers.append(ee.Reducer.median())
    if "min" in stats_set:
        reducers.append(ee.Reducer.min())
    if "max" in stats_set:
        reducers.append(ee.Reducer.max())
    if "std" in stats_set:
        reducers.append(ee.Reducer.stdDev())
    if percentiles:
        reducers.append(ee.Reducer.percentile(percentiles))

    if not reducers:
        raise ValueError("No statistics selected. Provide stats and/or percentiles.")

    reducer = reducers[0]
    for next_reducer in reducers[1:]:
        reducer = reducer.combine(next_reducer, sharedInputs=True)

    return reducer


def _validate_stats_params(
    stats: list[str] | None,
    percentiles: list[int] | None,
) -> list[str]:
    """Validate stats and percentiles parameters, return stats list."""
    if stats is None:
        stats = ["mean"]

    invalid_stats = [s for s in stats if s not in VALID_STATS]
    if invalid_stats:
        raise ValueError(f"Invalid stats: {invalid_stats}. Supported: {sorted(VALID_STATS)}")

    if percentiles:
        invalid_percentiles = [p for p in percentiles if not isinstance(p, int) or p < 0 or p > 100]
        if invalid_percentiles:
            raise ValueError(
                f"Percentiles must be integers between 0 and 100. Invalid: {invalid_percentiles}"
            )

    return stats


def _reduce_area_stats(
    image: ee.Image,
    geometry: ee.Geometry,
    scale: int,
    stats: list[str] | None = None,
    percentiles: list[int] | None = None,
) -> ee.Dictionary:
    """Extract summary statistics for area geometries."""
    if stats is None:
        stats = ["mean"]
    reducer = _build_reducer(stats, percentiles)
    return image.reduceRegion(
        reducer=reducer,
        geometry=geometry,
        scale=scale,
        maxPixels=1e12,
    )


def _format_stats(
    label: str,
    raw_stats: dict,
    is_point: bool,
    requested_stats: list[str] | None = None,
    percentiles: list[int] | None = None,
) -> dict:
    """Format stats into the nested structure."""
    if is_point:
        return {label: raw_stats.get(label)}

    if requested_stats is None:
        requested_stats = ["mean"]

    stat_mapping = {
        "mean": "mean",
        "min": "min",
        "max": "max",
        "median": "median",
        "std": "stdDev",
    }

    result: dict = {}

    # When there's only one stat and no percentiles, GEE returns just the band name
    # without a suffix (e.g., "eii" instead of "eii_mean")
    single_reducer = len(requested_stats) == 1 and not percentiles

    for stat in requested_stats:
        if single_reducer:
            result[stat] = raw_stats.get(label)
        else:
            suffix = stat_mapping.get(stat, stat)
            result[stat] = raw_stats.get(f"{label}_{suffix}")

    if percentiles:
        for p in percentiles:
            result[f"p{p}"] = raw_stats.get(f"{label}_p{p}")

    return {label: result}


def _to_geodataframe(
    geometry: ee.Geometry,
    properties: dict,
):
    try:
        import geopandas as gpd
        import pandas as pd
        from shapely.geometry import shape
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("geopandas and shapely are required for geodataframe output") from exc

    geom_info = geometry.getInfo()
    geom_shape = shape(geom_info)
    _ = properties.pop("geometry_type", None)
    data_rows = []
    columns = []
    for metric, metric_value in properties.items():
        if isinstance(metric_value, dict):
            for stat, value in metric_value.items():
                columns.append((metric, stat))
                data_rows.append(value)
        else:
            columns.append((metric, "value"))
            data_rows.append(metric_value)

    data = pd.DataFrame([data_rows], columns=pd.MultiIndex.from_tuples(columns))

    return gpd.GeoDataFrame(data, geometry=[geom_shape], crs="EPSG:4326")


def _available_memory_bytes() -> int | None:
    """Best-effort available memory for auto chunking decisions."""
    try:
        import psutil  # type: ignore

        return int(psutil.virtual_memory().available)
    except Exception:
        pass

    try:
        pages = os.sysconf("SC_AVPHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        if isinstance(pages, int) and isinstance(page_size, int):
            return int(pages * page_size)
    except Exception:
        pass

    return None


def _estimate_area_sq_km(geometry: ee.Geometry, max_error: float = 1) -> float:
    """Estimate AOI area using bounds to keep server-side geometry simple."""
    bounds = geometry.bounds(max_error)
    area_m2 = bounds.area(max_error)
    return float(ee.Number(area_m2).divide(1_000_000).getInfo())


def _estimate_bytes(
    area_sq_km: float,
    scale: int,
    band_count: int,
    dtype: str,
) -> float:
    bytes_per_value = {"float32": 4, "float64": 8}.get(dtype, 8)
    pixels = (area_sq_km * 1_000_000) / float(scale * scale)
    return pixels * band_count * bytes_per_value


def _download_image_to_geotiff(
    image: ee.Image,
    region: ee.Geometry,
    scale: int,
    crs: str,
    out_path: Path,
) -> Path:
    params = {
        "region": region,
        "scale": scale,
        "crs": crs,
        "fileFormat": "GEO_TIFF",
        "filePerBand": False,
    }
    url = image.getDownloadURL(params)

    tmp_dir = Path(tempfile.mkdtemp())
    try:
        tmp_download = tmp_dir / "download"
        with urllib.request.urlopen(url) as response, open(tmp_download, "wb") as handle:
            handle.write(response.read())

        if zipfile.is_zipfile(tmp_download):
            with zipfile.ZipFile(tmp_download) as zf:
                zf.extractall(tmp_dir)
            tif_paths = list(tmp_dir.glob("*.tif"))
            if not tif_paths:
                raise RuntimeError("Download completed but no GeoTIFF was found.")
            source_path = tif_paths[0]
        else:
            source_path = tmp_download
            if source_path.suffix.lower() != ".tif":
                source_path = source_path.with_suffix(".tif")
                tmp_download.rename(source_path)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source_path), out_path)
        return out_path
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _geotiff_to_xarray_dataset(
    path: Path,
    band_names: list[str],
    dtype: str,
):
    try:
        import numpy as np
        import rasterio
        import xarray as xr
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("xarray, rasterio, and numpy are required for raster downloads") from exc

    with rasterio.open(path) as src:
        data = src.read(masked=True)
        if dtype:
            data = data.astype(dtype)
        if isinstance(data, np.ma.MaskedArray):
            data = data.filled(np.nan)

        if data.shape[0] != len(band_names):
            band_names = [f"band_{idx + 1}" for idx in range(data.shape[0])]

        transform = src.transform
        width = src.width
        height = src.height
        x = transform.c + transform.a * (np.arange(width) + 0.5)
        y = transform.f + transform.e * (np.arange(height) + 0.5)

        dataset = xr.Dataset(coords={"x": x, "y": y})
        for idx, name in enumerate(band_names):
            dataset[name] = (("y", "x"), data[idx])

        dataset.attrs["crs"] = str(src.crs) if src.crs else None
        dataset.attrs["transform"] = tuple(transform)
        dataset.attrs["resolution"] = (transform.a, transform.e)
        return dataset


def _write_dataset_to_geotiff(
    dataset,
    path: Path,
    compression: str | None,
):
    try:
        import rasterio
        import rioxarray  # noqa: F401
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("rioxarray is required for GeoTIFF output") from exc

    crs = dataset.attrs.get("crs")
    transform = dataset.attrs.get("transform")
    if crs:
        dataset = dataset.rio.write_crs(crs, inplace=False)
    if transform:
        if isinstance(transform, list | tuple):
            if len(transform) == 6:
                transform = rasterio.Affine.from_gdal(*transform)
            else:
                transform = rasterio.Affine(*transform[:6])
        dataset = dataset.rio.write_transform(transform, inplace=False)

    kwargs = {"BIGTIFF": "YES"}
    if compression:
        kwargs["compress"] = compression
    dataset.rio.to_raster(path, **kwargs)

    band_names = list(dataset.data_vars)
    if band_names:
        try:
            with rasterio.open(path, "r+") as dst:
                if dst.count == len(band_names):
                    dst.descriptions = tuple(band_names)
        except Exception:
            pass


def _stream_tiles_to_geotiff(
    tile_paths: list[Path],
    out_path: Path,
    bounds: list[float],
    compression: str | None,
    band_names: list[str] | None = None,
):
    try:
        import rasterio
        from rasterio.windows import from_bounds
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("rasterio is required for GeoTIFF output") from exc

    if not tile_paths:
        raise ValueError("No tiles provided for streaming output.")

    with rasterio.open(tile_paths[0]) as first:
        res_x = first.transform.a
        res_y = first.transform.e
        count = first.count
        dtype = first.dtypes[0]
        crs = first.crs
        blockxsize = first.width
        blockysize = first.height

    min_lon, min_lat, max_lon, max_lat = bounds
    width = int(round((max_lon - min_lon) / res_x))
    height = int(round((max_lat - min_lat) / abs(res_y)))
    transform = rasterio.transform.from_bounds(min_lon, min_lat, max_lon, max_lat, width, height)

    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": count,
        "dtype": dtype,
        "crs": crs,
        "transform": transform,
        "tiled": True,
        "blockxsize": blockxsize,
        "blockysize": blockysize,
    }
    if compression:
        profile["compress"] = compression
        profile["BIGTIFF"] = "YES"
    else:
        profile["BIGTIFF"] = "YES"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as dst:
        for tile_path in tile_paths:
            with rasterio.open(tile_path) as src:
                src_bounds = src.bounds
                window = from_bounds(*src_bounds, transform=transform)
                data = src.read()
                dst.write(data, window=window)

        if band_names and dst.count == len(band_names):
            dst.descriptions = tuple(band_names)


def get_stats(
    geometry: ee.Geometry,
    aggregation_method: AggregationMethod = "min_fuzzy_logic",
    include_components: bool = True,
    scale: int = 300,
    stats: list[str] | None = None,
    percentiles: list[int] | None = None,
    compute_mode: ComputeMode = "precomputed",
    npp_year_range: list[str] = OBSERVED_NPP_YEAR_RANGE,
    include_seasonality: bool = True,
    natural_npp_asset_path: str = NATURAL_NPP_ASSET_PATH,
    observed_npp_asset_path: str = OBSERVED_NPP_ASSET_PATH,
    structural_asset_path: str | None = None,
    bii_year: int = DEFAULT_BII_YEAR,
    absolute_diff_percentile: str = "p95",
    output_format: OutputFormat = "dict",
):
    """
    Extract EII values or statistics for a given geometry.

    Args:
        geometry: AOI geometry (ee.Geometry/Feature/FeatureCollection, GeoPandas,
            shapely, or bbox tuple).
        aggregation_method: EII aggregation method. Defaults to "min_fuzzy_logic".
        include_components: If True, include statistics for functional,
            structural, and compositional integrity. Defaults to True.
        scale: Resolution in meters for the reduction. Defaults to 300.
        stats: Statistics to compute. Supported: "mean", "median", "min", "max",
            "std". Defaults to ["mean"].
        percentiles: Optional list of percentiles (0-100) to include.
        compute_mode: "precomputed" (default) or "on_the_fly".
        npp_year_range: Date range for observed NPP (on-the-fly mode).
        include_seasonality: Include seasonality in functional integrity.
        natural_npp_asset_path: Asset for pre-computed natural NPP.
        observed_npp_asset_path: Asset for observed NPP annual tiles.
        structural_asset_path: Asset for structural integrity core area.
        bii_year: BII year to select for compositional integrity.
        absolute_diff_percentile: Percentile key for NPP absolute difference penalty.
        output_format: "dict" (default) or "geodataframe".

    Returns:
        A structured dictionary or GeoDataFrame depending on output_format.
        Point geometries return single values; area geometries return summary
        statistics in nested format: {"eii": {"mean": ..., "min": ...}}.

    Example:
        >>> import ee
        >>> ee.Initialize()
        >>> polygon = ee.Geometry.Rectangle([-60, -10, -55, -5])
        >>> stats = get_stats(polygon, stats=["mean", "min", "max"])
        >>> print(stats)
    """
    stats = _validate_stats_params(stats, percentiles)
    geometry = normalize_client_input(geometry, target="geometry")
    layers = get_layers(
        layers="all" if include_components else "eii",
        aggregation_method=aggregation_method,
        compute_mode=compute_mode,
        geometry=geometry,
        npp_year_range=npp_year_range,
        include_seasonality=include_seasonality,
        natural_npp_asset_path=natural_npp_asset_path,
        observed_npp_asset_path=observed_npp_asset_path,
        structural_asset_path=structural_asset_path,
        bii_year=bii_year,
        absolute_diff_percentile=absolute_diff_percentile,
    )
    eii_image = layers["eii"]

    is_point = _is_point_geometry(geometry).getInfo()
    geometry_type = geometry.type().getInfo()

    if is_point:
        eii_stats = _reduce_point_values(eii_image, geometry, scale).getInfo()
        values = _format_stats("eii", eii_stats, True, stats, percentiles)
    else:
        eii_stats = _reduce_area_stats(eii_image, geometry, scale, stats, percentiles).getInfo()
        values = _format_stats("eii", eii_stats, False, stats, percentiles)

    if include_components:
        components = {k: v for k, v in layers.items() if k != "eii"}
        for name, image in components.items():
            component_label = name if name.endswith("_integrity") else f"{name}_integrity"
            if is_point:
                comp_stats = _reduce_point_values(image, geometry, scale).getInfo()
                values.update(_format_stats(component_label, comp_stats, True, stats, percentiles))
            else:
                comp_stats = _reduce_area_stats(
                    image, geometry, scale, stats, percentiles
                ).getInfo()
                values.update(_format_stats(component_label, comp_stats, False, stats, percentiles))

    result = {"geometry_type": geometry_type, "values": values}

    if output_format == "geodataframe":
        geo_props = {"geometry_type": result["geometry_type"]}
        geo_props.update(result["values"])
        return _to_geodataframe(geometry, geo_props)

    return result


def get_raster(
    geometry: ee.Geometry,
    aggregation_method: AggregationMethod = "min_fuzzy_logic",
    include_components: bool = True,
    scale: int = 300,
    compute_mode: ComputeMode = "precomputed",
    npp_year_range: list[str] = OBSERVED_NPP_YEAR_RANGE,
    include_seasonality: bool = True,
    natural_npp_asset_path: str = NATURAL_NPP_ASSET_PATH,
    observed_npp_asset_path: str = OBSERVED_NPP_ASSET_PATH,
    structural_asset_path: str | None = None,
    bii_year: int = DEFAULT_BII_YEAR,
    absolute_diff_percentile: str = "p95",
    max_area_sq_km: float | None = None,
    max_pixels: int | None = None,
    chunking: Literal["auto", "never", "always"] = "auto",
    tile_size_deg: float = 1.0,
    tmp_dir: str | None = None,
    out_path: str | None = None,
    output_format: Literal["memory", "geotiff"] = "memory",
    compression: str | None = None,
    max_memory_fraction: float = 0.25,
    dtype: Literal["float32", "float64"] = "float64",
    download_crs: str = "EPSG:4326",
):
    """
    Download EII rasters for an AOI as an xarray Dataset.

    Args:
        geometry: AOI geometry (ee.Geometry/Feature/FeatureCollection, GeoPandas,
            shapely, or bbox tuple).
        aggregation_method: EII aggregation method.
        include_components: If True, include component bands alongside EII.
        scale: Resolution in meters for export. Defaults to 300.
        compute_mode: "precomputed" (default) or "on_the_fly".
        npp_year_range: Date range for observed NPP (on-the-fly mode).
        include_seasonality: Include seasonality in functional integrity.
        natural_npp_asset_path: Asset for pre-computed natural NPP.
        observed_npp_asset_path: Asset for observed NPP annual tiles.
        structural_asset_path: Asset for structural integrity core area.
        bii_year: BII year to select for compositional integrity.
        absolute_diff_percentile: Percentile key for NPP absolute difference penalty.
        max_area_sq_km: Optional upper bound for AOI area before chunking.
        max_pixels: Optional upper bound for total pixels before chunking.
        chunking: "auto" (default), "never", or "always".
        tile_size_deg: Tile size for chunked downloads in degrees.
        tmp_dir: Optional directory for intermediate downloads.
        out_path: Optional output path. For chunked downloads, use a directory
            (multiple files) or a file path for single GeoTIFF streaming output.
        output_format: "memory" (default) or "geotiff".
        compression: Compression codec for GeoTIFF (default ZSTD).
        max_memory_fraction: Fraction of available RAM used for auto chunking.
        dtype: Output dtype for the xarray Dataset.
        download_crs: CRS for downloads (default EPSG:4326).

    Band order:
        "eii", "functional_integrity", "structural_integrity", "compositional_integrity".

    Returns:
        xarray.Dataset when output_format="memory"; otherwise a Path (single download)
        or list[Path] (chunked downloads) to the written files.
    """
    geometry = normalize_client_input(geometry, target="geometry")
    if output_format == "geotiff" and compression is None:
        compression = "ZSTD"
    layers = get_layers(
        layers="all" if include_components else "eii",
        aggregation_method=aggregation_method,
        compute_mode=compute_mode,
        geometry=geometry,
        npp_year_range=npp_year_range,
        include_seasonality=include_seasonality,
        natural_npp_asset_path=natural_npp_asset_path,
        observed_npp_asset_path=observed_npp_asset_path,
        structural_asset_path=structural_asset_path,
        bii_year=bii_year,
        absolute_diff_percentile=absolute_diff_percentile,
    )

    band_images = [layers["eii"].rename("eii")]
    band_names = ["eii"]
    if include_components:
        for name in ("functional", "structural", "compositional"):
            band = layers[name].rename(f"{name}_integrity")
            band_images.append(band)
            band_names.append(f"{name}_integrity")

    image = ee.Image.cat(band_images).clip(geometry)
    area_sq_km = _estimate_area_sq_km(geometry)
    pixel_count = (area_sq_km * 1_000_000) / float(scale * scale)
    estimated_bytes = _estimate_bytes(area_sq_km, scale, len(band_names), dtype)

    available_bytes = _available_memory_bytes()
    auto_max_area = None
    if available_bytes is not None:
        auto_max_area = (available_bytes * max_memory_fraction) / _estimate_bytes(
            1.0, scale, len(band_names), dtype
        )

    effective_max_area = max_area_sq_km if max_area_sq_km is not None else auto_max_area

    exceeds_area = effective_max_area is not None and area_sq_km > effective_max_area
    exceeds_pixels = max_pixels is not None and pixel_count > max_pixels
    exceeds_memory = available_bytes is not None and estimated_bytes > (
        available_bytes * max_memory_fraction
    )

    if chunking == "never" and (exceeds_area or exceeds_pixels or exceeds_memory):
        raise ValueError("AOI exceeds configured limits; enable chunking or increase limits.")

    should_chunk = chunking == "always" or (
        chunking == "auto" and (exceeds_area or exceeds_pixels or exceeds_memory)
    )

    if not should_chunk:
        output_dir = Path(tmp_dir) if tmp_dir else Path(tempfile.mkdtemp())
        output_path = Path(out_path) if out_path else output_dir / "eii_download.tif"
        if output_path.is_dir():
            output_path = output_path / "eii_download.tif"

        _download_image_to_geotiff(image, geometry, scale, download_crs, output_path)
        if output_format == "memory":
            try:
                dataset = _geotiff_to_xarray_dataset(output_path, band_names, dtype)
            finally:
                if out_path is None:
                    shutil.rmtree(output_dir, ignore_errors=True)
        else:
            dataset = None

        if output_format == "geotiff":
            final_path = Path(out_path) if out_path else output_path
            if compression:
                dataset = _geotiff_to_xarray_dataset(output_path, band_names, dtype)
                _write_dataset_to_geotiff(dataset, final_path, compression)
                if final_path != output_path and output_path.exists():
                    output_path.unlink(missing_ok=True)
            return final_path

        dataset.attrs.update(
            {
                "eii_method": aggregation_method,
                "scale_m": scale,
                "area_sq_km": area_sq_km,
            }
        )
        return dataset

    if out_path is not None and output_format != "geotiff" and not Path(out_path).is_dir():
        raise ValueError("Chunked downloads require out_path to be a directory.")

    stream_single_geotiff = (
        output_format == "geotiff" and out_path is not None and not Path(out_path).is_dir()
    )
    output_dir = (
        Path(out_path)
        if out_path and Path(out_path).is_dir()
        else Path(tmp_dir or tempfile.mkdtemp())
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    bounds = geometry.bounds(1).getInfo()["coordinates"][0]
    lons = [coord[0] for coord in bounds]
    lats = [coord[1] for coord in bounds]
    min_lon, max_lon = min(lons), max(lons)
    min_lat, max_lat = min(lats), max(lats)

    datasets = []
    output_paths: list[Path] = []
    tile_idx = 0
    lon = min_lon
    while lon < max_lon:
        next_lon = min(lon + tile_size_deg, max_lon)
        lat = min_lat
        while lat < max_lat:
            next_lat = min(lat + tile_size_deg, max_lat)
            tile_geom = ee.Geometry.Rectangle([lon, lat, next_lon, next_lat])
            tile_path = output_dir / f"eii_tile_{tile_idx}.tif"
            _download_image_to_geotiff(image, tile_geom, scale, download_crs, tile_path)
            if stream_single_geotiff:
                output_paths.append(tile_path)
            elif output_format == "geotiff":
                if compression:
                    dataset = _geotiff_to_xarray_dataset(tile_path, band_names, dtype)
                    compressed_path = output_dir / f"eii_tile_{tile_idx}_compressed.tif"
                    _write_dataset_to_geotiff(dataset, compressed_path, compression)
                    tile_path.unlink(missing_ok=True)
                    output_paths.append(compressed_path)
                else:
                    output_paths.append(tile_path)
            else:
                datasets.append(_geotiff_to_xarray_dataset(tile_path, band_names, dtype))
            tile_idx += 1
            lat = next_lat
        lon = next_lon

    if output_format == "geotiff" and stream_single_geotiff:
        assert out_path is not None
        final_path = Path(out_path)
        _stream_tiles_to_geotiff(
            output_paths,
            final_path,
            [min_lon, min_lat, max_lon, max_lat],
            compression,
            band_names,
        )
        if tmp_dir is None:
            shutil.rmtree(output_dir, ignore_errors=True)
        else:
            for tile_path in output_paths:
                tile_path.unlink(missing_ok=True)
        return final_path

    if output_format == "geotiff":
        return output_paths

    try:
        import xarray as xr
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("xarray is required for raster downloads") from exc

    dataset = xr.combine_by_coords(datasets, combine_attrs="override")
    dataset.attrs.update(
        {
            "eii_method": aggregation_method,
            "scale_m": scale,
            "area_sq_km": area_sq_km,
            "chunked": True,
            "tile_size_deg": tile_size_deg,
        }
    )

    if out_path is None and tmp_dir is None:
        shutil.rmtree(output_dir, ignore_errors=True)

    return dataset
