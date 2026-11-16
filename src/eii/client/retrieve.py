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

from eii.compute.settings import REFERENCE_YEARS

from .assets import get_asset_path, resolve_version
from .utils import normalize_client_input

OutputFormat = Literal["dict", "geodataframe"]
LayerSelection = Literal[
    "eii", "components", "all", "limiting_component", "wetland_flag", "fi_change"
]
CATEGORICAL_BANDS = ("limiting_component", "wetland_flag")
FI_CHANGE_BANDS = ("fi_state", "fi_change", "fi_change_sensors")
COMPONENT_BANDS = {
    "functional": "functional_integrity",
    "structural": "structural_integrity",
    "compositional": "compositional_integrity",
}


def get_layers(
    layers: LayerSelection = "all",
    version: str | None = None,
    reference_year: int | None = None,
) -> dict[str, ee.Image]:
    """
    Get pre-computed EII and/or component layers.

    Args:
        layers: Which layers to return: "eii", "components", "all",
            "limiting_component" or "wetland_flag" (both v2 only). "all" includes
            both in v2. "fi_change" (v2 only, not part of "all") returns the
            parts of functional integrity: "fi_state", "fi_change" and "fi_change_sensors".
        version: EII dataset version (e.g. "v1", "v2"). If None, uses the active
            session version (see eii.client.set_version), defaulting to config.yaml.
        reference_year: v2 reference year, one of REFERENCE_YEARS (observed NPP of the 3
            years ending in it and BII of that year). Defaults to the latest year; pass it
            explicitly for reproducible results. Not used for v1.

    Returns:
        Dictionary of Earth Engine images. Keys include "eii" and/or component
        names ("functional", "structural", "compositional"). Version v2 also
        exposes "limiting_component": 1 = functional, 2 = structural,
        3 = compositional. Ties prefer functional, then structural. It also exposes
        "wetland_flag": 0 = none, 1 = wetland, 2 = peatland, 3 = mangrove, where
        potential NPP and so functional integrity are less reliable.
        With layers="fi_change": "fi_state" (observed against potential NPP),
        "fi_change" (EVI change against the pixel's baseline, masked where no sensor
        serves) and "fi_change_sensors" (uint8, 0-2 sensors; 0 where functional
        integrity equals fi_state). "functional" is fi_state x fi_change.

    Example:
        >>> layers = get_layers(layers="all", version="v2", reference_year=2024)
        >>> eii = layers["eii"].clip(my_polygon)
    """
    version = resolve_version(version)
    if layers not in ("eii", "components", "all", "fi_change", *CATEGORICAL_BANDS):
        raise ValueError(
            "layers must be one of: 'eii', 'components', 'all', 'limiting_component', "
            "'wetland_flag', 'fi_change'"
        )
    if layers in (*CATEGORICAL_BANDS, "fi_change") and version != "v2":
        raise ValueError(f"{layers} is available only for version='v2'")
    image = _release_image(version, reference_year)

    result: dict[str, ee.Image] = {}
    if layers in ("eii", "all"):
        result["eii"] = image.select("eii")
    if layers in ("components", "all"):
        result.update({name: image.select(band) for name, band in COMPONENT_BANDS.items()})
    for band in CATEGORICAL_BANDS:
        if version == "v2" and layers in (band, "all"):
            result[band] = image.select(band)
    if layers == "fi_change":
        result.update({band: image.select(band) for band in FI_CHANGE_BANDS})
        # float32 in the tile (one GeoTIFF data type), a sensor count 0-2
        result["fi_change_sensors"] = result["fi_change_sensors"].toUint8()
    return result


def _release_image(version: str, reference_year: int | None) -> ee.Image:
    """The released multiband EII image: v1 is static, v2 has one mosaic per reference year."""
    if version != "v2":
        if reference_year is not None:
            raise ValueError("reference_year applies to version='v2' only")
        image: ee.Image = ee.Image(get_asset_path("eii", version))
        return image
    if reference_year is None:
        reference_year = REFERENCE_YEARS[-1]
    if reference_year not in REFERENCE_YEARS:
        raise ValueError(f"reference_year must be one of {list(REFERENCE_YEARS)}")
    collection = ee.ImageCollection(get_asset_path("eii", version))
    mosaic: ee.Image = collection.filter(ee.Filter.eq("reference_year", reference_year)).mosaic()
    return mosaic


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
    result = image.reduceRegion(
        reducer=reducer,
        geometry=geometry,
        scale=scale,
        maxPixels=1e12,
    )
    # sub-pixel polygons have no pixel centers inside, so reduceRegion returns
    # all nulls. Fall back to sampling at the centroid with the same reducer.
    has_data = result.values().reduce(ee.Reducer.count())
    centroid_result = image.reduceRegion(
        reducer=reducer,
        geometry=geometry.centroid(maxError=1),
        scale=scale,
        maxPixels=1e12,
    )
    return ee.Algorithms.If(has_data, result, centroid_result)


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
        with (
            urllib.request.urlopen(url) as response,
            open(tmp_download, "wb") as handle,
        ):
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
    include_components: bool = True,
    scale: int = 300,
    stats: list[str] | None = None,
    percentiles: list[int] | None = None,
    output_format: OutputFormat = "dict",
    version: str | None = None,
    reference_year: int | None = None,
    exclude_wetlands: bool = False,
):
    """
    Extract EII values or statistics for a given geometry.

    Args:
        geometry: AOI geometry (ee.Geometry/Feature/FeatureCollection, GeoPandas,
            shapely, or bbox tuple).
        include_components: If True, include statistics for functional,
            structural, and compositional integrity. Defaults to True.
        scale: Resolution in meters for the reduction. Defaults to 300.
        stats: Statistics to compute. Supported: "mean", "median", "min", "max",
            "std". Defaults to ["mean"].
        percentiles: Optional list of percentiles (0-100) to include.
        reference_year: v2 reference year, defaults to the latest; see ``get_layers``.
        output_format: "dict" (default) or "geodataframe".
        version: EII dataset version (e.g. "v1", "v2"). If None, uses the active
            session version (see eii.client.set_version), defaulting to config.yaml.
        exclude_wetlands: v2 only. Mask pixels with a nonzero ``wetland_flag``
            before the reduction. By default flagged pixels are only reported.

    Returns:
        A structured dictionary or GeoDataFrame depending on output_format.
        Point geometries return single values; area geometries return summary
        statistics in nested format: {"eii": {"mean": ..., "min": ...}}.
        Version v2 adds "wetland_flag": the flag code at a point, or
        {"share": ...}, the flagged share of valid EII pixels, for an area
        (computed before any exclusion).

    Example:
        >>> import ee
        >>> ee.Initialize()
        >>> polygon = ee.Geometry.Rectangle([-60, -10, -55, -5])
        >>> stats = get_stats(polygon, stats=["mean", "min", "max"])
        >>> print(stats)
    """
    stats = _validate_stats_params(stats, percentiles)
    version = resolve_version(version)
    if exclude_wetlands and version != "v2":
        raise ValueError("exclude_wetlands is available only for version='v2'")
    geometry = normalize_client_input(geometry, target="geometry")
    layers = get_layers(
        layers="all" if include_components else "eii",
        version=version,
        reference_year=reference_year,
    )
    flag = None
    if version == "v2":
        flag = get_layers("wetland_flag", version, reference_year)["wetland_flag"]
        flagged = flag.gt(0).updateMask(layers["eii"].mask()).rename("wetland_flag")
        if exclude_wetlands:
            layers = {name: image.updateMask(flag.eq(0)) for name, image in layers.items()}
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
        components = {name: layers[name] for name in ("functional", "structural", "compositional")}
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

    if flag is not None:
        if is_point:
            code = _reduce_point_values(flag, geometry, scale).getInfo() or {}
            values["wetland_flag"] = code.get("wetland_flag")
        else:
            share = _reduce_area_stats(flagged, geometry, scale).getInfo() or {}
            values["wetland_flag"] = {"share": share.get("wetland_flag")}

    result = {"geometry_type": geometry_type, "values": values}

    if output_format == "geodataframe":
        geo_props = {"geometry_type": result["geometry_type"]}
        geo_props.update(result["values"])
        return _to_geodataframe(geometry, geo_props)

    return result


def get_raster(
    geometry: ee.Geometry,
    include_components: bool = True,
    scale: int = 300,
    version: str | None = None,
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
    reference_year: int | None = None,
):
    """
    Download EII rasters for an AOI as an xarray Dataset.

    Args:
        geometry: AOI geometry (ee.Geometry/Feature/FeatureCollection, GeoPandas,
            shapely, or bbox tuple).
        include_components: If True, include component bands alongside EII.
        scale: Resolution in meters for export. Defaults to 300.
        reference_year: v2 reference year, defaults to the latest; see ``get_layers``.
        version: EII dataset version (e.g. "v1", "v2"). If None, uses the active
            session version (see eii.client.set_version), defaulting to config.yaml.
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
        Version v2 appends "limiting_component" and "wetland_flag" when
        include_components=True.
        These categorical codes are excluded from get_stats/get_zonal_stats reductions.

    Returns:
        xarray.Dataset when output_format="memory"; otherwise a Path (single download)
        or list[Path] (chunked downloads) to the written files.
    """
    geometry = normalize_client_input(geometry, target="geometry")
    if output_format == "geotiff" and compression is None:
        compression = "ZSTD"
    layers = get_layers(
        layers="all" if include_components else "eii",
        version=version,
        reference_year=reference_year,
    )

    band_images = [layers["eii"].rename("eii")]
    band_names = ["eii"]
    if include_components:
        for name in ("functional", "structural", "compositional"):
            band = layers[name].rename(f"{name}_integrity")
            band_images.append(band)
            band_names.append(f"{name}_integrity")
        for name in CATEGORICAL_BANDS:
            if name in layers:
                band_images.append(layers[name])
                band_names.append(name)

    image = ee.Image.cat(band_images).clip(geometry)
    if any(name in band_names for name in CATEGORICAL_BANDS):
        # A multiband GeoTIFF needs a common type for scores and category codes.
        image = image.toDouble() if dtype == "float64" else image.toFloat()
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
                "eii_method": "min_fuzzy_logic",
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
            "eii_method": "min_fuzzy_logic",
            "scale_m": scale,
            "area_sq_km": area_sq_km,
            "chunked": True,
            "tile_size_deg": tile_size_deg,
        }
    )

    if out_path is None and tmp_dir is None:
        shutil.rmtree(output_dir, ignore_errors=True)

    return dataset
