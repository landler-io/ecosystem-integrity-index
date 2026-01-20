"""
Client utilities for visualization.
"""

from __future__ import annotations

from collections.abc import Iterable
from contextlib import suppress
from math import floor, log10
from typing import Any, Literal, cast

import ee


def _normalize_geometry(geometry: ee.Geometry, max_error: float = 1) -> ee.Geometry:
    """Ensure geometry uses EPSG:4326 when a projection is defined."""
    crs = geometry.projection().crs()
    normalized = ee.Algorithms.If(
        ee.String(crs).compareTo("EPSG:4326").neq(0),
        geometry.transform("EPSG:4326", max_error),
        geometry,
    )
    return ee.Geometry(normalized)


def _is_bbox(value: object) -> bool:
    if not isinstance(value, list | tuple) or len(value) != 4:
        return False
    return all(isinstance(v, int | float) for v in value)


def _has_geo_interface(value: object) -> bool:
    return hasattr(value, "__geo_interface__")


def _to_geo_interface(value: object) -> dict[Any, Any]:
    if not _has_geo_interface(value):
        raise ValueError("Object does not provide __geo_interface__")
    return cast(dict[Any, Any], value.__geo_interface__)  # type: ignore[attr-defined]


def _normalize_feature_collection(
    features: ee.FeatureCollection, max_error: float = 1
) -> ee.FeatureCollection:
    """Normalize feature geometries to EPSG:4326."""

    def _normalize_feature(feature: ee.Feature) -> ee.Feature:
        return feature.setGeometry(_normalize_geometry(feature.geometry(), max_error))

    return features.map(_normalize_feature)


def normalize_client_input(
    value: object,
    target: Literal["geometry", "features"] = "geometry",
    max_error: float = 1,
) -> ee.Geometry | ee.FeatureCollection:
    """Normalize client geometry inputs to ee.Geometry or ee.FeatureCollection."""
    if target not in ("geometry", "features"):
        raise ValueError("target must be 'geometry' or 'features'")

    if target == "geometry":
        if isinstance(value, ee.Geometry):
            return _normalize_geometry(value, max_error)
        if isinstance(value, ee.Feature):
            return _normalize_geometry(value.geometry(), max_error)
        if isinstance(value, ee.FeatureCollection):
            return _normalize_geometry(value.geometry(), max_error)
        if _is_bbox(value):
            return _normalize_geometry(ee.Geometry.Rectangle(value), max_error)
        if _has_geo_interface(value):
            geo = _to_geo_interface(value)
            if geo.get("type") == "FeatureCollection":
                return _normalize_geometry(ee.FeatureCollection(geo).geometry(), max_error)
            if geo.get("type") == "Feature":
                return _normalize_geometry(ee.Feature(geo).geometry(), max_error)
            return _normalize_geometry(ee.Geometry(geo), max_error)
        raise ValueError(
            "Unsupported geometry input. Provide ee.Geometry, ee.Feature, "
            "ee.FeatureCollection, a shapely geometry, a GeoPandas object, or a bbox."
        )

    if isinstance(value, ee.FeatureCollection):
        return _normalize_feature_collection(value, max_error)
    if isinstance(value, ee.Feature):
        return _normalize_feature_collection(ee.FeatureCollection([value]), max_error)
    if isinstance(value, ee.Geometry):
        return _normalize_feature_collection(ee.FeatureCollection([ee.Feature(value)]), max_error)
    if _is_bbox(value):
        geometry = ee.Geometry.Rectangle(value)
        return _normalize_feature_collection(
            ee.FeatureCollection([ee.Feature(geometry)]), max_error
        )
    if _has_geo_interface(value):
        geo = _to_geo_interface(value)
        if geo.get("type") == "FeatureCollection":
            return _normalize_feature_collection(ee.FeatureCollection(geo), max_error)
        if geo.get("type") == "Feature":
            return _normalize_feature_collection(ee.FeatureCollection([ee.Feature(geo)]), max_error)
        return _normalize_feature_collection(
            ee.FeatureCollection([ee.Feature(ee.Geometry(geo))]), max_error
        )
    if isinstance(value, Iterable) and not isinstance(value, str | bytes):
        features = []
        for item in value:
            if isinstance(item, ee.Feature):
                features.append(item)
            elif isinstance(item, ee.Geometry):
                features.append(ee.Feature(item))
            elif _has_geo_interface(item):
                features.append(ee.Feature(_to_geo_interface(item)))
            else:
                raise ValueError("Unsupported item in feature iterable.")
        return _normalize_feature_collection(ee.FeatureCollection(features), max_error)
    raise ValueError(
        "Unsupported features input. Provide ee.FeatureCollection, ee.Feature, "
        "ee.Geometry, a GeoPandas object, a shapely geometry, or a bbox."
    )


def quicklook(
    geometry: ee.Geometry,
    min_bbox_km: float = 1.0,
    layer_color: str = "yellow",
    overlay_fields: list[str] | None = None,
):
    """
    Create a static plot of a geometry on a satellite basemap.

    Args:
        geometry: Earth Engine geometry or FeatureCollection, a GeoDataFrame,
            or a shapely geometry. Multiple features create a faceted plot.
        min_bbox_km: Minimum bounding box size (km). Ensures small geometries
            are shown with at least a 1km extent.
        layer_color: Color to use for the geometry outline.
        overlay_fields: Optional list of fields to display in an overlay box.
            If None, attempts to show EII mean and component means when present.

    Returns:
        A matplotlib Figure and Axes with the rendered map(s).
    """
    import contextily as ctx
    import geopandas as gpd
    import matplotlib.pyplot as plt
    from pyproj import Transformer
    from shapely.geometry import box, shape
    from shapely.ops import transform as shapely_transform

    properties: list[dict] = []

    if hasattr(geometry, "geometry") and hasattr(geometry, "to_crs"):
        gdf_input = geometry
        with suppress(Exception):
            gdf_input = gdf_input.to_crs("EPSG:4326")
        shapes = list(gdf_input.geometry)
        if hasattr(gdf_input, "iterrows"):
            for _, row in gdf_input.iterrows():
                properties.append(row.drop(labels=["geometry"]).to_dict())
    elif hasattr(geometry, "geom_type"):
        shapes = [geometry]
        properties = [{}]
    elif isinstance(geometry, ee.FeatureCollection):
        fc_info = geometry.getInfo()
        feats = fc_info.get("features", [])
        shapes = [shape(feat["geometry"]) for feat in feats]
        properties = [feat.get("properties", {}) for feat in feats]
    else:
        geometry = _normalize_geometry(geometry)
        shapes = [shape(geometry.getInfo())]
        properties = [{}]

    to_3857 = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True).transform
    to_4326 = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True).transform

    min_size_m = float(min_bbox_km) * 1000.0
    n = len(shapes)
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 6 * rows))
    axes_flat = list(axes.ravel()) if hasattr(axes, "ravel") else [axes]

    def _format_overlay_text(props: dict) -> str | None:
        if not props:
            return None

        metric_candidates = [
            "eii",
            "functional_integrity",
            "structural_integrity",
            "compositional_integrity",
        ]

        lines = []

        def _get_value(field: str | tuple):
            if field in props:
                return props.get(field)
            if isinstance(field, str) and ":" in field:
                metric, stat = field.split(":", 1)
                return props.get((metric, stat))
            return None

        if overlay_fields:
            for field in overlay_fields:
                value = _get_value(field)
                if value is not None:
                    if isinstance(field, tuple):
                        label = f"{field[0]} {field[1]}".replace("_", " ")
                    else:
                        label = field.replace("_", " ")
                    lines.append(
                        f"{label}: {value:.3f}"
                        if isinstance(value, int | float)
                        else f"{label}: {value}"
                    )
            return "\n".join(lines) if lines else None

        # Auto-detect mean values
        for metric in metric_candidates:
            for key in (f"{metric}_mean", metric, (metric, "mean")):
                if key in props:
                    value = props.get(key)
                    lines.append(
                        f"{metric.replace('_', ' ')}: {value:.3f}"
                        if isinstance(value, int | float)
                        else f"{metric.replace('_', ' ')}: {value}"
                    )
                    break

        return "\n".join(lines) if lines else None

    for idx, geom_shape in enumerate(shapes):
        ax = axes_flat[idx]
        geom_3857 = shapely_transform(to_3857, geom_shape)
        minx, miny, maxx, maxy = geom_3857.bounds
        width = maxx - minx
        height = maxy - miny
        expand_x = max(0.0, (min_size_m - width) / 2.0)
        expand_y = max(0.0, (min_size_m - height) / 2.0)

        bbox_3857 = box(minx - expand_x, miny - expand_y, maxx + expand_x, maxy + expand_y)
        bbox_4326 = shapely_transform(to_4326, bbox_3857)

        gdf_geom = gpd.GeoDataFrame(geometry=[geom_shape], crs="EPSG:4326").to_crs("EPSG:3857")
        gdf_bbox = gpd.GeoDataFrame(geometry=[bbox_4326], crs="EPSG:4326").to_crs("EPSG:3857")

        gdf_bbox.plot(ax=ax, facecolor="none", edgecolor="none")
        gdf_geom.plot(ax=ax, facecolor="none", edgecolor=layer_color, linewidth=2)
        ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, attribution=False)
        ax.set_xlim(bbox_3857.bounds[0], bbox_3857.bounds[2])
        ax.set_ylim(bbox_3857.bounds[1], bbox_3857.bounds[3])
        ax.text(
            0.995,
            0.005,
            "Esri World Imagery",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=6,
            color="white",
        )
        # Overlay stats box
        props = properties[idx] if idx < len(properties) else {}
        overlay_text = _format_overlay_text(props)
        if overlay_text:
            ax.text(
                0.01,
                0.99,
                overlay_text,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10,
                color="black",
                bbox={
                    "boxstyle": "round,pad=0.3",
                    "facecolor": layer_color,
                    "alpha": 0.7,
                    "edgecolor": "none",
                },
            )

        # Scale bar based on actual bbox width
        actual_width = bbox_3857.bounds[2] - bbox_3857.bounds[0]
        target_len = actual_width / 4
        magnitude = 10 ** floor(log10(target_len)) if target_len > 0 else 1000
        nice_values = [1, 2, 5, 10]
        bar_len = min(nice_values, key=lambda x: abs(x * magnitude - target_len)) * magnitude
        bar_len = int(bar_len)

        pad = actual_width * 0.02
        bar_x0 = bbox_3857.bounds[0] + pad
        bar_y0 = bbox_3857.bounds[1] + pad
        ax.plot([bar_x0, bar_x0 + bar_len], [bar_y0, bar_y0], color=layer_color, linewidth=2)

        bar_label = f"{int(bar_len / 1000)} km" if bar_len >= 1000 else f"{int(bar_len)} m"
        ax.text(
            bar_x0 + bar_len / 2.0,
            bar_y0 + pad * 0.6,
            bar_label,
            ha="center",
            va="bottom",
            fontsize=7,
            color=layer_color,
        )
        ax.set_axis_off()

    for ax in axes_flat[n:]:
        ax.set_axis_off()

    return fig, axes
