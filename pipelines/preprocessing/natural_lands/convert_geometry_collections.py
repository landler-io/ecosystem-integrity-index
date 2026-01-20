#!/usr/bin/env python3
"""
Convert GeometryCollections in GeoJSON files to Polygon/MultiPolygon only.

Drops LineStrings, Points, and other non-polygon geometry types.
Useful for cleaning up Earth Engine exports that contain mixed geometry types.

Usage:
    python convert_geometry_collections.py input.geojson
    python convert_geometry_collections.py input.geojson -o output.geojson
    python convert_geometry_collections.py *.geojson  # batch mode
"""

import argparse
import json
from pathlib import Path

from shapely.geometry import (
    GeometryCollection,
    MultiPolygon,
    Polygon,
    mapping,
    shape,
)


def extract_polygons(geometry):
    """
    Extract only Polygon geometries from any geometry type.

    Returns a list of Polygon objects, or empty list if none found.
    """
    geom = shape(geometry)

    if isinstance(geom, Polygon):
        return [geom] if geom.is_valid and not geom.is_empty else []

    if isinstance(geom, MultiPolygon):
        return [g for g in geom.geoms if g.is_valid and not g.is_empty]

    if isinstance(geom, GeometryCollection):
        polygons = []
        for g in geom.geoms:
            if isinstance(g, Polygon) and g.is_valid and not g.is_empty:
                polygons.append(g)
            elif isinstance(g, MultiPolygon):
                polygons.extend(p for p in g.geoms if p.is_valid and not p.is_empty)
        return polygons

    # Other geometry types (Point, LineString, etc.) - skip
    return []


def convert_feature(feature):
    """
    Convert a single GeoJSON feature to contain only polygon geometries.

    Returns the modified feature, or None if no polygons were found.
    """
    polygons = extract_polygons(feature["geometry"])

    if not polygons:
        return None

    # Use MultiPolygon if multiple polygons, otherwise single Polygon
    new_geom = polygons[0] if len(polygons) == 1 else MultiPolygon(polygons)

    feature = feature.copy()
    feature["geometry"] = mapping(new_geom)
    return feature


def convert_geojson(input_path, output_path=None):
    """
    Convert a GeoJSON file to contain only polygon geometries.

    Args:
        input_path: Path to input GeoJSON file
        output_path: Path to output file. If None, appends '_polygons_only' to input name.

    Returns:
        Tuple of (output_path, original_count, converted_count)
    """
    input_path = Path(input_path)

    if output_path is None:
        output_path = input_path.with_stem(input_path.stem + "_polygons_only")
    else:
        output_path = Path(output_path)

    with open(input_path) as f:
        data = json.load(f)

    original_count = len(data["features"])

    # Convert features, filtering out those with no polygons
    converted_features = []
    for feature in data["features"]:
        converted = convert_feature(feature)
        if converted is not None:
            converted_features.append(converted)

    data["features"] = converted_features

    with open(output_path, "w") as f:
        json.dump(data, f)

    return output_path, original_count, len(converted_features)


def main():
    parser = argparse.ArgumentParser(
        description="Convert GeometryCollections in GeoJSON to Polygon/MultiPolygon only.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        type=Path,
        help="Input GeoJSON file(s) to convert",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output file path (only valid with single input file)",
    )
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="Overwrite input files instead of creating new ones",
    )

    args = parser.parse_args()

    if args.output and len(args.input_files) > 1:
        parser.error("--output can only be used with a single input file")

    for input_file in args.input_files:
        if not input_file.exists():
            print(f"Warning: {input_file} does not exist, skipping")
            continue

        if args.inplace:
            output_file = input_file
        elif args.output:
            output_file = args.output
        else:
            output_file = None  # auto-generate

        try:
            out_path, orig_count, conv_count = convert_geojson(input_file, output_file)
            dropped = orig_count - conv_count
            print(
                f"{input_file.name}: {conv_count} polygon features retained, {dropped} dropped -> {out_path.name}"
            )
        except Exception as e:
            print(f"Error processing {input_file}: {e}")


if __name__ == "__main__":
    main()
