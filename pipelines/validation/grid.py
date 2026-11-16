"""The Eckert IV 1 km sampling grid of the external validation.

All layers are compared on Eckert IV (ESRI:54012, equal-area) at 1 km, the grid used by
Hill et al. (2026), section A2. Our layers (~300 m, EPSG:4326) are aggregated to it by the
area average of valid pixels. The pixel comparison with the WCMC EII runs on the WCMC
product grid instead (EPSG:4326, read from their GeoTIFF).
"""

from __future__ import annotations

from typing import Any

import numpy as np
from rasterio.crs import CRS
from rasterio.transform import Affine
from rasterio.warp import transform as warp_transform

NODATA = -9999.0


def eckert4_grid(res_m: float = 1000.0) -> tuple[CRS, Affine, int, int]:
    """Global Eckert IV (ESRI:54012) grid at ``res_m`` metres, snapped to whole pixels.

    Args:
        res_m: Cell size in metres.

    Returns:
        (crs, transform, width, height) of a grid that covers the full globe.
    """
    crs = CRS.from_string("ESRI:54012")
    wgs84 = CRS.from_string("EPSG:4326")
    xs, _ = warp_transform(wgs84, crs, [-180.0, 180.0], [0.0, 0.0])
    _, norths = warp_transform(wgs84, crs, [0.0], [90.0])
    _, souths = warp_transform(wgs84, crs, [0.0], [-90.0])
    west, east = min(xs), max(xs)
    south, north = min(souths), max(norths)
    width = int(np.ceil((east - west) / res_m))
    height = int(np.ceil((north - south) / res_m))
    return crs, Affine(res_m, 0.0, west, 0.0, -res_m, north), width, height


def eckert4_profile(res_m: float = 1000.0) -> dict[str, Any]:
    """rasterio profile of the Eckert IV sampling grid."""
    crs, transform, width, height = eckert4_grid(res_m)
    return {
        "driver": "GTiff",
        "dtype": "float32",
        "nodata": NODATA,
        "crs": crs,
        "width": width,
        "height": height,
        "transform": transform,
    }
