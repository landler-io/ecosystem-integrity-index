"""
Compositional Integrity Module
-----------------------------

This module handles functionality related to calculating compositional integrity
based on biodiversity metrics like the Biodiversity Intactness Index (BII).
"""

import ee

DEFAULT_BII_ASSET_PATH = "projects/ebx-data/assets/earthblox/IO/BIOINTACT"


def calculate_compositional_integrity(
    aoi: ee.Geometry | None = None,
    year: int = 2020,
    asset_path: str = DEFAULT_BII_ASSET_PATH,
) -> ee.Image:
    """
    Calculate compositional ecosystem integrity based on Biodiversity Intactness Index.

    Args:
        aoi: Area of interest (optional). If None, returns unclipped global image.
        year: Year to use for BII data.
        asset_path: Earth Engine ImageCollection path for BII data.

    Returns:
        Compositional integrity score (0-1 scale).
    """
    bii_collection = ee.ImageCollection(asset_path)
    start_date = f"{year}-01-01"
    end_date = f"{year}-12-31"

    year_filtered = bii_collection.filterDate(start_date, end_date)
    most_recent = bii_collection.sort("system:time_start", False).first()

    bii_image = ee.Image(
        ee.Algorithms.If(
            year_filtered.size().gt(0),
            year_filtered.mosaic(),
            most_recent,
        )
    )

    compositional_integrity = bii_image.rename("compositional_integrity")

    if aoi is not None:
        compositional_integrity = compositional_integrity.clip(aoi)

    return compositional_integrity
