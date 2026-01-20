"""
Statistical utilities for EII computation.
"""

from __future__ import annotations

import ee


def calculate_percentiles_from_collection(
    image_collection: ee.ImageCollection,
    band_name: str,
    percentiles: list[int] | None = None,
    samples_per_tile: int = 1000,
    scale: int = 300,
) -> dict:
    """
    Calculates percentiles for an image collection by sampling each image.
    This version pre-filters the collection to ensure only tiles with valid data are processed.
    """
    if percentiles is None:
        percentiles = [10, 20, 30, 40, 50, 60, 70, 80, 90]

    # Helper function to check if an image has any valid (unmasked) pixels.
    def check_has_data(image):
        # Use reduceRegion to get the pixel count.
        stats = image.select(band_name).reduceRegion(
            reducer=ee.Reducer.count(),
            geometry=image.geometry().bounds(),
            scale=scale * 100,  # Use a very coarse scale for a quick check
            bestEffort=True,  # Use bestEffort to avoid 'Too many pixels' error on large tiles.
        )
        # Set a property on the image. This will be null if there were no valid pixels.
        return image.set("has_data", stats.get(band_name))

    # Map the checker over the collection and filter for images that have the 'has_data' property.
    valid_collection = image_collection.map(check_has_data).filter(ee.Filter.notNull(["has_data"]))

    # --- Diagnostic Print ---
    # This will show how many tiles are being processed after filtering out empty ones.
    print("Initial tile count:", image_collection.size().getInfo())
    print("Tile count after filtering for valid data:", valid_collection.size().getInfo())
    # ---

    def sample_tile(image):
        """Samples a single image in the collection."""
        return image.select(band_name).sample(
            region=image.geometry().bounds(),  # Use bounds for robust geometry
            scale=scale,
            numPixels=samples_per_tile,
            seed=42,
            dropNulls=True,  # This is safe now because we know there's data
            tileScale=4,
        )

    # Map the sampling function over the *filtered* collection.
    samples = valid_collection.map(sample_tile).flatten()

    print(f"Number of samples: {samples.size().getInfo()}")

    # Final percentile calculation, now on a valid set of samples.
    percentile_reducer = ee.Reducer.percentile(percentiles)
    percentile_values = samples.reduceColumns(reducer=percentile_reducer, selectors=[band_name])

    return dict(percentile_values.getInfo())
