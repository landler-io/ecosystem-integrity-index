"""
Internal GEE utility functions.
"""

from __future__ import annotations

import time

import ee


def load_tiled_collection(asset_folder: str) -> ee.ImageCollection:
    """
    Load a folder of tiled images as an ImageCollection.

    Args:
        asset_folder: Path to GEE asset folder containing image tiles.

    Returns:
        ImageCollection of all images in the folder.
    """
    assets = ee.data.listAssets({"parent": asset_folder})
    asset_ids = [a["id"] for a in assets.get("assets", [])]

    if not asset_ids:
        raise ValueError(f"No assets found in folder: {asset_folder}")

    return ee.ImageCollection(asset_ids)


def mosaic_collection(collection: ee.ImageCollection) -> ee.Image:
    """
    Mosaic an ImageCollection into a single Image.

    Args:
        collection: ImageCollection to mosaic.

    Returns:
        Single mosaicked Image.
    """
    return collection.mosaic()


def wait_for_tasks(task_ids: list[str], poll_interval: int = 30) -> None:
    """
    Wait for Earth Engine tasks to complete.

    Args:
        task_ids: List of task IDs to monitor.
        poll_interval: Seconds between status checks.
    """
    while True:
        statuses = [ee.data.getTaskStatus(tid)[0] for tid in task_ids]
        active = sum(1 for s in statuses if s["state"] in ["READY", "RUNNING"])
        completed = sum(1 for s in statuses if s["state"] == "COMPLETED")
        failed = sum(1 for s in statuses if s["state"] in ["FAILED", "CANCELLED"])

        print(
            f"Tasks: {len(task_ids)} total, {active} active, {completed} completed, {failed} failed"
        )

        if active == 0:
            if failed > 0:
                print(f"Warning: {failed} tasks failed")
            else:
                print("All tasks completed successfully")
            break

        time.sleep(poll_interval)
