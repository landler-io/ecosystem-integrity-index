import time

import ee


def create_assets_folder(folder_path):
    """Recursively creates a folder structure in GEE."""
    try:
        ee.data.getAsset(folder_path)
    except ee.EEException:
        parts = folder_path.split("/")
        if len(parts) > 4:  # basic check to avoid going too high up 'projects/name/assets'
            parent = "/".join(parts[:-1])
            create_assets_folder(parent)

        try:
            ee.data.createAsset({"type": "Folder"}, folder_path)
            print(f"Created: {folder_path}")
        except ee.EEException as e:
            if "Already exists" not in str(e):
                print(f"Failed to create {folder_path}: {e}")


def get_status(task_ids):
    """Get the status of a specific export task."""
    task_info = ee.data.getTaskStatus(task_ids)
    return [s["state"] for s in task_info]


def wait_for_completion(export_tasks=None, id_list=None, wait=30):
    if id_list is None:
        id_list = [e["id"] for e in export_tasks]

    while True:
        status = get_status(id_list)
        failed = [s for s in status if s == "FAILED"]
        completed = [s for s in status if s == "COMPLETED"]
        active_states = ["RUNNING", "READY"]
        active = [s for s in status if s in active_states]

        print(
            f"Tasks: {len(status)} total, {len(active)} active, "
            f"{len(completed)} completed, {len(failed)} failed."
        )

        if not active:
            if failed:
                print(f"All active tasks have finished, but {len(failed)} tasks failed.")
            else:
                print("All tasks completed successfully.")
            return

        time.sleep(wait)
