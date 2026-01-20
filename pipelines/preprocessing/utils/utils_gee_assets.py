import contextlib

import ee


def make_asset_dir(folder_path: str) -> None:
    parts = folder_path.split("/")
    try:
        assets_idx = parts.index("assets")
    except ValueError:
        return

    for i in range(assets_idx + 1, len(parts)):
        partial_path = "/".join(parts[: i + 1])
        with contextlib.suppress(ee.EEException):
            ee.data.createAsset({"type": "Folder"}, partial_path)
