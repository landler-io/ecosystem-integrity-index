"""COG conversion of the v2 publication (pipelines/production/publish_v2.py)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipelines" / "production"))

import publish_v2  # noqa: E402


def test_to_cog_keeps_data_bands_and_tags(tmp_path: Path) -> None:
    data = np.random.default_rng(0).random((2, 64, 64), dtype=np.float32)
    data[:, :8, :8] = publish_v2.NODATA
    tile = tmp_path / "e010n40.tif"
    profile = {
        "driver": "GTiff",
        "width": 64,
        "height": 64,
        "count": 2,
        "dtype": "float32",
        "crs": "EPSG:4326",
        "transform": from_origin(10, 50, 1 / 336, 1 / 336),
        "nodata": publish_v2.NODATA,
    }
    with rasterio.open(tile, "w", **profile) as dst:
        dst.write(data)
        dst.descriptions = ("eii", "fi_change")
        dst.update_tags(reference_year="2024", fi_change_series="TV")

    cog = tmp_path / "cog" / "2024" / "e010n40.tif"
    cog.parent.mkdir(parents=True)
    record = publish_v2.to_cog(tile, cog, "2024/e010n40", {"reference_year": 2024})

    assert record is not None and record["grid"]["bands"] == ["eii", "fi_change"]
    with rasterio.open(cog) as src:
        assert np.array_equal(src.read(), data)
        assert src.descriptions == ("eii", "fi_change")
        assert src.tags()["reference_year"] == "2024"
        assert src.tags()["fi_change_series"] == "TV"
