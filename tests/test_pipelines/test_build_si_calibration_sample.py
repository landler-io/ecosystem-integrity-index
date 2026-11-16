"""Offline tests for the structural calibration sample: draw, land filter, class attachment."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipelines" / "training"))
import build_si_calibration_sample as si  # noqa: E402


def test_draw_is_seeded_equal_area_and_inside_the_box():
    lon, lat = si.draw_points()
    assert lon[0] == pytest.approx(98.62417748) and lat[0] == pytest.approx(-30.53964358)
    assert lon.min() >= -180 and lon.max() <= 180 and lat.min() >= -60 and lat.max() <= 75
    s = np.sin(np.radians(lat))
    edges = np.linspace(np.sin(np.radians(-60)), np.sin(np.radians(75)), 5)
    shares = np.histogram(s, edges)[0] / len(s)
    np.testing.assert_allclose(shares, 0.25, atol=0.005)


def test_land_filter_keeps_land_cell_centres_per_chunk(monkeypatch):
    responses = [
        {
            "features": [
                {"geometry": {"coordinates": [1.0, 2.0]}, "properties": {"land": 1}},
                {"geometry": {"coordinates": [3.0, 4.0]}, "properties": {"land": 0}},
            ]
        },
        {"features": [{"geometry": {"coordinates": [5.0, 6.0]}, "properties": {"land": 1}}]},
    ]
    ee = MagicMock()
    land = ee.ImageCollection.return_value.first.return_value.select.return_value.mask.return_value
    sample = land.And.return_value.rename.return_value.sampleRegions
    sample.return_value.getInfo.side_effect = responses
    monkeypatch.setitem(sys.modules, "ee", ee)

    out = si.land_points(np.array([0.1, 0.2, 0.3]), np.array([0.4, 0.5, 0.6]), chunk_size=2)

    assert out.to_numpy().tolist() == [[1.0, 2.0], [5.0, 6.0]]
    assert sample.call_count == 2
    assert all(c.kwargs["scale"] == 1000 and c.kwargs["geometries"] for c in sample.call_args_list)


def _raster(path: Path, data: np.ndarray, nodata: float | None) -> Path:
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype=data.dtype,
        crs="EPSG:4326",
        transform=from_origin(0, 2, 1, 1),
        nodata=nodata,
    ) as dst:
        dst.write(data, 1)
    return path


def test_classes_map_eco_id_to_biome_and_fill_nodata(tmp_path):
    eco = _raster(tmp_path / "eco.tif", np.array([[1, -9], [0, 400]], "int32"), -9)
    anthromes = _raster(tmp_path / "anthromes.tif", np.array([[61, 0], [32, 11]], "uint8"), None)
    points = pd.DataFrame({"lon": [0.5, 1.5, 0.5, 1.5, 5.0], "lat": [1.5, 1.5, 0.5, 0.5, 1.0]})

    out = si.attach_classes(points, eco, anthromes)

    lookup = pd.read_json(si.ECO_TO_BIOME, typ="series")
    assert out["biome"].tolist() == [lookup[1], 0, lookup[0], lookup[400], 0]
    assert out["anthrome"].tolist() == [61, 0, 32, 11, 0]


def test_complete_records_keeps_returned_rows_in_input_order(monkeypatch):
    image = MagicMock()
    image.sampleRegions.return_value.getInfo.side_effect = [
        {"features": [{"properties": {"row": 1}}, {"properties": {"row": 0}}]},
        {"features": []},
    ]
    monkeypatch.setattr(si, "record_image", lambda _ee: image)
    monkeypatch.setitem(sys.modules, "ee", MagicMock())
    points = pd.DataFrame({"lon": [1.0, 2.0, 3.0], "lat": [4.0, 5.0, 6.0]})

    out = si.complete_records(points, chunk_size=2)

    assert out["lon"].tolist() == [1.0, 2.0]
    assert image.sampleRegions.call_count == 2
