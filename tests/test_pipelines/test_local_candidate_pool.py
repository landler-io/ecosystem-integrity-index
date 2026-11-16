"""Offline tests for the local candidate pool: masks, tiles, selection."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import rasterio

from eii.training.sampling import (
    distance_to_mask_km,
    hmi_band,
    in_europe,
    inclusion_probabilities,
    natural_areas_arrays,
    one_per_cell,
)

PIPELINES = Path(__file__).resolve().parents[2] / "pipelines"
sys.path[:0] = [str(PIPELINES / "training"), str(PIPELINES / "preprocessing")]

import build_candidate_pool_local as pool  # noqa: E402
import local_sources as src  # noqa: E402


def test_natural_mask_buffers_high_pixels_by_a_ten_pixel_circle():
    series = np.zeros((2, 41, 41))
    series[1, 20, 20] = 0.5
    _, natural = natural_areas_arrays(series)
    # The 5x5 gap-fill spreads the high pixel to a 5x5 block (edge at row 18), then
    # the 10-pixel circle buffers that block.
    assert not natural[20, 20]
    assert not natural[20 - 12, 20]  # 10 px from the filled block edge
    assert natural[20 - 13, 20]
    assert not natural[18 - 7, 18 - 7]  # sqrt(98) <= 10
    assert natural[18 - 8, 18 - 8]  # sqrt(128) > 10


def test_natural_mask_gap_fill_and_nodata():
    series = np.full((1, 20, 20), np.nan)
    series[0, :, :10] = 0.05
    series[0, 5, 5] = np.nan  # interior gap: filled by the 5x5 focal max
    hmi_max, natural = natural_areas_arrays(series)
    assert hmi_max[5, 5] == pytest.approx(0.05)
    assert natural[5, 5]
    assert np.isnan(hmi_max[5, 15]) and not natural[5, 15]
    assert hmi_max[5, 11] == pytest.approx(0.05)  # filled two pixels into the gap


def test_hmi_max_takes_the_series_maximum():
    series = np.stack([np.full((7, 7), 0.02), np.full((7, 7), 0.15)])
    hmi_max, natural = natural_areas_arrays(series)
    assert np.allclose(hmi_max, 0.15)
    assert not natural.any()


def test_distance_to_mask_is_great_circle_and_capped():
    pixel = 0.01
    target = np.zeros((11, 201), dtype=bool)
    target[5, 100] = True
    lat = np.linspace(60.05, 59.95, 11)
    dist = distance_to_mask_km(target, lat, pixel, max_km=5.0)
    km_per_deg = np.pi / 180 * 6371.0088
    assert dist[5, 100] == 0
    assert dist[5, 108] == pytest.approx(8 * pixel * km_per_deg * np.cos(np.radians(60)), rel=1e-3)
    assert dist[1, 100] == pytest.approx(4 * pixel * km_per_deg, rel=1e-3)
    assert dist[5, 0] == 5.0
    assert (distance_to_mask_km(~target & False, lat, pixel, 5.0) == 5.0).all()


def test_hmi_bands():
    values = np.array([0.0, 0.049, 0.05, 0.1, 0.19, 0.2, 0.3, 0.31, np.nan])
    assert hmi_band(values).tolist() == [0, 0, 1, 2, 2, 3, 3, -1, -1]


def test_europe_box():
    assert in_europe(np.array([10.0, -30.0, 10.0]), np.array([50.0, 50.0, 30.0])).tolist() == [
        True,
        False,
        False,
    ]


def test_one_per_cell_keeps_the_smallest_key():
    cells = np.array([3, 1, 3, 1, 2])
    keys = np.array([0.5, 0.9, 0.1, 0.2, 0.7])
    assert one_per_cell(cells, keys).tolist() == [2, 3, 4]
    assert one_per_cell(np.array([], dtype=int), np.array([])).tolist() == []


def test_inclusion_probabilities_meet_every_target_in_expectation():
    rng = np.random.default_rng(0)
    biome = rng.choice([1, 2], size=20_000, p=[0.9, 0.1])
    europe = rng.random(20_000) < 0.05
    band = rng.choice([0, 1, 2, 3, -1], size=20_000, p=[0.6, 0.2, 0.1, 0.05, 0.05])
    p = inclusion_probabilities(biome, europe, band, {1: 900, 2: 500}, 300, 200)
    assert ((p >= 0) & (p <= 1)).all()
    assert p[biome == 2].sum() >= 500 - 1e-6
    assert p[biome == 1].sum() >= 900 - 1e-6
    assert p[europe].sum() >= 300 - 1e-6
    for b in range(4):
        assert p[band == b].sum() >= 200 - 1e-6
    assert p[biome == 2].min() == pytest.approx(500 / (biome == 2).sum())


def test_inclusion_probabilities_cap_at_one():
    p = inclusion_probabilities(np.array([1, 1]), np.zeros(2, bool), np.zeros(2, int), {1: 5})
    assert p.tolist() == [1.0, 1.0]


def test_tiles_partition_the_hmi_grid_without_gaps_or_overlap():
    rows_a, cols_a = pool.core_pixels((0, 45, 10, 55))
    rows_b, cols_b = pool.core_pixels((10, 45, 20, 55))
    rows_c, _ = pool.core_pixels((0, 35, 10, 45))
    assert cols_b[0] == cols_a[-1] + 1
    assert rows_a[-1] + 1 == rows_c[0]
    assert np.array_equal(rows_a, rows_b)
    s = src.GHM_PIXEL
    lon = src.GHM_X0 + (cols_a + 0.75) * s
    assert lon.min() >= 0 and lon.max() < 10


def test_tile_names_round_trip():
    tiles = pool.all_tiles()
    assert len(tiles) == 36 * 14
    assert all(pool.tile_name(*pool.tile_bounds(t)[:2]) == t for t in tiles)
    assert pool.tile_bounds("w010s20") == (-10, -20, 0, -10)


def test_source_tile_ids():
    assert src.worldcover_tiles((0, 45, 4, 48)) == ["N45E000", "N45E003"]


def test_ghm_columns_wrap_at_the_antimeridian():
    period = src.GHM_PERIOD
    assert src.wrapped_runs(10, 5, period) == [(0, 10, 5)]
    assert src.wrapped_runs(-3, 5, period) == [(0, period - 3, 3), (3, 0, 2)]
    assert src.wrapped_runs(period - 2, 5, period) == [(0, period - 2, 2), (2, 0, 3)]
    # One period west of the grid origin is 360 degrees, within half a pixel.
    assert abs(period * src.GHM_PIXEL - 360) < src.GHM_PIXEL / 2


def test_billing_project_comes_from_the_environment(monkeypatch):
    monkeypatch.delenv(src.BILLING_ENV, raising=False)
    with pytest.raises(RuntimeError, match=src.BILLING_ENV):
        src.bucket_url("datasets/x.tif")
    monkeypatch.setenv(src.BILLING_ENV, "my-project")
    assert src.bucket_url("datasets/x.tif").endswith("datasets/x.tif?userProject=my-project")


def test_wdpa_window_reads_a_local_raster_on_the_ghm_grid(tmp_path):
    from rasterio.transform import from_origin

    s = src.GHM_PIXEL
    data = np.array([[0, 1, 255], [1, 0, 1]], dtype="uint8")
    path = tmp_path / "wdpa.tif"
    with rasterio.open(path, "w", driver="GTiff", width=3, height=2, count=1, dtype="uint8",
                       crs="EPSG:4326", transform=from_origin(src.GHM_X0 + 10 * s,
                                                              src.GHM_Y0 - 20 * s, s, s),
                       nodata=255) as dst:  # fmt: skip
        dst.write(data, 1)
    window = src.read_wdpa_window(path, 19, 9, 4, 5)
    expected = np.zeros((4, 5), bool)
    expected[1, 2] = expected[2, 1] = expected[2, 3] = True
    assert np.array_equal(window, expected)
    with rasterio.open(tmp_path / "off.tif", "w", driver="GTiff", width=1, height=1, count=1,
                       dtype="uint8", crs="EPSG:4326",
                       transform=from_origin(src.GHM_X0 + s / 2, src.GHM_Y0, s, s)) as dst:  # fmt: skip
        dst.write(np.ones((1, 1), "uint8"), 1)
    with pytest.raises(ValueError, match="GHM grid"):
        src.read_wdpa_window(tmp_path / "off.tif", 0, 0, 1, 1)
