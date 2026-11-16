"""Tests for re-sampling the training response from the DMP V2.0.2 annual sums."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

rasterio = pytest.importorskip("rasterio")
from rasterio.transform import from_origin  # noqa: E402

REPO = Path(__file__).parents[2]
sys.path.insert(0, str(REPO / "pipelines" / "training"))
import training_table_dmp_v2 as table_dmp  # noqa: E402

PIXEL = 0.5
HEIGHT, WIDTH = 20, 1100  # spans several 512-pixel blocks
TRANSFORM = from_origin(-180, 10, PIXEL, PIXEL)


def write_annual(path: Path) -> np.ndarray:
    data = np.arange(HEIGHT * WIDTH, dtype="int32").reshape(HEIGHT, WIDTH)
    data[3, 7] = -1
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=WIDTH,
        height=HEIGHT,
        count=2,
        dtype="int32",
        crs="EPSG:4326",
        transform=TRANSFORM,
        nodata=-1,
        tiled=True,
        blockxsize=256,
        blockysize=256,
    ) as dst:
        dst.write(data, 1)
        dst.write(data * 2, 2)
    return data


def test_sampling_reads_the_pixel_containing_each_point(tmp_path):
    data = write_annual(tmp_path / "annual.tif")
    rng = np.random.default_rng(0)
    rows, cols = rng.integers(0, HEIGHT, 500), rng.integers(0, WIDTH, 500)
    # Anywhere inside the pixel, not only at its centre.
    lon = -180 + (cols + rng.uniform(0.01, 0.99, 500)) * PIXEL
    lat = 10 - (rows + rng.uniform(0.01, 0.99, 500)) * PIXEL
    out = table_dmp.sample_pixels(tmp_path / "annual.tif", lon, lat)
    assert np.array_equal(out[0], data[rows, cols])
    valid = data[rows, cols] >= 0
    assert np.array_equal(out[1][valid], data[rows, cols][valid] * 2)


def test_sampling_marks_outside_points_and_nodata(tmp_path):
    write_annual(tmp_path / "annual.tif")
    lon = np.array([-179.9, 400.0, -176.25])
    lat = np.array([50.0, 5.0, 8.25])  # above the raster, right of it, nodata pixel
    assert table_dmp.sample_pixels(tmp_path / "annual.tif", lon, lat)[0].tolist() == [-1] * 3


def make_table(n: int = 3) -> tuple[pd.DataFrame, pd.DataFrame]:
    years = range(2014, 2026)
    annual = pd.DataFrame({f"{b}_{y}": [1_000_000] * n for b in ("sum", "sd") for y in years})
    table = annual.copy() * 2
    table["longterm_avg_npp_sum"] = [2_000_000.0] * n
    table["selected"] = [True] * n
    table["reference"] = [True] * n
    table["response_source"] = ["old"] * n
    table["n_valid_years"] = [11] * n
    return table, annual


def test_rebuild_replaces_response_and_drops_rows_without_enough_years():
    table, annual = make_table()
    annual.loc[1, [f"sum_{y}" for y in range(2019, 2022)]] = -1  # 4 of 7 years left
    annual.loc[2, "sum_2018"] = -1  # outside the response window
    out, report = table_dmp.rebuild(table, annual)
    assert len(out) == 2
    assert (out["longterm_avg_npp_sum"] == 1_000_000).all()
    assert (out["n_valid_years"] == 7).all()
    assert out["sum_2018"].tolist() == [1_000_000, -1]
    assert (out["response_source"] == table_dmp.RESPONSE_SOURCE).all()
    assert report["dropped"] == report["dropped_selected"] == report["dropped_reference"] == 1
    assert report["response_median_ratio_new_over_old"] == pytest.approx(0.5)
