"""Offline tests of the external validation suite (``pipelines/validation``).

Synthetic rasters and arrays only: no network, no Earth Engine, no external data.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.transform import Affine, from_bounds

PIPELINES_DIR = Path(__file__).resolve().parents[2] / "pipelines"
if str(PIPELINES_DIR) not in sys.path:
    sys.path.insert(0, str(PIPELINES_DIR))

from validation import analysis, metrics, run_validation, wcmc  # noqa: E402
from validation.grid import eckert4_profile  # noqa: E402
from validation.warping import WarpJob, cos_lat_weights, read_full, run_warp  # noqa: E402


@pytest.fixture(scope="session", autouse=True)
def initialize_ee():
    """Override the Earth Engine session fixture: these tests run offline."""


def _write(path: Path, data: np.ndarray, transform: Affine, crs: str = "EPSG:4326", nodata=-9999):
    with rasterio.open(
        path, "w", driver="GTiff", width=data.shape[-1], height=data.shape[-2],
        count=1 if data.ndim == 2 else data.shape[0], dtype=data.dtype, crs=crs,
        transform=transform, nodata=nodata,
    ) as dst:  # fmt: skip
        dst.write(data if data.ndim == 3 else data[None])
    return path


def _quiet(_msg: str) -> None:
    """Discard warp progress messages."""


# --- warping ----------------------------------------------------------------------------


def test_warp_average_aggregates_3x3_and_skips_nodata(tmp_path):
    """300 m -> 900 m area average ignores nodata pixels inside the target cell."""
    src = np.arange(36, dtype="float32").reshape(6, 6)
    src[0, 0] = -9999
    path = _write(tmp_path / "src.tif", src, from_bounds(0, 0, 6, 6, 6, 6))
    grid = {"crs": CRS.from_epsg(4326), "transform": from_bounds(0, 0, 6, 6, 2, 2), "width": 2,
            "height": 2}  # fmt: skip
    out = tmp_path / "out.tif"
    run_warp(WarpJob(path, (1,), (out,), src_nodata=-9999), grid, workers=1, log=_quiet)
    got = read_full(out)
    block = src[:3, :3].ravel()[1:]  # nodata excluded
    assert got[0, 0] == pytest.approx(block.mean())
    assert got[1, 1] == pytest.approx(src[3:, 3:].mean())


def test_warp_multiband_scale_and_valid_range(tmp_path):
    data = np.stack([np.full((4, 4), 250, "int16"), np.full((4, 4), 9000, "int16")])
    path = _write(tmp_path / "src.tif", data, from_bounds(0, 0, 4, 4, 4, 4), nodata=-1)
    grid = {"crs": CRS.from_epsg(4326), "transform": from_bounds(0, 0, 4, 4, 2, 2), "width": 2,
            "height": 2}  # fmt: skip
    outs = (tmp_path / "a.tif", tmp_path / "b.tif")
    run_warp(WarpJob(path, (1, 2), outs, scale=0.01, valid=(0, 64)), grid, workers=1,
             log=_quiet)  # fmt: skip
    assert np.allclose(read_full(outs[0]), 2.5)
    assert np.isnan(read_full(outs[1])).all()  # 90 > 64 -> NaN


def test_warp_to_eckert_edge_blocks_do_not_fail(tmp_path):
    """Blocks straddling the Eckert IV domain edge warp (no GDAL error, no whole-source read)."""
    src = np.ones((180, 360), "float32")
    path = _write(tmp_path / "globe.tif", src, from_bounds(-180, -90, 180, 90, 360, 180))
    p = eckert4_profile(200_000.0)  # coarse global Eckert IV grid
    grid = {k: p[k] for k in ("crs", "transform", "width", "height")}
    out = tmp_path / "e.tif"
    run_warp(WarpJob(path, (1,), (out,), resampling=Resampling.nearest), grid, workers=1,
             block=16, log=_quiet)  # fmt: skip
    got = read_full(out)
    assert np.nanmax(got) == 1.0
    assert np.isnan(got[0, 0])  # map corner is off the globe
    assert np.isfinite(got[got.shape[0] // 2, got.shape[1] // 2])


def test_cos_lat_weights():
    grid = {"transform": from_bounds(-180, -90, 180, 90, 4, 2), "height": 2}
    assert np.allclose(cos_lat_weights(grid), np.cos(np.radians([45, -45])))


# --- analysis helpers -------------------------------------------------------------------


def test_influence_classes_per_layer():
    anth = np.array([61, 62, 63, 43, 53, 54, 11, 34, 0, np.nan], dtype="float32")
    assert analysis.influence_classes_anthromes(anth).tolist() == [0, 0, 0, 1, 1, 1, 2, 2, -1, -1]
    hf = np.array([0, 0.01, 3, 3.01, np.nan], dtype="float32")
    assert analysis.influence_classes_hf(hf).tolist() == [0, 1, 1, 2, -1]
    vl = np.array([0, 100, 0, 100, np.nan], dtype="float32")  # 0 = inside the area
    lo = np.array([0, 0, 100, 100, np.nan], dtype="float32")
    assert analysis.influence_classes_lia(vl, lo).tolist() == [0, 1, 0, 2, -1]  # very low wins


def test_class_histograms_and_auc():
    classes = np.array([0, 0, 2, 2, -1], "int8")
    values = np.array([0.95, 0.85, 0.15, 0.25, 0.5], "float32")
    hist = analysis.class_histograms(classes, values)
    assert hist.shape == (3, 100)
    assert hist[0, 95] == 1 and hist[0, 85] == 1 and hist[2].sum() == 2 and hist.sum() == 4
    assert analysis.auc_from_hist(hist[0], hist[2]) == pytest.approx(1.0)
    assert analysis.auc_from_hist(hist[2], hist[0]) == pytest.approx(0.0)
    assert analysis.auc_from_hist(hist[0], hist[0]) == pytest.approx(0.5)


def test_hill_mean_abs_r_reproduces_paper_values():
    """Hill's §A2.1 means follow from the Fig A4 pairwise r values."""
    means = analysis.hill_means(analysis.PAPER["A4_pairs"])
    for key, value in analysis.PAPER["A4_mean"].items():
        assert means[key] == pytest.approx(value, abs=0.005)  # paper rounds HFP 0.677 to 0.68


def test_unit_means_nan_aware():
    ids = np.array([1, 1, 2, 2, np.nan], "float32")
    vals = np.array([0.2, 0.4, np.nan, 1.0, 5.0], "float32")
    mean, count = analysis.unit_means(ids, vals, 3)
    assert mean[1] == pytest.approx(0.3) and mean[2] == pytest.approx(1.0)
    assert count.tolist() == [0, 2, 1]


def test_weighted_moments_match_numpy():
    rng = np.random.default_rng(0)
    x, y, w = rng.random(1000), rng.random(1000), rng.random(1000)
    y = 0.5 * x + 0.5 * y
    acc = np.zeros((8, 2))
    ids = (rng.random(1000) < 0.5).astype("int64")
    analysis.accumulate(acc, ids, x, y, w)
    res = analysis.finish(acc)
    for k in (0, 1):
        m = ids == k
        xm, ym, wm = x[m], y[m], w[m]
        assert res.loc[k, "bias"] == pytest.approx(np.average(ym - xm, weights=wm))
        assert res.loc[k, "rmse"] == pytest.approx(np.sqrt(np.average((ym - xm) ** 2, weights=wm)))
        c = np.cov(xm, ym, aweights=wm)
        assert res.loc[k, "pearson_r_w"] == pytest.approx(c[0, 1] / np.sqrt(c[0, 0] * c[1, 1]))


def test_draw_land_points_reproducible():
    p = eckert4_profile(500_000.0)
    grid = {k: p[k] for k in ("crs", "transform", "width", "height")}
    land = np.zeros((grid["height"], grid["width"]), bool)
    land[: grid["height"] // 2] = True  # northern half
    a = analysis.draw_land_points(land, grid, 200, seed=42)
    b = analysis.draw_land_points(land, grid, 200, seed=42)
    assert a.equals(b) and len(a) == 200
    assert land[a["row"], a["col"]].all()
    assert (a["lat"] > -5).all()


def test_stage_warp_rejects_a_filter_without_jobs(monkeypatch: pytest.MonkeyPatch) -> None:
    """A year skipped as incomplete must not pass as a finished warp."""
    monkeypatch.setattr(run_validation, "jobs", lambda _only: [])
    with pytest.raises(RuntimeError, match="eii_2025"):
        run_validation.stage_warp(["eii_2025"], workers=1)


def test_random_points_are_equal_area_and_reproducible():
    a = analysis.random_points(200_000, seed=42)
    assert np.array_equal(a, analysis.random_points(200_000, seed=42))
    sin_lat = np.sin(np.radians(a[:, 1]))
    assert abs(sin_lat.mean()) < 0.01 and abs(np.mean(np.abs(sin_lat)) - 0.5) < 0.01


def test_mosaic_vrt_places_tiles(tmp_path):
    left = _write(
        tmp_path / "a.tif", np.full((2, 2), 1.0, "float32"), from_bounds(0, 0, 2, 2, 2, 2)
    )
    right = _write(
        tmp_path / "b.tif", np.full((2, 2), 2.0, "float32"), from_bounds(2, 0, 4, 2, 2, 2)
    )
    vrt = run_validation.mosaic_vrt([left, right], tmp_path / "m.vrt")
    with rasterio.open(vrt) as src:
        assert (src.width, src.height) == (4, 2)
        assert src.read(1).tolist() == [[1, 1, 2, 2], [1, 1, 2, 2]]


# --- metrics ----------------------------------------------------------------------------


def test_pearson_spearman_perfect_nan_and_constant():
    x = np.linspace(0, 1, 100)
    y = 2 * x + 1
    y[0] = np.nan
    out = metrics.pearson_spearman(x, y)
    assert out["pearson_r"] == pytest.approx(1.0) and out["spearman_rho"] == pytest.approx(1.0)
    assert out["n"] == 99
    constant = metrics.pearson_spearman(np.ones(10), np.arange(10.0))
    assert np.isnan(constant["pearson_r"]) and constant["n"] == 10


def test_locke_conditions_mutually_exclusive():
    hf = np.array([1.0, 1.0, 9.0, 9.0])
    intensive = np.array([0.001, 0.6, 0.001, 0.3])
    out = metrics.locke_condition(hf, intensive)
    assert list(out) == ["large_wild_areas", "cities_farms", "shared_lands", "shared_lands"]


def test_parse_checksums_strips_binary_marker():
    text = "ABC123  *EII_1km/EII_2000_1km.tif\n\ndef456 README.txt\n"
    assert wcmc.parse_checksums(text) == {
        "EII_1km/EII_2000_1km.tif": "abc123",
        "README.txt": "def456",
    }
