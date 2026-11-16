"""Tests for the CLMS DMP to annual NPP aggregation (pipelines/preprocessing/npp)."""

import importlib.util
import sys
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import pytest

rasterio = pytest.importorskip("rasterio")
pytest.importorskip("numba")
from rasterio.transform import from_origin  # noqa: E402

NPP_DIR = Path(__file__).parents[2] / "pipelines" / "preprocessing" / "npp"
sys.path.insert(0, str(NPP_DIR.parent))
_spec = importlib.util.spec_from_file_location("annual_npp", NPP_DIR / "01_calculate_annual_npp.py")
npp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(npp)

SIZE = 8
TRANSFORM = from_origin(-180, 80, 1 / 336, 1 / 336)


def dekad_ends(year: int) -> list[date]:
    ends = []
    for month in range(1, 13):
        first_of_next = date(year + month // 12, month % 12 + 1, 1)
        ends += [date(year, month, 10), date(year, month, 20), first_of_next - timedelta(days=1)]
    return ends


def write(path: Path, data: np.ndarray, dtype: str, nodata: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    profile = {
        "driver": "GTiff",
        "width": SIZE,
        "height": SIZE,
        "count": 1,
        "dtype": dtype,
        "crs": "EPSG:4326",
        "transform": TRANSFORM,
        "nodata": nodata,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(data.astype(dtype), 1)


def make_year(
    root: Path, year: int, dmp, qflag=None, version: str = "V2.0.2", sensor: str = "OLCI"
) -> list[tuple[Path, Path]]:
    """Write one product per dekad. ``dmp``/``qflag`` map dekad index to an array."""
    pairs = []
    for i, end in enumerate(dekad_ends(year)):
        stamp = f"{end:%Y%m%d}0000"
        folder = root / f"{year}" / f"c_gls_DMP300-RT6_{stamp}_GLOBE_{sensor}_{version}_cog"
        name = f"RT6_{stamp}_GLOBE_{sensor}_{version}.tiff"
        dmp_path = folder / f"c_gls_DMP300-DMP-{name}"
        qflag_path = folder / f"c_gls_DMP300-QFLAG-{name}"
        write(dmp_path, dmp(i), "int16", -1)
        write(qflag_path, np.zeros((SIZE, SIZE)) if qflag is None else qflag(i), "uint8", 255)
        pairs.append((dmp_path, qflag_path))
    return pairs


def constant(_dekad: int) -> np.ndarray:
    return np.full((SIZE, SIZE), 1000)


@pytest.fixture
def output_dir(tmp_path, monkeypatch):
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr(npp, "NPP_OUTPUT_DIR", out)
    return out


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("c_gls_DMP300-DMP-RT6_201501100000_GLOBE_PROBAV_V2.0.2.tiff", 5.5),
        ("c_gls_DMP300-DMP-RT6_201501200000_GLOBE_PROBAV_V2.0.2.tiff", 15.5),
        ("c_gls_DMP300-DMP-RT6_201501310000_GLOBE_PROBAV_V2.0.2.tiff", 26.0),
        ("c_gls_DMP300-DMP-RT6_202102280000_GLOBE_OLCI_V2.0.2.tiff", 55.5),
        ("c_gls_DMP300-DMP-RT6_202012310000_GLOBE_OLCI_V2.0.2.tiff", 361.0),
    ],
)
def test_dekads_are_anchored_at_their_midpoint(name, expected):
    assert npp.get_day_of_year_from_filename(Path(name)) == expected


def test_rejects_dates_that_do_not_end_a_dekad():
    with pytest.raises(ValueError, match="dekad end"):
        npp.get_day_of_year_from_filename(Path("c_gls_DMP300-DMP-RT6_202102270000_X.tiff"))


def test_product_version_is_parsed():
    name = "c_gls_DMP300-DMP-RT6_201501100000_GLOBE_PROBAV_V2.0.2.tiff"
    assert npp.get_product_version(Path(name)) == ("PROBAV", "V2.0.2")


@pytest.mark.usefixtures("output_dir")
def test_constant_rate_gives_stored_physical_times_10000(tmp_path):
    pairs = make_year(tmp_path, 2021, constant)
    out = npp.process_year(2021, pairs)
    with rasterio.open(out) as src:
        total, std = src.read()
        assert src.scales == (1e-4, 1e-4)
        assert src.units == ("gC m-2 yr-1", "gC m-2 day-1")
        assert src.tags()["source"] == "CLMS DMP 300m 10-daily V2.0.2 RT6 (OLCI)"
    # 10 kg DM/ha/day (DN 1000) for 365 days = 182.5 gC m-2 yr-1 at 0.5 gC/gDM.
    assert np.all(total == 1_825_000)
    assert np.all(std == 0)


@pytest.mark.usefixtures("output_dir")
def test_water_missing_and_invalid_quality_are_masked(tmp_path):
    def dmp(i):
        values = np.full((SIZE, SIZE), 1000)
        values[0, 0] = -2  # water in every dekad
        values[0, 1] = -1  # missing in every dekad
        if i == 18:
            values[1, :] = 30000  # invalid outlier, flagged below
        return values

    def qflag(i):
        flags = np.zeros((SIZE, SIZE))
        flags[0, 2] = 255  # missing quality in every dekad
        if i == 18:
            flags[1, :4] = 1  # FAPAR invalid
            flags[1, 4:] = 2  # meteo invalid
        flags[2, :] = 8 | 64  # informational bits only: kept
        return flags

    out = npp.process_year(2021, make_year(tmp_path, 2021, dmp, qflag))
    with rasterio.open(out) as src:
        total = src.read(1)
    assert total[0, 0] == total[0, 1] == total[0, 2] == npp.NODATA_VALUE
    assert np.all(total[1] == 1_825_000)
    assert np.all(total[2] == 1_825_000)


@pytest.mark.usefixtures("output_dir")
def test_incomplete_year_is_rejected(tmp_path):
    pairs = make_year(tmp_path, 2021, constant)
    with pytest.raises(ValueError, match="35 dekads"):
        npp.process_year(2021, pairs[:-1])


@pytest.mark.usefixtures("output_dir")
def test_mixed_product_versions_are_rejected(tmp_path):
    pairs = make_year(tmp_path / "a", 2021, constant)
    other = make_year(tmp_path / "b", 2021, constant, version="V1.1.1")
    with pytest.raises(ValueError, match="Mixed product versions"):
        npp.process_year(2021, pairs[:18] + other[18:])
