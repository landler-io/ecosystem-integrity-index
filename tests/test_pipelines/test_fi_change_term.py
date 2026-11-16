"""Grid helpers of the FI change-term run (pipelines/production/fi_change_term.py)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipelines" / "production"))

import evi_sources as es  # noqa: E402
import fi_change_term as fct  # noqa: E402


def test_tile_box_covers_tile_and_pixels_map_back() -> None:
    h, v = 18, 4
    lon0, lat0, lon1, lat1 = es.tile_bounds_lonlat(h, v)
    r0, c0, hgt, w = fct.tile_box(h, v)
    tr = fct.box_transform(h, v)
    assert tr.c <= lon0 and tr.c + w * tr.a >= lon1
    assert tr.f >= lat1 and tr.f + hgt * tr.e <= lat0

    # Every native pixel of the tile lies in the tile and gets a cell inside the box.
    n = 48
    inside = np.zeros(hgt * w, bool)
    groups = {
        "biome": np.arange(hgt * w, dtype="int64"),
        "cls": np.zeros(hgt * w, "int8"),
        "inside": inside,
    }
    cell, _ = fct.pixel_groups(h, v, n, groups)
    assert cell.min() >= 0 and cell.max() < hgt * w
    rr, cc = np.indices((n, n))
    size = es.TILE_M / n
    x = es.X0 + h * es.TILE_M + (cc.ravel() + 0.5) * size
    y = es.Y0 - v * es.TILE_M - (rr.ravel() + 0.5) * size
    lat = np.degrees(y / es.EARTH_R)
    lon = np.degrees(x / (es.EARTH_R * np.cos(np.radians(lat))))
    ok = np.abs(lon) <= 180
    th, tv = es.tile_of(lon[ok], lat[ok])
    assert (th == h).all() and (tv == v).all()


def test_composite_dates() -> None:
    modis = es.modis_composite_dates((2001, 2001))
    assert len(modis) == 23 and modis[0] == "2001_01_01" and modis[1] == "2001_01_17"
    viirs = es.viirs_composite_dates((2012, 2013))
    assert viirs[0] == "2012_01_17" and len(viirs) == 45
    assert es.item_date("MOD13Q1.A2024353.h18v04.061.x").isoformat() == "2024-12-18"
