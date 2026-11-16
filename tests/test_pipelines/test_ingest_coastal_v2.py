"""Offline checks for the coastal predictor rasters: grid alignment and the side-table rules."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipelines/preprocessing"))
import ingest_coastal_v2 as ing  # noqa: E402


def test_grid_is_ee_300m_grid_and_global():
    assert ing.TRANSFORM.c / ing.PIXEL == -ing.COLS / 2  # aligned to (0, 0) like EE
    assert ing.TRANSFORM.c <= -180 and ing.TRANSFORM.c + ing.COLS * ing.PIXEL >= 180
    assert 89.99 < ing.TRANSFORM.f <= 90


def test_strip_matches_brute_force_and_side_table_rules(monkeypatch):
    monkeypatch.setattr(ing, "COLS", 256)
    row0 = ing.ROWS // 2 - 256
    lon, lat = ing.pixel_centres(np.array([row0 + 100]), np.array([20]))
    coast = ing.cKDTree(ing.unit_xyz(lon, lat))
    ocean = (ing.cKDTree(ing.unit_xyz(lon, lat)), np.array([100.0]))  # 100 cm amplitude
    salinity = (ocean[0], np.array([35.0]))
    dist, tidal, sss = ing.strip_values(row0, coast, [ocean] * 4, salinity)

    rows, cols = np.indices(dist.shape)
    xyz = ing.unit_xyz(*ing.pixel_centres(row0 + rows.ravel(), cols.ravel()))
    brute = ing.coastal_values(xyz, coast, [ocean] * 4, salinity)
    for raster, expected in zip((dist, tidal, sss), brute, strict=True):
        np.testing.assert_array_equal(raster.ravel(), expected)

    assert dist[100, 20] == 0 and dist[100, 250] == 50000
    km = ing.chord_to_km(np.linalg.norm(xyz - ing.unit_xyz(lon, lat), axis=1))
    near = km < ing.COAST_KM
    np.testing.assert_allclose(dist.ravel()[near], km[near] * 1000, atol=0.5)
    assert set(np.unique(tidal.ravel()[near])) == {8.0}  # 2 x 4 x 100 cm
    assert set(np.unique(sss.ravel()[near])) == {35.0}
    assert set(np.unique(tidal.ravel()[~near])) == {ing.NODATA}


def test_manifest_records_licence_and_pins_source():
    blob = type("Blob", (), {"bucket": type("B", (), {"name": "bucket"}), "name": "x.tif"})()
    blob.generation, blob.crc32c = 7, "crc"
    for name in ing.DATASETS:
        manifest = ing.manifest(name, blob, {"sha256": "sha"})
        assert manifest["name"].endswith(f"datasets/{name}/v1_0")
        assert {"license", "attribution", "author", "source"} <= manifest["properties"].keys()
        assert ("missingData" in manifest["bands"][0]) == (name != "osm_coast_distance")
    assert "ODbL" in ing.DATASETS["osm_coast_distance"]["properties"]["share_alike"]
