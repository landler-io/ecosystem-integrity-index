"""Tests for eii.client.retrieve module."""

import tempfile
from pathlib import Path

import ee
import pytest

from eii.client.retrieve import (
    _build_reducer,
    _validate_stats_params,
    get_layers,
    get_raster,
    get_stats,
)


class TestRetrieve:
    """Tests for retrieve module functions."""

    def test_stats_param_validation(self):
        """Test parameter validation (no GEE calls)."""
        # Default stats
        assert _validate_stats_params(None, None) == ["mean"]
        assert _validate_stats_params(["mean", "min", "max"], None) == ["mean", "min", "max"]
        assert _validate_stats_params(["mean"], [10, 50, 90]) == ["mean"]

        # Invalid stats
        with pytest.raises(ValueError, match="Invalid stats"):
            _validate_stats_params(["invalid"], None)

        # Invalid percentiles
        with pytest.raises(ValueError, match="Percentiles must be integers"):
            _validate_stats_params(["mean"], [150])
        with pytest.raises(ValueError, match="Percentiles must be integers"):
            _validate_stats_params(["mean"], [10.5])

        # Reducer creation
        assert isinstance(_build_reducer(["mean"]), ee.Reducer)
        assert isinstance(_build_reducer(["mean", "min", "max"], [10, 90]), ee.Reducer)
        with pytest.raises(ValueError, match="No statistics selected"):
            _build_reducer([])

    def test_get_layers(self):
        """Test layer selection and structure."""
        # 'all' includes EII and components
        layers = get_layers(layers="all")
        assert isinstance(layers, dict)
        assert all(key in layers for key in ["eii", "functional", "structural", "compositional"])
        for layer in layers.values():
            assert isinstance(layer, ee.Image)

        # 'eii' only
        layers = get_layers(layers="eii")
        assert "eii" in layers and "functional" not in layers
        assert "eii" in layers["eii"].bandNames().getInfo()

        # 'components' only
        layers = get_layers(layers="components")
        assert "eii" not in layers and "functional" in layers

        # Mode restrictions
        with pytest.raises(ValueError, match="pre-computed default method"):
            get_layers(aggregation_method="minimum", compute_mode="precomputed")
        with pytest.raises(ValueError):
            get_layers(compute_mode="invalid")

    def test_get_stats_polygon(self, amazon_small_polygon):
        """Test comprehensive stats extraction for polygon."""
        result = get_stats(
            amazon_small_polygon,
            stats=["mean", "min", "max", "std"],
            percentiles=[10, 90],
            include_components=True,
        )

        # Structure
        assert result["geometry_type"] in ("Polygon", "Rectangle")
        assert "values" in result

        # EII stats present
        eii = result["values"]["eii"]
        assert all(stat in eii for stat in ["mean", "min", "max", "std", "p10", "p90"])

        # Values in range
        for stat in ["mean", "min", "max"]:
            if eii[stat] is not None:
                assert 0 <= eii[stat] <= 1, f"EII {stat}={eii[stat]} out of range"

        # Components present and valid
        for comp in ["functional_integrity", "structural_integrity", "compositional_integrity"]:
            assert comp in result["values"]
            val = result["values"][comp]["mean"]
            if val is not None:
                assert 0 <= val <= 1, f"{comp}={val} out of range"

    def test_get_stats_variants(self, amazon_point, bbox_tuple, amazon_small_polygon):
        """Test point geometry, bbox input, and output formats."""
        # Point returns single value
        result = get_stats(amazon_point, include_components=False)
        assert result["geometry_type"] == "Point"
        assert isinstance(result["values"]["eii"], float | int | type(None))

        # Bbox tuple works
        result = get_stats(bbox_tuple, include_components=False)
        assert "eii" in result["values"]

        # Exclude components
        result = get_stats(amazon_small_polygon, include_components=False)
        assert "eii" in result["values"]
        assert "functional_integrity" not in result["values"]

        # GeoDataFrame output
        gpd = pytest.importorskip("geopandas")
        result = get_stats(amazon_small_polygon, output_format="geodataframe")
        assert isinstance(result, gpd.GeoDataFrame)

        # Validation errors
        with pytest.raises(ValueError, match="Invalid stats"):
            get_stats(amazon_small_polygon, stats=["invalid"])

    def test_get_raster_memory(self):
        """Test raster download to xarray Dataset."""
        xr = pytest.importorskip("xarray")
        np = pytest.importorskip("numpy")

        tiny_polygon = ee.Geometry.Rectangle([-62.6, -3.6, -62.5, -3.5])
        result = get_raster(tiny_polygon, scale=1000, include_components=True)

        # Structure
        assert isinstance(result, xr.Dataset)
        assert all(
            var in result.data_vars
            for var in [
                "eii",
                "functional_integrity",
                "structural_integrity",
                "compositional_integrity",
            ]
        )
        assert "x" in result.coords and "y" in result.coords
        assert "crs" in result.attrs

        # Values in range
        valid = result["eii"].values[~np.isnan(result["eii"].values)]
        if len(valid) > 0:
            assert valid.min() >= 0 and valid.max() <= 1

        # Without components
        result = get_raster(tiny_polygon, scale=1000, include_components=False)
        assert "eii" in result.data_vars
        assert "functional_integrity" not in result.data_vars

        # dtype parameter
        result = get_raster(tiny_polygon, scale=1000, dtype="float32", include_components=False)
        assert result["eii"].dtype == np.float32

    def test_get_raster_geotiff(self):
        """Test GeoTIFF output and chunking behavior."""
        pytest.importorskip("rasterio")
        pytest.importorskip("xarray")

        tiny_polygon = ee.Geometry.Rectangle([-62.6, -3.6, -62.5, -3.5])

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "test.tif"
            result = get_raster(
                tiny_polygon, scale=1000, output_format="geotiff", out_path=str(out_path)
            )
            assert isinstance(result, Path)
            assert result.exists()
            assert result.suffix == ".tif"

        # Chunking error with restrictive limits
        medium = ee.Geometry.Rectangle([-63.0, -4.0, -62.0, -3.0])
        with pytest.raises(ValueError, match="exceeds configured limits"):
            get_raster(medium, scale=100, chunking="never", max_area_sq_km=1)
