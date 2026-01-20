"""Tests for eii.client.analysis module."""

import ee
import pytest

from eii.client.analysis import compare_methods, get_zonal_stats


class TestAnalysis:
    """Tests for analysis module functions."""

    def test_zonal_stats_comprehensive(self, simple_feature_collection):
        """Test zonal stats with all features."""
        gpd = pytest.importorskip("geopandas")

        result = get_zonal_stats(
            simple_feature_collection,
            stats=["mean", "min", "max"],
            percentiles=[10, 90],
            include_components=True,
        )

        # Type and row count
        assert isinstance(result, gpd.GeoDataFrame)
        assert len(result) == 2
        assert result.geometry is not None

        # Columns present (MultiIndex)
        cols = result.columns.get_level_values(0)
        assert "eii" in cols
        # Components may have _integrity suffix or not
        assert any(c in cols for c in ["functional_integrity", "functional"])

        # Stats present
        stats_cols = result.columns.get_level_values(1)
        assert all(s in stats_cols for s in ["mean", "min", "max", "p10", "p90"])

        # Values in range
        eii_mean = result[("eii", "mean")]
        for val in eii_mean.dropna():
            assert 0 <= val <= 1, f"EII mean {val} out of range"

    def test_zonal_stats_inputs_and_options(self, simple_feature_collection, simple_feature):
        """Test various input types and options."""
        gpd = pytest.importorskip("geopandas")
        from shapely.geometry import box

        # Single feature
        fc = ee.FeatureCollection([simple_feature])
        result = get_zonal_stats(fc, include_components=False)
        assert len(result) == 1
        assert "eii" in result.columns.get_level_values(0)
        assert "functional_integrity" not in result.columns.get_level_values(0)

        # GeoDataFrame input
        gdf = gpd.GeoDataFrame(
            {"name": ["amazon", "europe"]},
            geometry=[box(-63, -4, -62.5, -3.5), box(10, 48, 10.5, 48.5)],
            crs="EPSG:4326",
        )
        result = get_zonal_stats(gdf, include_components=False)
        assert len(result) == 2

        # keep_columns preserves properties
        result = get_zonal_stats(
            simple_feature_collection, keep_columns=["name"], include_components=False
        )
        assert "name" in result.columns.get_level_values(0)

    @pytest.mark.slow
    def test_compare_methods_on_the_fly(self):
        """Test compare_methods with on_the_fly compute (slow)."""
        tiny = ee.Geometry.Rectangle([-62.55, -3.55, -62.5, -3.5])
        result = compare_methods(tiny, compute_mode="on_the_fly", scale=1000)

        assert isinstance(result, ee.Dictionary)

        info = result.getInfo()
        expected = ["minimum", "product", "min_fuzzy_logic", "geometric_mean"]
        for method in expected:
            assert method in info, f"Missing method: {method}"
            val = info[method]
            assert isinstance(val, int | float | type(None)), f"{method} not numeric"
            if val is not None:
                assert 0 <= val <= 1, f"{method}={val} out of range"

        # Methods should produce different values
        values = [v for v in info.values() if v is not None]
        if len(values) >= 2:
            assert {round(v, 4) for v in values}.__len__() >= 2, "Methods should differ"

    def test_compare_methods_precomputed_and_errors(self):
        """Test precomputed mode behavior and error handling."""
        tiny = ee.Geometry.Rectangle([-62.55, -3.55, -62.5, -3.5])

        # Precomputed only supports default method, so compare_methods fails
        with pytest.raises(ValueError, match="pre-computed default method"):
            compare_methods(tiny, compute_mode="precomputed").getInfo()

        # Invalid compute_mode
        with pytest.raises(ValueError):
            compare_methods(tiny, compute_mode="invalid").getInfo()

        # But single method in precomputed works via get_stats
        from eii.client import get_stats

        result = get_stats(tiny, include_components=False, scale=1000)
        val = result["values"]["eii"]["mean"]
        assert val is None or 0 <= val <= 1
