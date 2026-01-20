"""Tests for eii.client.utils module."""

import ee
import pytest

from eii.client.utils import normalize_client_input


class TestNormalizeClientInput:
    """Tests for normalize_client_input function."""

    def test_normalize_to_geometry(
        self,
        amazon_small_polygon,
        amazon_point,
        simple_feature,
        simple_feature_collection,
        bbox_tuple,
    ):
        """Test all input types normalize to ee.Geometry."""
        # ee.Geometry passthrough
        result = normalize_client_input(amazon_small_polygon, target="geometry")
        assert isinstance(result, ee.Geometry)
        info = result.getInfo()
        assert info["type"] in ("Polygon", "Rectangle")

        # Point geometry
        result = normalize_client_input(amazon_point, target="geometry")
        assert result.getInfo()["type"] == "Point"

        # ee.Feature extracts geometry
        result = normalize_client_input(simple_feature, target="geometry")
        assert isinstance(result, ee.Geometry)

        # ee.FeatureCollection dissolves to geometry
        result = normalize_client_input(simple_feature_collection, target="geometry")
        assert isinstance(result, ee.Geometry)
        assert result.getInfo()["type"] in (
            "Polygon",
            "MultiPolygon",
            "GeometryCollection",
        )

        # Bbox tuple creates Rectangle
        result = normalize_client_input(bbox_tuple, target="geometry")
        assert isinstance(result, ee.Geometry)
        assert result.getInfo()["type"] in ("Polygon", "Rectangle")

    def test_normalize_to_features(
        self,
        amazon_small_polygon,
        simple_feature,
        simple_feature_collection,
        bbox_tuple,
    ):
        """Test all input types normalize to ee.FeatureCollection."""
        # ee.FeatureCollection passthrough
        result = normalize_client_input(simple_feature_collection, target="features")
        assert isinstance(result, ee.FeatureCollection)
        assert result.size().getInfo() == 2

        # ee.Feature wraps in collection
        result = normalize_client_input(simple_feature, target="features")
        assert isinstance(result, ee.FeatureCollection)
        assert result.size().getInfo() == 1

        # ee.Geometry wraps in Feature then Collection
        result = normalize_client_input(amazon_small_polygon, target="features")
        assert isinstance(result, ee.FeatureCollection)
        assert result.size().getInfo() == 1

        # Bbox tuple creates single-feature collection
        result = normalize_client_input(bbox_tuple, target="features")
        assert isinstance(result, ee.FeatureCollection)
        assert result.size().getInfo() == 1

    def test_geo_interface_inputs(self):
        """Test shapely and geopandas inputs via __geo_interface__."""
        pytest.importorskip("shapely")
        from shapely.geometry import box

        # Shapely geometry to ee.Geometry
        shapely_poly = box(-63.0, -4.0, -62.5, -3.5)
        result = normalize_client_input(shapely_poly, target="geometry")
        assert isinstance(result, ee.Geometry)
        assert result.getInfo()["type"] == "Polygon"

        # Shapely to FeatureCollection
        result = normalize_client_input(shapely_poly, target="features")
        assert isinstance(result, ee.FeatureCollection)
        assert result.size().getInfo() == 1

        # GeoDataFrame support
        gpd = pytest.importorskip("geopandas")
        gdf = gpd.GeoDataFrame(geometry=[shapely_poly], crs="EPSG:4326")

        result = normalize_client_input(gdf, target="geometry")
        assert isinstance(result, ee.Geometry)

        result = normalize_client_input(gdf, target="features")
        assert isinstance(result, ee.FeatureCollection)

    def test_crs_handling(self, amazon_small_polygon):
        """Test CRS handling and reprojection."""
        # EPSG:4326 coords stay valid
        result = normalize_client_input(amazon_small_polygon, target="geometry")
        bounds = result.bounds().getInfo()
        coords = bounds["coordinates"][0]
        lons = [c[0] for c in coords]
        lats = [c[1] for c in coords]
        assert all(-180 <= lon <= 180 for lon in lons)
        assert all(-90 <= lat <= 90 for lat in lats)

        # Projected geometry reprojected to WGS84
        coords = [
            [500000, 5500000],
            [510000, 5500000],
            [510000, 5510000],
            [500000, 5510000],
        ]
        projected_geom = ee.Geometry.Polygon([coords], proj="EPSG:32632", geodesic=False)
        result = normalize_client_input(projected_geom, target="geometry")
        bounds = result.bounds().getInfo()
        bound_coords = bounds["coordinates"][0]
        lons = [c[0] for c in bound_coords]
        lats = [c[1] for c in bound_coords]
        # Should be valid WGS84 in central Europe
        assert all(-180 <= lon <= 180 for lon in lons)
        assert all(-90 <= lat <= 90 for lat in lats)
        assert any(5 <= lon <= 15 for lon in lons)
        assert any(45 <= lat <= 55 for lat in lats)

    def test_invalid_inputs_raise(self, amazon_small_polygon):
        """Test error handling for invalid inputs."""
        # Invalid target
        with pytest.raises(ValueError, match="target must be"):
            normalize_client_input(amazon_small_polygon, target="invalid")

        # Unsupported input type
        with pytest.raises(ValueError, match="Unsupported"):
            normalize_client_input("not a geometry", target="geometry")

        with pytest.raises(ValueError, match="Unsupported"):
            normalize_client_input("not a feature", target="features")
