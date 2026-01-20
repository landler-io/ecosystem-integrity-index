"""
Shared fixtures for EII tests.

These tests use real Google Earth Engine calls. Run with authenticated credentials:
    ee.Authenticate()  # one-time setup
    pytest tests/
"""

from __future__ import annotations

import ee
import pytest


@pytest.fixture(scope="session", autouse=True)
def initialize_ee():
    """Initialize Earth Engine once per test session."""
    try:
        ee.Initialize(opt_url="https://earthengine-highvolume.googleapis.com")
    except Exception:
        ee.Initialize()


# Small test geometries to keep GEE calls fast


@pytest.fixture
def amazon_point() -> ee.Geometry:
    """A point in the Amazon rainforest (high integrity expected)."""
    return ee.Geometry.Point([-62.5, -3.5])


@pytest.fixture
def amazon_small_polygon() -> ee.Geometry:
    """A small polygon in the Amazon (~50km x 50km)."""
    return ee.Geometry.Rectangle([-63.0, -4.0, -62.5, -3.5])


@pytest.fixture
def europe_small_polygon() -> ee.Geometry:
    """A small polygon in central Europe (lower integrity expected)."""
    return ee.Geometry.Rectangle([10.0, 48.0, 10.5, 48.5])


@pytest.fixture
def bbox_tuple() -> tuple[float, float, float, float]:
    """A bbox tuple for testing input normalization."""
    return (-63.0, -4.0, -62.5, -3.5)


@pytest.fixture
def simple_feature(amazon_small_polygon) -> ee.Feature:
    """A simple ee.Feature wrapping a polygon."""
    return ee.Feature(amazon_small_polygon, {"name": "test_feature"})


@pytest.fixture
def simple_feature_collection(amazon_small_polygon, europe_small_polygon) -> ee.FeatureCollection:
    """A FeatureCollection with two features for zonal stats testing."""
    return ee.FeatureCollection(
        [
            ee.Feature(amazon_small_polygon, {"name": "amazon", "id": 1}),
            ee.Feature(europe_small_polygon, {"name": "europe", "id": 2}),
        ]
    )
