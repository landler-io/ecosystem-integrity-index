"""Fixtures for the pipeline tests: offline, on small synthetic inputs."""

from __future__ import annotations

import pytest


@pytest.fixture(scope="session", autouse=True)
def initialize_ee():
    """Override the repository-wide Earth Engine fixture; pipeline tests are offline."""
    yield
