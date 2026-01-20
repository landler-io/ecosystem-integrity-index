"""
Ecosystem Integrity Index (EII)

A Python package for accessing and computing the Ecosystem Integrity Index,
a global measure of ecosystem condition based on deviation from natural
potential Net Primary Productivity (NPP).

Submodules:
    - eii.client: Lightweight retrieval of pre-computed EII from GEE assets
    - eii.compute: Full EII calculation from scratch
    - eii.training: Model training utilities (optional dependencies)

Quick Start:
    >>> import ee
    >>> from eii.client import get_eii_for_polygon
    >>>
    >>> ee.Initialize()
    >>> polygon = ee.Geometry.Rectangle([-60, -10, -55, -5])
    >>> stats = get_eii_for_polygon(polygon)
"""

try:
    from ._version import __version__
except ImportError:
    __version__ = "0.0.0.dev0"

__all__ = ["__version__"]
