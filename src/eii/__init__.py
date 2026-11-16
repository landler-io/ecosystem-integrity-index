"""
Ecosystem Integrity Index (EII)

A Python package for accessing and computing the Ecosystem Integrity Index, a global 0-1
measure of ecosystem condition that combines functional, structural and compositional
integrity (1 = fully intact).

Submodules:
    - eii.client: Retrieval of the published EII from Earth Engine assets
    - eii.compute: The v2 calculation as numpy functions (optional dependencies)
    - eii.training: Potential-NPP model training utilities (optional dependencies)
    - eii.modulation: Natural Capital modulation of EII

Quick Start:
    >>> import ee
    >>> from eii.client import get_stats
    >>>
    >>> ee.Initialize()
    >>> polygon = ee.Geometry.Rectangle([-60, -10, -55, -5])
    >>> stats = get_stats(polygon)
"""

try:
    from ._version import __version__
except ImportError:
    __version__ = "0.0.0.dev0"

__all__ = ["__version__"]
