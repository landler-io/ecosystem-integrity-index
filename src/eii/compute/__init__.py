"""
EII v2 calculation on local rasters.

Pure numpy functions for the three integrity components and their combination. The
global production run (10 degree tiles, inputs, publication) is
``pipelines/production/compute_eii_v2.py``; results are retrieved with ``eii.client``.

Requires additional dependencies: pip install "ecosystem-integrity-index[compute]"

Example:
    >>> import numpy as np
    >>> from eii.compute import combine_components
    >>> float(combine_components(np.array(0.5), np.array(1.0), np.array(1.0)))
    0.5
"""

from importlib import import_module
from typing import Any

_MODULES = {
    "observed_npp": "functional",
    "observed_npp_window": "functional",
    "deviation_tolerance_integrity": "functional",
    "recentring_sums": "functional",
    "recentring_factors": "functional",
    "recentre_potential": "functional",
    "wetland_flag": "functional",
    "structural_integrity": "structural",
    "structural_halo_m": "structural",
    "habitat_quality": "structural",
    "compositional_integrity": "integrity",
    "combine_components": "integrity",
    "limiting_component": "integrity",
    "encode_modis_evi": "change",
    "encode_viirs_evi": "change",
    "annual_evi": "change",
    "complete_years": "change",
    "natural_class": "change",
    "correct_drift": "change",
    "chain_factor": "change",
    "chained": "change",
    "baseline_median": "change",
    "change_ratio": "change",
    "combine_ratios": "change",
    "tolerance_h": "change",
    "fi_change": "change",
    "functional_with_change": "change",
    "series_for": "change",
    "constants_sha256": "change",
}

__all__ = list(_MODULES)


def __getattr__(name: str) -> Any:
    """Import numpy/scipy code on first use, so ``eii.compute.settings`` stays light."""
    if name not in _MODULES:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f"{__name__}.{_MODULES[name]}"), name)
    globals()[name] = value
    return value
