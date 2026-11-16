"""
Google Earth Engine asset paths and metadata for EII data products.

This module contains the paths to pre-computed EII data assets in Google Earth
Engine. Asset paths are versioned: every product lives under a dataset version
(e.g. ``v1``, ``v2``). The version can be chosen in three ways, in order of
precedence:

1. Per call - pass ``version=...`` to a retrieval function or ``get_asset_path``.
2. Per session - ``set_version("v1")`` or the ``use_version("v1")`` context manager.
3. Project default - the ``gee.version`` field in ``config.yaml``.

v2 products live under ``V2_PRODUCTS_BASE``. ``gee.asset_root`` in ``config.yaml`` only sets the
root of the v1 asset paths.
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Any

from ..config import EII_ASSET_ROOT
from ..config import VERSION as DEFAULT_VERSION

# Root of the v1 asset paths (``gee.asset_root``); v2 paths do not depend on it
EII_ASSET_BASE = EII_ASSET_ROOT

# v2 products are cloud-optimized GeoTIFFs in GCS, registered as Earth Engine image collections
V2_PRODUCTS_BASE = "projects/landler-open-data/assets/products/eii/v2_0"
V2_RELEASE = f"{V2_PRODUCTS_BASE}/release/global"

# Asset definitions. Each ``template`` carries a ``{version}`` placeholder that is
# filled in by :func:`get_asset_path` using the requested or active version. A key
# named after a version (e.g. ``"v2"``) overrides the template for that version; v2
# release assets are image collections with one image per 10 degree tile and
# reference year (property ``reference_year``).
_ASSET_DEFS: dict[str, dict[str, Any]] = {
    "eii": {
        "template": f"{EII_ASSET_BASE}/global/eii_global_{{version}}",
        "v2": V2_RELEASE,
        "description": "Combined Ecosystem Integrity Index",
        "resolution": 300,
        "bands": ["eii"],
        "default_band": "eii",
    },
    "components": {
        "template": f"{EII_ASSET_BASE}/global/eii_global_{{version}}",
        "v2": V2_RELEASE,
        "description": "EII components from the global multiband asset",
        "resolution": 300,
        "bands": [
            "functional_integrity",
            "structural_integrity",
            "compositional_integrity",
        ],
    },
    "wetland_flag": {
        "template": f"{EII_ASSET_BASE}/global/eii_global_{{version}}",
        "v2": V2_RELEASE,
        "description": (
            "v2 wetland flag from GLWD v2 (0 none, 1 wetland, 2 peatland, 3 mangrove); "
            "an annotation where functional integrity is less reliable, not a score"
        ),
        "resolution": 300,
        "bands": ["wetland_flag"],
        "versions": ["v2"],
    },
    "natural_npp": {
        "template": (
            f"{EII_ASSET_BASE}/intermediate/functional/predictions/natural_npp_{{version}}"
        ),
        "v2": f"{V2_PRODUCTS_BASE}/intermediate/functional/potential_npp_olci_gapfill",
        "description": (
            "Pre-computed natural (potential) NPP: v1 mean/std; v2 potential_npp, "
            "lower_68, upper_68, lower_90, upper_90 in gC m-2 yr-1"
        ),
    },
}

# Active dataset version for the current session. Defaults to config.yaml.
_active_version: str = DEFAULT_VERSION


def get_version() -> str:
    """Return the dataset version currently used for asset resolution."""
    return _active_version


def set_version(version: str) -> None:
    """
    Set the active dataset version for the current session.

    Affects every subsequent retrieval that does not pass an explicit ``version``.

    Args:
        version: Dataset version string, e.g. "v1" or "v2".
    """
    global _active_version
    if not isinstance(version, str) or not version:
        raise ValueError(f"version must be a non-empty string, got {version!r}")
    _active_version = version


@contextmanager
def use_version(version: str):
    """
    Temporarily set the active dataset version within a ``with`` block.

    Example:
        >>> with use_version("v1"):
        ...     stats = get_stats(polygon)  # reads v1 assets
    """
    previous = _active_version
    set_version(version)
    try:
        yield version
    finally:
        set_version(previous)


def resolve_version(version: str | None = None) -> str:
    """Return ``version`` if given, otherwise the active session version."""
    return version if version is not None else _active_version


def get_asset_path(asset_name: str, version: str | None = None) -> str:
    """
    Get the Earth Engine asset path for a given asset name and dataset version.

    Args:
        asset_name: Name of the asset (e.g., 'eii', 'components', 'natural_npp').
        version: Dataset version. If None, uses the active session version
            (see :func:`set_version` / :func:`use_version`), which defaults to
            ``config.yaml``.

    Returns:
        The full Earth Engine asset path.

    Raises:
        KeyError: If the asset name is not found.

    Example:
        >>> get_asset_path("eii", version="v2")
        'projects/landler-open-data/assets/products/eii/v2_0/release/global'
        >>> get_asset_path("eii", version="v1")
        'projects/landler-open-data/assets/eii/global/eii_global_v1'
    """
    if asset_name not in _ASSET_DEFS:
        available = ", ".join(_ASSET_DEFS.keys())
        raise KeyError(f"Unknown asset '{asset_name}'. Available: {available}")
    version = resolve_version(version)
    definition = _ASSET_DEFS[asset_name]
    return str(definition.get(version) or definition["template"].format(version=version))


def get_asset_info(asset_name: str, version: str | None = None) -> dict[str, Any]:
    """
    Get full metadata for an asset, with its path resolved for the given version.

    Args:
        asset_name: Name of the asset.
        version: Dataset version. If None, uses the active session version.

    Returns:
        Dictionary containing asset metadata (path, description, bands, etc.).
    """
    if asset_name not in _ASSET_DEFS:
        available = ", ".join(_ASSET_DEFS.keys())
        raise KeyError(f"Unknown asset '{asset_name}'. Available: {available}")
    info = {
        k: v
        for k, v in _ASSET_DEFS[asset_name].items()
        if k != "template" and not (k[0] == "v" and k[1:].isdigit())
    }
    info["path"] = get_asset_path(asset_name, version)
    return info


# Backward-compatible static registry, resolved at the default (config) version.
# Prefer get_asset_path()/get_asset_info() for version-aware access.
ASSETS: dict[str, dict[str, Any]] = {name: get_asset_info(name) for name in _ASSET_DEFS}
