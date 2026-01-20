"""
Google Earth Engine asset paths and metadata for EII data products.

This module contains the paths to pre-computed EII data assets in Google Earth Engine.
"""

from typing import Any

from ..config import EII_ASSET_ROOT

# Base path for EII assets
EII_ASSET_BASE = EII_ASSET_ROOT

# Base path for source datasets
DATASETS_BASE = "projects/landler-open-data/assets/datasets"

# Asset paths and metadata
ASSETS: dict[str, dict[str, Any]] = {
    "eii": {
        "path": f"{EII_ASSET_BASE}/global/eii_global_v1",
        "description": "Combined Ecosystem Integrity Index",
        "resolution": 300,
        "bands": ["eii"],
        "default_band": "eii",
    },
    "components": {
        "path": f"{EII_ASSET_BASE}/global/eii_global_v1",
        "description": "EII components from the global multiband asset",
        "resolution": 300,
        "bands": [
            "functional_integrity",
            "structural_integrity",
            "compositional_integrity",
        ],
    },
    "structural_core_area": {
        "path": f"{EII_ASSET_BASE}/products/v1/structural_integrity/core_area",
        "description": "Structural integrity (quality-weighted core area)",
        "resolution": 300,
        "bands": ["structural_integrity"],
    },
    "npp_predictions": {
        "path": f"{EII_ASSET_BASE}/predictions/npp",
        "description": "NPP model predictions",
        "resolution": 300,
        "bands": ["potential_npp", "actual_npp", "relative_npp", "npp_difference"],
    },
    "npp_model": {
        "path": f"{EII_ASSET_BASE}/models/potential_npp_classifier",
        "description": "Trained NPP potential model",
        "type": "classifier",
    },
}

# Source dataset paths (for reference)
SOURCE_DATASETS = {
    "chelsa_npp": f"{DATASETS_BASE}/chelsa/npp/chelsa_npp_1981_2010_v2-1",
    "clms_npp": f"{DATASETS_BASE}/clms/npp/annual",
    "soilgrids_sand": f"{DATASETS_BASE}/soilgrids/sand/sand_15-30cm_mean_gapfilled",
    "soilgrids_clay": f"{DATASETS_BASE}/soilgrids/clay/clay_15-30cm_mean_gapfilled",
    "soilgrids_phh2o": f"{DATASETS_BASE}/soilgrids/phh2o/phh2o_15-30cm_mean_gapfilled",
    "hmi_masks": f"{DATASETS_BASE}/natural_lands/hmi_masks/v1",
    "ecoregions": "RESOLVE/ECOREGIONS/2017",
    "wdpa": "WCMC/WDPA/current/polygons",
}

# Reference dates for data products
REFERENCE_DATES = {
    "eii": "2024-12-31",
    "landcover": "2023-01-01",
}


def get_asset_path(asset_name: str) -> str:
    """
    Get the Earth Engine asset path for a given asset name.

    Args:
        asset_name: Name of the asset (e.g., 'eii', 'components', 'npp_predictions')

    Returns:
        The full Earth Engine asset path.

    Raises:
        KeyError: If the asset name is not found.

    Example:
        >>> path = get_asset_path("eii")
        >>> print(path)
        'projects/landler-open-data/assets/eii/scores/eii'
    """
    if asset_name not in ASSETS:
        available = ", ".join(ASSETS.keys())
        raise KeyError(f"Unknown asset '{asset_name}'. Available: {available}")
    return str(ASSETS[asset_name]["path"])


def get_asset_info(asset_name: str) -> dict[str, Any]:
    """
    Get full metadata for an asset.

    Args:
        asset_name: Name of the asset

    Returns:
        Dictionary containing asset metadata (path, description, bands, etc.)
    """
    if asset_name not in ASSETS:
        available = ", ".join(ASSETS.keys())
        raise KeyError(f"Unknown asset '{asset_name}'. Available: {available}")
    return ASSETS[asset_name].copy()
