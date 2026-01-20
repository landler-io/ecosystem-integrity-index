"""
Structural integrity calculation using quality-weighted core area metrics.

Structural integrity measures landscape fragmentation AND habitat quality by:
1. Identifying core habitat (interior areas unaffected by edge effects)
2. Weighting core pixels by habitat quality class based on HMI

Quality Classes:
    HMI < 0.1: Pristine (weight 4)
    HMI 0.1-0.2: Low-impact (weight 3)
    HMI 0.2-0.3: Moderate (weight 2)
    HMI 0.3-0.4: Semi-natural (weight 1)
    HMI >= 0.4: Modified (weight 0)

Score interpretation:
    1.0 = All pristine core habitat
    0.25 = All semi-natural core habitat
    0.0 = No core habitat (fragmented or modified)
"""

from __future__ import annotations

import ee

HMI_ASSET = "projects/sat-io/open-datasets/GHM/HM_2022_300M"
HMI_THRESHOLD = 0.4
DEFAULT_EDGE_DEPTH_M = 300
DEFAULT_NEIGHBORHOOD_M = 5000
DEFAULT_SCALE_M = 300

QUALITY_THRESHOLDS = [0.1, 0.2, 0.3, 0.4]
QUALITY_WEIGHTS = [4, 3, 2, 1]
MAX_QUALITY_WEIGHT = 4


def _create_quality_class(hmi: ee.Image) -> ee.Image:
    """
    Create quality class image from HMI.

    Higher values indicate higher habitat quality (less human modification).

    Args:
        hmi: Human Modification Index image (0-1 scale).

    Returns:
        Quality class image (0-4) where 4=pristine, 1=semi-natural, 0=modified.
    """
    quality = ee.Image(0)
    quality = quality.where(hmi.lt(0.4), 1)  # Semi-natural
    quality = quality.where(hmi.lt(0.3), 2)  # Moderate
    quality = quality.where(hmi.lt(0.2), 3)  # Low-impact
    quality = quality.where(hmi.lt(0.1), 4)  # Pristine
    return quality.rename("quality_class")


def calculate_structural_integrity(
    aoi: ee.Geometry | None = None,
    edge_depth_m: int = DEFAULT_EDGE_DEPTH_M,
    neighborhood_m: int = DEFAULT_NEIGHBORHOOD_M,
    scale_m: int = DEFAULT_SCALE_M,
) -> ee.Image:
    """
    Calculate structural integrity using quality-weighted core area.

    Core area is habitat that survives erosion by the edge depth. Each core
    pixel is weighted by its habitat quality class (pristine=4, semi-natural=1).
    The final score reflects both fragmentation AND habitat quality.

    Args:
        aoi: Area of interest geometry. If None, returns global unclipped image
            suitable for tiled export pipelines.
        edge_depth_m: Edge effect penetration depth in meters. Areas within
            this distance of habitat boundaries are considered edge-affected.
            Default 300m based on literature (typical range 100-500m).
        neighborhood_m: Landscape analysis radius in meters. Default 5km
            balances local actionability with landscape-scale processes.
        scale_m: Output scale in meters for reprojection. Default 300m matches
            the native HMI resolution.

    Returns:
        Structural integrity score (0-1) where:
        - 1.0 = all pristine core habitat
        - 0.25 = all semi-natural core habitat
        - 0.0 = no core habitat (highly fragmented or modified)
    """
    hmi = _load_hmi()

    habitat_binary = hmi.lt(HMI_THRESHOLD).rename("habitat")

    core_habitat = habitat_binary.focal_min(
        radius=edge_depth_m,
        kernelType="circle",
        units="meters",
    ).rename("core")

    quality_class = _create_quality_class(hmi)

    # Weighted core is 0 for non-habitat and edge pixels, 1-4 for core based on quality
    # Do NOT mask - we want 0s to contribute to the neighborhood mean
    weighted_core = core_habitat.multiply(quality_class)

    structural_integrity = (
        weighted_core.reduceNeighborhood(
            reducer=ee.Reducer.mean(),
            kernel=ee.Kernel.circle(radius=neighborhood_m, units="meters"),
        )
        .divide(MAX_QUALITY_WEIGHT)
        .unmask(0)
        .reproject(crs="EPSG:4326", scale=scale_m)
        .rename("structural_integrity")
    )

    if aoi is not None:
        structural_integrity = structural_integrity.clip(aoi)

    return structural_integrity


def _load_hmi() -> ee.Image:
    """Load Global Human Modification Index (2022, all threats combined)."""
    return ee.ImageCollection(HMI_ASSET).filter(ee.Filter.eq("threat_code", "AA")).first()
