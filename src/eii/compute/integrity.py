"""
Combination of the three integrity components into the Ecosystem Integrity Index.
"""

from __future__ import annotations

import numpy as np

from .settings import DEFAULT_AGGREGATION_METHOD, AggregationMethod


def compositional_integrity(bii: np.ndarray, land: np.ndarray) -> np.ndarray:
    """Biodiversity Intactness clamped to [0, 1], NaN off ``land``.

    BII v1.1 has no nodata value: sea, inland water and source gaps are exactly 0. The
    global driver sets those zeros to NaN before averaging onto the grid; ``land`` still
    masks the rest of the sea.
    """
    return np.where(land, np.clip(bii, 0, 1), np.nan).astype("float32")


def combine_components(
    functional: np.ndarray,
    structural: np.ndarray,
    compositional: np.ndarray,
    method: AggregationMethod = DEFAULT_AGGREGATION_METHOD,
) -> np.ndarray:
    """Combine integrity components into EII; NaN wherever any component is NaN.

    Args:
        functional: Functional integrity (0-1).
        structural: Structural integrity (0-1).
        compositional: Compositional integrity (0-1).
        method: ``minimum``; ``product``; ``geometric_mean``; or ``min_fuzzy_logic``
            (default), the minimum times the fuzzy sum ``y + z - yz`` of the other two.

    Returns:
        float32 EII (0-1).
    """
    components = np.stack([functional, structural, compositional]).astype("float64")
    if method == "minimum":
        eii = components.min(axis=0)
    elif method == "product":
        eii = components.prod(axis=0)
    elif method == "geometric_mean":
        eii = np.cbrt(components.prod(axis=0))
    elif method == "min_fuzzy_logic":
        low, mid, high = np.sort(components, axis=0)
        eii = low * (mid + high - mid * high)
    else:
        raise ValueError(f"Unknown method '{method}'")
    # np.sort moves NaN last, so mask explicitly.
    return np.where(np.isfinite(components).all(axis=0), eii, np.nan).astype("float32")


def limiting_component(
    functional: np.ndarray, structural: np.ndarray, compositional: np.ndarray
) -> np.ndarray:
    """Smallest component: 1 = functional, 2 = structural, 3 = compositional, 0 if any is NaN.

    Ties prefer functional, then structural, then compositional.
    """
    components = np.stack([functional, structural, compositional])
    valid = np.isfinite(components).all(axis=0)
    code = np.argmin(np.where(valid, components, 0), axis=0) + 1
    return np.where(valid, code, 0).astype("uint8")
