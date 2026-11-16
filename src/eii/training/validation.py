"""
Spatial cross-validation folds of the potential-NPP model.

Fold generators return one integer fold id per row, so the same folds can be stored as a
column and reused. ``spatial_buffer_mask`` removes training rows near the test blocks.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from sklearn.model_selection import GroupKFold

N_FOLDS = 5
SPATIAL_BUFFER_DEG = 1.0


def spatial_block_ids(lat: Any, lon: Any, block_deg: float) -> np.ndarray:
    """Integer id of the ``block_deg`` lat/lon block of each point (wraps at the antimeridian)."""
    n_cols = int(np.ceil(360.0 / block_deg))
    row = np.floor((np.asarray(lat, dtype="float64") + 90.0) / block_deg).astype("int64")
    col = np.floor((np.asarray(lon, dtype="float64") + 180.0) / block_deg).astype("int64")
    return np.asarray(row * n_cols + np.mod(col, n_cols))


def group_kfold_folds(groups: Any, n_folds: int = N_FOLDS) -> np.ndarray:
    """Assign whole groups to ``n_folds`` folds of similar size (sklearn ``GroupKFold``)."""
    groups = np.asarray(groups)
    folds = np.full(len(groups), -1, dtype="int64")
    splitter = GroupKFold(n_splits=n_folds)
    for fold, (_, test_index) in enumerate(splitter.split(np.zeros(len(groups)), groups=groups)):
        folds[test_index] = fold
    return folds


def spatial_buffer_mask(
    lat: Any,
    lon: Any,
    test_mask: Any,
    block_deg: float,
    buffer_deg: float = SPATIAL_BUFFER_DEG,
) -> np.ndarray:
    """
    Rows lying within ``buffer_deg`` (Chebyshev, in degrees) of a block that holds test rows.

    Dropping these rows from the training side only keeps every test
    point while removing near-duplicates across the block edge. Test rows are
    always inside their own block and are therefore flagged as well; callers use
    the mask as ``train = ~test & ~mask``.
    """
    lat = np.asarray(lat, dtype="float64")
    lon = np.asarray(lon, dtype="float64")
    test_blocks = np.unique(spatial_block_ids(lat[test_mask], lon[test_mask], block_deg))
    n_cols = int(np.ceil(360.0 / block_deg))
    row_lo = np.floor((lat - buffer_deg + 90.0) / block_deg).astype("int64")
    row_hi = np.floor((lat + buffer_deg + 90.0) / block_deg).astype("int64")
    col_lo = np.floor((lon - buffer_deg + 180.0) / block_deg).astype("int64")
    col_hi = np.floor((lon + buffer_deg + 180.0) / block_deg).astype("int64")
    span = int(np.ceil(2 * buffer_deg / block_deg)) + 1

    near = np.zeros(len(lat), dtype=bool)
    for d_row in range(span):
        row = row_lo + d_row
        for d_col in range(span):
            col = col_lo + d_col
            valid = (row <= row_hi) & (col <= col_hi)
            ids = row * n_cols + np.mod(col, n_cols)
            near |= valid & np.isin(ids, test_blocks)
    return near
