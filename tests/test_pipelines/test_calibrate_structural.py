"""Structural calibration: point scores equal the raster function, and the pick rule."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipelines" / "training"))
import calibrate_structural as cs  # noqa: E402

from eii.compute.structural import landscape_quality  # noqa: E402


def test_point_scores_match_landscape_quality():
    rng = np.random.default_rng(0)
    res, lat = 0.0026949458523585, 48.0
    half_r, half_c = 105, 160
    hmi = np.clip(gaussian_filter(rng.random((2 * half_r + 1, 2 * half_c + 1)), 6) * 3 - 1, 0, 1)
    hmi[:10, :20] = np.nan
    top_lat = lat + half_r * res + res / 2
    nat = 1 - hmi
    row = cs.point_scores(hmi.astype("float32"), nat.astype("float32"), lat, res, (0.3, 0.4))
    for mid in (0.3, 0.4):
        ref = landscape_quality(hmi, top_lat, res, midpoint=mid)
        # The raster uses one kernel per 0.5 deg strip, the point its own latitude.
        assert abs(row[f"s_{mid}"] - ref[half_r, half_c]) < 2e-3
    assert 0 <= row["t"] <= 1 and 0 <= row["nmean"] <= 1


def test_criteria_pick_gentlest_passing_midpoint():
    n = 200
    frame = pd.DataFrame({
        "biome": np.repeat([1, 4], n // 2),
        "anthrome": np.tile([61, 21], n // 2),
        "t": np.tile([1.0, 0.4], n // 2),
        "nmean": np.tile([1.0, 0.4], n // 2),
        "s_0.3": np.tile([0.99, 0.30], n // 2),
        "s_0.4": np.tile([0.99, 0.44], n // 2),
        "s_0.5": np.tile([1.00, 0.60], n // 2),
    })  # fmt: skip
    table = cs.criteria(frame, (0.3, 0.4, 0.5))
    assert table["c1"].tolist() == [True, True, False]
    assert table["c2"].all()
    assert cs.pick(table) == 0.4
