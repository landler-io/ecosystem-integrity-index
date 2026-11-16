"""Offline tests for the Earth Engine candidate-pool sampler: join, masks, chunk cache, merge."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipelines" / "training"))
import sample_candidate_pool_ee as pool  # noqa: E402

PREDICTORS = dict.fromkeys(pool.PREDICTORS, 1.0)
ANNUAL = dict.fromkeys(pool.ANNUAL_COLUMNS, 7)


def _image(rows: list[dict]) -> MagicMock:
    image = MagicMock()
    sampled = image.unmask.return_value.sampleRegions.return_value
    sampled.getInfo.return_value = {"features": [{"properties": r} for r in rows]}
    image.sampleRegions.return_value = sampled
    return image


def _points(n: int) -> pd.DataFrame:
    frame = pd.DataFrame({c: np.arange(n, dtype=float) for c in pool.POINT_COLUMNS})
    frame["tile"] = "e010n40"
    return frame.rename_axis("point_id")


def test_chunk_joins_on_point_id_and_masks_predictors(monkeypatch):
    stack = _image(
        [
            {"point_id": 11, **PREDICTORS, "sand": pool.MASKED},
            {"point_id": 10, **PREDICTORS},
        ]
    )
    annual = _image([{"point_id": 10, **ANNUAL}])
    monkeypatch.setattr(pool, "predictor_stack", lambda _ee: stack)
    monkeypatch.setattr(pool, "annual_npp", lambda _ee: (annual, "projection"))
    chunk = _points(3).set_axis([10, 11, 12]).rename_axis("point_id")

    out = pool.sample_chunk(MagicMock(), chunk)

    assert out.index.tolist() == [10, 11, 12]
    assert list(out.columns) == [*pool.PREDICTORS, *pool.ANNUAL_COLUMNS]
    assert out.loc[10, "sand"] == 1.0 and np.isnan(out.loc[11, "sand"])
    assert out.loc[10, "sum_2014"] == 7 and np.isnan(out.loc[11, "sum_2014"])
    assert out.loc[12].isna().all()
    stack.unmask.assert_called_once_with(pool.MASKED)


def test_extract_resumes_from_cache_and_merge_orders_columns(tmp_path, monkeypatch):
    points = _points(5)
    calls = []

    def fake_chunk(_ee, chunk):
        calls.append(chunk.index.tolist())
        values = pd.DataFrame({**PREDICTORS, **ANNUAL}, index=chunk.index)
        values.loc[chunk.index[0], "sum_2020"] = np.nan
        return values

    monkeypatch.setattr(pool, "sample_chunk", fake_chunk)
    pool.extract(None, points, tmp_path, chunk_size=2, workers=1, limit=1)
    assert calls == [[0, 1]]
    assert pool.merge(points, tmp_path) is None

    pool.extract(None, points, tmp_path, chunk_size=2, workers=1, limit=None)
    assert calls == [[0, 1], [2, 3], [4]]
    merged = pool.merge(points, tmp_path)

    assert list(merged.columns) == [
        *pool.POINT_COLUMNS,
        *pool.PREDICTORS,
        *pool.ANNUAL_COLUMNS,
        "tile",
    ]
    assert merged.index.tolist() == [0, 1, 2, 3, 4]
    assert merged["sum_2020"].dtype == "int64"
    assert merged.loc[[0, 2, 4], "sum_2020"].eq(pool.MASKED).all()
    assert merged.loc[1, "sum_2020"] == 7
