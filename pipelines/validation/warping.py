"""Parallel blockwise warping of global rasters onto the validation grids.

Warping goes through ``rasterio.warp.reproject``, so no GDAL command-line tools are needed:
the destination grid is split into blocks, each worker process reads only the source window
under its block (computed from a lattice of transformed points, so blocks at the edge of the
Eckert IV domain cannot trigger a whole-globe read), warps it (GDAL resampling, nodata
excluded from averages) and the parent writes it into a tiled, compressed float32 GeoTIFF
with NaN nodata. Blocks whose source window would exceed ``MAX_SRC_VALUES`` are split in
four.

Several bands of one pixel-interleaved source are warped in a single read, so our
nine-band EII tiles are decompressed once per run, not once per band.
"""

from __future__ import annotations

import math
import os
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pyproj
import rasterio
import rasterio.errors
import rasterio.warp
import rasterio.windows
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.warp import reproject
from rasterio.windows import Window

OUT_PROFILE = {
    "driver": "GTiff",
    "dtype": "float32",
    "nodata": np.nan,
    "tiled": True,
    "blockxsize": 512,
    "blockysize": 512,
    "compress": "zstd",
    "predictor": 3,
    "bigtiff": "if_safer",
}


@dataclass(frozen=True)
class WarpJob:
    """One source raster warped onto one grid.

    Attributes:
        src: Source raster (GeoTIFF or VRT).
        bands: Source band indexes to warp.
        outs: One output path per band (``None`` skips writing that band).
        resampling: GDAL resampling (``average`` for continuous, ``nearest``/``mode``
            for categorical layers).
        src_crs: CRS override for sources without a tagged CRS.
        src_nodata: Nodata override (also used when the source declares none).
        scale: Multiplier applied after warping (e.g. 0.01 for HII stored x100).
        valid: Inclusive (min, max) range kept after scaling; outside becomes NaN.
    """

    src: Path
    bands: tuple[int, ...]
    outs: tuple[Path | None, ...]
    resampling: Resampling = Resampling.average
    src_crs: str | None = None
    src_nodata: float | None = None
    scale: float = 1.0
    valid: tuple[float, float] | None = None
    tags: dict[str, str] = field(default_factory=dict)

    def done(self) -> bool:
        return all(o is None or o.exists() for o in self.outs)


MAX_SRC_VALUES = 40_000_000  # per block read (bands x rows x cols), 160 MB as float32


def _src_window(
    src: rasterio.DatasetReader, src_crs: CRS, grid: dict[str, Any], window: Window, n: int = 33
) -> Window | None:
    """Source window covering a destination block, or ``None`` if it is off the globe.

    The block is sampled on an ``n`` x ``n`` lattice and transformed to the source CRS;
    points outside the projection domain (non-finite) are dropped. The bounds are padded
    by 5 % plus two source pixels, so GDAL never falls back to reading the whole source.
    """
    t = rasterio.windows.transform(window, grid["transform"])
    cc, rr = np.meshgrid(np.linspace(0, window.width, n), np.linspace(0, window.height, n))
    xs, ys = t * (cc.ravel(), rr.ravel())
    sx, sy = _to_src(grid["crs"], src_crs, xs, ys)
    ok = np.isfinite(sx) & np.isfinite(sy) & (np.abs(sx) < 1e15) & (np.abs(sy) < 1e15)
    if src_crs.is_geographic:
        ok &= (np.abs(sx) <= 180.5) & (np.abs(sy) <= 90.5)
    if not ok.any():
        return None
    sx, sy = sx[ok], sy[ok]
    padx = 0.05 * (sx.max() - sx.min()) + 2 * abs(src.transform.a)
    pady = 0.05 * (sy.max() - sy.min()) + 2 * abs(src.transform.e)
    win = rasterio.windows.from_bounds(
        sx.min() - padx, sy.min() - pady, sx.max() + padx, sy.max() + pady, src.transform
    )
    c0, r0 = math.floor(win.col_off), math.floor(win.row_off)
    win = Window(
        c0, r0, math.ceil(win.col_off + win.width) - c0, math.ceil(win.row_off + win.height) - r0
    )
    try:
        win = win.intersection(Window(0, 0, src.width, src.height))
    except rasterio.errors.WindowError:
        return None
    return win if win.width > 0 and win.height > 0 else None


def _to_src(dst_crs: CRS, src_crs: CRS, xs: Any, ys: Any) -> tuple[np.ndarray, np.ndarray]:
    """Transform points with pyproj; points outside the domain come back non-finite.

    ``rasterio.warp.transform`` raises on the first such point instead.
    """
    tr = pyproj.Transformer.from_crs(dst_crs.to_wkt(), src_crs.to_wkt(), always_xy=True)
    with np.errstate(all="ignore"):
        sx, sy = tr.transform(np.asarray(xs, "float64"), np.asarray(ys, "float64"), errcheck=False)
    return np.asarray(sx), np.asarray(sy)


EDGE_BLOCK = 32  # px; edge blocks at most this wide use nearest-cell mapping


def _split(job: WarpJob, grid: dict[str, Any], window: Window, out: np.ndarray) -> np.ndarray:
    """Warp ``window`` as four quadrants into ``out``."""
    half_w, half_h = int(window.width) // 2, int(window.height) // 2
    for c0, r0, w, h in (
        (0, 0, half_w, half_h),
        (half_w, 0, int(window.width) - half_w, half_h),
        (0, half_h, half_w, int(window.height) - half_h),
        (half_w, half_h, int(window.width) - half_w, int(window.height) - half_h),
    ):
        sub = Window(window.col_off + c0, window.row_off + r0, w, h)
        _, data = _warp_block(job, grid, sub)
        out[:, r0 : r0 + h, c0 : c0 + w] = data
    return out


def _nearest_cells(
    data: np.ndarray, src_transform: Any, src_crs: CRS, grid: dict[str, Any], window: Window
) -> np.ndarray:
    """Nearest source cell under each destination cell centre (off-domain cells NaN)."""
    out = np.full((data.shape[0], int(window.height), int(window.width)), np.nan, "float32")
    t = rasterio.windows.transform(window, grid["transform"])
    cc, rr = np.meshgrid(np.arange(window.width) + 0.5, np.arange(window.height) + 0.5)
    xs, ys = t * (cc.ravel(), rr.ravel())
    sx, sy = _to_src(grid["crs"], src_crs, xs, ys)
    ok = np.isfinite(sx) & np.isfinite(sy)
    col = np.full(sx.shape, -1, dtype="int64")
    row = np.full(sy.shape, -1, dtype="int64")
    inv = ~src_transform
    c, r = inv * (sx[ok], sy[ok])
    col[ok], row[ok] = np.floor(c).astype("int64"), np.floor(r).astype("int64")
    ok &= (col >= 0) & (col < data.shape[2]) & (row >= 0) & (row < data.shape[1])
    flat = out.reshape(data.shape[0], -1)
    flat[:, ok] = data[:, row[ok], col[ok]]
    return out


def _warp_block(job: WarpJob, grid: dict[str, Any], window: Window) -> tuple[Window, np.ndarray]:
    """Warp one destination block from an explicitly bounded source read."""
    out = np.full((len(job.bands), int(window.height), int(window.width)), np.nan, "float32")
    with rasterio.open(job.src) as src:
        src_crs = CRS.from_string(job.src_crs) if job.src_crs else src.crs
        win = _src_window(src, src_crs, grid, window)
        if win is None:
            return window, out
        if len(job.bands) * win.width * win.height > MAX_SRC_VALUES and window.width > 64:
            return window, _split(job, grid, window, out)
        raw = src.read(list(job.bands), window=win)
        nodata = job.src_nodata if job.src_nodata is not None else src.nodata
        data = raw.astype("float32")
        if nodata is not None:
            data[raw == nodata] = np.nan
        if np.issubdtype(raw.dtype, np.floating):
            data[~np.isfinite(raw)] = np.nan
        del raw
        src_transform = src.window_transform(win)
        try:
            reproject(
                data,
                out,
                src_transform=src_transform,
                src_crs=src_crs,
                src_nodata=np.nan,
                dst_transform=rasterio.windows.transform(window, grid["transform"]),
                dst_crs=grid["crs"],
                dst_nodata=np.nan,
                resampling=job.resampling,
                num_threads=1,
            )
        except Exception:  # GDAL raises CPLE_AppDefinedError (not a RasterioError)
            # The block straddles the edge of the projection domain: split it until the
            # pieces are inside (normal warp) or tiny (nearest-cell mapping below).
            if window.width > EDGE_BLOCK:
                del data
                return window, _split(job, grid, window, out)
            out = _nearest_cells(data, src_transform, src_crs, grid, window)
    if job.scale != 1.0:
        out *= np.float32(job.scale)
    if job.valid is not None:
        out[(out < job.valid[0]) | (out > job.valid[1])] = np.nan
    return window, out


def blocks(width: int, height: int, size: int = 1024) -> list[Window]:
    """Destination windows covering a ``width`` x ``height`` grid."""
    return [
        Window(c, r, min(size, width - c), min(size, height - r))
        for r in range(0, height, size)
        for c in range(0, width, size)
    ]


def run_warp(
    job: WarpJob,
    grid: dict[str, Any],
    workers: int = 8,
    block: int = 1024,
    log: Callable[[str], None] = print,
) -> list[Path]:
    """Warp ``job`` onto ``grid`` (keys crs, transform, width, height); skip if done."""
    if job.done():
        return [o for o in job.outs if o is not None]
    profile = OUT_PROFILE | {
        "crs": grid["crs"],
        "transform": grid["transform"],
        "width": grid["width"],
        "height": grid["height"],
        "count": 1,
    }
    tmp = {o: o.with_suffix(".partial.tif") for o in job.outs if o is not None}
    for path in tmp.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    dsts = {o: rasterio.open(p, "w", **profile) for o, p in tmp.items()}
    windows = blocks(grid["width"], grid["height"], block)
    env = {"GDAL_CACHEMAX": "256", "GDAL_NUM_THREADS": "1"}
    os.environ.update(env)
    try:
        with ProcessPoolExecutor(workers) as pool:
            pending = {pool.submit(_warp_block, job, grid, w) for w in windows}
            for i, fut in enumerate(as_completed(pending), 1):
                window, data = fut.result()
                pending.discard(fut)  # release the block once written
                for k, out in enumerate(job.outs):
                    if out is not None:
                        dsts[out].write(data[k], 1, window=window)
                if i % max(1, len(windows) // 10) == 0:
                    log(f"  {job.src.name}: {i}/{len(windows)} blocks")
    finally:
        for dst in dsts.values():
            dst.update_tags(**job.tags, source=str(job.src))
            dst.close()
    for out, path in tmp.items():
        path.rename(out)
    return list(tmp)


def read_full(path: Path) -> np.ndarray:
    """Whole band 1 as float32 with NaN nodata."""
    with rasterio.open(path) as src:
        data = src.read(1).astype("float32", copy=False)
        if src.nodata is not None and np.isfinite(src.nodata):
            data[data == src.nodata] = np.nan
    return data


def sample_cells(path: Path, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    """Values of band 1 at integer (row, col) cells of a grid-aligned raster."""
    return read_full(path)[rows, cols]


def grid_of(path: Path) -> dict[str, Any]:
    """Grid dict (crs, transform, width, height) of an existing raster."""
    with rasterio.open(path) as src:
        return {
            "crs": src.crs,
            "transform": src.transform,
            "width": src.width,
            "height": src.height,
        }


def cos_lat_weights(grid: dict[str, Any]) -> np.ndarray:
    """Per-row cell-area weights (cos latitude of row centres) of an EPSG:4326 grid."""
    t = grid["transform"]
    lat = t.f + t.e * (np.arange(grid["height"]) + 0.5)
    return np.cos(np.radians(lat)).astype("float64")


__all__: Sequence[str] = [
    "WarpJob",
    "blocks",
    "cos_lat_weights",
    "grid_of",
    "read_full",
    "run_warp",
    "sample_cells",
]
