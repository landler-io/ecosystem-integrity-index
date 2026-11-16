#!/usr/bin/env python3
"""FI change term inputs: per-sensor EVI ratio r on the 300 m EII grid.

The change term compares recent vegetation greenness with a fixed baseline. For each sensor,
annual mean EVI is computed from the 16-day composites (snow composites take the pixel's 5th
percentile of valid EVI), corrected for sensor drift with the stored per-biome and per-land-cover
factors, and divided by the median over the sensor's baseline years:
r = mean of the 3-year window ending in the reference year / baseline median. The numpy rules and
stored constants are ``eii.compute.change``; the readers are ``evi_sources.py``.

Reference years 2021-2025 are served by the TV series (MODIS Terra and VIIRS Suomi NPP), so
``--sensor T`` and ``--sensor V`` are all a run of v2.0 needs. The TNc and Nc series, with
NOAA-20 VIIRS chained to Suomi NPP by the median annual ratio over the overlap years (Nc), take
over once the Suomi NPP record (then MODIS Terra) no longer covers a reference-year window.
They are not wired up yet: ``--sensor VN`` (V plus NOAA-20) stores the chain factor ``f.tif``
and the Nc baseline median, but no stage writes ``r_Nc.vrt``. ``compute_eii_v2.py`` therefore
stops with an error for a reference year that ``series_for`` assigns to TNc or Nc (possible
only after a VN run); without N, a year Suomi NPP does not cover gets FI = FI_state.

Per sinusoidal tile and sensor (``sensor`` stage):

- read every 16-day composite of the sensor's record: T MOD13Q1 061 2001 onward from Planetary
  Computer, LP DAAC HDF4 for composites it lacks; V VNP13A1 002 2012 onward and N VJ113A1 002
  2018 onward from LP DAAC HDF5. N runs only together with V (``--sensor VN``, which reruns V
  of the tile): Nc needs both;
- annual means and snow p5 (``annual_evi``), drift correction with the biome and WorldCover
  class of the 300 m cell holding the pixel centre, baseline median, r per reference year (T, V);
- store per pixel, deflate GeoTIFFs on the sinusoidal grid: ``annual.tif`` (int16 EVI x 1e4,
  baseline start to the last year), ``p5.tif``, ``baseline_median.tif`` (int16 x 1e4), for N
  also ``f.tif`` (float32 chain factor) and the Nc baseline median;
- r per reference year, area-averaged onto the 300 m grid over the tile's lon/lat box, kept only
  in cells whose centre lies in the tile: ``r_300m.tif`` (float32, one band per reference year).

``tiles`` writes the land tiles of MODIS (``tiles_T.txt``) and VIIRS (``tiles_V.txt``, for
``V`` and ``VN``): those with a granule on one composite date in NASA CMR, inside the EII grid's
latitudes (60 S to 80 N).
``sensor`` skips finished tiles (``done.json``), retries a failed composite read and a failed
tile, and records a tile that still fails in ``<out>/failed/<sensor>_<tile>.json``.

``vrt`` mosaics the ``r_300m.tif`` of T and V into ``r_T.vrt`` and ``r_V.vrt`` and writes
``completeness.json`` (complete sensor years from the composite dates read, also for N).
``compute_eii_v2.py run --fi-change <out>`` combines them per reference year (``series_for``).

Layout: ``<out>/<sensor>/<tile>/``, scratch (composite codes, deleted after the tile) in
``<out>/scratch``. Credentials: see ``evi_sources.py``.

Example:
    python pipelines/production/fi_change_term.py tiles fi_change
    python pipelines/production/fi_change_term.py sensor fi_change --sensor T \\
        --tiles @fi_change/tiles_T.txt --ecoregions ECOREGIONS.tif
    python pipelines/production/fi_change_term.py sensor fi_change --sensor V \\
        --tiles @fi_change/tiles_V.txt --ecoregions ECOREGIONS.tif
    python pipelines/production/fi_change_term.py vrt fi_change
"""

from __future__ import annotations

import argparse
import json
import resource
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import evi_sources as es  # noqa: E402

from eii.compute import (  # noqa: E402
    annual_evi,
    baseline_median,
    chain_factor,
    chained,
    change_ratio,
    complete_years,
    correct_drift,
    encode_modis_evi,
    natural_class,
)
from eii.compute.settings import FI_CHANGE_BASELINES, REFERENCE_YEARS  # noqa: E402

LAST_YEAR = max(REFERENCE_YEARS)
GRID_PER_DEG = 336  # EII grid: 1/336 degree, cell centres on multiples of 1/336
GRID_ORIGIN = (-180 - 0.5 / GRID_PER_DEG, 80 + 0.5 / GRID_PER_DEG)  # west, north edges
SINU = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
RECORD = {"T": (2001, LAST_YEAR), "V": (2012, LAST_YEAR), "N": (2018, LAST_YEAR)}
"""Composite years read per sensor (the record behind snow p5)."""
STORE_FROM = {"T": FI_CHANGE_BASELINES["T"][0], "V": FI_CHANGE_BASELINES["V"][0], "N": 2018}
INT16_NODATA = -32768
EII_LATITUDES = (-60.0, 80.0)
TILE_DATE = "2025-01-01"  # CMR granules of this composite give the land tiles
RETRIES = 3  # attempts per composite read and per tile

_SHARED: dict = {}  # memmaps shared with forked workers


def log(msg: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def peak_gb() -> float:
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 2)


def cpu_s() -> float:
    """User + system CPU seconds of this process and its finished children (annual workers)."""
    return sum(r.ru_utime + r.ru_stime for r in (resource.getrusage(resource.RUSAGE_SELF),
                                                  resource.getrusage(resource.RUSAGE_CHILDREN)))  # fmt: skip


def retried(fn, what: str):  # noqa: ANN001, ANN201
    """``fn()``, retried ``RETRIES`` times with a growing pause (LP DAAC, PC, CMR errors)."""
    for attempt in range(RETRIES):
        try:
            return fn()
        except Exception as e:  # noqa: BLE001 - transient server errors; re-raised at the end
            if attempt == RETRIES - 1:
                raise
            log(f"retry {attempt + 1} {what}: {type(e).__name__}: {str(e)[:200]}")
            time.sleep(30 * (attempt + 1))


def parse_tile(tile: str) -> tuple[int, int]:
    return int(tile[1:3]), int(tile[4:6])


# --------------------------------------------------------------------------------------
# Grids
# --------------------------------------------------------------------------------------


def tile_box(h: int, v: int) -> tuple[int, int, int, int]:
    """(row0, col0, height, width) of the sinusoidal tile's lon/lat box on the 300 m grid."""
    lon0, lat0, lon1, lat1 = es.tile_bounds_lonlat(h, v)
    (x0, y0), g = GRID_ORIGIN, GRID_PER_DEG
    c0, r0 = int(np.floor((lon0 - x0) * g)), int(np.floor((y0 - lat1) * g))
    return r0, c0, int(np.ceil((y0 - lat0) * g)) - r0, int(np.ceil((lon1 - x0) * g)) - c0


def box_transform(h: int, v: int):  # noqa: ANN201
    from rasterio.transform import from_origin

    r0, c0, _, _ = tile_box(h, v)
    (x0, y0), g = GRID_ORIGIN, GRID_PER_DEG
    return from_origin(x0 + c0 / g, y0 - r0 / g, 1 / g, 1 / g)


def sinu_transform(h: int, v: int, n: int):  # noqa: ANN201
    from rasterio.transform import from_origin

    return from_origin(es.X0 + h * es.TILE_M, es.Y0 - v * es.TILE_M, es.TILE_M / n, es.TILE_M / n)


def box_groups(h: int, v: int, out: Path, ecoregions: Path) -> dict[str, np.ndarray]:
    """Biome, WorldCover natural class and in-tile mask of each 300 m cell of the tile box.

    Biome: RESOLVE ecoregion raster, nearest neighbour (``compute_eii_v2.tile_biome``);
    WorldCover: class at the cell centre from the 32x overview. Cached in ``<out>/grid``.
    """
    path = out / "grid" / f"h{h:02d}v{v:02d}.npz"
    if path.exists():
        with np.load(path) as z:
            return {k: z[k] for k in z.files}
    import compute_eii_v2 as CE
    import rasterio
    from rasterio.warp import Resampling, reproject

    _, _, hgt, w = tile_box(h, v)
    tr = box_transform(h, v)
    rows, cols = np.indices((hgt, w))
    lon = tr.c + (cols.ravel() + 0.5) * tr.a
    lat = tr.f + (rows.ravel() + 0.5) * tr.e
    th, tv = es.tile_of(lon, lat)
    inside = (th == h) & (tv == v)
    eco = np.zeros((hgt, w), "int64")
    with rasterio.open(ecoregions) as src:
        reproject(rasterio.band(src, 1), eco, dst_transform=tr, dst_crs="EPSG:4326",
                  resampling=Resampling.nearest)  # fmt: skip
    lookup = CE.eco_to_biome_lookup()
    known = (eco >= 0) & (eco < len(lookup))
    biome = np.where(known, lookup[np.where(known, eco, 0)], 0).ravel()
    wc = np.zeros(hgt * w, "int16")
    wc[inside] = es.worldcover_codes(lon[inside], lat[inside])
    groups = {"biome": np.where(biome > 0, biome, -1).astype("int16"),
              "cls": natural_class(wc).astype("int8"), "inside": inside}  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **groups)
    return groups


def pixel_groups(h: int, v: int, n: int, groups: dict) -> tuple[np.ndarray, np.ndarray]:
    """Biome and natural class of each native pixel: those of the 300 m cell at its centre."""
    size = es.TILE_M / n
    rr, cc = np.indices((n, n))
    x = es.X0 + h * es.TILE_M + (cc.ravel() + 0.5) * size
    y = es.Y0 - v * es.TILE_M - (rr.ravel() + 0.5) * size
    lat = np.degrees(y / es.EARTH_R)
    with np.errstate(all="ignore"):
        lon = np.degrees(x / (es.EARTH_R * np.cos(np.radians(lat))))
    tr = box_transform(h, v)
    _, _, hgt, w = tile_box(h, v)
    gr = np.clip(np.floor((lat - tr.f) / tr.e).astype(int), 0, hgt - 1)
    gc = np.clip(np.floor(np.nan_to_num(lon - tr.c) / tr.a).astype(int), 0, w - 1)
    k = gr * w + gc
    return groups["biome"][k].astype(int), groups["cls"][k].astype(int)


def to_box(arr: np.ndarray, h: int, v: int, n: int, inside: np.ndarray, threads: int) -> np.ndarray:
    """Area average of (bands, n, n) sinusoidal values onto the tile box; NaN outside the tile."""
    from rasterio.warp import Resampling, reproject

    _, _, hgt, w = tile_box(h, v)
    out = np.full((arr.shape[0], hgt, w), np.nan, "float32")
    for b in range(arr.shape[0]):
        reproject(np.asarray(arr[b], "float32"), out[b], src_transform=sinu_transform(h, v, n),
                  src_crs=SINU, dst_transform=box_transform(h, v), dst_crs="EPSG:4326",
                  resampling=Resampling.average, src_nodata=np.nan, dst_nodata=np.nan,
                  num_threads=threads)  # fmt: skip
    out[:, ~inside.reshape(hgt, w)] = np.nan
    return out


def write_tif(path: Path, arr: np.ndarray, nodata: float, transform, crs: str,
              descriptions: list[str] | None = None) -> float:  # fmt: skip
    """Deflate tiled GeoTIFF of (bands, h, w); returns its size in MB."""
    import rasterio

    float_ = arr.dtype.kind == "f"
    profile = {"driver": "GTiff", "tiled": True, "blockxsize": 512, "blockysize": 512,
               "compress": "DEFLATE", "predictor": 3 if float_ else 2, "zlevel": 6,
               "bigtiff": "IF_SAFER", "width": arr.shape[2], "height": arr.shape[1],
               "count": arr.shape[0], "dtype": arr.dtype.name, "nodata": nodata,
               "transform": transform, "crs": crs}  # fmt: skip
    part = path.with_name(path.name + ".part")
    with rasterio.open(part, "w", **profile) as dst:
        dst.write(arr)
        if descriptions:
            dst.descriptions = tuple(descriptions)
    part.rename(path)
    return path.stat().st_size / 1e6


def int16(x: np.ndarray) -> np.ndarray:
    """EVI (or median) x 1e4 as int16, ``INT16_NODATA`` where missing."""
    with np.errstate(invalid="ignore"):
        return np.where(np.isfinite(x), np.round(x * 1e4), INT16_NODATA).astype("int16")


# --------------------------------------------------------------------------------------
# Read
# --------------------------------------------------------------------------------------


def read_composites(sensor: str, h: int, v: int, folder: Path, threads: int) -> dict:
    """Every composite of the sensor's record over the tile -> ``codes.npy`` (int16) + dates."""
    import rasterio

    n = es.SENSORS[sensor][2]
    stats = {"bytes": 0.0, "granules": 0, "read_s": 0.0, "chunks": "", "lock": threading.Lock()}
    jobs: list[tuple[str, object]] = []
    pc_items = 0
    if sensor == "T":
        dates = es.modis_composite_dates(RECORD["T"])
        signer = es.Signer()
        items = es.wanted(es.stac_items(signer, h, v, RECORD["T"]), dates)
        pc_items = len(items)

        def pc(item: dict) -> np.ndarray:
            bands = []
            for key in ("250m_16_days_EVI", "250m_16_days_pixel_reliability"):
                url = signer.sign(item["assets"][key]["href"])
                with rasterio.Env(**es.GDAL_ENV), rasterio.open(url) as src:
                    bands.append(src.read(1))
                    bh, bw = src.block_shapes[0]
                    size = sum(src.block_size(1, i, j) for i in range(-(-src.height // bh))
                               for j in range(-(-src.width // bw)))  # fmt: skip
                with stats["lock"]:
                    stats["bytes"] += size
            return encode_modis_evi(*bands)

        have = {es.item_date(it["id"]).strftime("%Y_%m_%d"): it for it in items}
        jobs += [(day, lambda it=it: pc(it)) for day, it in have.items()]
        missing = [day for day in dates if day not in have]
        lp = es.granules("T", h, v, missing) if missing else []
        jobs += [(day, lambda u=u: es.read_hdf4(u, stats, folder / "tmp")) for day, u in lp]
    else:
        lp = es.granules(sensor, h, v, es.viirs_composite_dates(RECORD[sensor]))
        jobs += [(day, lambda u=u: es.read_viirs(u, stats)) for day, u in lp]
    jobs.sort(key=lambda j: j[0])
    codes = np.lib.format.open_memmap(folder / "codes.npy", "w+", "int16", (len(jobs), n, n))
    t0 = time.time()

    def run(k: int) -> None:
        codes[k] = retried(jobs[k][1], f"{sensor} h{h:02d}v{v:02d} {jobs[k][0]}")

    with ThreadPoolExecutor(threads) as pool:
        list(pool.map(run, range(len(jobs))))
    codes.flush()
    dates = [day for day, _ in jobs]
    (folder / "dates.json").write_text(json.dumps(dates))
    return {"composites": len(jobs), "pc_items": pc_items, "lp_granules": stats["granules"],
            "read_GB": round(stats["bytes"] / 1e9, 2), "read_s": round(time.time() - t0, 1)}  # fmt: skip


# --------------------------------------------------------------------------------------
# Annual means
# --------------------------------------------------------------------------------------


def _annual_strip(args: tuple[int, int]) -> tuple[float, float]:
    """Annual means and snow p5 of rows r0:r1 (forked worker): (seconds, worker peak RSS GB)."""
    t0 = time.time()
    r0, r1 = args
    codes, cy, years = _SHARED["codes"], _SHARED["composite_years"], _SHARED["years"]
    n = codes.shape[2]
    block = np.asarray(codes[:, r0:r1, :]).reshape(codes.shape[0], -1).T
    mean, p5 = annual_evi(block, cy, years)
    _SHARED["annual"][:, r0:r1, :] = mean.T.reshape(len(years), r1 - r0, n)
    _SHARED["p5"][r0:r1, :] = p5.reshape(r1 - r0, n)
    return time.time() - t0, peak_gb()


def annual_means(folder: Path, workers: int, strip: int) -> tuple[np.ndarray, np.ndarray, dict]:
    """(years, n, n) annual means and (n, n) p5 from ``codes.npy`` (fork pool over row strips)."""
    import multiprocessing as mp

    dates = json.loads((folder / "dates.json").read_text())
    codes = np.load(folder / "codes.npy", mmap_mode="r")
    n = codes.shape[1]
    years = list(range(int(dates[0][:4]), LAST_YEAR + 1))
    annual = np.lib.format.open_memmap(folder / "annual.npy", "w+", "float32", (len(years), n, n))
    p5 = np.lib.format.open_memmap(folder / "p5.npy", "w+", "float32", (n, n))
    _SHARED.update(codes=codes, composite_years=[int(x[:4]) for x in dates], years=years,
                   annual=annual, p5=p5)  # fmt: skip
    t0 = time.time()
    strips = [(r, min(r + strip, n)) for r in range(0, n, strip)]
    with ProcessPoolExecutor(workers, mp_context=mp.get_context("fork")) as pool:
        done = list(pool.map(_annual_strip, strips))
    annual.flush()
    p5.flush()
    timing = {"seconds": round(time.time() - t0, 1), "workers": workers, "strip_rows": strip,
              "process_hours": round(sum(s for s, _ in done) / 3600, 3),
              "peak_rss_worker_gb": max(m for _, m in done)}  # fmt: skip
    return np.load(folder / "annual.npy"), np.load(folder / "p5.npy"), timing | {"years": years}


# --------------------------------------------------------------------------------------
# One tile and sensor
# --------------------------------------------------------------------------------------


def score_sensor(sensor: str, x: np.ndarray, years: list[int], biome: np.ndarray,
                 cls: np.ndarray) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """Drift-corrected baseline median and r per reference year (pixels) of one sensor."""
    xc = correct_drift(x, sensor, years, biome, cls)
    median = baseline_median(xc, years, sensor)
    r = np.stack([change_ratio(xc, years, y, median) for y in REFERENCE_YEARS])
    return median, r.astype("float32")


def padded(x: np.ndarray, years: list[int], full: list[int]) -> np.ndarray:
    """(pixels, ``full``) series from (pixels, ``years``); NaN in the other years."""
    out = np.full((len(x), len(full)), np.nan, "float32")
    out[:, [full.index(y) for y in years]] = x
    return out


def run_tile(sensor: str, tile: str, out: Path, args: argparse.Namespace) -> dict:
    h, v = parse_tile(tile)
    groups = box_groups(h, v, out, args.ecoregions)
    sensors = ["V", "N"] if sensor == "VN" else [sensor]
    timing: dict = {"tile": tile}
    if (out / sensors[-1] / tile / "done.json").exists():  # VN: N needs V's float annual means
        return timing
    annual: dict[str, tuple[np.ndarray, list[int]]] = {}
    for s in sensors:
        folder = out / s / tile
        scratch = out / "scratch" / s / tile
        scratch.mkdir(parents=True, exist_ok=True)
        folder.mkdir(parents=True, exist_ok=True)
        t: dict = {"sensor": s}
        t0_wall, t0_cpu = time.time(), cpu_s()
        if not (scratch / "dates.json").exists():
            t["read"] = read_composites(s, h, v, scratch, args.threads) | {"peak_gb": peak_gb()}
            log(f"{tile} {s} read: {t['read']}")
        a, p5, t["annual"] = annual_means(scratch, args.workers, args.strip)
        log(f"{tile} {s} annual: {t['annual']}")
        dates = json.loads((scratch / "dates.json").read_text())
        n = a.shape[1]
        years = t["annual"].pop("years")
        t0 = time.time()
        biome, cls = pixel_groups(h, v, n, groups)
        x = a.reshape(len(years), -1).T
        annual[s] = (x, years)
        crs, tr = SINU, sinu_transform(h, v, n)
        stored = {"p5": write_tif(folder / "p5.tif", int16(p5)[None], INT16_NODATA, tr, crs)}
        j0 = years.index(STORE_FROM[s])
        stored["annual"] = write_tif(
            folder / "annual.tif", int16(a[j0:]), INT16_NODATA, tr, crs,
            [str(y) for y in years[j0:]])  # fmt: skip
        if s in ("T", "V"):
            median, r = score_sensor(s, x, years, biome, cls)
            stored["baseline_median"] = write_tif(
                folder / "baseline_median.tif", int16(median.reshape(1, n, n)), INT16_NODATA,
                tr, crs)  # fmt: skip
            r_box = to_box(r.reshape(-1, n, n), h, v, n, groups["inside"], args.threads)
            stored["r_300m"] = write_tif(
                folder / "r_300m.tif", r_box, np.nan, box_transform(h, v), "EPSG:4326",
                [str(y) for y in REFERENCE_YEARS])  # fmt: skip
            t["share_pixels_with_r"] = {y: float(np.isfinite(r[j]).mean())
                                        for j, y in enumerate(REFERENCE_YEARS)}  # fmt: skip
        else:  # N: chain factor and the Nc baseline (Nc = V before 2021, N x f from 2021)
            full = list(range(RECORD["V"][0], LAST_YEAR + 1))
            vx, nx = padded(*annual["V"], full), padded(x, years, full)
            f = chain_factor(vx, nx, full)
            median, _ = score_sensor("Nc", chained(vx, nx, full, f), full, biome, cls)
            stored["f"] = write_tif(folder / "f.tif", f.reshape(1, n, n).astype("float32"),
                                    np.nan, tr, crs)  # fmt: skip
            stored["baseline_median_Nc"] = write_tif(
                folder / "baseline_median_Nc.tif", int16(median.reshape(1, n, n)), INT16_NODATA,
                tr, crs)  # fmt: skip
            t["share_pixels_with_f"] = float(np.isfinite(f).mean())
        t["score_store_s"] = round(time.time() - t0, 1)
        t["stored_MB"] = {k: round(val, 1) for k, val in stored.items()}
        t["peak_gb"] = peak_gb()
        t["wall_s"] = round(time.time() - t0_wall, 1)
        t["cpu_h"] = round((cpu_s() - t0_cpu) / 3600, 3)
        (folder / "dates.json").write_text(json.dumps(dates))
        (folder / "done.json").write_text(json.dumps(t, indent=1, default=float) + "\n")
        if not args.keep_scratch:
            shutil.rmtree(scratch)
        log(f"{tile} {s} done: {json.dumps(t, default=float)}")
        timing[s] = t
    return timing


def land_tiles(sensor: str) -> list[str]:
    """Tiles with a granule on ``TILE_DATE`` in CMR whose latitudes overlap the EII grid."""
    import requests

    short, version, _ = es.SENSORS[sensor]
    params = {"short_name": short, "version": version, "page_size": 2000,
              "temporal": f"{TILE_DATE}T00:00:00Z,{TILE_DATE}T23:59:59Z"}  # fmt: skip

    def query() -> list[dict]:
        r = requests.get(es.CMR, params=params, timeout=120)
        r.raise_for_status()
        entries: list[dict] = r.json()["feed"]["entry"]
        return entries

    tiles = set()
    for e in retried(query, "CMR tiles"):
        h, v = (int(x) for x in e["title"].split(".")[2][1:].split("v"))
        _, lat0, _, lat1 = es.tile_bounds_lonlat(h, v)
        if lat1 > EII_LATITUDES[0] + 1e-6 and lat0 < EII_LATITUDES[1] - 1e-6:
            tiles.add(f"h{h:02d}v{v:02d}")
    return sorted(tiles)


def stage_tiles(out: Path) -> None:
    """``tiles_T.txt`` (MOD13Q1) and ``tiles_V.txt`` (VNP13A1) in ``out``."""
    for sensor in ("T", "V"):
        tiles = land_tiles(sensor)
        (out / f"tiles_{sensor}.txt").write_text("\n".join(tiles) + "\n")
        log(f"{sensor}: {len(tiles)} tiles")


def run_tile_retried(sensor: str, tile: str, out: Path, args: argparse.Namespace) -> None:
    """``run_tile`` with ``RETRIES`` attempts; a tile that still fails is recorded, not raised."""
    failed = out / "failed" / f"{sensor}_{tile}.json"
    for attempt in range(RETRIES):
        try:
            timing = run_tile(sensor, tile, out, args)
            failed.unlink(missing_ok=True)
            log(f"{tile} {sensor}: {json.dumps(timing, default=float)[:400]}")
            return
        except Exception as e:  # noqa: BLE001 - one tile must not stop the run
            log(f"{tile} {sensor} attempt {attempt + 1} failed: {type(e).__name__}: {e}")
            failed.parent.mkdir(parents=True, exist_ok=True)
            failed.write_text(json.dumps({"tile": tile, "sensor": sensor, "attempts":
                                          attempt + 1, "error": f"{type(e).__name__}: {e}"}) + "\n")  # fmt: skip
            time.sleep(60)


# --------------------------------------------------------------------------------------
# Mosaic and completeness
# --------------------------------------------------------------------------------------


def stage_vrt(out: Path) -> None:
    """``r_<sensor>.vrt`` per sensor and ``completeness.json`` (complete years per sensor)."""
    record: dict = {}
    for s in ("T", "V", "N"):
        tiles = sorted(p.parent for p in (out / s).glob("*/done.json"))
        if not tiles:
            continue
        dates: dict[str, set[str]] = {}
        for t in tiles:
            for day in json.loads((t / "dates.json").read_text()):
                dates.setdefault(day[:4], set()).add(day)
        record[s] = {"tiles": len(tiles),
                     "complete_years": sorted(complete_years([x for v in dates.values()
                                                              for x in v])),
                     "composites_per_year": {y: len(v) for y, v in sorted(dates.items())}}  # fmt: skip
        if s in ("T", "V"):
            sources = [str(t / "r_300m.tif") for t in tiles]
            listing = out / f"r_{s}.txt"
            listing.write_text("\n".join(sources) + "\n")
            subprocess.run(["gdalbuildvrt", "-q", "-input_file_list", str(listing),
                            str(out / f"r_{s}.vrt")], check=True)  # fmt: skip
    (out / "completeness.json").write_text(json.dumps(record, indent=1) + "\n")
    log(json.dumps({s: v["complete_years"] for s, v in record.items()}))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("stage", choices=["tiles", "sensor", "vrt"])
    parser.add_argument("out", type=Path, help="work and output folder")
    parser.add_argument("--sensor", choices=["T", "V", "VN"],
                        help="T MODIS Terra, V VIIRS Suomi NPP, VN V plus NOAA-20 (Nc inputs)")  # fmt: skip
    parser.add_argument("--tiles", nargs="*", default=[],
                        help="sinusoidal tiles, e.g. h18v04, or @file (one per line)")  # fmt: skip
    parser.add_argument("--ecoregions", type=Path,
                        help="RESOLVE ecoregion raster (ECO_ID), as for compute_eii_v2.py")  # fmt: skip
    parser.add_argument("--threads", type=int, default=16, help="read threads")
    parser.add_argument("--workers", type=int, default=16, help="annual-mean processes")
    parser.add_argument("--strip", type=int, default=10, help="rows per annual-mean task")
    parser.add_argument("--keep-scratch", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.stage == "vrt":
        stage_vrt(args.out)
        return
    if args.stage == "tiles":
        stage_tiles(args.out)
        return
    if not args.sensor or not args.tiles or args.ecoregions is None:
        parser.error("sensor needs --sensor, --tiles and --ecoregions")
    tiles = [t for x in args.tiles for t in (Path(x[1:]).read_text().split() if x.startswith("@")
                                             else [x])]  # fmt: skip
    for tile in tiles:
        run_tile_retried(args.sensor, tile, args.out, args)


if __name__ == "__main__":
    main()
