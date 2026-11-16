#!/usr/bin/env python3
"""Build the near-pristine control table of the FI change term (``fi_change_controls.csv.gz``).

``derive_fi_change_constants.py`` derives the stored change-term constants from this table.

Stages (``<work>`` holds every intermediate and the result):

- ``select``: equal-area random draws (seed ``SEED``, ``N_DRAWS`` draws between 60 S and 80 N)
  snapped to the 300 m EII grid (1/336 deg; one candidate per cell). A cell is a candidate if
  it has a RESOLVE biome, GHM 2022 HMI < ``HMI_MAX`` at its centre and no exclusion point within
  ``EXCLUDE_KM``; the first ``N_CONTROLS // 14 * OVERSAMPLE`` candidates per biome (draw order)
  are kept and get a random ``draw_key``. A candidate is pure if ESA WorldCover 2021 tree,
  shrubland or grassland covers at least ``PURE`` of the cell (10 m pixels, majority class).
  The controls are the pure candidates with the smallest draw keys per biome, ``N_CONTROLS`` in
  total, split equally over the biomes (a biome with too few gives its share to the others),
  less those with HMI >= ``HMI_MAX`` at their stored position (see below)
  -> ``controls.csv.gz`` (``point_id`` in draw-key order, cell, centre, biome, WorldCover).
- ``extract --sensor T|V|N``: the 16-day EVI composites of the sensor's record at the pixel
  holding each control, encoded as in production (``fi_change_term.py``), per sinusoidal tile
  -> ``evi/<sensor>/<tile>.csv.gz``. Finished tiles are skipped.
- ``assemble``: annual mean EVI per sensor (``eii.compute.annual_evi``; snow takes the pixel's
  5th percentile over the record) on the composite dates of all tiles together: a date read in
  some tile but not in another counts there as an invalid composite, which short-gap filling
  may interpolate. Then the 15 potential-NPP predictors at the cell centres and ``fi_state``
  per reference year from a ``compute_eii_v2.py run`` -> ``fi_change_controls.csv.gz``, the
  layout ``derive_fi_change_constants.py`` reads.

Sources:

- RESOLVE Ecoregions 2017 raster (``ECO_ID``; ``datasets/ecoregions/v1_0`` synced by
  ``sync_predictor_sources.py``, requester pays: ``EII_BILLING_PROJECT``), mapped to biomes by
  ``pipelines/preprocessing/reference/eco_id_to_biome.json``;
- GHM 2022 (``HMv20240801_2022s_AA_300.tif``, Zenodo 10.5281/zenodo.14502573), downloaded to
  ``--hmi`` if missing;
- ESA WorldCover 2021 v200 (public AWS bucket, no credentials);
- MOD13Q1 061 (Planetary Computer; LP DAAC for ``T_FILL``), VNP13A1 002 and VJ113A1 002 (LP
  DAAC): the readers and credentials of ``evi_sources.py`` (Planetary Computer key, Earthdata
  token in ``~/.edl_token``, ``gdal_translate`` with the HDF4 driver);
- the predictor stack of ``build_predictor_stack.py`` (``stack.vrt``);
- the tiles of a ``compute_eii_v2.py run`` (band ``fi_state``).

``fi_state`` is the functional state before the change term (observed against potential NPP).
The run computes it independently of the change-term constants, and the selection uses no FI,
so a run made with the stored constants is a valid source. Only k_c depends on it.

Reproducing the published table:

- The published selection stored candidate positions and draw keys with 7 significant digits
  and read WorldCover and the final HMI check there; it read the Terra composites at those
  positions rounded again to 6 significant digits (up to ~0.001 deg off the cell centre where
  |lon| >= 100), and the VIIRS composites at the cell centres with 7 decimals. ``stored``
  repeats this: many centres lie on a 500 m VIIRS pixel edge, where the last digit decides.
- Cells within ``EXCLUDE_KM`` of earlier internal evaluation points were left out
  (``--exclude``, CSV with ``lon`` and ``lat``). That list is not published; without it the
  candidate pool, and the controls, differ.
- The Terra series holds the composites Planetary Computer had in November 2026 plus
  ``T_FILL`` (Earth Engine copy, ``encode_ee_modis``); it lacks 2023-02-18 and 2024-08-12
  everywhere, 2024-12-18 and 2025-04-07 in most or some tiles. Planetary Computer items added
  since then change the result.
- The published ``fi_state`` comes from an earlier run on the potential NPP before gap
  filling (same FI_state rule); 18 % of the controls differ from the released run by more
  than 0.001. With the released run's ``fi_state``, ``derive_fi_change_constants.py`` gives the
  same drift factors, h models and k_c; only the p10 values in ``k_c.json`` change.
- The published ``point_id`` values are row ids of a larger internal table; here they are 0 to
  n - 1. Rows follow the draw key; 11 pairs share a stored key and are in arbitrary order there,
  by ``cell_id`` here. LightGBM bagging depends on row order, so only the published row order
  gives the stored h models bit for bit (drift factors are unaffected, k_c stays 1.0).

Example:
    python pipelines/training/build_fi_change_controls.py select work \\
        --ecoregions ECOREGIONS.tif --hmi HMv20240801_2022s_AA_300.tif --exclude exclusion.csv
    python pipelines/training/build_fi_change_controls.py extract work --sensor T
    python pipelines/training/build_fi_change_controls.py extract work --sensor V
    python pipelines/training/build_fi_change_controls.py extract work --sensor N
    python pipelines/training/build_fi_change_controls.py assemble work \\
        --run eii_v2 --predictor-stack predictor_stack/stack.vrt
"""

from __future__ import annotations

import argparse
import logging
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent / "production"), str(HERE.parent / "preprocessing")]
import fi_change_term as fct  # noqa: E402
from derive_fi_change_constants import YEARS  # noqa: E402
from local_sources import GHM_2022, ZENODO_URL, download  # noqa: E402

from eii.compute import annual_evi, encode_modis_evi  # noqa: E402
from eii.compute.settings import BIOME_NAMES, EVI_INVALID_CODE, REFERENCE_YEARS  # noqa: E402

es = fct.es
SEED = 55
N_DRAWS = 8_000_000
N_CONTROLS = 20_000
OVERSAMPLE = 3
HMI_MAX = 0.05
PURE = 0.75
NATURAL_WORLDCOVER = (10, 20, 30)
EXCLUDE_KM = 5.0
EARTH_KM = 6371.0
RES = 1 / 336
NCOLS = 120960
SENSORS = ("T", "V", "N")
T_FILL = tuple(f"2025_{d}" for d in ("07_12", "07_28", "08_13", "08_29", "09_14", "09_30",
                                     "10_16", "11_01"))  # fmt: skip
"""MOD13Q1 composites missing on Planetary Computer that the published series took from the Earth
Engine copy (here LP DAAC, ``encode_ee_modis``). ``fi_change_term.py`` fills every gap."""

log = logging.getLogger("build_fi_change_controls")


# --------------------------------------------------------------------------------------
# Sampling
# --------------------------------------------------------------------------------------


def stored(x: np.ndarray, spec: str) -> np.ndarray:
    """``x`` as written to text with the format ``spec`` (e.g. ``.7g``) and read back."""
    return np.array([float(format(v, spec)) for v in np.asarray(x, "float64")])


def sample(path: Path | str, lon: np.ndarray, lat: np.ndarray, bands: list[int],
           threads: int = 16) -> np.ndarray:  # fmt: skip
    """Bands at the points (pixel holding the point, nodata NaN), threads over latitude chunks."""
    import rasterio

    lon, lat = np.asarray(lon, "float64"), np.asarray(lat, "float64")
    out = np.full((len(lon), len(bands)), np.nan)

    def read(idx: np.ndarray) -> np.ndarray:
        with rasterio.open(path) as src:
            v = np.array(list(src.sample(zip(lon[idx], lat[idx], strict=True), indexes=bands)),
                         "float64").reshape(len(idx), len(bands))  # fmt: skip
            if src.nodata is not None:
                v[v == src.nodata] = np.nan
        return v

    chunks = [c for c in np.array_split(np.argsort(-lat, kind="stable"), threads) if len(c)]
    with ThreadPoolExecutor(max(len(chunks), 1)) as pool:
        for idx, v in zip(chunks, pool.map(read, chunks), strict=True):
            out[idx] = v
    return out


def sample_tiles(folder: Path, lon: np.ndarray, lat: np.ndarray, band: str) -> np.ndarray:
    """Band ``band`` (by description) of the tiles ``folder/*.tif`` at the points."""
    import rasterio

    out = np.full(len(lon), np.nan)
    for path in sorted(folder.glob("*.tif")):
        with rasterio.open(path) as src:
            b, index = src.bounds, src.descriptions.index(band) + 1
        idx = np.flatnonzero((lon >= b.left) & (lon < b.right) & (lat > b.bottom) & (lat <= b.top))
        if len(idx):
            out[idx] = sample(path, lon[idx], lat[idx], [index], threads=1)[:, 0]
    return out


# --------------------------------------------------------------------------------------
# Selection
# --------------------------------------------------------------------------------------


def draw_cells(n: int, rng: np.random.Generator) -> pd.DataFrame:
    """Grid cells under ``n`` equal-area random points between 60 S and 80 N, in draw order."""
    lat = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    lon = rng.uniform(-180, 180, n)
    keep = (lat < 80 - RES) & (lat > -60)
    col = np.rint((lon[keep] + 180.0) / RES).astype("int64") % NCOLS
    row = np.rint((80.0 - lat[keep]) / RES).astype("int64")
    cells = pd.DataFrame({"cell_id": row * NCOLS + col, "lon": -180.0 + col * RES,
                          "lat": 80.0 - row * RES})  # fmt: skip
    return cells.drop_duplicates("cell_id").reset_index(drop=True)


def within_km(lon: np.ndarray, lat: np.ndarray, ref: pd.DataFrame, km: float) -> np.ndarray:
    """True where a point lies within ``km`` (great circle) of a point of ``ref`` (lon, lat)."""
    from scipy.spatial import cKDTree

    def xyz(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        x, y = np.radians(np.asarray(x, "float64")), np.radians(np.asarray(y, "float64"))
        return np.c_[np.cos(y) * np.cos(x), np.cos(y) * np.sin(x), np.sin(y)]

    if ref.empty:
        return np.zeros(len(lon), bool)
    chord = 2 * np.sin(km / EARTH_KM / 2)
    dist, _ = cKDTree(xyz(ref["lon"], ref["lat"])).query(xyz(lon, lat), distance_upper_bound=chord)
    return np.isfinite(dist)


def worldcover_cells(tile: str, cells: pd.DataFrame) -> pd.DataFrame:
    """Majority WorldCover class and its share of the 10 m pixels in each cell of one 3 deg tile."""
    import rasterio
    from rasterio.windows import from_bounds

    rows = []
    try:
        src = rasterio.open(es.WORLDCOVER.format(tile=tile))
    except rasterio.RasterioIOError:  # no tile (ocean)
        return pd.DataFrame({"cell_id": cells["cell_id"], "worldcover": 0, "wc_share": np.nan})
    with src:
        for cid, x, y in cells[["cell_id", "lon", "lat"]].itertuples(index=False):
            win = from_bounds(x - RES / 2, y - RES / 2, x + RES / 2, y + RES / 2, src.transform)
            v = src.read(1, window=win.round_offsets().round_lengths(), boundless=True)
            counts = np.bincount(v.ravel(), minlength=101)
            code = int(counts.argmax())
            rows.append((cid, code, counts[code] / max(v.size, 1)))
    return pd.DataFrame(rows, columns=["cell_id", "worldcover", "wc_share"])


def worldcover_shares(cells: pd.DataFrame, cache: Path, threads: int) -> pd.DataFrame:
    """``worldcover_cells`` for every cell, cached in ``cache`` (only new cells are read)."""
    import rasterio

    done = pd.read_csv(cache) if cache.exists() else pd.DataFrame(columns=["cell_id"])
    todo = cells[~cells["cell_id"].isin(done["cell_id"])]
    if len(todo):
        lat3 = (np.floor(todo["lat"] / 3) * 3).astype(int)
        lon3 = (np.floor(todo["lon"] / 3) * 3).astype(int)
        tile = [f"{'N' if a >= 0 else 'S'}{abs(a):02d}{'E' if b >= 0 else 'W'}{abs(b):03d}"
                for a, b in zip(lat3, lon3, strict=True)]  # fmt: skip
        groups = list(todo.groupby(np.array(tile)))
        log.info("worldcover: %d cells in %d tiles", len(todo), len(groups))
        env = rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR")
        with env, ThreadPoolExecutor(threads) as pool:
            parts = list(pool.map(lambda g: worldcover_cells(*g), groups))
        done = pd.concat([done, *parts], ignore_index=True)
        done.to_csv(cache, index=False, float_format="%.6g")
    return cells[["cell_id"]].merge(done, on="cell_id", how="left")


def water_fill(available: pd.Series, total: int) -> pd.Series:
    """Equal quota per group, capped at what each group has; the rest go to the others."""
    quota = pd.Series(0, index=available.index)
    left = total
    while left > 0:
        open_ = quota.index[quota < available]
        q = left // len(open_) if len(open_) else 0
        if q == 0:
            break
        add = np.minimum(available[open_] - quota[open_], q)
        quota[open_] += add
        left -= int(add.sum())
    return quota


def pick_controls(cands: pd.DataFrame, n: int = N_CONTROLS) -> pd.DataFrame:
    """The pure candidates with the smallest draw keys, ``water_fill`` quota per biome.

    Args:
        cands: Candidates with ``biome``, ``draw_key``, ``worldcover`` and ``wc_share``.
        n: Controls in total.

    Returns:
        The controls in draw-key order (ties: ``cell_id``).
    """
    pure = cands["worldcover"].isin(NATURAL_WORLDCOVER) & (cands["wc_share"] >= PURE)
    pts = cands[pure].sort_values(["draw_key", "cell_id"], kind="stable")
    quota = water_fill(pts.groupby("biome").size(), n)
    return pts[pts.groupby("biome").cumcount() < pts["biome"].map(quota)].reset_index(drop=True)


def stage_select(args: argparse.Namespace) -> None:
    import compute_eii_v2 as ce

    if not args.hmi.exists():
        record, name, md5 = GHM_2022
        download(ZENODO_URL.format(record=record, name=name), args.hmi, md5)
    rng = np.random.default_rng(SEED)
    cells = draw_cells(N_DRAWS, rng)
    lon, lat = cells["lon"].to_numpy(), cells["lat"].to_numpy()
    eco = sample(args.ecoregions, lon, lat, [1], args.threads)[:, 0]
    lookup = ce.eco_to_biome_lookup()
    known = np.isfinite(eco) & (eco >= 0) & (eco < len(lookup))
    cells["biome"] = np.where(known, lookup[np.where(known, eco, 0).astype("int64")], 0)
    cells = cells[cells["biome"] > 0].reset_index(drop=True)
    lon, lat = cells["lon"].to_numpy(), cells["lat"].to_numpy()
    log.info("select: %d land cells of %d draws", len(cells), N_DRAWS)
    hmi = sample(args.hmi, lon, lat, [1], args.threads)[:, 0]
    ref = pd.read_csv(args.exclude, usecols=["lon", "lat"]) if args.exclude else pd.DataFrame()
    cands = cells[(hmi < HMI_MAX) & ~within_km(lon, lat, ref, EXCLUDE_KM)]
    cands = cands.groupby("biome").head(N_CONTROLS // len(BIOME_NAMES) * OVERSAMPLE)
    # the draw keys continue the generator after the draws, as in the published selection
    cands = cands.assign(draw_key=stored(rng.random(len(cands)), ".7g"),
                         lon=stored(cands["lon"], ".7g"),
                         lat=stored(cands["lat"], ".7g"))  # fmt: skip
    log.info("select: %d candidates", len(cands))
    args.work.mkdir(parents=True, exist_ok=True)
    wc = worldcover_shares(cands, args.work / "worldcover.csv.gz", args.threads)
    controls = pick_controls(cands.merge(wc, on="cell_id"))
    hmi = sample(args.hmi, controls["lon"], controls["lat"], [1], args.threads)[:, 0]
    controls = controls[hmi < HMI_MAX].reset_index(drop=True)
    controls.insert(0, "point_id", np.arange(len(controls)))
    row, col = np.divmod(controls["cell_id"].to_numpy(), NCOLS)
    controls["lon"], controls["lat"] = -180.0 + col / 336, 80.0 - row / 336
    cols = ["point_id", "cell_id", "lon", "lat", "biome", "worldcover"]
    controls[cols].to_csv(args.work / "controls.csv.gz", index=False, float_format="%.10g")
    log.info("select: %d controls -> %s", len(controls), args.work / "controls.csv.gz")


# --------------------------------------------------------------------------------------
# EVI
# --------------------------------------------------------------------------------------


def pixel_of(lon: np.ndarray, lat: np.ndarray, h: int, v: int, n: int) -> tuple[np.ndarray, ...]:
    """Row and column of the pixel holding each point in sinusoidal tile (h, v), n x n pixels."""
    x, y = es.sinusoidal(lon, lat)
    size = es.TILE_M / n
    col = np.floor((x - (es.X0 + h * es.TILE_M)) / size).astype(int)
    row = np.floor(((es.Y0 - v * es.TILE_M) - y) / size).astype(int)
    return np.clip(row, 0, n - 1), np.clip(col, 0, n - 1)


def modis_points(signer: es.Signer, item: dict, row: np.ndarray, col: np.ndarray) -> np.ndarray:
    """Encoded MOD13Q1 codes of one Planetary Computer item at (row, col)."""
    import rasterio
    from rasterio.windows import Window

    bands = []
    for key in ("250m_16_days_EVI", "250m_16_days_pixel_reliability"):
        url = signer.sign(item["assets"][key]["href"])
        with rasterio.Env(**es.GDAL_ENV), rasterio.open(url) as src:
            bands.append([src.read(1, window=Window(c, r, 1, 1))[0, 0]
                          for r, c in zip(row, col, strict=True)])  # fmt: skip
    return encode_modis_evi(*bands)


def encode_ee_modis(evi: np.ndarray, reliability: np.ndarray) -> np.ndarray:
    """``encode_modis_evi`` as applied to the Earth Engine copy, where reliability -1 (no
    reliability given) keeps the EVI."""
    return encode_modis_evi(evi, np.where(np.asarray(reliability) == -1, 0, reliability))


def tile_codes(sensor: str, h: int, v: int, row: np.ndarray, col: np.ndarray, threads: int,
               tmp: Path) -> tuple[np.ndarray, list[str]]:  # fmt: skip
    """Encoded codes (points, composites) of a sensor's record in one tile, and their dates."""
    stats = {"bytes": 0.0, "granules": 0, "read_s": 0.0, "chunks": "", "lock": threading.Lock()}
    jobs = {}
    if sensor == "T":
        dates = es.modis_composite_dates(fct.RECORD["T"])
        signer = es.Signer()
        for item in es.wanted(es.stac_items(signer, h, v, fct.RECORD["T"]), dates):
            day = es.item_date(item["id"]).strftime("%Y_%m_%d")
            jobs[day] = lambda item=item: modis_points(signer, item, row, col)
        missing = [d for d in T_FILL if d not in jobs]
        for day, url in es.granules("T", h, v, missing) if missing else []:
            jobs[day] = lambda url=url: es.read_hdf4(url, stats, tmp, encode_ee_modis)[row, col]
    else:
        for day, url in es.granules(sensor, h, v, es.viirs_composite_dates(fct.RECORD[sensor])):
            jobs[day] = lambda url=url: es.read_viirs(url, stats)[row, col]
    dates = sorted(jobs)
    if not dates:
        return np.empty((len(row), 0), "int16"), dates
    with ThreadPoolExecutor(threads) as pool:
        cols = list(pool.map(lambda d: fct.retried(jobs[d], f"{sensor} h{h:02d}v{v:02d} {d}"),
                             dates))  # fmt: skip
    return np.stack(cols, 1), dates


def stage_extract(args: argparse.Namespace) -> None:
    controls = pd.read_csv(args.work / "controls.csv.gz")
    lon, lat = controls["lon"].to_numpy(), controls["lat"].to_numpy()
    th, tv = es.tile_of(lon, lat)
    names = [f"h{a:02d}v{b:02d}" for a, b in zip(th, tv, strict=True)]
    groups = pd.Series(names).groupby(names).indices
    folder = args.work / "evi" / args.sensor
    folder.mkdir(parents=True, exist_ok=True)
    n = es.SENSORS[args.sensor][2]
    for k, (tile, idx) in enumerate(sorted(groups.items()), 1):
        path = folder / f"{tile}.csv.gz"
        if (args.tiles and tile not in args.tiles) or path.exists():
            continue
        h, v = fct.parse_tile(tile)
        x, y = lon[idx], lat[idx]
        if args.sensor == "T":  # the published positions (module docstring)
            x, y = stored(stored(x, ".7g"), ".6g"), stored(stored(y, ".7g"), ".6g")
        else:
            x, y = stored(x, ".7f"), stored(y, ".7f")
        row, col = pixel_of(x, y, h, v, n)
        codes, dates = tile_codes(args.sensor, h, v, row, col, args.threads, args.work / "tmp")
        frame = pd.DataFrame(codes, columns=dates)
        frame.insert(0, "point_id", controls["point_id"].to_numpy()[idx])
        frame.to_csv(path, index=False)
        log.info("extract %s %d/%d %s: %d points, %d composites", args.sensor, k, len(groups),
                 tile, len(idx), len(dates))  # fmt: skip


# --------------------------------------------------------------------------------------
# Table
# --------------------------------------------------------------------------------------


def annual_series(codes: pd.DataFrame, sensor: str) -> pd.DataFrame:
    """Annual mean EVI (5 significant digits) from composite codes on the union of all dates.

    Args:
        codes: Codes indexed by ``point_id``, one ``YYYY_MM_DD`` column per composite; NaN
            where a point's tile has no such composite (counted as invalid).
        sensor: Sensor letter of the column names.

    Returns:
        ``evi_<sensor>_<year>`` columns for ``YEARS``, indexed by ``point_id``.
    """
    codes = codes.reindex(columns=sorted(codes.columns))
    x = codes.fillna(EVI_INVALID_CODE).to_numpy("int64")
    mean, _ = annual_evi(x, [int(c[:4]) for c in codes.columns], YEARS)
    out = pd.DataFrame(stored(mean.ravel(), ".5g").reshape(mean.shape), index=codes.index,
                       columns=[f"evi_{sensor}_{y}" for y in YEARS])  # fmt: skip
    return out


def control_table(controls: pd.DataFrame, predictors: pd.DataFrame,
                  fi_state: dict[int, np.ndarray], evi: pd.DataFrame) -> pd.DataFrame:  # fmt: skip
    """The table of ``derive_fi_change_constants.py``: ids, groups, predictors, FI_state, EVI.

    Args:
        controls: ``point_id``, ``biome`` and ``worldcover`` per control.
        predictors: The predictor columns, in ``controls`` order.
        fi_state: FI_state per reference year, in ``controls`` order.
        evi: ``evi_<sensor>_<year>`` columns indexed by ``point_id`` (missing points: NaN).

    Returns:
        One row per control, columns in the published order.
    """
    from eii.training.settings import PREDICTOR_VARIABLES

    out = controls[["point_id", "biome", "worldcover"]].reset_index(drop=True)
    out[list(PREDICTOR_VARIABLES)] = predictors[list(PREDICTOR_VARIABLES)].to_numpy("float64")
    for y in REFERENCE_YEARS:
        out[f"fi_state_{y}"] = np.asarray(fi_state[y], "float64")
    cols = [f"evi_{s}_{y}" for s in SENSORS for y in YEARS]
    series = evi.reindex(index=out["point_id"], columns=cols).to_numpy("float64")
    return pd.concat([out, pd.DataFrame(series, columns=cols)], axis=1)


def stage_assemble(args: argparse.Namespace) -> None:
    import rasterio

    from eii.training.settings import PREDICTOR_VARIABLES

    controls = pd.read_csv(args.work / "controls.csv.gz")
    lon, lat = controls["lon"].to_numpy(), controls["lat"].to_numpy()
    th, tv = es.tile_of(lon, lat)
    tiles = {f"h{a:02d}v{b:02d}" for a, b in zip(th, tv, strict=True)}
    parts = []
    for s in SENSORS:
        missing = sorted(t for t in tiles if not (args.work / "evi" / s / f"{t}.csv.gz").exists())
        if missing:
            raise SystemExit(f"extract --sensor {s} is missing {len(missing)} tiles: {missing}")
        files = sorted((args.work / "evi" / s).glob("*.csv.gz"))
        codes = pd.concat([pd.read_csv(f) for f in files]).set_index("point_id")
        parts.append(annual_series(codes, s))
    with rasterio.open(args.predictor_stack) as src:
        bands = [src.descriptions.index(p) + 1 for p in PREDICTOR_VARIABLES]
    pred = sample(args.predictor_stack, lon, lat, bands, args.threads)
    state = {y: sample_tiles(args.run / str(y), lon, lat, "fi_state") for y in REFERENCE_YEARS}
    table = control_table(controls, pd.DataFrame(pred, columns=list(PREDICTOR_VARIABLES)), state,
                          pd.concat(parts, axis=1))  # fmt: skip
    path = args.work / "fi_change_controls.csv.gz"
    table.to_csv(path, index=False, float_format="%.6g")
    log.info("assemble: %d controls -> %s", len(table), path)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("stage", choices=["select", "extract", "assemble"])
    parser.add_argument("work", type=Path, help="work and output folder")
    parser.add_argument("--ecoregions", type=Path, help="select: RESOLVE ecoregion raster")
    parser.add_argument("--hmi", type=Path, help="select: GHM 2022 raster (downloaded if missing)")
    parser.add_argument("--exclude", type=Path, help="select: CSV of points (lon, lat) to avoid")
    parser.add_argument("--sensor", choices=SENSORS,
                        help="extract: T Terra, V Suomi NPP, N NOAA-20")  # fmt: skip
    parser.add_argument("--tiles", nargs="*", default=[], help="extract: only these tiles")
    parser.add_argument("--run", type=Path, help="assemble: compute_eii_v2.py run output")
    parser.add_argument("--predictor-stack", type=Path, help="assemble: stack.vrt")
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()
    need = {"select": ("ecoregions", "hmi"), "extract": ("sensor",),
            "assemble": ("run", "predictor_stack")}[args.stage]  # fmt: skip
    if missing := [f"--{n.replace('_', '-')}" for n in need if getattr(args, n) is None]:
        parser.error(f"{args.stage} needs {', '.join(missing)}")
    {"select": stage_select, "extract": stage_extract, "assemble": stage_assemble}[args.stage](args)


if __name__ == "__main__":
    main()
