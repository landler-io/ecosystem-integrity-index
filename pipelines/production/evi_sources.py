"""Readers of the EVI composites behind the FI change term, and the grids they sit on.

- MODIS Terra MOD13Q1 061 (250 m, 4800 x 4800 per tile) from Microsoft Planetary Computer;
  composites missing there come from the LP DAAC HDF4 granules.
- VIIRS VNP13A1 002 (Suomi NPP) and VJ113A1 002 (NOAA-20), 500 m, 2400 x 2400 per tile, from
  the LP DAAC HDF5 granules (NASA CMR search, Earthdata login).
- ESA WorldCover 2021 v200 from the public AWS bucket.

Credentials: the Planetary Computer subscription key is ``PC_SDK_SUBSCRIPTION_KEY`` or the line
``PC_SDK_SUBSCRIPTION_KEY=...`` in ``PC_SETTINGS_FILE`` (default
``~/.planetarycomputer/settings.env``); without a key, access is anonymous and rate limited.
The Earthdata bearer token is read from ``EDL_TOKEN_FILE`` (default ``~/.edl_token``), else the
``~/.netrc`` login applies on the redirect.

The LP DAAC HDF4 granules are converted with the ``gdal_translate`` command line tool, which
needs GDAL's HDF4 driver; set ``GDAL_TRANSLATE`` to a build that has it if the one on ``PATH``
does not.
"""

from __future__ import annotations

import io
import json
import os
import re
import subprocess
import tempfile
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import date, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from eii.compute import encode_modis_evi, encode_viirs_evi

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
SAS = "https://planetarycomputer.microsoft.com/api/sas/v1/token/{account}/{container}"
CMR = "https://cmr.earthdata.nasa.gov/search/granules.json"
EARTH_R = 6371007.181
TILE_M = 1111950.5197665
X0, Y0 = -20015109.354, 10007554.677
"""MODIS sinusoidal grid: sphere radius, tile size and upper-left corner (m)."""
PC_COLLECTION, PC_PREFIX = "modis-13Q1-061", "MOD13Q1"
SENSORS = {
    "T": ("MOD13Q1", "061", 4800),
    "V": ("VNP13A1", "002", 2400),
    "N": ("VJ113A1", "002", 2400),
}
"""Sensor -> (product, collection version, pixels per tile side)."""
GDAL_ENV = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "GDAL_HTTP_MULTIPLEX": "YES",
    "GDAL_HTTP_MAX_RETRY": "5",
    "GDAL_HTTP_RETRY_DELAY": "2",
    "GDAL_CACHEMAX": 64,
    "CPL_VSIL_CURL_CHUNK_SIZE": "65536",
}
MODIS_SDS = {"evi": "250m 16 days EVI", "rel": "250m 16 days pixel reliability"}
MODIS_GRID = "MODIS_Grid_16DAY_250m_500m_VI"
VIIRS_DATASETS = {"evi": "16 days EVI", "rel": "16 days pixel reliability",
                  "qa": "16 days VI Quality"}  # fmt: skip
VIIRS_FIRST_COMPOSITE = date(2012, 1, 17)
WORLDCOVER = (
    "/vsicurl/https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
    "ESA_WorldCover_10m_2021_v200_{tile}_Map.tif"
)
# HDF4 driver required; the GDAL on PATH (e.g. of a virtualenv) may lack it.
GDAL_TRANSLATE = os.environ.get("GDAL_TRANSLATE", "gdal_translate")

_local = threading.local()


# --------------------------------------------------------------------------------------
# Sinusoidal grid
# --------------------------------------------------------------------------------------


def sinusoidal(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sinusoidal x and y (m) of lon/lat degrees."""
    lat_r = np.radians(lat)
    return EARTH_R * np.radians(lon) * np.cos(lat_r), EARTH_R * lat_r


def tile_of(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Horizontal and vertical sinusoidal tile number of lon/lat points."""
    x, y = sinusoidal(lon, lat)
    return np.floor((x - X0) / TILE_M).astype(int), np.floor((Y0 - y) / TILE_M).astype(int)


def tile_bounds_lonlat(h: int, v: int) -> tuple[float, float, float, float]:
    """lon/lat box of a sinusoidal tile (|lon| grows with |lat|, so the corners bound it)."""
    lat0 = np.degrees((Y0 - (v + 1) * TILE_M) / EARTH_R)
    lat1 = np.degrees((Y0 - v * TILE_M) / EARTH_R)
    x = np.array([X0 + h * TILE_M, X0 + (h + 1) * TILE_M])
    lon = np.degrees(x[:, None] / (EARTH_R * np.cos(np.radians([lat0, lat1]))[None, :]))
    return max(-180.0, lon.min()), lat0, min(180.0, lon.max()), lat1


def modis_composite_dates(years: tuple[int, int]) -> list[str]:
    """``YYYY_MM_DD`` of the MOD13Q1 16-day composites (DOY 1 + 16 k) in ``years``."""
    out = []
    for y in range(years[0], years[1] + 1):
        for k in range(23):
            out.append((date(y, 1, 1) + timedelta(days=16 * k)).strftime("%Y_%m_%d"))
    return out


def viirs_composite_dates(years: tuple[int, int]) -> list[str]:
    """``YYYY_MM_DD`` of the VIIRS 16-day composites (DOY 1 + 16 k, from 17 January 2012)."""
    return [d for d in modis_composite_dates(years)
            if date(*map(int, d.split("_"))) >= VIIRS_FIRST_COMPOSITE]  # fmt: skip


def item_date(item_id: str) -> date:
    """Composite date from the ``.AYYYYDDD.`` part of a MODIS or VIIRS granule id."""
    match = re.search(r"\.A(\d{4})(\d{3})\.", item_id)
    if match is None:
        raise ValueError(f"no .AYYYYDDD. date in {item_id}")
    y, doy = match.groups()
    return date(int(y), 1, 1) + timedelta(days=int(doy) - 1)


# --------------------------------------------------------------------------------------
# Planetary Computer (MOD13Q1)
# --------------------------------------------------------------------------------------


def subscription_key() -> str | None:
    """Planetary Computer key from the environment or ``PC_SETTINGS_FILE``; None if absent."""
    if key := os.environ.get("PC_SDK_SUBSCRIPTION_KEY"):
        return key
    path = Path(os.environ.get("PC_SETTINGS_FILE", Path.home() / ".planetarycomputer/settings.env"))
    if path.exists():
        for line in path.read_text().splitlines():
            if line.startswith("PC_SDK_SUBSCRIPTION_KEY="):
                return line.split("=", 1)[1].strip().strip("\"'")
    return None


class Signer:
    """SAS tokens per storage container, refreshed 5 minutes before expiry."""

    def __init__(self) -> None:
        self.key = subscription_key()
        self.tokens: dict[str, tuple[str, float]] = {}
        self.lock = threading.Lock()

    def headers(self) -> dict[str, str]:
        return {"Ocp-Apim-Subscription-Key": self.key} if self.key else {}

    def sign(self, href: str) -> str:
        account = href.split("//")[1].split(".")[0]
        container = href.split(".net/")[1].split("/")[0]
        name = f"{account}/{container}"
        with self.lock:
            token, expiry = self.tokens.get(name, ("", 0.0))
            if time.time() > expiry - 300:
                req = urllib.request.Request(
                    SAS.format(account=account, container=container), headers=self.headers()
                )
                body = json.load(urllib.request.urlopen(req, timeout=60))
                token, expiry = body["token"], pd.Timestamp(body["msft:expiry"]).timestamp()
                self.tokens[name] = (token, expiry)
        return f"{href}?{token}"


def stac_items(signer: Signer, h: int, v: int, years: tuple[int, int]) -> list[dict]:
    """All MOD13Q1 STAC items of one tile and year range (paged)."""
    body: dict[str, Any] | None = {
        "collections": [PC_COLLECTION],
        "datetime": f"{years[0]}-01-01T00:00:00Z/{years[1]}-12-31T23:59:59Z",
        "query": {"modis:horizontal-tile": {"eq": h}, "modis:vertical-tile": {"eq": v}},
        "limit": 1000,
    }
    items: list[dict] = []
    url = STAC
    while body is not None:
        for attempt in range(6):
            try:
                req = urllib.request.Request(
                    url, data=json.dumps(body).encode(),
                    headers={"Content-Type": "application/json", **signer.headers()},
                )  # fmt: skip
                page = json.load(urllib.request.urlopen(req, timeout=120))
                break
            except Exception:  # noqa: BLE001 - transient STAC API errors; retried
                if attempt == 5:
                    raise
                time.sleep(10 * (attempt + 1))
        items += [f for f in page["features"] if f["id"].startswith(PC_PREFIX)]
        nxt = [link for link in page.get("links", []) if link.get("rel") == "next"]
        body = nxt[0].get("body") if nxt else None
        url = nxt[0]["href"] if nxt else STAC
    return items


def wanted(items: list[dict], dates: list[str]) -> list[dict]:
    """Items on ``dates`` (``YYYY_MM_DD``), one per date, sorted by date."""
    keep = set(dates)
    by_date = {}
    for it in items:
        d = item_date(it["id"]).strftime("%Y_%m_%d")
        if d in keep:
            by_date[d] = it
    return [by_date[d] for d in sorted(by_date)]


# --------------------------------------------------------------------------------------
# LP DAAC (MOD13Q1 fill, VNP13A1, VJ113A1)
# --------------------------------------------------------------------------------------


def session() -> Any:
    """Per-thread ``requests`` session with the Earthdata auth.

    ``requests`` drops the header on the redirect to S3, so the token only goes to the LP DAAC
    host.
    """
    import requests

    if not hasattr(_local, "session"):
        _local.session = requests.Session()
        token = Path(os.environ.get("EDL_TOKEN_FILE", Path.home() / ".edl_token"))
        if token.exists():
            _local.session.headers["Authorization"] = f"Bearer {token.read_text().strip()}"
    return _local.session


def download(url: str) -> bytes:
    """Whole granule in memory (one streamed GET, retried).

    HTTP byte ranges into the HDF5 files are slow and the chunk index is scattered, so one
    40-90 MB download beats many range requests.
    """
    for attempt in range(5):
        try:
            with session().get(url, timeout=300) as r:
                r.raise_for_status()
                return r.content
        except Exception:  # noqa: BLE001 - transient HTTP errors; retried
            if attempt == 4:
                raise
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(url)


def granules(sensor: str, h: int, v: int, dates: list[str]) -> list[tuple[str, str]]:
    """(``YYYY_MM_DD``, HTTPS URL) of a sensor's LP DAAC granules of one tile on ``dates``."""
    import requests

    short, version, _ = SENSORS[sensor]
    suffix = ".hdf" if sensor == "T" else ".h5"  # MODIS HDF4, VIIRS HDF5
    keep = set(dates)
    years = sorted({d[:4] for d in dates})
    params = {"short_name": short, "version": version, "page_size": 2000,
              "temporal": f"{years[0]}-01-01T00:00:00Z,{years[-1]}-12-31T23:59:59Z",
              "attribute[]": [f"int,HORIZONTALTILENUMBER,{h}", f"int,VERTICALTILENUMBER,{v}"]}  # fmt: skip
    out, headers = {}, {}
    while True:
        r = requests.get(CMR, params=params, headers=headers, timeout=120)
        r.raise_for_status()
        entries = r.json()["feed"]["entry"]
        for e in entries:
            d = item_date(e["title"]).strftime("%Y_%m_%d")
            href = [x["href"] for x in e["links"] if x["href"].startswith("https")
                    and x["href"].endswith(suffix)]  # fmt: skip
            if d in keep and href:
                out[d] = href[0]  # the latest production of a date wins
        if "CMR-Search-After" not in r.headers or not entries:
            break
        headers = {"CMR-Search-After": r.headers["CMR-Search-After"]}
    return sorted(out.items())


def read_hdf4(url: str, stats: dict, tmp: Path, encode: Any = encode_modis_evi) -> np.ndarray:
    """Codes of a whole MOD13Q1 HDF4 granule (downloaded, converted in ``tmp``) by ``encode``."""
    import rasterio

    t0 = time.time()
    data = download(url)
    tmp.mkdir(parents=True, exist_ok=True)
    vals = {}
    with tempfile.TemporaryDirectory(dir=tmp) as d:
        src = Path(d) / "granule.hdf"
        src.write_bytes(data)
        for key, sds in MODIS_SDS.items():
            tif = Path(d) / f"{key}.tif"
            name = f'HDF4_EOS:EOS_GRID:"{src}":{MODIS_GRID}:{sds}'
            subprocess.run([GDAL_TRANSLATE, "-q", name, str(tif)], check=True)
            with rasterio.open(tif) as r:
                vals[key] = r.read(1)
    with stats["lock"]:
        stats["bytes"] += len(data)
        stats["granules"] += 1
        stats["read_s"] += time.time() - t0
    return encode(vals["evi"], vals["rel"])


def h5_datasets(f: Any) -> dict[str, Any]:
    """The EVI, reliability and VI quality datasets of a VIIRS granule."""
    import h5py

    found: dict[str, Any] = {}

    def visit(name: str, obj: Any) -> None:
        for key, suffix in VIIRS_DATASETS.items():
            if isinstance(obj, h5py.Dataset) and name.endswith(suffix):
                found[key] = obj

    f.visititems(visit)
    if set(found) != set(VIIRS_DATASETS):
        raise KeyError(f"VIIRS datasets not found, have {sorted(found)}")
    return found


def read_viirs(url: str, stats: dict) -> np.ndarray:
    """Encoded codes of a whole VNP13A1 or VJ113A1 granule (downloaded whole)."""
    import h5py

    t0 = time.time()
    data = download(url)
    with h5py.File(io.BytesIO(data), "r") as f:
        ds = h5_datasets(f)
        out = encode_viirs_evi(ds["evi"][:], ds["rel"][:], ds["qa"][:])
        chunks = str(ds["evi"].chunks)
    with stats["lock"]:
        stats["bytes"] += len(data)
        stats["granules"] += 1
        stats["chunks"] = chunks
        stats["read_s"] += time.time() - t0
    return out


# --------------------------------------------------------------------------------------
# WorldCover
# --------------------------------------------------------------------------------------


def worldcover_codes(lon: np.ndarray, lat: np.ndarray, threads: int = 12) -> np.ndarray:
    """WorldCover 2021 class code at each point from each 3-degree tile's 32x overview.

    Each tile's overview is read whole; points in tiles that do not exist (open ocean) get 0.
    """
    import rasterio

    out = np.zeros(len(lon), "int16")
    lat0, lon0 = np.floor(lat / 3).astype(int) * 3, np.floor(lon / 3).astype(int) * 3
    tiles = pd.DataFrame({"a": lat0, "o": lon0}).groupby(["a", "o"]).indices

    def read(key: tuple[int, int]) -> None:
        a, o = key
        name = f"{'N' if a >= 0 else 'S'}{abs(a):02d}{'E' if o >= 0 else 'W'}{abs(o):03d}"
        idx = tiles[key]
        env = rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", AWS_NO_SIGN_REQUEST="YES")
        try:
            with env, rasterio.open(WORLDCOVER.format(tile=name), overview_level=4) as src:
                arr = src.read(1)
                row, col = rasterio.transform.rowcol(src.transform, lon[idx], lat[idx])
        except rasterio.errors.RasterioIOError:
            return
        row = np.clip(np.asarray(row), 0, arr.shape[0] - 1)
        col = np.clip(np.asarray(col), 0, arr.shape[1] - 1)
        out[idx] = arr[row, col]

    with ThreadPoolExecutor(threads) as pool:
        list(pool.map(read, tiles))
    return out
