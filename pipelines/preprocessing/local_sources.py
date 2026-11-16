"""
Input rasters of the potential-NPP candidate pool: GHM, WorldCover, WDPA and ecoregions.

GHM and WorldCover are read from their public distributions (Zenodo, AWS) and the RESOLVE
ecoregion raster from the requester-pays GCS bucket, through GDAL ``/vsicurl/``, in windows or
by block. Reads from the bucket bill the Google Cloud project ``billing_project()`` (env
``EII_BILLING_PROJECT``). WDPA is not distributed: ``read_wdpa_window`` reads a local raster
built from your own WDPA download. Each entry records the Earth Engine asset it mirrors.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import subprocess
import urllib.request
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data" / "candidate_pool_local"
SOURCE_DIR = DATA_DIR / "sources"

BUCKET = "landler-opendata-gee-assets"
BILLING_ENV = "EII_BILLING_PROJECT"
GCLOUD = os.environ.get("GCLOUD", "gcloud")

# Global Human Modification v3 (Theobald et al. 2025, Sci Data 12:489), the files behind
# projects/sat-io/open-datasets/GHM/HM_1990_2020_OVERALL_300M and HM_2022_300M (AA).
# Zenodo 10.5281/zenodo.14449495 (1990-2020, v1.0.0) and 10.5281/zenodo.14502573 (2022).
GHM_SERIES = {
    1990: ("14449495", "HMv20240801_1990c_AA.tif", "c0062e42cacc017d74eb564a9738aacb"),
    1995: ("14449495", "HMv20240801_1995c_AA.tif", "a4de37fac993955681d5342cab9d5abe"),
    2000: ("14449495", "HMv20240801_2000c_AA.tif", "f9b68bc9e4943063389c3de708c1d3fb"),
    2005: ("14449495", "HMv20240801_2005c_AA.tif", "d58f015826018abb0d53ec23642a39da"),
    2010: ("14449495", "HMv20240801_2010c_AA.tif", "afb6d3f231f82cb815b62a67254fcb3d"),
    2015: ("14449495", "HMv20240801_2015c_AA.tif", "57591b61b5d5c359cf2ebad615ce8e6d"),
    2020: ("14449495", "HMv20240801_2020c_AA.tif", "faa86714ef11b1f2a481324b141522ff"),
}
GHM_2022 = ("14502573", "HMv20240801_2022s_AA_300.tif", "41719a65898c22b5ab74529b2c1052f2")
ZENODO_URL = "https://zenodo.org/api/records/{record}/files/{name}/content"
# Grid of every GHM file and of the WDPA raster: EE "300 m" in EPSG:4326.
GHM_PIXEL = 0.0026949458523585646
GHM_X0, GHM_Y0 = -180.00082337073326, 90.00041168536663
GHM_PERIOD = round(360 / GHM_PIXEL)  # columns per 360 degrees of longitude

# ESA WorldCover 10 m 2021 v200 (ESA/WorldCover/v200), 3x3-degree COGs on AWS open data.
WORLDCOVER_URL = (
    "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
    "ESA_WorldCover_10m_2021_v200_{name}_Map.tif"
)
WORLDCOVER_PIXEL = 1 / 12000

# The bucket's GCS COG copies of the EE annotation assets.
BUCKET_PREFIXES = {
    "ecoregion": "datasets/ecoregions/v1_0/resolve_ecoregions_rasterized",
}
# WDPA high protection was removed from the bucket and EE (WDPA forbids redistribution); the
# raster (IUCN Ia, Ib, II, -1000 m buffer, 1 = protected) must be rebuilt from WDPA on the
# GHM grid and passed as a local file (``read_wdpa_window``).

# Global grid of each sharded source: affine (a, b, c, d, e, f), size, band count, nodata.
SHARD_GRIDS = {
    "ecoregion": {
        "transform": (
            0.04491576420597607, 0.0, -180.0223829375521,
            0.0, -0.04491576420597607, 90.01119146877605,
        ),
        "height": 3340, "width": 8016, "bands": 1, "nodata": -2147483648,
    },
}  # fmt: skip

_GDAL_ENV = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "GDAL_HTTP_MAX_RETRY": "5",
    "GDAL_HTTP_RETRY_DELAY": "2",
    "VSI_CACHE": "TRUE",
}


@contextmanager
def gdal_env(gcs: bool = True) -> Iterator[None]:
    """
    GDAL options as process environment variables, visible to worker threads.

    With ``gcs`` a fresh bearer token (valid ~1 h) is sent with every HTTP request, so
    keep Zenodo/S3 reads outside such a block.
    """
    options = dict(_GDAL_ENV)
    if gcs:
        token = subprocess.run(
            [GCLOUD, "auth", "print-access-token"], capture_output=True, text=True, check=True
        ).stdout.strip()
        options["GDAL_HTTP_HEADERS"] = f"Authorization: Bearer {token}"
    else:
        # Large chunks for whole-window reads from Zenodo/S3; GCS block reads keep the
        # 16 KB default so sparse point reads do not over-fetch (paid egress).
        options["CPL_VSIL_CURL_CHUNK_SIZE"] = str(1 << 20)
    previous = {k: os.environ.get(k) for k in options}
    os.environ.update(options)
    try:
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def billing_project() -> str:
    """Google Cloud project billed for requester-pays reads of ``BUCKET`` (``BILLING_ENV``)."""
    project = os.environ.get(BILLING_ENV)
    if not project:
        raise RuntimeError(
            f"gs://{BUCKET} is requester pays: set {BILLING_ENV} (or --billing-project) to a "
            "Google Cloud project you can bill"
        )
    return project


def bucket_url(name: str) -> str:
    return (
        f"/vsicurl/https://storage.googleapis.com/{BUCKET}/{name}?userProject={billing_project()}"
    )


def zenodo_url(record: str, name: str) -> str:
    return "/vsicurl/" + ZENODO_URL.format(record=record, name=name)


# --------------------------------------------------------------------------------------
# Grids and reads
# --------------------------------------------------------------------------------------


class Shards:
    """A tiled EE export in the bucket: shards named by their pixel offsets on one grid."""

    def __init__(self, key: str):
        self.key = key
        prefix = BUCKET_PREFIXES[key]
        grid = SHARD_GRIDS[key]
        self.transform = grid["transform"]
        self.height, self.width = grid["height"], grid["width"]
        self.band_count = grid["bands"]
        self.nodata = grid["nodata"]
        names = _bucket_listing(prefix)
        offsets = []
        for name in names:
            match = re.search(r"(\d{10})-(\d{10})\.tif$", name)
            offsets.append((int(match[1]), int(match[2])) if match else (0, 0))
        rows = sorted({r for r, _ in offsets}) + [self.height]
        cols = sorted({c for _, c in offsets}) + [self.width]
        self.shards = [
            (name, r, c, rows[rows.index(r) + 1] - r, cols[cols.index(c) + 1] - c)
            for name, (r, c) in zip(names, offsets, strict=True)
        ]

    def rowcol(self, lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        a = self.transform
        return (
            np.floor((np.asarray(lat) - a[5]) / a[4]).astype("int64"),
            np.floor((np.asarray(lon) - a[2]) / a[0]).astype("int64"),
        )

    def sample(self, lon: np.ndarray, lat: np.ndarray, workers: int = 16) -> np.ndarray:
        """Nearest-pixel values ``(n, bands)`` by block reads; NaN for nodata or outside."""
        rows, cols = self.rowcol(lon, lat)
        out = np.full((len(rows), self.band_count), np.nan)
        tasks = []
        for name, r0, c0, h, w in self.shards:
            inside = np.flatnonzero((rows >= r0) & (rows < r0 + h) & (cols >= c0) & (cols < c0 + w))
            if len(inside):
                tasks.append((name, inside, rows[inside] - r0, cols[inside] - c0))
        blocks = []
        for name, idx, r, c in tasks:
            block_id = (r // 512) * 100_000 + c // 512
            for b in np.unique(block_id):
                sel = block_id == b
                blocks.append((name, idx[sel], r[sel], c[sel]))

        def read(task: tuple) -> None:
            name, idx, r, c = task
            win = Window(c.min(), r.min(), c.max() - c.min() + 1, r.max() - r.min() + 1)
            with rasterio.open(bucket_url(name)) as ds:
                data = ds.read(window=win).astype("float64")
                if ds.nodata is not None:
                    data[data == ds.nodata] = np.nan
            out[idx] = data[:, r - r.min(), c - c.min()].T

        with ThreadPoolExecutor(workers) as pool:
            list(pool.map(read, blocks))
        return out

    def read_window(self, row0: int, col0: int, height: int, width: int) -> np.ndarray:
        """Band 1 over a global-grid window, mosaicked from the shards (NaN = nodata)."""
        out = np.full((height, width), np.nan, dtype="float32")
        for name, r0, c0, h, w in self.shards:
            top, left = max(row0, r0), max(col0, c0)
            bottom, right = min(row0 + height, r0 + h), min(col0 + width, c0 + w)
            if top >= bottom or left >= right:
                continue
            with rasterio.open(bucket_url(name)) as ds:
                data = ds.read(1, window=Window(left - c0, top - r0, right - left, bottom - top))
                data = data.astype("float32")
                if ds.nodata is not None:
                    data[data == ds.nodata] = np.nan
            out[top - row0 : bottom - row0, left - col0 : right - col0] = data
        return out


def _bucket_listing(prefix: str) -> list[str]:
    """Object names under ``prefix`` (cached locally; listing is a class-A operation)."""
    cache_file = SOURCE_DIR / "listings" / (prefix.replace("/", "__") + ".json")
    if cache_file.exists():
        return json.loads(cache_file.read_text())
    lines = subprocess.run(
        [
            GCLOUD,
            "storage",
            "ls",
            f"--billing-project={billing_project()}",
            f"gs://{BUCKET}/{prefix}/",
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    names = sorted(line.split(f"gs://{BUCKET}/", 1)[1] for line in lines if line.endswith(".tif"))
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    cache_file.write_text(json.dumps(names))
    return names


def read_ghm_window(
    record: str, name: str, row0: int, col0: int, height: int, width: int
) -> np.ndarray:
    """
    GHM window on the GHM grid (float32, NaN outside land), wrapped at the antimeridian.

    Rows beyond the grid are NaN. Columns wrap with ``GHM_PERIOD``; the grid is not exactly
    360 degrees wide, so wrapped columns sit up to half a pixel off (halo only).
    """
    out = np.full((height, width), np.nan, dtype="float32")
    local = SOURCE_DIR / "ghm" / name
    with rasterio.open(local if local.exists() else zenodo_url(record, name)) as ds:
        top, bottom = max(row0, 0), min(row0 + height, ds.height)
        for offset, col, n in wrapped_runs(col0, width, GHM_PERIOD):
            n = min(n, ds.width - col)
            if n > 0 and bottom > top:
                win = Window(col, top, n, bottom - top)
                out[top - row0 : bottom - row0, offset : offset + n] = ds.read(1, window=win)
    return out


def read_wdpa_window(path: Path, row0: int, col0: int, height: int, width: int) -> np.ndarray:
    """
    WDPA high protection (bool) over a GHM-grid window, from a local raster on the GHM grid.

    The raster (> 0 = protected) may cover any part of the grid; pixels outside it, nodata and
    0 are False.
    """
    with rasterio.open(path) as ds:
        t = ds.transform
        col, row = (t.c - GHM_X0) / GHM_PIXEL, (GHM_Y0 - t.f) / GHM_PIXEL
        if not (
            math.isclose(t.a, GHM_PIXEL, rel_tol=1e-6)
            and math.isclose(-t.e, GHM_PIXEL, rel_tol=1e-6)
            and abs(col - round(col)) < 1e-3
            and abs(row - round(row)) < 1e-3
        ):
            raise ValueError(f"{path} is not on the GHM grid")
        window = Window(col0 - round(col), row0 - round(row), width, height)
        data = ds.read(1, window=window, boundless=True, fill_value=0)
        valid = data != ds.nodata if ds.nodata is not None else np.ones(data.shape, bool)
    return valid & (data > 0)


def wrapped_runs(col0: int, width: int, period: int) -> list[tuple[int, int, int]]:
    """``(offset, column, length)`` runs covering columns ``col0 .. col0 + width`` mod ``period``."""
    runs, done = [], 0
    while done < width:
        col = (col0 + done) % period
        n = min(width - done, period - col)
        runs.append((done, col, n))
        done += n
    return runs


def download_ghm() -> None:
    """Download the eight GHM files once (~76 GB; whole files stream ~17x faster than windows)."""
    for record, name, md5 in [*GHM_SERIES.values(), GHM_2022]:
        download(ZENODO_URL.format(record=record, name=name), SOURCE_DIR / "ghm" / name, md5)


# --------------------------------------------------------------------------------------
# WorldCover
# --------------------------------------------------------------------------------------


def worldcover_tiles(bounds: tuple[float, float, float, float]) -> list[str]:
    """Names (lower-left corner, e.g. N48E009) of the 3-degree tiles overlapping ``bounds``."""
    west, south, east, north = bounds
    names = []
    for lat in range(math.floor(south / 3) * 3, math.ceil(north / 3) * 3, 3):
        for lon in range(math.floor(west / 3) * 3, math.ceil(east / 3) * 3, 3):
            ns, ew = ("N" if lat >= 0 else "S"), ("E" if lon >= 0 else "W")
            names.append(f"{ns}{abs(lat):02d}{ew}{abs(lon):03d}")
    return names


def _url_exists(url: str) -> bool:
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=60):
            return True
    except urllib.error.HTTPError as error:
        if error.code in (403, 404):
            return False
        raise


def worldcover_vrt(bounds: tuple[float, float, float, float], path: Path) -> Path:
    """VRT over the existing WorldCover COGs that overlap ``bounds`` (ocean tiles are absent)."""
    sources = []
    west, north = math.floor(bounds[0] / 3) * 3, math.ceil(bounds[3] / 3) * 3
    for name in worldcover_tiles(bounds):
        url = WORLDCOVER_URL.format(name=name)
        if not _url_exists(url):
            continue
        lat = int(name[1:3]) * (1 if name[0] == "N" else -1)
        lon = int(name[4:7]) * (1 if name[3] == "E" else -1)
        x_off, y_off = round((lon - west) * 12000), round((north - lat - 3) * 12000)
        sources.append(
            f'<SimpleSource><SourceFilename relativeToVRT="0">/vsicurl/{url}</SourceFilename>'
            f'<SourceBand>1</SourceBand><SrcRect xOff="0" yOff="0" xSize="36000" ySize="36000"/>'
            f'<DstRect xOff="{x_off}" yOff="{y_off}" xSize="36000" ySize="36000"/></SimpleSource>'
        )
    cols = round((math.ceil(bounds[2] / 3) * 3 - west) * 12000)
    rows = round((north - math.floor(bounds[1] / 3) * 3) * 12000)
    path.write_text(
        f'<VRTDataset rasterXSize="{cols}" rasterYSize="{rows}"><SRS>EPSG:4326</SRS>'
        f"<GeoTransform>{west}, {WORLDCOVER_PIXEL}, 0, {north}, 0, {-WORLDCOVER_PIXEL}"
        '</GeoTransform><VRTRasterBand dataType="Byte" band="1"><NoDataValue>0</NoDataValue>'
        + "".join(sources)
        + "</VRTRasterBand></VRTDataset>"
    )
    return path


def worldcover_mode(vrt: Path, transform: rasterio.Affine, shape: tuple[int, int]) -> np.ndarray:
    """Modal WorldCover class on a destination grid (0 = no data)."""
    from rasterio.warp import Resampling, reproject

    out = np.zeros(shape, dtype="uint8")
    with rasterio.open(vrt) as ds:
        reproject(
            rasterio.band(ds, 1),
            out,
            src_nodata=0,
            dst_transform=transform,
            dst_crs="EPSG:4326",
            dst_nodata=0,
            resampling=Resampling.mode,
            warp_mem_limit=256,
        )
    return out


# --------------------------------------------------------------------------------------
# Downloads
# --------------------------------------------------------------------------------------


def download(url: str, path: Path, md5: str | None = None) -> Path:
    """``url`` to ``path`` once (via ``.part``), checked against ``md5`` if given."""
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp, digest = path.with_suffix(path.suffix + ".part"), hashlib.md5()
    with urllib.request.urlopen(url, timeout=600) as response, tmp.open("wb") as handle:
        while chunk := response.read(1 << 22):
            digest.update(chunk)
            handle.write(chunk)
    if md5 and digest.hexdigest() != md5:
        raise ValueError(f"md5 mismatch for {url}: {digest.hexdigest()} != {md5}")
    tmp.rename(path)
    return path
