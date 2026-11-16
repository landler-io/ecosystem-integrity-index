"""Fetch and verify the published WCMC EII V1 rasters (Hill et al. 2026, CC BY 4.0).

The portal item's download (``https://wcmc.io/eii_download``) resolves to the public S3
object ``EII_V1.zip``. It carries ``EII_V1/EII_1km/EII_<year>_1km.tif`` for 2000, 2005,
2010, 2015 and 2020 plus ``EII_V1/CHECKSUMS.sha256`` and a README. Only the composite EII
is published; the component layers are not distributed as rasters.
"""

from __future__ import annotations

import hashlib
import json
import zipfile
from collections.abc import Iterator
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import rasterio

CHUNK_BYTES = 1 << 20
ZIP_URL = "https://datadownload-production.s3.us-east-1.amazonaws.com/EII_V1.zip"
ZIP_BYTES = 2_419_027_098
YEARS = (2000, 2005, 2010, 2015, 2020)
# Per-cell mean EII over valid land cells and the 2020 land fraction, from the WCMC README.
README_MEAN_EII = {2000: 0.7846, 2005: 0.7877, 2010: 0.7853, 2015: 0.7798, 2020: 0.7898}
README_LAND_FRACTION_2020 = 0.2317


def parse_checksums(text: str) -> dict[str, str]:
    """Parse a ``sha256sum``-format CHECKSUMS.sha256 into {member: hexdigest}."""
    out: dict[str, str] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        digest, _, name = line.partition(" ")
        name = name.strip().lstrip("*")
        out[name] = digest.lower()
    return out


def sha256_file(path: Path, chunk: int = CHUNK_BYTES) -> str:
    """Streaming SHA256 of a local file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_zip_members(
    zip_path: Path,
    checksums: dict[str, str],
    members: list[str] | None = None,
) -> dict[str, str]:
    """Verify zip members against ``checksums``; returns {member: sha256}.

    Raises ``ValueError`` on the first mismatch or missing member. ``members``
    restricts verification to a subset (e.g. only the years we extract).
    """
    wanted = set(members) if members else set(checksums)
    verified: dict[str, str] = {}
    with zipfile.ZipFile(zip_path) as zf:
        for name in sorted(wanted):
            if name not in checksums:
                raise ValueError(f"no checksum recorded for {name}")
            with zf.open(name) as handle:
                digest = hashlib.sha256()
                for block in iter(lambda: handle.read(CHUNK_BYTES), b""):
                    digest.update(block)
            actual = digest.hexdigest()
            if actual != checksums[name]:
                raise ValueError(f"checksum mismatch for {name}: {actual} != {checksums[name]}")
            verified[name] = actual
    return verified


def download(url: str, dest: Path, expected_bytes: int | None = None) -> Path:
    """Stream ``url`` to ``dest`` (no auth; the WCMC zip is public). Returns ``dest``."""
    import urllib.request

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as response, tmp.open("wb") as out:
        while True:
            block = response.read(CHUNK_BYTES)
            if not block:
                break
            out.write(block)
    if expected_bytes is not None and tmp.stat().st_size != expected_bytes:
        tmp.unlink(missing_ok=True)
        raise ValueError(f"downloaded {tmp.stat().st_size} bytes, expected {expected_bytes}")
    tmp.replace(dest)
    return dest


def extract_years(
    zip_path: Path, years: list[int], dest_dir: Path, prefix: str = "EII_V1"
) -> dict[int, Path]:
    """Extract ``EII_<year>_1km.tif`` members into ``dest_dir``; returns {year: path}."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    out: dict[int, Path] = {}
    with zipfile.ZipFile(zip_path) as zf:
        names = set(zf.namelist())
        for year in years:
            member = f"{prefix}/EII_1km/EII_{year}_1km.tif"
            if member not in names:
                raise ValueError(f"{member} missing from {zip_path}")
            target = dest_dir / f"EII_{year}_1km.tif"
            with zf.open(member) as src, target.open("wb") as dst:
                for block in iter(lambda: src.read(CHUNK_BYTES), b""):
                    dst.write(block)
            out[year] = target
    return out


def iter_valid_stats(tif_path: Path, block_rows: int = 2048) -> Iterator[tuple[float, int, int]]:
    """Yield (sum of valid values, valid px, total px) per block; nodata/NaN excluded."""
    with rasterio.open(tif_path) as src:
        nodata = src.nodata
        for row in range(0, src.height, block_rows):
            window = rasterio.windows.Window(0, row, src.width, min(block_rows, src.height - row))
            data = src.read(1, window=window).astype("float64")
            valid = np.isfinite(data)
            if nodata is not None and np.isfinite(nodata):
                valid &= data != nodata
            yield float(np.sum(data[valid])), int(valid.sum()), int(data.size)


def sanity_stats(tif_path: Path) -> dict[str, float]:
    """Global valid-cell mean and valid (land) fraction of a published EII raster."""
    total_sum = valid_px = total_px = 0
    for block_sum, block_valid, block_total in iter_valid_stats(tif_path):
        total_sum += block_sum
        valid_px += block_valid
        total_px += block_total
    return {
        "mean_eii": total_sum / valid_px if valid_px else float("nan"),
        "valid_fraction": valid_px / total_px if total_px else float("nan"),
        "valid_px": valid_px,
        "total_px": total_px,
    }


def fetch_and_verify(dest_dir: Path) -> dict[str, Any]:
    """Download the zip, verify the yearly rasters against its checksums and extract them.

    Writes ``EII_V1.zip`` and ``EII_<year>_1km.tif`` into ``dest_dir`` and records the
    global mean and land fraction of each year next to the values in the WCMC README.

    Args:
        dest_dir: Target directory (``inputs/wcmc`` under the data root).

    Returns:
        The manifest, also written to ``dest_dir/wcmc_manifest.json``.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    zip_path = dest_dir / "EII_V1.zip"
    if not zip_path.exists():
        download(ZIP_URL, zip_path, expected_bytes=ZIP_BYTES)
    with zipfile.ZipFile(zip_path) as zf:
        raw_checksums = parse_checksums(zf.read("EII_V1/CHECKSUMS.sha256").decode())
        names = set(zf.namelist())
    # CHECKSUMS.sha256 paths are relative to EII_V1/ (e.g. "*EII_1km/EII_2000_1km.tif")
    checksums = {n if n in names else f"EII_V1/{n}": d for n, d in raw_checksums.items()}
    members = [f"EII_V1/EII_1km/EII_{y}_1km.tif" for y in YEARS]
    verified = verify_zip_members(zip_path, checksums, members)
    extracted = extract_years(zip_path, list(YEARS), dest_dir)
    sanity: dict[str, dict[str, float]] = {
        str(year): sanity_stats(path) | {"readme_mean_eii": README_MEAN_EII[year]}
        for year, path in extracted.items()
    }
    sanity["2020"]["readme_land_fraction"] = README_LAND_FRACTION_2020
    manifest = {
        "zip": str(zip_path),
        "zip_sha256": sha256_file(zip_path),
        "zip_bytes": zip_path.stat().st_size,
        "source_url": ZIP_URL,
        "license": "CC BY 4.0",
        "citation": "Hill et al. 2026, Ecological Indicators 190, 115314",
        "fetched_at": datetime.now(timezone.utc).isoformat(),
        "verified_members": verified,
        "extracted": {str(y): str(p) for y, p in extracted.items()},
        "sanity": sanity,
    }
    (dest_dir / "wcmc_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return manifest
