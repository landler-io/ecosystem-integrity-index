#!/usr/bin/env python3
"""Global EII v2 on the local stack: one GeoTIFF per 10 deg land tile and reference year.

The grid is that of the potential-NPP tiles (``predict_potential_npp_v2.py``): CLMS NPP,
EPSG:4326, 1/336 deg. Inputs:

- potential NPP and its 68 % interval (``<potential>/tiles``, ``potential_npp.vrt``, gC);
- the natural inter-annual CV of NPP (``--natural-cv``, ``natural_cv.vrt`` of
  ``predict_potential_npp_v2.py --natural-cv``);
- annual CLMS DMP V2.0.2 NPP sums from ``--clms-dir`` (``pipelines/preprocessing/npp``, band 1,
  stored units);
- BII v1.1, read from the requester-pays bucket with GDAL (``gcs_env``: the project billed is
  ``EII_BILLING_PROJECT`` or ``--billing-project``);
- GHM 2022 (Zenodo 10.5281/zenodo.14502573), downloaded once to ``--hmi`` if missing;
- GLWD v2 dominant class and cover (``--glwd-class``, ``--glwd-cover``; Lehner et al. 2025,
  figshare 10.6084/m9.figshare.28519994, ``GLWD_v2_0_combined_classes_tif.zip``);
- the RESOLVE ecoregion raster (``--ecoregions``, ``datasets/ecoregions/v1_0`` as synced by
  ``sync_predictor_sources.py``), mapped to biomes by
  ``pipelines/preprocessing/reference/eco_id_to_biome.json``;
- the per-sensor EVI ratios of ``fi_change_term.py`` (``--fi-change``) and the predictor stack
  (``--predictor-stack``, band descriptions = predictor names) for the FI change term.

Per tile and year, with ``eii.compute`` and the defaults of ``eii.compute.settings``:

- functional state: observed NPP (3-year window, >= 2 valid years) against potential by
  ``deviation_tolerance``: deviations beyond the tolerance ``k sqrt(u^2 + cv^2)`` (u the relative
  68 % interval half-width on the side of the deviation, cv the natural CV), deficit k
  ``TOLERANCE_K``, surplus k ``TOLERANCE_K_SURPLUS``, surplus reaching 0 at ``SURPLUS_ZERO``.
  Potential NPP and its interval are recentred per biome and observed window: the factor is the
  ratio of mean observed to mean potential NPP on reference pixels (HMI <
  ``RECENTRE_REFERENCE_HMI_MAX`` and BII >= ``RECENTRE_REFERENCE_BII_MIN`` of the inputs);
  biomes with <= ``RECENTRE_MIN_REFERENCE_PIXELS`` reference pixels keep factor 1;
- FI change term: the serving sensor series of each reference year (``series_for`` on the
  complete sensor years), its sensors' r from ``r_<sensor>.vrt`` combined per cell, the
  tolerance h from the stored models on the predictor stack, FI_change;
  ``functional_integrity`` = FI_state x FI_change (FI_state where there is no change value);
- structural: v2 landscape quality of the logistic habitat quality of HMI
  (``QUALITY_MIDPOINT``), HMI read on its native grid with a halo of the kernel radius
  (26.5 km plus the 300 m opening), then bilinear onto the grid;
- compositional: BII area-averaged onto the grid (its zeros, i.e. water and source gaps,
  as missing), clamped to [0, 1], masked outside land. Land is where potential NPP exists
  (MERIT DEM land);
- EII by ``min_fuzzy_logic``, masked where any component is missing, plus the limiting
  component (1 F, 2 S, 3 C), the wetland flag (0 none, 1 wetland, 2 peatland, 3 mangrove) and
  the bands ``fi_state``, ``fi_change`` and ``fi_change_sensors`` (0-2; 0 where FI = FI_state).
  The serving series and the change-term constants' SHA-256 go into the tags.

Stages:

- ``factors``: first pass of recentring. Per-biome sums of observed and potential NPP on
  reference pixels, streamed over all tiles, written per year to
  ``<out>/recentring/<year>.json`` (factors, pixel counts, sums) and ``<year>_sums.npy``
  (the ratio histogram of ``recentring_sums``). ``run`` runs it for missing years.
- ``run``: ``<out>/<year>/<tile>.tif`` with float32 bands ``BANDS`` and the three change-term
  bands, nodata -9999, and one ``qa.csv`` row per tile and year; ``<out>/run_manifest.json``
  records the configuration and the recentring factors. Existing outputs are skipped.

Example:
    python pipelines/production/compute_eii_v2.py run potential_npp_v2 eii_v2 \\
        --clms-dir npp_annual_sums_v2 --natural-cv potential_npp_v2/natural_cv.vrt \\
        --hmi eii_inputs/HMv20240801_2022s_AA_300.tif \\
        --glwd-class eii_inputs/GLWD_v2_0_main_class.tif \\
        --glwd-cover eii_inputs/GLWD_v2_0_area_pct.tif \\
        --ecoregions predictor_stack/ecoregions/resolve_ecoregions_rasterized/resolve_ecoregions_rasterized_0000.tif \\
        --fi-change fi_change --predictor-stack predictor_stack/stack.vrt
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window

PREPROCESSING = Path(__file__).resolve().parents[1] / "preprocessing"
sys.path.insert(0, str(PREPROCESSING))
from build_predictor_stack import NODATA, PROFILE, warp  # noqa: E402
from local_sources import GHM_2022, ZENODO_URL, billing_project, download  # noqa: E402

from eii.compute import (  # noqa: E402
    combine_components,
    combine_ratios,
    compositional_integrity,
    constants_sha256,
    deviation_tolerance_integrity,
    fi_change,
    functional_with_change,
    limiting_component,
    observed_npp,
    observed_npp_window,
    recentre_potential,
    recentring_factors,
    recentring_sums,
    series_for,
    structural_halo_m,
    structural_integrity,
    tolerance_h,
    wetland_flag,
)
from eii.compute.settings import (  # noqa: E402
    BIOME_NAMES,
    FI_CHANGE_K,
    FI_CHANGE_SERIES,
    METERS_PER_DEGREE,
    N_BIOMES,
    POTENTIAL_NPP_FLOOR,
    QUALITY_MIDPOINT,
    RECENTRE_MIN_REFERENCE_PIXELS,
    RECENTRE_RATIO_BINS,
    RECENTRE_REFERENCE_BII_MIN,
    RECENTRE_REFERENCE_HMI_MAX,
    REFERENCE_YEARS,
    SURPLUS_ZERO,
    TOLERANCE_K,
    TOLERANCE_K_SURPLUS,
)

BUCKET = "landler-opendata-gee-assets"
NPP_SOURCE = "CLMS DMP 300m V2.0.2 RT6 annual sums"
BII = f"/vsigs/{BUCKET}/datasets/bii/v1_1/bii_{{year}}.tif"
ECO_TO_BIOME = PREPROCESSING / "reference" / "eco_id_to_biome.json"
BANDS = [
    "eii",
    "functional_integrity",
    "structural_integrity",
    "compositional_integrity",
    "limiting_component",
    "wetland_flag",
]
MARGIN = 1  # pixels read around each tile; padded arrays are cropped to the tile


def gcs_env() -> None:
    """GDAL /vsigs/ settings for the requester-pays BII reads, inherited by the workers.

    GDAL finds the credentials itself (``GOOGLE_APPLICATION_CREDENTIALS``, ``~/.boto``, the
    Compute Engine metadata server); if ``GOOGLE_APPLICATION_CREDENTIALS`` is unset, the gcloud
    application-default credentials file (user or service account) is used when it exists.
    """
    adc = Path.home() / ".config/gcloud/application_default_credentials.json"
    if "GOOGLE_APPLICATION_CREDENTIALS" not in os.environ and adc.exists():
        os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = str(adc)
    os.environ.update(
        GS_USER_PROJECT=billing_project(),
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        GDAL_HTTP_MAX_RETRY="5",
        GDAL_HTTP_RETRY_DELAY="2",
    )


def clms_vrt(year: int, out: Path, root: Path, band: int = 1) -> Path:
    """One-band VRT (1 annual sum, 2 intra-annual daily SD) of a year's annual NPP GeoTIFF."""
    path = root / f"annual_npp_{year}.tif"
    with rasterio.open(path) as src:
        width, height, t = src.width, src.height, src.transform
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        f'<VRTDataset rasterXSize="{width}" rasterYSize="{height}"><SRS>EPSG:4326</SRS>'
        f"<GeoTransform>{t.c}, {t.a}, 0, {t.f}, 0, {t.e}</GeoTransform>"
        '<VRTRasterBand dataType="Int32" band="1"><NoDataValue>-1</NoDataValue>'
        f'<SimpleSource><SourceFilename relativeToVRT="0">{escape(str(path))}</SourceFilename>'
        f'<SourceBand>{band}</SourceBand><SrcRect xOff="0" yOff="0" xSize="{width}" '
        f'ySize="{height}"/><DstRect xOff="0" yOff="0" xSize="{width}" ySize="{height}"/>'
        "</SimpleSource></VRTRasterBand></VRTDataset>\n"
    )
    return out


def read_padded(src: rasterio.DatasetReader, window: Window, bands: list[int]) -> np.ndarray:
    """float32 bands over ``window``; nodata and pixels outside the dataset are NaN."""
    out = np.full((len(bands), int(window.height), int(window.width)), np.nan, dtype="float32")
    c0, r0 = int(window.col_off), int(window.row_off)
    left, top = max(c0, 0), max(r0, 0)
    right = min(c0 + int(window.width), src.width)
    bottom = min(r0 + int(window.height), src.height)
    if right > left and bottom > top:
        data = src.read(bands, window=Window(left, top, right - left, bottom - top))
        data = data.astype("float32")
        if src.nodata is not None:
            data[data == src.nodata] = np.nan
        out[:, top - r0 : bottom - r0, left - c0 : right - c0] = data
    return out


def read_aligned(path: Path | str, transform, size: int, bands: list[int]) -> np.ndarray:
    """Bands of a raster on the tile's grid, with ``MARGIN`` pixels around the tile."""
    with rasterio.open(path) as src:
        col = round((transform.c - src.transform.c) / src.transform.a) - MARGIN
        row = round((transform.f - src.transform.f) / src.transform.e) - MARGIN
        return read_padded(src, Window(col, row, size + 2 * MARGIN, size + 2 * MARGIN), bands)


def read_warped(
    path: Path | str, transform, size: int, halo_deg: tuple[float, float] = (0.0, 0.0)
) -> tuple[np.ndarray, rasterio.Affine]:
    """Band 1 of another grid around the tile (plus ``halo_deg`` x/y) and its transform."""
    west, south, east, north = rasterio.transform.array_bounds(size, size, transform)
    with rasterio.open(path) as src:
        pad = 2 * abs(src.transform.a)
        box = rasterio.windows.from_bounds(
            west - halo_deg[0] - pad,
            south - halo_deg[1] - pad,
            east + halo_deg[0] + pad,
            north + halo_deg[1] + pad,
            transform=src.transform,
        )
        col, row = math.floor(box.col_off), math.floor(box.row_off)
        width = math.ceil(box.col_off + box.width) - col
        window = Window(col, row, width, math.ceil(box.row_off + box.height) - row)
        return read_padded(src, window, [1])[0], src.window_transform(window)


def to_tile(data: np.ndarray, src_transform, transform, size: int, resampling: Resampling):
    return warp(data[None], src_transform, transform, size, resampling)[0]


def tile_structural(hmi_path: Path, transform, size: int) -> np.ndarray:
    """Structural integrity on the tile grid, from HMI on its native grid with a halo."""
    radius_deg = structural_halo_m() / METERS_PER_DEGREE
    west, south, east, north = rasterio.transform.array_bounds(size, size, transform)
    max_lat = min(max(abs(south), abs(north)) + radius_deg, 89.0)
    halo = (radius_deg / math.cos(math.radians(max_lat)), radius_deg)
    # No wrap at the antimeridian: halo pixels beyond +-180 count as no data.
    hmi, hmi_transform = read_warped(hmi_path, transform, size, halo)
    hmi[(hmi < 0) | (hmi > 1)] = np.nan
    s = structural_integrity(hmi, hmi_transform.f, hmi_transform.a)
    return to_tile(s, hmi_transform, transform, size, Resampling.bilinear)


def tile_bii(year: int, transform, size: int) -> np.ndarray:
    """BII of ``year`` area-averaged onto the tile; its zeros (water, source gaps) are NaN."""
    bii, bii_transform = read_warped(BII.format(year=year), transform, size)
    # BII is exactly 0 on sea, inland water and source gaps (east of 174 E), never on land.
    bii[bii == 0] = np.nan
    return to_tile(bii, bii_transform, transform, size, Resampling.average)


def tile_biome(cfg: dict, transform, size: int) -> np.ndarray:
    """RESOLVE biome number (0 unknown) per tile pixel, with ``MARGIN`` pixels of 0 around."""
    eco = to_tile(*read_warped(cfg["ecoregions"], transform, size), transform, size,
                  Resampling.nearest)  # fmt: skip
    lookup = cfg["eco_to_biome"]
    index = np.nan_to_num(eco, nan=-1).astype("int64")
    known = (index >= 0) & (index < len(lookup))
    biome = np.where(known, lookup[np.where(known, index, 0)], 0)
    return np.pad(biome, MARGIN)


def eco_to_biome_lookup() -> np.ndarray:
    """Biome number per RESOLVE ``ECO_ID`` (0 where the ID is not in the table)."""
    table = {int(k): int(v) for k, v in json.loads(ECO_TO_BIOME.read_text()).items()}
    lookup = np.zeros(max(table) + 1, dtype="int64")
    for eco_id, biome in table.items():
        lookup[eco_id] = biome
    return lookup


def reference_pixels(cfg: dict, transform, size: int, year: int) -> np.ndarray:
    """Recentring reference pixels of one tile (padded like the potential arrays)."""
    hmi, hmi_transform = read_warped(cfg["hmi"], transform, size)
    hmi[(hmi < 0) | (hmi > 1)] = np.nan
    hmi = to_tile(hmi, hmi_transform, transform, size, Resampling.average)
    bii = tile_bii(year, transform, size)
    with np.errstate(invalid="ignore"):
        reference = (hmi < RECENTRE_REFERENCE_HMI_MAX) & (bii >= RECENTRE_REFERENCE_BII_MIN)
    return np.pad(reference, MARGIN)


def factors_tile(tile: str, years: list[int], cfg: dict) -> dict[int, np.ndarray]:
    """First recentring pass: ``recentring_sums`` of one tile for each year."""
    with rasterio.open(cfg["potential"] / "tiles" / f"{tile}.tif") as src:
        transform, size = src.transform, src.width
    potential = read_aligned(cfg["potential"] / "potential_npp.vrt", transform, size, [1])[0]
    if not np.isfinite(potential).any():
        return {}
    biome = tile_biome(cfg, transform, size)
    annual: dict[int, np.ndarray] = {}
    sums = {}
    for year in years:
        for y in observed_npp_window(year):
            if y not in annual:
                annual[y] = read_aligned(cfg["clms"][y], transform, size, [1])[0]
        actual = observed_npp(np.stack([annual[y] for y in observed_npp_window(year)]))
        reference = reference_pixels(cfg, transform, size, year)
        inner = np.zeros_like(reference)
        inner[MARGIN:-MARGIN, MARGIN:-MARGIN] = True  # margins belong to the neighbours
        sums[year] = recentring_sums(actual, potential, biome, reference & inner)
    return sums


def factors_path(out: Path, year: int) -> Path:
    return out / "recentring" / f"{year}.json"


def compute_factors(years: list[int], tiles: list[str], cfg: dict, workers: int) -> None:
    """Accumulate the per-biome sums over all tiles and write each year's factors."""
    total = {year: np.zeros((N_BIOMES, RECENTRE_RATIO_BINS + 1, 3)) for year in years}
    t0 = time.perf_counter()
    with ProcessPoolExecutor(workers) as pool:
        futures = {pool.submit(factors_tile, t, years, cfg): t for t in tiles}
        for i, future in enumerate(as_completed(futures), start=1):
            for year, sums in future.result().items():
                total[year] = total[year] + sums
            print(f"factors {i}/{len(tiles)} {futures[future]}", flush=True)
    for year in years:
        sums = total[year]
        factors = recentring_factors(sums)
        record = {
            "year": year,
            "window": observed_npp_window(year),
            **recentring_config(cfg),
            "tiles": tiles,
            "seconds": round(time.perf_counter() - t0, 1),
            "biomes": {
                str(b): {
                    "name": BIOME_NAMES[b],
                    "reference_pixels": int(sums[b, :, 0].sum()),
                    "sum_observed": float(sums[b, :, 1].sum()),
                    "sum_potential": float(sums[b, :, 2].sum()),
                    "factor": float(factors[b]),
                }
                for b in BIOME_NAMES
            },
        }
        path = factors_path(cfg["out"], year)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.save(path.with_name(f"{year}_sums.npy"), sums)
        path.write_text(json.dumps(record, indent=1) + "\n")
        print(json.dumps({"year": year, "factors": factor_summary(record)}), flush=True)


def recentring_config(cfg: dict) -> dict:
    """The settings that define a set of recentring factors."""
    return {
        "reference": "hmi_bii",
        "min_reference_pixels": RECENTRE_MIN_REFERENCE_PIXELS,
        "potential": str(cfg["potential"]),
        "hmi": Path(cfg["hmi"]).name,
        "hmi_max": RECENTRE_REFERENCE_HMI_MAX,
        "bii": "datasets/bii/v1_1",
        "bii_min": RECENTRE_REFERENCE_BII_MIN,
    }


def factor_summary(record: dict) -> dict[str, float]:
    return {b: round(v["factor"], 4) for b, v in record["biomes"].items()}


def load_factors(cfg: dict, year: int, tiles: list[str]) -> tuple[np.ndarray, dict]:
    """Factors of one year; refuses factors from other settings or another tile set."""
    record = json.loads(factors_path(cfg["out"], year).read_text())
    expected = recentring_config(cfg)
    stale = {k for k, v in expected.items() if record.get(k) != v}
    if stale or not set(tiles) <= set(record["tiles"]):
        raise SystemExit(
            f"{factors_path(cfg['out'], year)} was computed with other "
            f"{sorted(stale) or 'tiles'}; delete it (and the year's outputs) to recompute"
        )
    factors = np.ones(max(BIOME_NAMES) + 1)
    for b, v in record["biomes"].items():
        factors[int(b)] = v["factor"]
    return factors, record


def write(path: Path, bands: dict[str, np.ndarray], transform, tags: dict) -> None:
    size = bands["eii"].shape[0]
    image = np.stack(list(bands.values())).astype("float32")
    image[~np.isfinite(image)] = NODATA
    part = path.with_name(path.name + ".part")
    with rasterio.open(
        part, "w", width=size, height=size, count=len(bands), transform=transform, **PROFILE
    ) as dst:
        dst.write(image)
        dst.descriptions = tuple(bands)
        dst.update_tags(**{k: str(v) for k, v in tags.items()})
    part.rename(path)


def fi_change_config(root: Path, stack: Path, years: list[int], workers: int) -> dict:
    """Serving series per reference year, predictor bands and constants of the change term."""
    from eii.training.settings import PREDICTOR_VARIABLES

    record = json.loads((root / "completeness.json").read_text())
    complete = {s: set(v["complete_years"]) for s, v in record.items()}
    series = {y: series_for(y, complete) for y in years}
    for y, c in series.items():
        missing = [s for s in FI_CHANGE_SERIES.get(c, ()) if not (root / f"r_{s}.vrt").exists()]
        if missing:
            raise SystemExit(f"{y} is served by {c}, but {root} has no r for {missing}")
    with rasterio.open(stack) as src:
        bands = [src.descriptions.index(p) + 1 for p in PREDICTOR_VARIABLES]
    return {"dir": root, "series": series, "stack": stack, "bands": bands,
            "sha256": constants_sha256(), "threads": max(1, (os.cpu_count() or 1) // workers)}  # fmt: skip


def tile_fi_change(
    fc: dict, year: int, transform, size: int, land: np.ndarray, predictors: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """FI_change and the number of contributing sensors of one tile and reference year."""
    inner = (slice(MARGIN, -MARGIN), slice(MARGIN, -MARGIN))
    change = np.full(land.shape, np.nan)
    sensors = np.where(land, 0.0, np.nan)
    series = fc["series"][year]
    if series is None:
        return change, sensors
    band = REFERENCE_YEARS.index(year) + 1
    members = FI_CHANGE_SERIES[series]
    rs = [read_aligned(fc["dir"] / f"r_{s}.vrt", transform, size, [band])[0][inner]
          for s in members]  # fmt: skip
    r, count = combine_ratios(rs)
    sel = land & (count > 0)
    h = tolerance_h(series, predictors[:, sel].T, None if len(members) == 1 else count[sel],
                    threads=fc["threads"])  # fmt: skip
    score = fi_change(r[sel], h, FI_CHANGE_K[series])
    change[sel] = score
    sensors[sel] = np.where(np.isfinite(score), count[sel], 0)
    return change, sensors


def run_tile(tile: str, years: list[int], cfg: dict) -> list[dict]:
    with rasterio.open(cfg["potential"] / "tiles" / f"{tile}.tif") as src:
        transform, size = src.transform, src.width
    inner = (slice(MARGIN, -MARGIN), slice(MARGIN, -MARGIN))
    potential, lower, upper = read_aligned(
        cfg["potential"] / "potential_npp.vrt", transform, size, [1, 2, 3]
    )
    natural_cv = read_aligned(cfg["natural_cv"], transform, size, [1])[0]
    land = np.isfinite(potential[inner])
    if not land.any():
        return []
    biome = tile_biome(cfg, transform, size)
    structural = tile_structural(cfg["hmi"], transform, size)
    structural = np.where(land, structural, np.nan)
    glwd = [
        to_tile(*read_warped(cfg[k], transform, size), transform, size, Resampling.nearest)
        for k in ("glwd_class", "glwd_cover")
    ]
    flag = np.where(land, wetland_flag(*glwd), np.nan)

    fc = cfg["fi_change"]
    predictors = read_aligned(fc["stack"], transform, size, fc["bands"])[(slice(None), *inner)]

    annual: dict[int, np.ndarray] = {}
    rows = []
    for year in years:
        out = cfg["out"] / str(year) / f"{tile}.tif"
        if out.exists():
            continue
        t0 = time.perf_counter()
        for y in observed_npp_window(year):
            if y not in annual:
                annual[y] = read_aligned(cfg["clms"][y], transform, size, [1])[0]
        actual = observed_npp(np.stack([annual[y] for y in observed_npp_window(year)]))
        interval = recentre_potential(potential, lower, upper, biome, cfg["factors"][year])
        fi_state = deviation_tolerance_integrity(
            actual,
            *interval,
            k=cfg["tolerance_k"],
            natural_cv=natural_cv,
            k_surplus=cfg["tolerance_k_surplus"],
            surplus_zero=cfg["surplus_zero"],
        )[inner]
        change, sensors = tile_fi_change(fc, year, transform, size, land, predictors)
        functional = functional_with_change(fi_state, change)
        bii = tile_bii(year, transform, size)
        compositional = compositional_integrity(bii, land)
        eii = combine_components(functional, structural, compositional)
        limiting = limiting_component(functional, structural, compositional).astype("float32")
        limiting[limiting == 0] = np.nan
        bands = {
            "eii": eii,
            "functional_integrity": functional,
            "structural_integrity": structural,
            "compositional_integrity": compositional,
            "limiting_component": limiting,
            "wetland_flag": flag,
            "fi_state": fi_state,
            "fi_change": change,
            "fi_change_sensors": sensors,
        }
        out.parent.mkdir(parents=True, exist_ok=True)
        write(out, bands, transform, cfg["tags"] | cfg["year_tags"].get(year, {})
              | {"reference_year": year})  # fmt: skip
        rows.append(
            {
                "tile": tile,
                "year": year,
                "land_px": int(land.sum()),
                **{f"valid_{k}": int(np.isfinite(v).sum()) for k, v in bands.items()},
                **{
                    f"mean_{k}": float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
                    for k, v in bands.items()
                    if k.endswith(("eii", "integrity", "fi_state", "fi_change"))
                },
                "seconds": round(time.perf_counter() - t0, 1),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["factors", "run"])
    parser.add_argument("potential", type=Path, help="Root of predict_potential_npp_v2.py")
    parser.add_argument("out", type=Path)
    parser.add_argument(
        "--clms-dir", type=Path, required=True, help="annual_npp_YEAR.tif directory"
    )
    parser.add_argument(
        "--hmi", type=Path, required=True, help="GHM 2022 GeoTIFF (downloaded if missing)"
    )
    parser.add_argument(
        "--ecoregions", type=Path, required=True, help="RESOLVE ecoregion raster (ECO_ID)"
    )
    parser.add_argument("--glwd-class", type=Path)
    parser.add_argument("--glwd-cover", type=Path)
    parser.add_argument("--natural-cv", type=Path, help="natural_cv.vrt (inter-annual CV)")
    parser.add_argument(
        "--deficit-k",
        dest="tolerance_k",
        type=float,
        default=TOLERANCE_K,
        help="Deficit tolerance multiplier",
    )
    parser.add_argument(
        "--surplus-k",
        dest="tolerance_k_surplus",
        type=float,
        default=TOLERANCE_K_SURPLUS,
        help="Surplus tolerance multiplier",
    )
    parser.add_argument(
        "--surplus-zero",
        type=float,
        default=SURPLUS_ZERO,
        help="Relative surplus (A - P) / max(P, floor) at which functional integrity reaches 0",
    )
    parser.add_argument("--years", type=int, nargs="*", default=list(REFERENCE_YEARS))
    parser.add_argument("--tiles", nargs="*", help="Only these tile names (e.g. e000n40)")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--fi-change", type=Path, help="Output directory of fi_change_term.py")
    parser.add_argument("--predictor-stack", type=Path, help="stack.vrt of the predictor stack")
    parser.add_argument(
        "--billing-project",
        help="Google Cloud project billed for the requester-pays BII reads "
        "(default: EII_BILLING_PROJECT)",
    )
    args = parser.parse_args()
    if args.stage == "run":
        required = ["glwd_class", "glwd_cover", "natural_cv", "fi_change", "predictor_stack"]
        missing = [f"--{k.replace('_', '-')}" for k in required if getattr(args, k) is None]
        if missing:
            parser.error(f"run needs {' '.join(missing)}")
    if args.billing_project:
        os.environ["EII_BILLING_PROJECT"] = args.billing_project

    gcs_env()
    args.out.mkdir(parents=True, exist_ok=True)
    if not args.hmi.exists():
        record, name, md5 = GHM_2022
        download(ZENODO_URL.format(record=record, name=name), args.hmi, md5)
    needed = sorted({y for year in args.years for y in observed_npp_window(year)})
    cfg = {
        "potential": args.potential,
        "out": args.out,
        "clms": {
            y: clms_vrt(y, args.out / "vrt" / f"clms_npp_{y}.vrt", args.clms_dir) for y in needed
        },
        "natural_cv": args.natural_cv,
        "tolerance_k": args.tolerance_k,
        "tolerance_k_surplus": args.tolerance_k_surplus,
        "surplus_zero": args.surplus_zero,
        "glwd_class": args.glwd_class,
        "glwd_cover": args.glwd_cover,
        "hmi": args.hmi,
        "ecoregions": args.ecoregions,
        "eco_to_biome": eco_to_biome_lookup(),
        "factors": {},
        "year_tags": {},
    }
    tiles = sorted(p.stem for p in (args.potential / "tiles").glob("*.tif"))
    if args.tiles:
        tiles = [t for t in tiles if t in args.tiles]

    missing = [y for y in args.years if not factors_path(args.out, y).exists()]
    if args.stage == "factors":
        missing = list(args.years)
    if missing:
        print(json.dumps({"factors": missing, "tiles": len(tiles)}), flush=True)
        compute_factors(missing, tiles, cfg, args.workers)
    if args.stage == "factors":
        return
    records = {}
    for year in args.years:
        cfg["factors"][year], records[year] = load_factors(cfg, year, tiles)
        cfg["year_tags"][year] = {"recentring_factors": json.dumps(factor_summary(records[year]))}

    cfg["tags"] = {
        "hmi": GHM_2022[1],
        "structural_quality_midpoint": QUALITY_MIDPOINT,
        "bii": "datasets/bii/v1_1",
        "npp": NPP_SOURCE,
        "tolerance_k": args.tolerance_k,
        "tolerance_k_surplus": args.tolerance_k_surplus,
        "surplus_zero": args.surplus_zero,
        "potential_npp_floor": POTENTIAL_NPP_FLOOR,
        "natural_cv": "interannual_cv model",
        "recentring": "hmi_bii",
        "recentring_min_reference_pixels": RECENTRE_MIN_REFERENCE_PIXELS,
        "glwd": "GLWD v2.0",
    }
    fc = fi_change_config(args.fi_change, args.predictor_stack, args.years, args.workers)
    cfg["fi_change"] = fc
    cfg["tags"]["fi_change_constants_sha256"] = json.dumps(fc["sha256"])
    for year, series in fc["series"].items():
        cfg["year_tags"].setdefault(year, {}).update(
            fi_change_series=series or "none",
            fi_change_k=FI_CHANGE_K[series] if series else "none",
        )
    manifest = {
        "tags": cfg["tags"],
        "inputs": {
            "potential": str(args.potential),
            "clms_annual": str(args.clms_dir),
            "natural_cv": str(args.natural_cv),
            "hmi": str(args.hmi),
            "glwd_class": str(args.glwd_class),
            "glwd_cover": str(args.glwd_cover),
            "ecoregions": str(args.ecoregions),
        },
        "years": args.years,
        "tiles": len(tiles),
        "fi_change": {
            "dir": str(args.fi_change),
            "predictor_stack": str(args.predictor_stack),
            "series": {str(y): c for y, c in fc["series"].items()},
            "constants_sha256": fc["sha256"],
        },
        "recentring": {
            str(y): {"file": str(factors_path(args.out, y)), "factors": factor_summary(r)}
            for y, r in records.items()
        },
    }
    manifest_path = args.out / "run_manifest.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        manifest["years"] = sorted(set(previous.get("years", [])) | set(args.years))
        manifest["recentring"] = previous.get("recentring", {}) | manifest["recentring"]
    manifest_path.write_text(json.dumps(manifest, indent=1) + "\n")

    todo = [
        t for t in tiles if any(not (args.out / str(y) / f"{t}.tif").exists() for y in args.years)
    ]
    print(json.dumps({"tiles": len(tiles), "to_run": len(todo), "years": args.years}), flush=True)
    qa_path = args.out / "qa.csv"
    with ProcessPoolExecutor(args.workers) as pool:
        futures = {pool.submit(run_tile, t, args.years, cfg): t for t in todo}
        for i, future in enumerate(as_completed(futures), start=1):
            for row in future.result():
                new = not qa_path.exists()
                with qa_path.open("a", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=list(row))
                    if new:
                        writer.writeheader()
                    writer.writerow(row)
            print(f"{i}/{len(todo)} {futures[future]}", flush=True)


if __name__ == "__main__":
    main()
