"""External validation of EII v2 against independent layers and the published WCMC EII.

Stages (``python pipelines/validation/run_validation.py <stage>``):

- ``fetch-wcmc``: download the published WCMC EII (``EII_V1.zip``), verify its checksums and
  extract the yearly 1 km rasters into ``inputs/wcmc``.
- ``warp``: every external layer onto the Eckert IV 1 km sampling grid (``E``) and our EII
  onto the published WCMC EPSG:4326 1 km grid (``W``).
- ``points``: 500,000 equal-area random land points (seed 42) sampled from all ``E`` rasters.
- ``analyse``: correlations and area summaries on ``E`` and the pixel comparison with the
  WCMC EII on ``W``; CSV tables to ``<run dir>/outputs/tables``.
- ``figures``: PNG figures from the tables to ``<run dir>/outputs/figures``.
- ``all``: ``warp``, ``points``, ``analyse`` and ``figures`` in order.

Paths come from environment variables (relative defaults resolve against the working
directory):

- ``EII_VALIDATION_DATA_ROOT`` (default ``./validation_data``): external inputs under
  ``inputs/`` and the external layers warped to both grids under ``work/``, shared by runs.
- ``EII_VALIDATION_EII_ROOT`` (default ``./eii_v2``): our EII v2 product as written by
  ``pipelines/production/compute_eii_v2.py`` (``<year>/<tile>.tif``).
- ``EII_VALIDATION_ECOREGIONS`` (default ``<data root>/inputs/layers/ecoregions/
  resolve_ecoregions.tif``): RESOLVE ecoregions rasterized to ecoregion ids.
- ``EII_VALIDATION_RUN_DIR`` (default ``<data root>/runs/<name of the EII root>``): our EII on
  both grids, the point table, the log and the outputs of this run.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

import rasterio
from rasterio.enums import Resampling

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from validation.grid import eckert4_profile  # noqa: E402
from validation.warping import WarpJob, grid_of, run_warp  # noqa: E402

DATA = Path(os.environ.get("EII_VALIDATION_DATA_ROOT", "validation_data"))
EII_ROOT = Path(os.environ.get("EII_VALIDATION_EII_ROOT", "eii_v2"))
ECOREGIONS = Path(
    os.environ.get(
        "EII_VALIDATION_ECOREGIONS", DATA / "inputs/layers/ecoregions/resolve_ecoregions.tif"
    )
)
RUN_DIR = Path(os.environ.get("EII_VALIDATION_RUN_DIR", DATA / "runs" / EII_ROOT.name))
LAYERS = DATA / "inputs/layers"
WCMC_DIR = DATA / "inputs/wcmc"
E_LAYER_DIR = DATA / "work/e_grid/layers"
W_LAYER_DIR = DATA / "work/w_grid/layers"
OUT = RUN_DIR / "outputs"
LOG = RUN_DIR / "validation.log"
SEED = 42
PRODUCT = "EII v2.0.0"
COLLECTION = "projects/landler-open-data/assets/products/eii/v2_0/release/global"
N_POINTS = 500_000

# Our side: 2021 is the year nearest to the latest WCMC map (2020), 2025 the latest year.
YEARS = (2021, 2025)
COMPARISON_YEAR = 2021
RUNS: dict[str, int] = {f"eii_{year}": year for year in YEARS}
MIN_TILES = 342  # land tiles per reference year
BANDS = ("eii", "fi", "si", "ci")  # first four bands of every product tile

CK = LAYERS / "chaplin-kramer/extracted"
CGLS = LAYERS / "cgls_lc100_2015/PROBAV_LC100_global_v3.0.1_2015-base"

# name -> source, resampling and value handling on the E grid.
E_LAYERS: dict[str, dict[str, Any]] = {
    "wcmc_eii_2020": {"src": WCMC_DIR / "EII_2020_1km.tif"},
    "hii_2013": {
        "src": LAYERS / "hfi_2013/hii_2013-01-01.tif",
        "scale": 0.01,
        "valid": (0.0, 64.0),
        "note": "WCS Human Influence Index 2013 (stored x100), the human footprint layer",
    },
    "anthromes_2000": {
        "src": LAYERS / "anthromes/anthro2_a2000.tif",
        "resampling": Resampling.nearest,
        "src_nodata": 0,
    },
    "lia_very_low": {
        "src": LAYERS / "low_impact_areas/Very_Low_Impact.tif",
        "resampling": Resampling.nearest,
        "src_nodata": -128,
    },
    "lia_low": {
        "src": LAYERS / "low_impact_areas/Low_Impact.tif",
        "resampling": Resampling.nearest,
        "src_nodata": -128,
    },
    "bhi_hab_2015": {
        "src": LAYERS / "bhi/data/BHI_v2_Proportion_of_Habitat/BILBI_P_BHIv2_Habitat_2015.tif"
    },
    "bhi_spe_2015": {
        "src": LAYERS / "bhi/data/BHI_v2_Proportion_of_Species/BILBI_P_BHIv2_Species_2015.tif"
    },
    "msa_2015": {"src": LAYERS / "msa/TerrestrialMSA_2015_World.tif", "valid": (0.0, 1.0)},
    "intactness_2013": {
        "src": LAYERS / "ecoregion_intactness/Intactness_2013.tif",
        "scale": 0.001,
        "valid": (0.0, 1.0),
    },
    "pollination_potential": {
        "src": CK
        / "pollination_potential_10s_cur_compressed_md5_20e27265b1627b6fbdb10b5fa504e745.tif"
    },
    "nutrient_potential": {
        "src": CK / "nutrient_potential_10s_cur_compressed_md5_9e0ae4df4e399350087c1c2d00d0a1bb.tif"
    },
    "n_retention": {
        "src": CK / "worldclim_esa_2015_n_retention_md5_d10a396b8f0fee70dd3bbd3524a6a97c.tif"
    },
    "ndr_2015": {
        "src": CK / "2015_aligned_ndr_compressed_md5_52442a88444f8f7fa648f0014b9d3a63.tif"
    },
    "gsoc": {"src": LAYERS / "gsoc/GSOCSEQ.GSOCMAP1-5-0.tif", "valid": (0.0, 1e4)},
    "defaunation_2019": {
        "src": LAYERS / "defaunation/DefInd.tif",
        "src_crs": "ESRI:54009",
        "src_nodata": -1.7e308,
        "valid": (0.0, 1.0),
    },
    "crops_2015": {
        "src": Path(f"{CGLS}_Crops-CoverFraction-layer_EPSG-4326.tif"),
        "valid": (0.0, 100.0),
    },
    "builtup_2015": {
        "src": Path(f"{CGLS}_BuiltUp-CoverFraction-layer_EPSG-4326.tif"),
        "valid": (0.0, 100.0),
    },
    "ecoregion": {
        "src": ECOREGIONS,
        "resampling": Resampling.nearest,
    },
}
W_LAYERS: dict[str, dict[str, Any]] = {
    "ecoregion": E_LAYERS["ecoregion"],
    "sayre": {
        "src": LAYERS / "sayre_wte/unpacked/commondata/raster_data/WorldEcosystem.tif",
        "resampling": Resampling.mode,
    },
}


def log(msg: str) -> None:
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as f:
        f.write(line + "\n")


def e_grid() -> dict[str, Any]:
    p = eckert4_profile(1000.0)
    return {k: p[k] for k in ("crs", "transform", "width", "height")}


def w_grid() -> dict[str, Any]:
    return grid_of(WCMC_DIR / "EII_2020_1km.tif")


def our_e(name: str, band: str) -> Path:
    return RUN_DIR / "e_grid" / name / f"{band}.tif"


def our_w(name: str) -> Path:
    return RUN_DIR / "w_grid" / f"{name}_eii.tif"


def tiles(year: int) -> list[Path]:
    return sorted((EII_ROOT / str(year)).glob("*.tif"))


def mosaic_vrt(files: list[Path], out: Path) -> Path:
    """GDAL VRT mosaic of same-grid, same-band GeoTIFF tiles (no copying)."""
    metas = []
    for path in files:
        with rasterio.open(path) as src:
            metas.append((path, src.transform, src.width, src.height))
            count, dtype, nodata, crs = src.count, src.dtypes[0], src.nodata, src.crs
    pixel = metas[0][1].a
    left = min(t.c for _, t, *_ in metas)
    top = max(t.f for _, t, *_ in metas)
    right = max(t.c + w * pixel for _, t, w, _ in metas)
    bottom = min(t.f - h * pixel for _, t, _, h in metas)
    width, height = round((right - left) / pixel), round((top - bottom) / pixel)
    gdal_type = {"float32": "Float32", "uint16": "UInt16"}[dtype]
    nd = f"<NoDataValue>{nodata}</NoDataValue>" if nodata is not None else ""
    bands = []
    for b in range(1, count + 1):
        sources = "".join(
            f'<SimpleSource><SourceFilename relativeToVRT="0">{escape(str(path))}'
            f"</SourceFilename><SourceBand>{b}</SourceBand>"
            f'<SrcRect xOff="0" yOff="0" xSize="{w}" ySize="{h}"/>'
            f'<DstRect xOff="{round((t.c - left) / pixel)}" yOff="{round((top - t.f) / pixel)}" '
            f'xSize="{w}" ySize="{h}"/></SimpleSource>'
            for path, t, w, h in metas
        )
        bands.append(
            f'<VRTRasterBand dataType="{gdal_type}" band="{b}">{nd}{sources}</VRTRasterBand>'
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        f'<VRTDataset rasterXSize="{width}" rasterYSize="{height}">'
        f"<SRS>{escape(crs.to_wkt())}</SRS>"
        f"<GeoTransform>{left}, {pixel}, 0, {top}, 0, {-pixel}</GeoTransform>"
        f"{''.join(bands)}</VRTDataset>\n"
    )
    return out


def run_vrt(name: str) -> Path:
    vrt = RUN_DIR / "our_eii" / f"{name}.vrt"
    if not vrt.exists():
        mosaic_vrt(tiles(RUNS[name]), vrt)
    return vrt


def complete(year: int) -> bool:
    """All land tiles of a reference year exist and none is still being written."""
    found = tiles(year)
    return len(found) >= MIN_TILES and all(settled(t) for t in found)


def settled(path: Path, age_s: float = 120.0) -> bool:
    """File exists and has not been modified for ``age_s``."""
    return path.exists() and time.time() - path.stat().st_mtime > age_s


def layer_job(name: str, spec: dict[str, Any], out_dir: Path) -> WarpJob:
    return WarpJob(
        src=spec["src"],
        bands=(1,),
        outs=(out_dir / f"{name}.tif",),
        resampling=spec.get("resampling", Resampling.average),
        src_crs=spec.get("src_crs"),
        src_nodata=spec.get("src_nodata"),
        scale=spec.get("scale", 1.0),
        valid=spec.get("valid"),
    )


def jobs(only: list[str] | None = None) -> list[tuple[str, WarpJob, dict[str, Any]]]:
    """All warp jobs as (label, job, grid)."""
    eg, wg = e_grid(), w_grid()
    out: list[tuple[str, WarpJob, dict[str, Any]]] = []
    for name, year in RUNS.items():
        if not complete(year):
            continue
        tags = {"run": name, "reference_year": str(year), "product": PRODUCT}
        e_outs = tuple(our_e(name, band) for band in BANDS)
        out.append((f"E:{name}", WarpJob(run_vrt(name), (1, 2, 3, 4), e_outs, src_nodata=-9999.0,
                                         tags=tags), eg))  # fmt: skip
        out.append((f"W:{name}", WarpJob(run_vrt(name), (1,), (our_w(name),), src_nodata=-9999.0,
                                         tags=tags), wg))  # fmt: skip
    for name, spec in E_LAYERS.items():
        if settled(spec["src"]):
            out.append((f"E:{name}", layer_job(name, spec, E_LAYER_DIR), eg))
    for name, spec in W_LAYERS.items():
        if settled(spec["src"]):
            out.append((f"W:{name}", layer_job(name, spec, W_LAYER_DIR), wg))
    if only:
        out = [j for j in out if any(o in j[0] for o in only)]
    return out


def stage_warp(only: list[str] | None, workers: int) -> None:
    selected = jobs(only)
    unmatched = [o for o in only or [] if not any(o in label for label, _, _ in selected)]
    if unmatched:
        raise RuntimeError(f"no warp job for {unmatched}: inputs missing or still being written")
    for label, job, grid in selected:
        if job.done():
            continue
        t0 = time.time()
        log(f"warp {label} <- {job.src}")
        run_warp(job, grid, workers=workers, log=log)
        log(f"warp {label} done in {time.time() - t0:.0f} s")


def write_manifest() -> Path:
    """Record the product, years and settings of this run next to its outputs."""
    with rasterio.open(tiles(COMPARISON_YEAR)[0]) as src:
        skip = ("AREA_OR_POINT", "reference_year")
        # Absolute paths are local to the machine that wrote the tiles.
        tags = {k: v for k, v in src.tags().items() if k not in skip and not v.startswith("/")}
    manifest = {
        "product": PRODUCT,
        "collection": COLLECTION,
        "collection_note": "one image per 10 degree tile and reference year "
        "(property reference_year)",
        "reference_years": list(YEARS),
        "comparison_year": COMPARISON_YEAR,
        "wcmc_year": 2020,
        "seed": SEED,
        "n_points": N_POINTS,
        "sampling_grid": "Eckert IV (ESRI:54012) 1 km, area average of valid source pixels",
        "comparison_grid": "WCMC EII_2020_1km.tif grid (EPSG:4326), cos(latitude) area weights",
        "product_tags": tags,
    }
    path = OUT / "run_manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=1) + "\n")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "stage", choices=["fetch-wcmc", "warp", "points", "analyse", "figures", "all", "list"]
    )
    parser.add_argument("--only", nargs="*", help="substring filter on warp job labels")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--parts", nargs="*", default=["points", "a5", "a6", "wcmc"],
                        help="analyse parts: points, a5, a6, wcmc")  # fmt: skip
    args = parser.parse_args()
    stages = ["warp", "points", "analyse", "figures"] if args.stage == "all" else [args.stage]
    for stage in stages:
        if stage == "list":
            for label, job, _ in jobs(args.only):
                print(("done " if job.done() else "todo ") + label, job.src)
        elif stage == "fetch-wcmc":
            from validation.wcmc import fetch_and_verify

            print(json.dumps(fetch_and_verify(WCMC_DIR), indent=1))
        elif stage == "warp":
            stage_warp(args.only, args.workers)
        elif stage == "points":
            from validation.analysis import stage_points

            stage_points()
        elif stage == "analyse":
            from validation.analysis import stage_analyse

            write_manifest()
            stage_analyse(args.parts)
        elif stage == "figures":
            from validation.figures import main as figures

            figures()
        log(f"stage {stage} ok")


if __name__ == "__main__":
    main()
