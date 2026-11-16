#!/usr/bin/env python3
"""Global v2 potential NPP with prediction intervals, per 10 deg land tile of the predictor stack.

Model and calibration are those of ``fit_potential_npp_v2.py`` (``model_dir``):

- potential: LightGBM Tweedie mean model times the global bias correction ``c``;
- intervals: asymmetric QRF-scaled design ``norm_asym``,
  ``[P_c - k_b (q50 - q16), P_c + k_b (q84 - q50)]``, where k_b is the per-biome split-conformal
  factor for the 68 % or 90 % level (pooled factor where the biome is unknown) and the lower
  bound is clipped at 0.

Biome is the RESOLVE ecoregion raster that the training rows were annotated with
(``datasets/ecoregions/v1_0``, ~5 km, as synced by ``sync_predictor_sources.py``), remapped by
``pipelines/preprocessing/reference/eco_id_to_biome.json`` and taken by nearest neighbour.
Pixels with any predictor missing stay nodata. Predictors are read from the stack by band name in
model-card order, so the stack may hold extra bands.

``--natural-cv`` predicts the natural inter-annual CV of NPP instead (``interannual_cv_model`` in
the model card: LightGBM times its bias correction) into ``<out_root>/natural_cv/<tile>.tif``
(one float32 band ``natural_cv``) and ``natural_cv.vrt``.

Output: one GeoTIFF per tile on the CLMS NPP grid with float32 bands ``potential_npp``,
``lower_68``, ``upper_68``, ``lower_90``, ``upper_90`` in gC m-2 yr-1, nodata -9999, plus
``potential_npp.vrt``. Existing tiles are skipped; ``qa.csv`` gets one row per tile.

Example:
    python pipelines/training/predict_potential_npp_v2.py predictor_stack \\
        models/potential_npp_v2 potential_npp_v2
    python pipelines/training/predict_potential_npp_v2.py predictor_stack \\
        models/potential_npp_v2 potential_npp_v2 --natural-cv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import joblib
import lightgbm as lgb
import numpy as np
import rasterio
from rasterio.enums import Resampling

PREPROCESSING = Path(__file__).resolve().parents[1] / "preprocessing"
sys.path.insert(0, str(PREPROCESSING))
from build_predictor_stack import NODATA, PROFILE, mosaic_vrt, read_window, warp  # noqa: E402

from eii.compute.settings import BIOME_NAMES  # noqa: E402

ECO_TO_BIOME = PREPROCESSING / "reference" / "eco_id_to_biome.json"
BANDS = ["potential_npp", "lower_68", "upper_68", "lower_90", "upper_90"]
DESIGN = "norm_asym"


def biome_k_lookup(k_by_name: dict[str, float]) -> np.ndarray:
    """Factor per biome number 0-255; unknown biomes take the pooled factor."""
    lookup = np.full(256, k_by_name["__pooled__"], dtype="float32")
    for number, name in BIOME_NAMES.items():
        if name in k_by_name:
            lookup[number] = k_by_name[name]
    return lookup


def band_indexes(stack_root: Path, predictors: list[str]) -> list[int]:
    """1-based stack band index of each predictor, in model-card order."""
    with rasterio.open(stack_root / "stack.vrt") as src:
        names = list(src.descriptions)
    missing = [p for p in predictors if p not in names]
    if missing:
        raise ValueError(f"predictors missing from the stack: {missing}")
    return [names.index(p) + 1 for p in predictors]


def predict_tile(tile_path: Path, out: Path, models: dict, chunk: int) -> dict:
    with rasterio.open(tile_path) as src:
        stack = src.read(models["bands"]).astype("float64")
        transform, size = src.transform, src.width
        valid = (stack != src.nodata).all(axis=0)
    x = stack[:, valid].T
    n = len(x)
    p = np.empty(n)
    q = np.empty((n, 3))
    for start in range(0, n, chunk):
        part = x[start : start + chunk]
        p[start : start + chunk] = models["lgbm"].predict(part)
        q[start : start + chunk] = models["qrf"].predict(part, quantiles=[0.16, 0.5, 0.84])
    p_c = models["c"] * p
    lo_s = np.maximum(q[:, 1] - q[:, 0], models["eps"])
    hi_s = np.maximum(q[:, 2] - q[:, 1], models["eps"])

    bounds = rasterio.transform.array_bounds(size, size, transform)
    eco, eco_transform, _ = read_window(models["ecoregion"], bounds)
    eco = warp(eco, eco_transform, transform, size, Resampling.nearest)[0]
    eco = np.nan_to_num(eco, nan=0).astype("int64")[valid]
    known = (eco >= 0) & (eco < len(models["eco_to_biome"]))
    biome = np.where(known, models["eco_to_biome"][np.where(known, eco, 0)], 0)

    bands = {"potential_npp": p_c}
    for level in ("68", "90"):
        k = models["k"][level][biome]
        bands[f"lower_{level}"] = np.maximum(p_c - k * lo_s, 0)
        bands[f"upper_{level}"] = p_c + k * hi_s
    image = np.full((len(BANDS), size, size), NODATA, dtype="float32")
    for i, name in enumerate(BANDS):
        image[i][valid] = bands[name]

    part = out.with_name(out.name + ".part")
    with rasterio.open(
        part, "w", width=size, height=size, count=len(BANDS), transform=transform, **PROFILE
    ) as dst:
        dst.write(image)
        dst.descriptions = tuple(BANDS)
        dst.update_tags(units="gC m-2 yr-1", design=DESIGN, c=str(models["c"]))
    part.rename(out)
    return {
        "tile": tile_path.stem,
        "pixels": n,
        "potential_mean": float(p_c.mean()) if n else np.nan,
        "rel_half_width_68_median": float(
            np.median((bands["upper_68"] - bands["lower_68"])[p_c >= 50] / 2 / p_c[p_c >= 50])
        )
        if (p_c >= 50).any()
        else np.nan,
        "unknown_biome_share": float((biome == 0).mean()) if n else np.nan,
    }


def predict_cv_tile(
    tile_path: Path, out: Path, model: lgb.Booster, c: float, bands: list[int], chunk: int
) -> dict:
    with rasterio.open(tile_path) as src:
        stack = src.read(bands).astype("float64")
        transform, size = src.transform, src.width
        valid = (stack != src.nodata).all(axis=0)
    x = stack[:, valid].T
    cv = np.empty(len(x))
    for start in range(0, len(x), chunk):
        cv[start : start + chunk] = c * model.predict(x[start : start + chunk])
    image = np.full((1, size, size), NODATA, dtype="float32")
    image[0][valid] = np.maximum(cv, 0)
    part = out.with_name(out.name + ".part")
    with rasterio.open(
        part, "w", width=size, height=size, count=1, transform=transform, **PROFILE
    ) as dst:
        dst.write(image)
        dst.descriptions = ("natural_cv",)
        dst.update_tags(units="1", model="interannual_cv_model", c=str(c))
    part.rename(out)
    return {
        "tile": tile_path.stem,
        "pixels": int(valid.sum()),
        "cv_median": float(np.median(cv)) if len(cv) else np.nan,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stack_root", type=Path)
    parser.add_argument("model_dir", type=Path)
    parser.add_argument("out_root", type=Path)
    parser.add_argument("--chunk", type=int, default=250_000)
    parser.add_argument("--tiles", nargs="*", help="Only these tile names (e.g. e000n40)")
    parser.add_argument("--natural-cv", action="store_true", help="Predict natural_cv only")
    args = parser.parse_args()

    card = json.loads((args.model_dir / "model_card.json").read_text())
    if args.natural_cv:
        run_natural_cv(args, card)
        return
    calibration = card["calibration"]
    eco_lookup = {int(k): v for k, v in json.loads(ECO_TO_BIOME.read_text()).items()}
    eco_to_biome = np.zeros(max(eco_lookup) + 1, dtype="int64")
    for eco_id, biome_number in eco_lookup.items():
        eco_to_biome[eco_id] = biome_number
    ecoregion = sorted((args.stack_root / "ecoregions").rglob("*.tif"))
    if len(ecoregion) != 1:
        raise FileNotFoundError(f"expected one ecoregion raster under {args.stack_root}/ecoregions")
    models = {
        "lgbm": lgb.Booster(model_file=str(args.model_dir / card["mean_model"]["file"])),
        "qrf": joblib.load(args.model_dir / card["interval_model"]["file"]),
        "c": float(calibration["c"]),
        "eps": float(calibration["eps"]),
        "k": {lvl: biome_k_lookup(k) for lvl, k in calibration["designs"][DESIGN].items()},
        "eco_to_biome": eco_to_biome,
        "ecoregion": ecoregion[0],
        "bands": band_indexes(args.stack_root, card["predictors"]),
    }

    tiles_dir = args.out_root / "tiles"
    tiles_dir.mkdir(parents=True, exist_ok=True)
    tiles = sorted((args.stack_root / "stack").glob("*.tif"))
    if args.tiles:
        tiles = [t for t in tiles if t.stem in args.tiles]
    todo = [t for t in tiles if not (tiles_dir / t.name).exists()]
    print(json.dumps({"tiles": len(tiles), "to_predict": len(todo)}), flush=True)
    qa_path = args.out_root / "qa.csv"
    for i, tile in enumerate(todo, start=1):
        t0 = time.perf_counter()
        qa = predict_tile(tile, tiles_dir / tile.name, models, args.chunk)
        qa["seconds"] = round(time.perf_counter() - t0, 1)
        new = not qa_path.exists()
        with qa_path.open("a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(qa))
            if new:
                writer.writeheader()
            writer.writerow(qa)
        print(f"{i}/{len(todo)} {json.dumps(qa)}", flush=True)
    outputs = sorted(tiles_dir.glob("*.tif"))
    if outputs:
        mosaic_vrt(outputs, args.out_root / "potential_npp.vrt")
        print(f"potential_npp.vrt: {len(outputs)} tiles", flush=True)


def run_natural_cv(args: argparse.Namespace, card: dict) -> None:
    spec = card["interannual_cv_model"]
    bands = band_indexes(args.stack_root, spec["predictors"])
    model = lgb.Booster(model_file=str(args.model_dir / spec["file"]))
    out_dir = args.out_root / "natural_cv"
    out_dir.mkdir(parents=True, exist_ok=True)
    tiles = sorted((args.stack_root / "stack").glob("*.tif"))
    if args.tiles:
        tiles = [t for t in tiles if t.stem in args.tiles]
    todo = [t for t in tiles if not (out_dir / t.name).exists()]
    print(json.dumps({"tiles": len(tiles), "to_predict": len(todo)}), flush=True)
    for i, tile in enumerate(todo, start=1):
        qa = predict_cv_tile(tile, out_dir / tile.name, model, float(spec["c"]), bands, args.chunk)
        print(f"{i}/{len(todo)} {json.dumps(qa)}", flush=True)
    outputs = sorted(out_dir.glob("*.tif"))
    if outputs:
        mosaic_vrt(outputs, args.out_root / "natural_cv.vrt")
        print(f"natural_cv.vrt: {len(outputs)} tiles", flush=True)


if __name__ == "__main__":
    main()
