"""Analyses of the external validation: correlations, area summaries and the WCMC comparison.

Inputs are the rasters written by ``run_validation.py warp``: every layer on the Eckert IV
1 km grid (``E``) and our EII on the published WCMC 1 km EPSG:4326 grid (``W``). The
analyses follow the cross-indicator validation in the supplement of Hill et al. (2026),
section A2, so our numbers sit next to theirs for the same layers.

- Correlations (their Figs A3, A4, A8, A9) use 500,000 equal-area random land points (seed
  42) on ``E``, Pearson r with p plus Spearman rho, pairwise-complete per pair.
- Area summaries (Figs A5, A6, Locke et al. conditions) use the full ``E`` grid (1 km^2 cells).
- The comparison with the WCMC EII uses the full ``W`` grid with cos(latitude) cell weights.
"""

from __future__ import annotations

import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyproj
import rasterio
from rasterio.windows import Window
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from validation.metrics import locke_condition, pearson_spearman  # noqa: E402
from validation.run_validation import (  # noqa: E402
    COMPARISON_YEAR,
    E_LAYER_DIR,
    N_POINTS,
    OUT,
    RUN_DIR,
    RUNS,
    SEED,
    W_LAYER_DIR,
    WCMC_DIR,
    e_grid,
    log,
    our_e,
    our_w,
)
from validation.warping import cos_lat_weights, grid_of, read_full  # noqa: E402

TABLES = OUT / "tables"
FIGS = OUT / "figures"
POINTS = RUN_DIR / "points_e.csv.gz"
ECO_TO_BIOME = Path(__file__).resolve().parents[1] / "preprocessing/reference/eco_id_to_biome.json"
COMPARISON_RUN = f"eii_{COMPARISON_YEAR}"

BIOME_NAMES: dict[int, str] = {
    1: "Tropical moist broadleaf",
    2: "Tropical dry broadleaf",
    3: "Tropical conifer",
    4: "Temperate broadleaf",
    5: "Temperate conifer",
    6: "Boreal",
    7: "Savanna",
    8: "Temperate grassland",
    9: "Flooded grassland",
    10: "Montane grassland",
    11: "Tundra",
    12: "Mediterranean",
    13: "Desert and xeric",
    14: "Mangrove",
}

# Hill et al. (2026) supplement, values read from Figs A3-A9 and the section A2.1 text.
PAPER = {
    "A3": {("S", "C"): 0.538, ("S", "F"): 0.099, ("C", "F"): 0.090},
    "A4_pairs": {
        ("hfp", "bhi_hab"): -0.658,
        ("hfp", "bhi_spe"): -0.667,
        ("hfp", "msa"): -0.589,
        ("hfp", "eii"): -0.793,
        ("bhi_hab", "bhi_spe"): 0.993,
        ("bhi_hab", "msa"): 0.691,
        ("bhi_hab", "eii"): 0.747,
        ("bhi_spe", "msa"): 0.677,
        ("bhi_spe", "eii"): 0.761,
        ("msa", "eii"): 0.657,
    },
    "A4_mean": {"hfp": 0.68, "bhi_hab": 0.699, "bhi_spe": 0.702, "eii": 0.7395, "msa": 0.6535},
    "A5": 0.835,
    "A8": {"pollination": -0.010, "nutrient": -0.058, "soc": 0.010, "n_retention": 0.008},
    # Hill's "BHI" is the species layer (by its value range). They invert the defaunation
    # index (1 - x) so that high = intact; DefInd.tif already stores that orientation
    # (retained fauna: mean 0.87 over the tropics, r = -0.71 with the human footprint), so
    # their value sits on the raw layer.
    "A9": {"bhi_spe": 0.517, "msa": 0.436, "retained_fauna": 0.233},
    "sayre": {"mean": 0.68, "sd": 0.14, "n": 422, "le_0.5": 87, "le_0.25": 25},
}

# Indicator columns of the point table (E grid layer name -> column).
A4_INDICATORS = {"hfp": "hii_2013", "bhi_hab": "bhi_hab_2015", "bhi_spe": "bhi_spe_2015",
                 "msa": "msa_2015"}  # fmt: skip
A8_LAYERS = {"pollination": "pollination_potential", "nutrient": "nutrient_potential",
             "soc": "gsoc", "n_retention": "n_retention", "ndr": "ndr_2015"}  # fmt: skip


def available_runs() -> list[str]:
    return [n for n in RUNS if our_e(n, "eii").exists()]


def random_points(n: int, seed: int) -> np.ndarray:
    """Equal-area random lon/lat pairs on the sphere (uniform in lon and sin(lat))."""
    rng = np.random.default_rng(seed)
    lon = rng.uniform(-180.0, 180.0, n)
    lat = np.degrees(np.arcsin(rng.uniform(-1.0, 1.0, n)))
    return np.column_stack([lon, lat])


# --------------------------------------------------------------------------------------
# Points
# --------------------------------------------------------------------------------------


def lonlat_to_cells(lon: np.ndarray, lat: np.ndarray, grid: dict[str, Any]) -> tuple[Any, Any]:
    """Row/col of lon/lat points on a projected grid."""
    tr = pyproj.Transformer.from_crs("EPSG:4326", grid["crs"].to_wkt(), always_xy=True)
    x, y = tr.transform(lon, lat)
    col, row = ~grid["transform"] * (np.asarray(x), np.asarray(y))
    return np.floor(row).astype("int64"), np.floor(col).astype("int64")


def draw_land_points(land: np.ndarray, grid: dict[str, Any], n: int, seed: int) -> pd.DataFrame:
    """First ``n`` equal-area random points (seed fixed) that fall on ``land`` cells."""
    pts = random_points(8 * n, seed)
    row, col = lonlat_to_cells(pts[:, 0], pts[:, 1], grid)
    inside = (row >= 0) & (row < land.shape[0]) & (col >= 0) & (col < land.shape[1])
    keep = np.zeros(len(pts), bool)
    keep[inside] = land[row[inside], col[inside]]
    idx = np.flatnonzero(keep)[:n]
    if len(idx) < n:
        raise RuntimeError(f"only {len(idx)} land points among {len(pts)} candidates")
    return pd.DataFrame({"lon": pts[idx, 0], "lat": pts[idx, 1], "row": row[idx], "col": col[idx]})


def stage_points() -> Path:
    """Sample every ``E`` raster at 500k equal-area land points."""
    grid = e_grid()
    land = np.isfinite(read_full(E_LAYER_DIR / "wcmc_eii_2020.tif"))
    land |= np.isfinite(read_full(our_e(COMPARISON_RUN, "si")))
    frame = draw_land_points(land, grid, N_POINTS, SEED)
    del land
    rows, cols = frame["row"].to_numpy(), frame["col"].to_numpy()
    for path in sorted(E_LAYER_DIR.glob("*.tif")):
        if ".partial" not in path.name:
            frame[path.stem] = read_full(path)[rows, cols]
    for name in available_runs():
        for band in ("eii", "fi", "si", "ci"):
            if our_e(name, band).exists():
                frame[f"{name}__{band}"] = read_full(our_e(name, band))[rows, cols]
    # RESOLVE id 0 (rock and ice) is also the raster's ocean fill, so it counts as missing.
    frame.loc[frame["ecoregion"] == 0, "ecoregion"] = np.nan
    to_biome = {int(k): int(v) for k, v in json.loads(ECO_TO_BIOME.read_text()).items()}
    frame["biome"] = frame["ecoregion"].map(lambda e: to_biome.get(int(e)) if e == e else None)
    POINTS.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(POINTS, index=False, compression="gzip")
    log(f"points: {len(frame)} rows, {frame.shape[1]} columns -> {POINTS}")
    return POINTS


# --------------------------------------------------------------------------------------
# Correlation helpers
# --------------------------------------------------------------------------------------


def corr_row(x: Any, y: Any, **labels: Any) -> dict[str, Any]:
    """Pearson/Spearman with p and n over pairwise-complete rows, plus labels."""
    return labels | pearson_spearman(x, y)


def mean_abs_r(r: dict[tuple[str, str], float], layer: str, exclude: Sequence[str] = ()) -> float:
    """Hill's 'mean correlation across all other indicators' (absolute r)."""
    vals = [
        v for (a, b), v in r.items() if layer in (a, b) and not ({a, b} - {layer}) & set(exclude)
    ]
    return float(np.mean(np.abs(vals))) if vals else float("nan")


def hill_means(r: dict[tuple[str, str], float]) -> dict[str, float]:
    """The five mean |r| values of section A2.1 (BHI-hab excludes BHI-spe and vice versa)."""
    return {
        "hfp": mean_abs_r(r, "hfp"),
        "bhi_hab": mean_abs_r(r, "bhi_hab", exclude=["bhi_spe"]),
        "bhi_spe": mean_abs_r(r, "bhi_spe", exclude=["bhi_hab"]),
        "msa": mean_abs_r(r, "msa"),
        "eii": mean_abs_r(r, "eii"),
    }


# --------------------------------------------------------------------------------------
# Fig A3: components
# --------------------------------------------------------------------------------------


def a3_components(pts: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name in available_runs():
        cols = {"S": f"{name}__si", "C": f"{name}__ci", "F": f"{name}__fi"}
        if not all(c in pts for c in cols.values()):
            continue
        for a, b in (("S", "C"), ("S", "F"), ("C", "F")):
            row = corr_row(pts[cols[a]], pts[cols[b]], run=name, pair=f"{a}-{b}")
            row["paper_r"] = PAPER["A3"][(a, b)]
            rows.append(row)
        # limiting component at 1 km (argmin of the aggregated components)
        trio = pts[[cols["F"], cols["S"], cols["C"]]].to_numpy()
        ok = np.isfinite(trio).all(axis=1)
        lim = np.argmin(trio[ok], axis=1)
        for k, comp_name in enumerate(("F", "S", "C")):
            rows.append(
                {"run": name, "pair": f"limiting={comp_name}", "share": float(np.mean(lim == k))}
            )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# Fig A4: EII vs four indicators
# --------------------------------------------------------------------------------------


def eii_versions() -> dict[str, str]:
    out = {"wcmc_2020": "wcmc_eii_2020"}
    out |= {name: f"{name}__eii" for name in available_runs()}
    return out


def a4_indicators(pts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pairwise r (Fig A4 layout) and Hill's mean |r| per EII version.

    Computed twice: on pairwise-complete rows (as GGally) and on the rows valid in all
    four indicators, WCMC and every run ("common"), so versions share one sample.
    """
    versions = eii_versions()
    common = pts[list(A4_INDICATORS.values()) + list(versions.values())].notna().all(axis=1)
    pair_rows, mean_rows = [], []
    for sample, mask in (("pairwise", pd.Series(True, index=pts.index)), ("common", common)):
        sub = pts[mask]
        ind_r: dict[tuple[str, str], float] = {}
        keys = list(A4_INDICATORS)
        for i, a in enumerate(keys):
            for b in keys[i + 1 :]:
                row = corr_row(sub[A4_INDICATORS[a]], sub[A4_INDICATORS[b]], sample=sample,
                               version="-", a=a, b=b)  # fmt: skip
                row["paper_r"] = PAPER["A4_pairs"][(a, b)]
                ind_r[(a, b)] = row["pearson_r"]
                pair_rows.append(row)
        for version, col in versions.items():
            r = dict(ind_r)
            for a in keys:
                row = corr_row(sub[A4_INDICATORS[a]], sub[col], sample=sample, version=version,
                               a=a, b="eii")  # fmt: skip
                row["paper_r"] = PAPER["A4_pairs"][(a, "eii")]
                r[(a, "eii")] = row["pearson_r"]
                pair_rows.append(row)
            means = hill_means(r)
            mean_rows.append(
                {"sample": sample, "version": version, "n_common": int(mask.sum())}
                | {f"mean_abs_r_{k}": v for k, v in means.items()}
                | {"eii_rank": int(1 + sum(means[k] > means["eii"] for k in means if k != "eii"))}
            )
    return pd.DataFrame(pair_rows), pd.DataFrame(mean_rows)


# --------------------------------------------------------------------------------------
# Fig A5: ecoregion means vs Beyer intactness
# --------------------------------------------------------------------------------------


def unit_means(ids: np.ndarray, values: np.ndarray, n_units: int) -> tuple[np.ndarray, np.ndarray]:
    """Mean of ``values`` and valid-cell count per integer unit id (NaN-aware)."""
    ok = np.isfinite(ids) & np.isfinite(values)
    idx = ids[ok].astype("int64")
    count = np.bincount(idx, minlength=n_units).astype("float64")
    total = np.bincount(idx, weights=values[ok].astype("float64"), minlength=n_units)
    with np.errstate(invalid="ignore", divide="ignore"):
        return total / count, count


def a5_ecoregions(min_cells: int = 100) -> tuple[pd.DataFrame, pd.DataFrame]:
    eco = read_full(E_LAYER_DIR / "ecoregion.tif")
    eco[eco == 0] = np.nan  # ocean fill (and rock and ice)
    n_units = int(np.nanmax(eco)) + 1
    intact, n_int = unit_means(eco, read_full(E_LAYER_DIR / "intactness_2013.tif"), n_units)
    frame = pd.DataFrame({"ecoregion": np.arange(n_units), "intactness": intact, "n_int": n_int})
    sources = {"wcmc_2020": E_LAYER_DIR / "wcmc_eii_2020.tif"}
    sources |= {name: our_e(name, "eii") for name in available_runs()}
    for version, path in sources.items():
        mean, count = unit_means(eco, read_full(path), n_units)
        frame[version] = np.where(count >= min_cells, mean, np.nan)
    frame = frame[frame["n_int"] >= min_cells]
    rows = [
        corr_row(frame["intactness"], frame[v], version=v, paper_r=PAPER["A5"], min_cells=min_cells)
        for v in sources
    ]
    return frame, pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# Fig A6: influence classes and Locke conditions (full E grid)
# --------------------------------------------------------------------------------------

CLASS_CODES = {"very_low": 0, "low": 1, "high": 2}
LOCKE_CODES = {"large_wild_areas": 0, "shared_lands": 1, "cities_farms": 2}


def influence_classes_anthromes(code: np.ndarray) -> np.ndarray:
    """Table A10, anthromes column: 61/62/63 very low, 43/53/54 low, other land high."""
    out = np.full(code.shape, -1, "int8")
    land = np.isfinite(code) & (code > 0) & (code < 70)
    out[land] = 2
    out[np.isin(code, (43, 53, 54))] = 1
    out[np.isin(code, (61, 62, 63))] = 0
    return out


def influence_classes_hf(hf: np.ndarray) -> np.ndarray:
    """Table A10, human footprint column: HF = 0 very low, 0 < HF <= 3 low, else high."""
    out = np.full(hf.shape, -1, "int8")
    ok = np.isfinite(hf)
    out[ok] = 2
    out[ok & (hf > 0) & (hf <= 3)] = 1
    out[ok & (hf == 0)] = 0
    return out


def influence_classes_lia(very_low: np.ndarray, low: np.ndarray) -> np.ndarray:
    """Table A10, Low Impact Areas column: very low wins over low (no overlap).

    The Jacobson et al. (2019) rasters code cells inside the (very) low impact area as 0 and
    outside as 100 (-128 is ocean): 0 in ``Low_Impact.tif`` covers 56 % of land, the share
    the paper reports for Low Impact Areas.
    """
    out = np.full(very_low.shape, -1, "int8")
    ok = np.isfinite(very_low) | np.isfinite(low)
    out[ok] = 2
    out[low == 0] = 1
    out[very_low == 0] = 0
    return out


def class_histograms(classes: np.ndarray, values: np.ndarray, n_classes: int = 3) -> np.ndarray:
    """Cell counts per class and 0.01 value bin, shape (n_classes, 100)."""
    ok = (classes >= 0) & np.isfinite(values)
    bins = np.clip((values[ok] * 100).astype("int64"), 0, 99)
    flat = np.bincount(classes[ok].astype("int64") * 100 + bins, minlength=n_classes * 100)
    return flat.reshape(n_classes, 100)


def auc_from_hist(high: np.ndarray, low: np.ndarray) -> float:
    """P(value in class ``high`` > value in class ``low``) from binned counts (ties 0.5)."""
    if high.sum() == 0 or low.sum() == 0:
        return float("nan")
    ph, pl = high / high.sum(), low / low.sum()
    below = np.concatenate([[0.0], np.cumsum(pl)[:-1]])
    return float(np.sum(ph * (below + 0.5 * pl)))


def hist_summary(hist: np.ndarray, names: Sequence[str]) -> list[dict[str, Any]]:
    centres = (np.arange(100) + 0.5) / 100
    rows = []
    for k, name in enumerate(names):
        h = hist[k]
        n = h.sum()
        cdf = np.cumsum(h) / n if n else np.full(100, np.nan)

        def q(p: float, cdf: np.ndarray = cdf, n: int = n) -> float:
            return float(centres[min(np.searchsorted(cdf, p), 99)]) if n else float("nan")

        rows.append(
            {
                "class": name,
                "area_km2": int(n),
                "mean": float((h * centres).sum() / n) if n else float("nan"),
                "p10": q(0.1),
                "median": q(0.5),
                "p90": q(0.9),
                "share_ge_0.9": float(h[90:].sum() / n) if n else float("nan"),
            }
        )
    return rows


def a6_influence() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Histograms + summaries of S and EII per influence class / Locke condition."""
    classes = {
        "anthromes": influence_classes_anthromes(read_full(E_LAYER_DIR / "anthromes_2000.tif")),
        "human_footprint": influence_classes_hf(read_full(E_LAYER_DIR / "hii_2013.tif")),
        "low_impact_areas": influence_classes_lia(
            read_full(E_LAYER_DIR / "lia_very_low.tif"), read_full(E_LAYER_DIR / "lia_low.tif")
        ),
    }
    hf = read_full(E_LAYER_DIR / "hii_2013.tif")
    intensive = (
        np.nan_to_num(read_full(E_LAYER_DIR / "crops_2015.tif"))
        + np.nan_to_num(read_full(E_LAYER_DIR / "builtup_2015.tif"))
    ) / 100.0
    locke = locke_condition(hf, intensive)
    locke_codes = np.full(hf.shape, -1, "int8")
    for name, code in LOCKE_CODES.items():
        locke_codes[(locke == name) & np.isfinite(hf)] = code
    del locke, hf, intensive
    classes["locke"] = locke_codes
    values = {
        "wcmc_eii_2020": E_LAYER_DIR / "wcmc_eii_2020.tif",
    }
    for name in available_runs():
        values[f"S_{RUNS[name]}"] = our_e(name, "si")
        values[name] = our_e(name, "eii")
    hist_rows, summ_rows = [], []
    for vname, path in values.items():
        v = read_full(path)
        for lname, cls in classes.items():
            names = list(LOCKE_CODES) if lname == "locke" else list(CLASS_CODES)
            hist = class_histograms(cls, v)
            for k, cname in enumerate(names):
                for b in range(100):
                    hist_rows.append((vname, lname, cname, b / 100, int(hist[k, b])))
            for row in hist_summary(hist, names):
                summ_rows.append({"value": vname, "layer": lname} | row)
            summ_rows.append(
                {
                    "value": vname,
                    "layer": lname,
                    "class": "AUC",
                    "auc_vl_high": auc_from_hist(hist[0], hist[2]),
                    "auc_low_high": auc_from_hist(hist[1], hist[2]),
                    "auc_vl_low": auc_from_hist(hist[0], hist[1]),
                }
            )
        del v
    hist_df = pd.DataFrame(hist_rows, columns=["value", "layer", "class", "bin_lo", "km2"])
    return hist_df, pd.DataFrame(summ_rows)


# --------------------------------------------------------------------------------------
# Figs A8 and A9
# --------------------------------------------------------------------------------------


def a8_functional(pts: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name in available_runs():
        for key, col in A8_LAYERS.items():
            if col in pts:
                row = corr_row(pts[f"{name}__fi"], pts[col], run=name, layer=key)
                row["paper_r"] = PAPER["A8"].get(key, np.nan)
                rows.append(row)
    return pd.DataFrame(rows)


def a9_compositional(pts: pd.DataFrame) -> pd.DataFrame:
    rows = []
    layers = {
        "bhi_spe": pts.get("bhi_spe_2015"),
        "bhi_hab": pts.get("bhi_hab_2015"),
        "msa": pts.get("msa_2015"),
        "retained_fauna": pts.get("defaunation_2019"),
        "defaunation": 1.0 - pts["defaunation_2019"] if "defaunation_2019" in pts else None,
    }
    for comp in available_runs():
        if f"{comp}__ci" not in pts:
            continue
        for key, series in layers.items():
            if series is not None:
                row = corr_row(pts[f"{comp}__ci"], series, run=comp, layer=key)
                row["paper_r"] = PAPER["A9"].get(key, np.nan)
                rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
# Our EII vs the published WCMC EII 2020 on the WCMC grid
# --------------------------------------------------------------------------------------

STATS = ("w", "x", "y", "xx", "yy", "xy", "ad", "d2")


def accumulate(acc: np.ndarray, ids: np.ndarray, x: np.ndarray, y: np.ndarray, w: np.ndarray):
    """Add weighted moment sums of (x = WCMC, y = ours) per unit id into ``acc``."""
    n = acc.shape[1]
    d = y - x
    for k, term in enumerate((1.0, x, y, x * x, y * y, x * y, np.abs(d), d * d)):
        acc[k] += np.bincount(ids, weights=w * term, minlength=n)


def finish(acc: np.ndarray) -> pd.DataFrame:
    w, x, y, xx, yy, xy, ad, d2 = acc
    with np.errstate(invalid="ignore", divide="ignore"):
        mx, my = x / w, y / w
        cov = xy / w - mx * my
        r = cov / np.sqrt((xx / w - mx**2) * (yy / w - my**2))
        return pd.DataFrame(
            {
                "area_weight": w,
                "wcmc_mean": mx,
                "ours_mean": my,
                "bias": my - mx,
                "mae": ad / w,
                "rmse": np.sqrt(d2 / w),
                "pearson_r_w": r,
            }
        )


def wcmc_compare(chunk: int = 1500) -> dict[str, Any]:
    """Shared-mask metrics per run, plus biome/ecoregion/Sayre breakdowns and own-mask means."""
    wcmc_path = WCMC_DIR / "EII_2020_1km.tif"
    grid = grid_of(wcmc_path)
    weights = cos_lat_weights(grid)
    eco_path, sayre_path = W_LAYER_DIR / "ecoregion.tif", W_LAYER_DIR / "sayre.tif"
    to_biome = np.zeros(2000, "int64")
    for k, v in json.loads(ECO_TO_BIOME.read_text()).items():
        to_biome[int(k)] = int(v)
    runs = [n for n in RUNS if our_w(n).exists()]
    n_eco, n_sayre = 2000, 1000
    acc = {
        r: {
            "global": np.zeros((8, 1)),
            "biome": np.zeros((8, 16)),
            "ecoregion": np.zeros((8, n_eco)),
            "sayre": np.zeros((8, n_sayre)),
        }
        for r in runs
    }
    # own-mask unit means (Hill's Fig 6 style) for WCMC years and our runs
    own = {k: np.zeros((2, n_sayre)) for k in ["wcmc_2020", *runs]}
    own_biome = {k: np.zeros((2, 16)) for k in ["wcmc_2000", "wcmc_2005", "wcmc_2010",
                                                "wcmc_2015", "wcmc_2020", *runs]}  # fmt: skip
    coverage = {r: {"shared": 0.0, "wcmc_only": 0.0, "ours_only": 0.0} for r in runs}
    diffs: dict[str, list[np.ndarray]] = {r: [] for r in runs}
    rng = np.random.default_rng(SEED)
    files = {"wcmc_2020": wcmc_path, "eco": eco_path, "sayre": sayre_path}
    files |= {f"wcmc_{y}": WCMC_DIR / f"EII_{y}_1km.tif" for y in (2000, 2005, 2010, 2015)}
    files |= {r: our_w(r) for r in runs}
    srcs = {k: rasterio.open(p) for k, p in files.items()}
    try:
        for r0 in range(0, grid["height"], chunk):
            h = min(chunk, grid["height"] - r0)
            win = Window(0, r0, grid["width"], h)
            read = {k: s.read(1, window=win).astype("float32") for k, s in srcs.items()}
            w = np.broadcast_to(weights[r0 : r0 + h, None], (h, grid["width"]))
            eco = np.nan_to_num(read["eco"], nan=0).astype("int64").clip(0, n_eco - 1)
            sayre_raw = read["sayre"]
            sayre = (
                np.where(np.isfinite(sayre_raw), sayre_raw, 0).astype("int64").clip(0, n_sayre - 1)
            )
            biome = to_biome[eco]
            x = read["wcmc_2020"]
            for key in own_biome:
                v = read[key]
                ok = np.isfinite(v)
                own_biome[key][0] += np.bincount(biome[ok], weights=w[ok], minlength=16)
                own_biome[key][1] += np.bincount(biome[ok], weights=w[ok] * v[ok], minlength=16)
            for key in own:
                v = read[key]
                ok = np.isfinite(v) & (sayre > 0)
                own[key][0] += np.bincount(sayre[ok], weights=w[ok], minlength=n_sayre)
                own[key][1] += np.bincount(sayre[ok], weights=w[ok] * v[ok], minlength=n_sayre)
            for run in runs:
                y = read[run]
                fx, fy = np.isfinite(x), np.isfinite(y)
                ok = fx & fy
                coverage[run]["shared"] += float(w[ok].sum())
                coverage[run]["wcmc_only"] += float(w[fx & ~fy].sum())
                coverage[run]["ours_only"] += float(w[fy & ~fx].sum())
                xo, yo, wo = x[ok].astype("float64"), y[ok].astype("float64"), w[ok]
                a = acc[run]
                accumulate(a["global"], np.zeros(len(xo), "int64"), xo, yo, wo)
                accumulate(a["biome"], biome[ok], xo, yo, wo)
                accumulate(a["ecoregion"], eco[ok], xo, yo, wo)
                accumulate(a["sayre"], sayre[ok], xo, yo, wo)
                if len(xo):  # area-proportional subsample for difference quantiles
                    take = rng.random(len(xo)) < 0.02 * wo
                    diffs[run].append((yo - xo)[take].astype("float32"))
            log(f"wcmc comparison rows {r0 + h}/{grid['height']}")
    finally:
        for s in srcs.values():
            s.close()
    tables: dict[str, Any] = {}
    glob_rows = []
    for run in runs:
        g = finish(acc[run]["global"]).iloc[0].to_dict()
        d = np.concatenate(diffs[run]) if diffs[run] else np.array([np.nan])
        cov = coverage[run]
        total = sum(cov.values())
        g |= {
            "run": run,
            "diff_q01": float(np.quantile(d, 0.01)),
            "diff_q99": float(np.quantile(d, 0.99)),
            "share_abs_diff_gt_0.2": float(np.mean(np.abs(d) > 0.2)),
            "shared_area_frac": cov["shared"] / total,
            "wcmc_only_area_frac": cov["wcmc_only"] / total,
            "ours_only_area_frac": cov["ours_only"] / total,
        }
        glob_rows.append(g)
        for unit in ("biome", "ecoregion", "sayre"):
            t = finish(acc[run][unit])
            t.insert(0, unit, np.arange(len(t)))
            t.insert(0, "run", run)
            t = t[(t["area_weight"] > 0) & (t[unit] > 0)]
            tables.setdefault(unit, []).append(t)
        np.save(RUN_DIR / f"wcmc_diff_sample_{run}.npy", d)
    tables["global"] = pd.DataFrame(glob_rows)
    for unit in ("biome", "ecoregion", "sayre"):
        tables[unit] = pd.concat(tables[unit], ignore_index=True)
    tables["biome"]["biome_name"] = tables["biome"]["biome"].map(BIOME_NAMES)
    # Hill-style per-ecosystem summary on each product's own mask
    sayre_rows = []
    for key, (wsum, vsum) in own.items():
        with np.errstate(invalid="ignore", divide="ignore"):
            m = vsum / wsum
        m = m[(wsum > 0) & (np.arange(n_sayre) > 0)]
        sayre_rows.append(
            {
                "version": key,
                "n_ecosystems": int(len(m)),
                "mean_of_means": float(np.mean(m)),
                "sd_of_means": float(np.std(m, ddof=1)),
                "n_le_0.5": int(np.sum(m <= 0.5)),
                "n_le_0.25": int(np.sum(m <= 0.25)),
            }
        )
    tables["sayre_summary"] = pd.DataFrame(sayre_rows)
    biome_rows = []
    for key, (wsum, vsum) in own_biome.items():
        for b in range(1, 15):
            if wsum[b] > 0:
                biome_rows.append({"version": key, "biome": b, "biome_name": BIOME_NAMES[b],
                                   "mean_own_mask": vsum[b] / wsum[b]})  # fmt: skip
    tables["biome_own_mask"] = pd.DataFrame(biome_rows)
    return tables


def wcmc_points_corr(pts: pd.DataFrame) -> pd.DataFrame:
    """Spearman (and unweighted Pearson) of our EII vs WCMC on the equal-area points."""
    return pd.DataFrame(
        [corr_row(pts["wcmc_eii_2020"], pts[f"{n}__eii"], run=n) for n in available_runs()]
    )


def diff_map(run: str, factor: int = 10) -> tuple[np.ndarray, list[float]]:
    """Block-mean (ours - WCMC) on the W grid, decimated by ``factor`` for plotting."""
    with (
        rasterio.open(WCMC_DIR / "EII_2020_1km.tif") as a,
        rasterio.open(our_w(run)) as b,
    ):
        h, w = a.height // factor, a.width // factor
        out = np.full((h, w), np.nan, "float32")
        step = 100 * factor
        for r0 in range(0, h * factor, step):
            hh = min(step, h * factor - r0)
            win = Window(0, r0, w * factor, hh)
            d = b.read(1, window=win) - a.read(1, window=win)
            d = d.reshape(hh // factor, factor, w, factor)
            with np.errstate(invalid="ignore"):
                out[r0 // factor : (r0 + hh) // factor] = np.nanmean(d, axis=(1, 3))
        bounds = list(a.bounds)
    return out, bounds


# --------------------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------------------


def write(frame: pd.DataFrame, name: str) -> Path:
    TABLES.mkdir(parents=True, exist_ok=True)
    path = TABLES / f"{name}.csv"
    frame.to_csv(path, index=False, float_format="%.6g")
    log(f"wrote {path} ({len(frame)} rows)")
    return path


def stage_analyse(parts: Sequence[str] = ("points", "a5", "a6", "wcmc")) -> None:
    import warnings

    warnings.filterwarnings("ignore", category=stats.ConstantInputWarning)
    if "points" in parts:
        pts = pd.read_csv(POINTS)
        write(a3_components(pts), "a3_components")
        pairs, means = a4_indicators(pts)
        write(pairs, "a4_pairs")
        write(means, "a4_mean_abs_r")
        write(a8_functional(pts), "a8_functional")
        write(a9_compositional(pts), "a9_compositional")
        write(wcmc_points_corr(pts), "wcmc_points_corr")
    if "a5" in parts:
        eco, a5 = a5_ecoregions()
        write(eco, "a5_ecoregion_means")
        write(a5, "a5_ecoregion_corr")
    if "a6" in parts:
        hist, summ = a6_influence()
        write(hist, "a6_histograms")
        write(summ, "a6_summary")
    if "wcmc" in parts:
        tables = wcmc_compare()
        for name, frame in tables.items():
            write(frame, f"wcmc_{name}")
