"""Figures of the methodology report (docs/report), from the published EII v2 tiles.

Three figures, written to ``docs/report/figures``:

- ``eii_map_<year>``: global EII at 0.05 deg in the Equal Earth projection (EPSG:8857).
- ``eii_by_landcover_<year>``: EII and its components by ESA WorldCover 2021 class, at an
  equal-area random sample of land points. The sample is stored next to the figure
  (``landcover_sample_<year>.csv.gz``) and reused when present.
- ``eii_by_influence_class_<year>``: EII by human-influence class, from the area histograms
  of the validation suite (``outputs/tables/a6_histograms.csv`` and ``a6_summary.csv``).

``--release`` is a directory with one subdirectory per reference year holding the 10 deg
tiles of the release (``<year>/<tile>.tif``, bands as published), e.g. a local copy of
``gs://landler-opendata-gee-assets/products/eii/v2_0/release/global``. WorldCover is read
from the public ``esa-worldcover`` bucket on AWS (no credentials needed).

Usage::

    python pipelines/validation/report_figures.py --release <dir> [--year 2025]
"""

from __future__ import annotations

import argparse
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, ListedColormap  # noqa: E402
from pyproj import Transformer  # noqa: E402
from rasterio.enums import Resampling  # noqa: E402
from rasterio.transform import from_origin  # noqa: E402
from rasterio.warp import reproject  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "report" / "figures"
TABLES = Path(__file__).resolve().parent / "outputs" / "tables"
WORLDCOVER = "/vsis3/esa-worldcover/v200/2021/map"
SEED = 42
N_SPHERE = 120_000
MAP_RES = 0.05
EII_PIXEL = 1 / 336

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
GREENS = LinearSegmentedColormap.from_list(
    "eii",
    [(0, "#f4efd0"), (0.4, "#d3e6a6"), (0.7, "#8cc47a"), (0.85, "#3f9450"), (1, "#0b4524")],
)
"""Single-hue sequential ramp with more lightness steps above 0.7, where most land lies."""
NO_EII = "#d9d8d3"

BANDS = {1: "eii", 2: "functional", 3: "structural", 4: "compositional"}
BAND_LABEL = {
    "eii": "EII",
    "functional": "Functional",
    "structural": "Structural",
    "compositional": "Compositional",
}
WORLDCOVER_CLASSES = {
    10: "Tree cover",
    20: "Shrubland",
    30: "Grassland",
    90: "Herbaceous wetland",
    95: "Mangroves",
    100: "Moss and lichen",
    60: "Bare / sparse vegetation",
    70: "Snow and ice",
    40: "Cropland",
    50: "Built-up",
}


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 8,
            "text.color": INK,
            "axes.labelcolor": INK2,
            "axes.edgecolor": MUTED,
            "axes.linewidth": 0.6,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.color": INK2,
            "ytick.color": INK2,
            "axes.titlesize": 8.5,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "legend.frameon": False,
            "pdf.fonttype": 42,
        }
    )


def save(fig: plt.Figure, name: str, ext: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.{ext}"
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", path, f"{path.stat().st_size / 1e6:.2f} MB")


def tiles(release: Path, year: int) -> dict[tuple[int, int], Path]:
    """Tile path by (west, south) corner in whole degrees."""
    index = {}
    for path in sorted((release / str(year)).glob("*.tif")):
        with rasterio.open(path) as ds:
            b = ds.bounds
        index[(round(b.left / 10) * 10, round(b.bottom / 10) * 10)] = path
    return index


def global_grid(release: Path, year: int) -> tuple[np.ndarray, np.ndarray]:
    """EII and a land indicator (structural integrity valid) at MAP_RES, by tile averages."""
    cache = Path(tempfile.gettempdir()) / f"eii_report_grid_{year}_{MAP_RES}.npz"
    if cache.exists():
        z = np.load(cache)
        return z["eii"], z["land"]
    n = round(10 / MAP_RES)
    eii = np.full((round(180 / MAP_RES), round(360 / MAP_RES)), np.nan, "float32")
    land = np.zeros(eii.shape, bool)

    def one(item: tuple[tuple[int, int], Path]) -> None:
        (west, south), path = item
        with rasterio.open(path) as ds:
            a = ds.read([1, 3], out_shape=(2, n, n), resampling=Resampling.average, masked=True)
        r0, c0 = round((90 - south - 10) / MAP_RES), round((west + 180) / MAP_RES)
        eii[r0 : r0 + n, c0 : c0 + n] = a[0].filled(np.nan)
        land[r0 : r0 + n, c0 : c0 + n] = ~np.ma.getmaskarray(a[1])

    with ThreadPoolExecutor(8) as pool:
        list(pool.map(one, tiles(release, year).items()))
    np.savez_compressed(cache, eii=eii, land=land)
    return eii, land


def fig_map(release: Path, year: int) -> None:
    eii, land = global_grid(release, year)
    src = from_origin(-180, 90, MAP_RES, MAP_RES)
    fwd = Transformer.from_crs("EPSG:4326", "EPSG:8857", always_xy=True)
    inv = Transformer.from_crs("EPSG:8857", "EPSG:4326", always_xy=True)
    south, north = -60.0, 80.0
    xmax = fwd.transform(180, 0)[0]
    ymin, ymax = fwd.transform(0, south)[1], fwd.transform(0, north)[1]
    width = 2400
    res = 2 * xmax / width
    height = round((ymax - ymin) / res)
    dst = from_origin(-xmax, ymax, res, res)
    yc = ymax - (np.arange(height) + 0.5) * res
    xc = -xmax + (np.arange(width) + 0.5) * res
    edge = fwd.transform(np.full(height, 180.0), inv.transform(np.zeros(height), yc)[1])[0]
    outside = np.abs(xc)[None, :] > edge[:, None]

    def warp(a: np.ndarray, method: Resampling) -> np.ndarray:
        out = np.full((height, width), np.nan, "float32")
        reproject(a.astype("float32"), out, src_transform=src, src_crs="EPSG:4326",
                  dst_transform=dst, dst_crs="EPSG:8857", src_nodata=np.nan, dst_nodata=np.nan,
                  resampling=method)  # fmt: skip
        out[outside] = np.nan
        return out

    e = warp(eii, Resampling.average)
    lnd = warp(np.where(land, 1.0, np.nan), Resampling.nearest)
    extent = (-xmax, xmax, ymin, ymax)

    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    no_eii = np.where(np.isfinite(lnd) & ~np.isfinite(e), 0.0, np.nan)
    ax.imshow(no_eii, extent=extent, cmap=ListedColormap([NO_EII]), interpolation="none")
    im = ax.imshow(e, extent=extent, cmap=GREENS, vmin=0, vmax=1, interpolation="none")
    lat = np.linspace(south, north, 141)
    lon = np.linspace(-180, 180, 361)
    frame = np.concatenate(
        [
            np.c_[np.full_like(lat, -180), lat],
            np.c_[lon, np.full_like(lon, north)],
            np.c_[np.full_like(lat, 180), lat[::-1]],
            np.c_[lon[::-1], np.full_like(lon, south)],
        ]
    )
    ax.plot(*fwd.transform(frame[:, 0], frame[:, 1]), color=MUTED, lw=0.5)
    ax.set_axis_off()
    ax.set_xlim(-xmax * 1.005, xmax * 1.005)
    ax.set_ylim(ymin - 5 * res, ymax + 5 * res)
    cax = ax.inset_axes((0.3, -0.07, 0.4, 0.03))
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label(f"Ecosystem Integrity Index, {year}", color=INK2)
    cb.outline.set_linewidth(0.4)
    cax.tick_params(length=2)
    save(fig, f"eii_map_{year}", "png")


def worldcover_class(item: tuple[str, pd.DataFrame]) -> pd.Series:
    """Dominant WorldCover class in the EII pixel around each point (0 if unavailable)."""
    path, pts = item
    out = pd.Series(0, index=pts.index, dtype="int64")
    try:
        ds = rasterio.open(path)
    except rasterio.errors.RasterioIOError:
        return out
    with ds:
        for i, p in pts.iterrows():
            left = -180 + np.floor((p.lon + 180) / EII_PIXEL) * EII_PIXEL
            top = 90 - np.floor((90 - p.lat) / EII_PIXEL) * EII_PIXEL
            r0, c0 = ds.index(left, top)
            r1, c1 = ds.index(left + EII_PIXEL, top - EII_PIXEL)
            r0, c0 = max(r0, 0), max(c0, 0)
            r1, c1 = min(r1, ds.height), min(c1, ds.width)
            if r1 <= r0 or c1 <= c0:
                continue
            a = ds.read(1, window=((r0, r1), (c0, c1))).ravel()
            a = a[a > 0]
            if a.size:
                out[i] = np.bincount(a).argmax()
    return out


def worldcover_tile(lon: float, lat: float, root: str) -> str:
    s, w = int(np.floor(lat / 3) * 3), int(np.floor(lon / 3) * 3)
    ns, ew = ("N" if s >= 0 else "S"), ("E" if w >= 0 else "W")
    return f"{root}/ESA_WorldCover_10m_2021_v200_{ns}{abs(s):02d}{ew}{abs(w):03d}_Map.tif"


def landcover_sample(release: Path, year: int, worldcover: str) -> pd.DataFrame:
    path = OUT / f"landcover_sample_{year}.csv.gz"
    if path.exists():
        return pd.read_csv(path)
    rng = np.random.default_rng(SEED)
    pts = pd.DataFrame(
        {
            "lon": rng.uniform(-180, 180, N_SPHERE),
            "lat": np.degrees(np.arcsin(rng.uniform(-1, 1, N_SPHERE))),
        }
    )
    index = tiles(release, year)
    key = list(zip(np.floor(pts.lon / 10) * 10, np.floor(pts.lat / 10) * 10, strict=True))
    pts["tile"] = [index.get((int(w), int(s))) for w, s in key]
    pts = pts.dropna(subset=["tile"])
    for name in BANDS.values():
        pts[name] = np.nan
    for tile, sub in pts.groupby("tile"):
        with rasterio.open(tile) as ds:
            xy = zip(sub.lon, sub.lat, strict=True)
            vals = np.array(list(ds.sample(xy, indexes=list(BANDS))))
            vals[vals == ds.nodata] = np.nan
        pts.loc[sub.index, list(BANDS.values())] = vals
    pts = pts.dropna(subset=["eii"]).drop(columns="tile")
    pts["wc_tile"] = [
        worldcover_tile(x, y, worldcover) for x, y in zip(pts.lon, pts.lat, strict=True)
    ]
    os.environ.setdefault("AWS_NO_SIGN_REQUEST", "YES")
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"), ThreadPoolExecutor(16) as pool:
        pts["worldcover"] = pd.concat(pool.map(worldcover_class, pts.groupby("wc_tile")))
    pts = pts.drop(columns="wc_tile").round(5)
    OUT.mkdir(parents=True, exist_ok=True)
    pts.to_csv(path, index=False)
    print("wrote", path, len(pts), "points")
    return pts


def fig_landcover(pts: pd.DataFrame, year: int) -> None:
    counts = pts["worldcover"].value_counts()
    classes = [c for c in WORLDCOVER_CLASSES if counts.get(c, 0) >= 100]
    fig, axes = plt.subplots(1, 4, figsize=(6.6, 2.8), sharey=True, layout="constrained")
    y = np.arange(len(classes))
    for ax, name in zip(axes, BANDS.values(), strict=True):
        q = np.array(
            [pts.loc[pts.worldcover == c, name].dropna().quantile([0.1, 0.25, 0.5, 0.75, 0.9])
             for c in classes]
        )  # fmt: skip
        ax.hlines(y, q[:, 0], q[:, 4], color=SERIES[0], lw=0.8)
        ax.hlines(y, q[:, 1], q[:, 3], color=SERIES[0], lw=4.0, alpha=0.45)
        ax.plot(q[:, 2], y, "o", ms=3.5, color=SERIES[0], mec="white", mew=0.5)
        ax.set_xlim(-0.03, 1.03)
        ax.set_xticks([0, 0.5, 1])
        ax.grid(axis="y", visible=False)
        ax.set_title(BAND_LABEL[name])
    axes[0].set_yticks(y, [f"{WORLDCOVER_CLASSES[c]} ({counts[c]:,})" for c in classes])
    axes[0].invert_yaxis()
    axes[0].tick_params(axis="y", length=0)
    fig.supxlabel(f"Value, {year}", fontsize=8, color=INK2)
    save(fig, f"eii_by_landcover_{year}", "pdf")


def fig_influence(year: int) -> None:
    h = pd.read_csv(TABLES / "a6_histograms.csv")
    h = h[h["value"] == f"eii_{year}"]
    auc = pd.read_csv(TABLES / "a6_summary.csv")
    auc = auc[(auc["value"] == f"eii_{year}") & (auc["class"] == "AUC")].set_index("layer")
    layers = {
        "anthromes": "Anthromes",
        "human_footprint": "Human footprint",
        "low_impact_areas": "Low Impact Areas",
        "locke": "Locke et al. conditions",
    }
    labels = {
        "very_low": "Very low",
        "low": "Low",
        "high": "High",
        "large_wild_areas": "Large wild areas",
        "shared_lands": "Shared lands",
        "cities_farms": "Cities and farms",
    }
    fig, axes = plt.subplots(1, 4, figsize=(6.6, 2.3), sharey=True, layout="constrained")
    edges = np.arange(0, 1.0001, 0.05)
    for ax, (key, title) in zip(axes, layers.items(), strict=True):
        sub = h[h["layer"] == key]
        for k, cls in enumerate(dict.fromkeys(sub["class"])):
            c = sub[sub["class"] == cls]
            area, _ = np.histogram(c["bin_lo"] + 0.005, bins=edges, weights=c["km2"])
            ax.stairs(100 * area / area.sum(), edges, color=SERIES[k], lw=1.2, label=labels[cls])
        ax.set_title(title)
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1])
        ax.text(0.03, 0.62, f"AUC {auc.loc[key, 'auc_vl_high']:.2f}", transform=ax.transAxes,
                color=INK2)  # fmt: skip
        ax.legend(loc="upper left", fontsize=6.5, handlelength=1.2)
    axes[0].set_ylabel("Share of class area (%)")
    fig.supxlabel(f"EII, {year}", fontsize=8, color=INK2)
    save(fig, f"eii_by_influence_class_{year}", "pdf")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--release", type=Path, required=True, help="dir with <year>/<tile>.tif")
    parser.add_argument("--year", type=int, default=2025)
    parser.add_argument("--worldcover", default=WORLDCOVER, help="WorldCover 2021 map tiles")
    args = parser.parse_args()
    style()
    fig_map(args.release, args.year)
    fig_landcover(landcover_sample(args.release, args.year, args.worldcover), args.year)
    fig_influence(args.year)


if __name__ == "__main__":
    main()
