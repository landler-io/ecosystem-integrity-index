"""Figures of the external validation, drawn from the ``analysis`` tables.

Static PNGs. Palette: three categorical colours that stay distinguishable under
colour-vision deficiency, a single-hue blue ramp for point density and a blue/red diverging
ramp with a grey midpoint for differences.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, PowerNorm  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from validation.analysis import (  # noqa: E402
    COMPARISON_RUN,
    FIGS,
    PAPER,
    POINTS,
    TABLES,
    available_runs,
    diff_map,
)
from validation.run_validation import RUNS  # noqa: E402

SURFACE, INK, INK2, MUTED, GRID, AXIS = (
    "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
)  # fmt: skip
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
SEQ = LinearSegmentedColormap.from_list(
    "seq", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"]
)
DIV = LinearSegmentedColormap.from_list(
    "div", ["#0d366b", "#3987e5", "#b7d3f6", "#f0efec", "#f3b9b8", "#e34948", "#8f1f1f"]
)
WCMC_COLOUR = SERIES[1]
LABEL = {"wcmc_2020": "WCMC EII 2020", "wcmc_eii_2020": "WCMC EII 2020"}
LABEL |= {name: f"EII v2 {year}" for name, year in RUNS.items()}
LABEL |= {f"S_{year}": f"Structural {year}" for year in RUNS.values()}
LAYER_LABEL = {
    "pollination": "Pollination",
    "nutrient": "Nutrient retention",
    "soc": "SOC",
    "n_retention": "N retention",
    "ndr": "NDR",
    "bhi_spe": "BHI species",
    "bhi_hab": "BHI habitat",
    "msa": "MSA",
    "retained_fauna": "Retained fauna",
}


def colour(version: str) -> str:
    if version.startswith("wcmc"):
        return WCMC_COLOUR
    return SERIES[0] if version == COMPARISON_RUN else SERIES[2]


def style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "font.family": "sans-serif",
            "font.size": 9,
            "text.color": INK,
            "axes.labelcolor": INK2,
            "axes.edgecolor": AXIS,
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "axes.titlesize": 10,
            "axes.titleweight": "bold",
            "legend.frameon": False,
        }
    )


def save(fig: plt.Figure, name: str) -> Path:
    FIGS.mkdir(parents=True, exist_ok=True)
    path = FIGS / f"{name}.png"
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)
    return path


def hexbin(ax: plt.Axes, x: pd.Series, y: pd.Series, extent=(0, 1, 0, 1)) -> None:
    ok = np.isfinite(x) & np.isfinite(y)
    ax.hexbin(x[ok], y[ok], gridsize=60, cmap=SEQ, mincnt=1, extent=extent, linewidths=0,
              norm=PowerNorm(0.35))  # fmt: skip
    ax.set_xlim(extent[:2])
    ax.set_ylim(extent[2:])
    ax.grid(False)


def fig_a3(pts: pd.DataFrame, run: str = COMPARISON_RUN) -> None:
    t = pd.read_csv(TABLES / "a3_components.csv")
    t = t[(t["run"] == run) & t["pearson_r"].notna()].set_index("pair")
    cols = {"S": f"{run}__si", "C": f"{run}__ci", "F": f"{run}__fi"}
    names = {"S": "Structural", "C": "Compositional", "F": "Functional"}
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.3))
    for ax, (a, b) in zip(axes, (("S", "C"), ("S", "F"), ("C", "F")), strict=True):
        hexbin(ax, pts[cols[a]], pts[cols[b]])
        row = t.loc[f"{a}-{b}"]
        ax.set_title(f"{names[a]} vs {names[b]}", loc="left")
        ax.text(0.03, 0.04, f"r = {row.pearson_r:.3f}\nHill et al. r = {row.paper_r:.3f}",
                transform=ax.transAxes, color=INK, fontsize=8.5)  # fmt: skip
        ax.set_xlabel(names[a])
        ax.set_ylabel(names[b])
    n_lo, n_hi = int(t["n"].min()), int(t["n"].max())
    fig.suptitle(f"Component correlations, {LABEL[run]} ({n_lo:,}-{n_hi:,} equal-area points "
                 "per pair, colour = point density, power scale)", x=0.01, ha="left",
                 fontsize=10)  # fmt: skip
    fig.tight_layout()
    save(fig, "fig_a3_components")


def fig_a4() -> None:
    pairs = pd.read_csv(TABLES / "a4_pairs.csv")
    means = pd.read_csv(TABLES / "a4_mean_abs_r.csv")
    pairs = pairs[(pairs["sample"] == "common") & (pairs["b"] == "eii")]
    means = means[means["sample"] == "common"].set_index("version")
    versions = ["wcmc_2020", *available_runs()]
    titles = {"hfp": "Human footprint (|r|)", "bhi_hab": "BHI habitat", "bhi_spe": "BHI species",
              "msa": "MSA"}  # fmt: skip
    fig, axes = plt.subplots(1, 5, figsize=(13, 2.6), sharey=True)
    y = np.arange(len(versions))
    colors = [colour(v) for v in versions]
    for ax, key in zip(axes[:4], titles, strict=True):
        sub = pairs[pairs["a"] == key].set_index("version")
        vals = [abs(sub.loc[v, "pearson_r"]) for v in versions]
        ax.barh(y, vals, color=colors, height=0.6)
        ax.axvline(abs(PAPER["A4_pairs"][(key, "eii")]), color=INK, lw=1, ls="--")
        for yy, v in zip(y, vals, strict=True):
            ax.text(v + 0.01, yy, f"{v:.2f}", va="center", fontsize=7, color=INK2)
        ax.set_title(titles[key], loc="left")
        ax.set_xlim(0, 1)
    ax = axes[4]
    vals = [means.loc[v, "mean_abs_r_eii"] for v in versions]
    ax.barh(y, vals, color=colors, height=0.6)
    ax.axvline(PAPER["A4_mean"]["eii"], color=INK, lw=1, ls="--")
    for yy, v in zip(y, vals, strict=True):
        ax.text(v + 0.01, yy, f"{v:.2f}", va="center", fontsize=7, color=INK2)
    ax.set_title("Mean |r| across indicators", loc="left")
    ax.set_xlim(0, 1)
    axes[0].set_yticks(y, [LABEL[v] for v in versions])
    axes[0].invert_yaxis()
    n_common = int(means["n_common"].iloc[0])
    fig.text(0.01, -0.08, f"Pearson |r| with each indicator on the {n_common:,} equal-area points "
             "valid in all layers. Dashed line: value for the WCMC EII in Hill et al. (2026), "
             "Fig A4.", fontsize=8, color=INK2)  # fmt: skip
    fig.suptitle("EII against four biodiversity indicators", x=0.01, ha="left")
    fig.tight_layout()
    save(fig, "fig_a4_indicators")


def fig_a5() -> None:
    eco = pd.read_csv(TABLES / "a5_ecoregion_means.csv")
    corr = pd.read_csv(TABLES / "a5_ecoregion_corr.csv").set_index("version")
    versions = [v for v in ("wcmc_2020", *available_runs()) if v in eco]
    fig, axes = plt.subplots(1, len(versions), figsize=(3.6 * len(versions), 3.4), sharey=True)
    for ax, v in zip(np.atleast_1d(axes), versions, strict=True):
        ax.scatter(eco[v], eco["intactness"], s=9, color=colour(v), edgecolors=SURFACE,
                   linewidths=0.4)  # fmt: skip
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(LABEL[v], loc="left")
        ax.text(0.03, 0.9, f"r = {corr.loc[v, 'pearson_r']:.3f}\n"
                f"n = {int(corr.loc[v, 'n'])} ecoregions", transform=ax.transAxes, fontsize=8.5)  # fmt: skip
        ax.set_xlabel("Mean EII per ecoregion")
    np.atleast_1d(axes)[0].set_ylabel("Ecoregion intactness (Beyer et al. 2020)")
    fig.suptitle("Ecoregion means against ecoregion intactness", x=0.01, ha="left")
    fig.tight_layout()
    save(fig, "fig_a5_ecoregions")


def fig_a6(value: str, name: str) -> None:
    h = pd.read_csv(TABLES / "a6_histograms.csv")
    h = h[h["value"] == value]
    layers = {"anthromes": "Anthromes", "human_footprint": "Human footprint (HII 2013)",
              "low_impact_areas": "Low Impact Areas", "locke": "Locke et al. conditions"}  # fmt: skip
    labels = {"very_low": "Very low", "low": "Low", "high": "High",
              "large_wild_areas": "Large wild areas", "shared_lands": "Shared lands",
              "cities_farms": "Cities and farms"}  # fmt: skip
    fig, axes = plt.subplots(2, 4, figsize=(13, 5.2), sharex=True)
    for j, (lkey, ltitle) in enumerate(layers.items()):
        sub = h[h["layer"] == lkey]
        classes = list(dict.fromkeys(sub["class"]))
        piv = sub.pivot(index="bin_lo", columns="class", values="km2")[classes]
        x = piv.index.to_numpy() + 0.005
        tot = piv.sum(axis=1).replace(0, np.nan)
        for k, c in enumerate(classes):
            axes[0, j].plot(x, piv[c] / 1e6, color=SERIES[k], lw=1.6, label=labels[c])
            axes[1, j].plot(x, piv[c] / tot, color=SERIES[k], lw=1.6, label=labels[c])
        axes[0, j].set_title(ltitle, loc="left")
        axes[1, j].set_ylim(0, 1)
        axes[1, j].set_xlabel(LABEL[value])
        axes[0, j].legend(fontsize=7.5, loc="upper left")
    axes[0, 0].set_ylabel("Area (million km$^2$) per 0.01 bin")
    axes[1, 0].set_ylabel("Proportion of area in bin")
    fig.suptitle(f"{LABEL[value]} by human-influence class (full Eckert IV 1 km grid)", x=0.01,
                 ha="left")  # fmt: skip
    fig.tight_layout()
    save(fig, name)


def fig_a6_auc() -> None:
    s = pd.read_csv(TABLES / "a6_summary.csv")
    s = s[s["class"] == "AUC"]
    values = ["wcmc_eii_2020", *available_runs(), *[f"S_{RUNS[n]}" for n in available_runs()]]
    layers = {"anthromes": "Anthromes", "human_footprint": "Human footprint",
              "low_impact_areas": "Low Impact Areas", "locke": "Locke conditions"}  # fmt: skip
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.0), sharey=True)
    y = np.arange(len(values))
    colors = [MUTED if v.startswith("S_") else colour(v) for v in values]
    for ax, (lkey, ltitle) in zip(axes, layers.items(), strict=True):
        sub = s[s["layer"] == lkey].set_index("value")
        vals = [sub.loc[v, "auc_vl_high"] if v in sub.index else np.nan for v in values]
        ax.barh(y, vals, color=colors, height=0.6)
        ax.axvline(0.5, color=MUTED, lw=0.8)
        ax.set_xlim(0.4, 1)
        ax.set_title(ltitle, loc="left")
        for yy, v in zip(y, vals, strict=True):
            if np.isfinite(v):
                ax.text(v + 0.005, yy, f"{v:.2f}", va="center", fontsize=7, color=INK2)
    axes[0].set_yticks(y, [LABEL[v] for v in values])
    axes[0].invert_yaxis()
    fig.text(0.01, -0.05, "AUC = probability that a random km$^2$ of the least-influenced class "
             "(very low / large wild areas) scores higher than one of the most-influenced class "
             "(high / cities and farms); 0.5 = no separation.", fontsize=8, color=INK2)  # fmt: skip
    fig.suptitle("Separation of human-influence classes", x=0.01, ha="left")
    fig.tight_layout()
    save(fig, "fig_a6_auc")


def heat(ax: plt.Axes, table: pd.DataFrame, title: str):
    im = ax.imshow(table.to_numpy(float), cmap=DIV, vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(table.shape[1]), table.columns, rotation=30, ha="right")
    ax.set_yticks(range(table.shape[0]), table.index)
    ax.grid(False)
    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            v = table.iat[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.3f}", ha="center", va="center", fontsize=7.5,
                        color="#ffffff" if abs(v) > 0.4 else INK)  # fmt: skip
    ax.set_title(title, loc="left")
    return im


def fig_a8_a9() -> None:
    a8 = pd.read_csv(TABLES / "a8_functional.csv")
    a9 = pd.read_csv(TABLES / "a9_compositional.csv")
    t8 = a8.pivot(index="run", columns="layer", values="pearson_r")
    cols8 = [c for c in ("pollination", "nutrient", "soc", "n_retention", "ndr") if c in t8]
    t8 = t8.reindex(available_runs())[cols8]
    t8.loc["Hill et al. (2026)"] = [PAPER["A8"].get(c, np.nan) for c in cols8]
    t9 = a9.pivot(index="run", columns="layer", values="pearson_r")
    cols9 = [c for c in ("bhi_spe", "bhi_hab", "msa", "retained_fauna") if c in t9]
    t9 = t9.reindex(available_runs())[cols9]
    t9.loc["Hill et al. (2026)"] = [PAPER["A9"].get(c, np.nan) for c in cols9]
    for t in (t8, t9):
        t.index = [LABEL.get(i, i) for i in t.index]
        t.columns = [LAYER_LABEL[c] for c in t.columns]
    fig, axes = plt.subplots(
        2, 1, figsize=(7, 5.5), gridspec_kw={"height_ratios": [1, 1], "hspace": 0.6}
    )
    heat(axes[0], t8, "Functional integrity vs nature's contributions and soil carbon")
    im = heat(axes[1], t9, "Compositional integrity vs biodiversity layers")
    fig.colorbar(im, ax=axes, shrink=0.5, label="Pearson r")
    save(fig, "fig_a8_a9")


def fig_wcmc_scatter(pts: pd.DataFrame) -> None:
    runs = [r for r in available_runs() if f"{r}__eii" in pts]
    g = pd.read_csv(TABLES / "wcmc_global.csv").set_index("run")
    fig, axes = plt.subplots(1, len(runs), figsize=(3.6 * len(runs), 3.5), sharey=True)
    for ax, run in zip(np.atleast_1d(axes), runs, strict=True):
        hexbin(ax, pts["wcmc_eii_2020"], pts[f"{run}__eii"])
        ax.plot([0, 1], [0, 1], color=MUTED, lw=0.8)
        ax.set_title(LABEL[run], loc="left")
        r = g.loc[run]
        ax.text(0.03, 0.80, f"r = {r.pearson_r_w:.3f}\nbias = {r.bias:+.3f}\nRMSE = {r.rmse:.3f}",
                transform=ax.transAxes, fontsize=8)  # fmt: skip
        ax.set_xlabel("WCMC EII 2020")
    np.atleast_1d(axes)[0].set_ylabel("EII v2")
    fig.suptitle("EII v2 against the WCMC EII 2020 (density of the equal-area points; "
                 "statistics from the full 1 km grid, area-weighted)", x=0.01, ha="left",
                 fontsize=9)  # fmt: skip
    fig.tight_layout()
    save(fig, "fig_wcmc_scatter")


def fig_wcmc_diffmap(run: str = COMPARISON_RUN) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.6))
    d, (w, s, e, n) = diff_map(run)
    im = ax.imshow(d, cmap=DIV, vmin=-0.5, vmax=0.5, extent=(w, e, s, n), interpolation="none")
    ax.set_ylim(-58, 84)
    ax.grid(False)
    ax.set_title(f"{LABEL[run]} minus WCMC EII 2020", loc="left")
    fig.colorbar(im, ax=ax, shrink=0.75, label="EII difference (EII v2 - WCMC)")
    save(fig, "fig_wcmc_diffmap")


def fig_wcmc_biome(run: str = COMPARISON_RUN) -> None:
    b = pd.read_csv(TABLES / "wcmc_biome.csv")
    sub = b[b["run"] == run].set_index("biome_name").sort_values("wcmc_mean")
    fig, ax = plt.subplots(figsize=(8, 5))
    y = np.arange(len(sub))
    ax.hlines(y, sub["wcmc_mean"], sub["ours_mean"], color=AXIS, lw=1.2, zorder=2)
    for col, label, c in (("wcmc_mean", "WCMC EII 2020", WCMC_COLOUR),
                          ("ours_mean", LABEL[run], colour(run))):  # fmt: skip
        ax.scatter(sub[col], y, color=c, s=40, label=label, zorder=3, edgecolors=SURFACE,
                   linewidths=1.5)  # fmt: skip
    ax.set_yticks(y, sub.index)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Area-weighted mean EII on the shared mask")
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title(f"Mean EII per WWF biome: {LABEL[run]} and WCMC EII 2020", loc="left")
    save(fig, "fig_wcmc_biome")


def main() -> None:
    style()
    pts = pd.read_csv(POINTS)
    fig_a3(pts)
    fig_a4()
    fig_a5()
    fig_a6(f"S_{RUNS[COMPARISON_RUN]}", "fig_a6_structural")
    fig_a6("wcmc_eii_2020", "fig_a6_wcmc_eii")
    fig_a6(COMPARISON_RUN, "fig_a6_eii")
    fig_a6_auc()
    fig_a8_a9()
    fig_wcmc_scatter(pts)
    fig_wcmc_diffmap()
    fig_wcmc_biome()


if __name__ == "__main__":
    main()
