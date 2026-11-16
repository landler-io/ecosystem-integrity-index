"""
Training sample of the potential-NPP model: thinning, quotas and provenance on the candidate pool.

The candidate pool comes from ``build_candidate_pool_local.py`` (``points`` plus predictors and
annual CLMS NPP sampled at each point; see its docstring). ``select`` re-checks eligibility at the
point (WorldCover, valid NPP years, pristine mask, biome), thins by distance, fits climate
clusters, draws per-biome quotas (``eii.training.sampling.draw_design``: area-proportional with a
floor, a cap per biome and 4-degree block, a minimum per climate cluster), adds provenance
columns and writes a sampling report.

The published table ``npp_training_samples_v2_wp11_blockcap_seed42_20261027.csv.gz`` was built
with the defaults and ``--keep-relaxed`` (every thinned row that is eligible apart from the
pristine mask, flagged ``reference`` and ``selected``)::

    python pipelines/training/sample_training_v2.py select \
        --pool npp_candidate_pool_v2_local_seed42_20261027.csv.gz --keep-relaxed \
        --output-csv npp_training_samples_v2_wp11_blockcap_seed42_20261027.csv.gz \
        --report training_selection.md
"""

from __future__ import annotations

import argparse
import math
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from eii.compute.settings import BIOME_NAMES, NPP_STORED_SCALE
from eii.training.sampling import (
    EXCLUDED_WORLDCOVER,
    HMI_MASK_NORTH_EDGE,
    MAX_PER_BLOCK,
    MIN_CLUSTER_SHARE,
    draw_design,
    thin_by_distance,
)

BIOME_LAND_SHARES = Path(__file__).resolve().parent / "biome_land_shares.csv"
"""Biome land shares from a global random sample of 41,907 land points."""
POOL_TOTAL = 200_000
POOL_MIN = 6_000
SCREEN_YEARS = (2014, 2024)
SCREEN_MIN_VALID_YEARS = 8
"""Eligibility screen on the pool's annual NPP: valid years within ``SCREEN_YEARS``."""
RESPONSE_SOURCE = "projects/landler-open-data/assets/datasets/clms_npp/v1_1/annual"
"""Annual CLMS NPP sampled into the published candidate pool."""
BANDS = {"sum": "longterm_avg_npp_sum", "sd": "longterm_avg_npp_sd"}
WORLDCOVER_CLASSES = {
    10: "tree",
    20: "shrub",
    30: "grass",
    40: "cropland",
    50: "built",
    60: "bare",
    70: "snow_ice",
    80: "water",
    90: "wetland",
    95: "mangrove",
    100: "moss_lichen",
}
ENV_CLUSTER_COLUMNS = (
    "mean_annual_temp",
    "annual_precip",
    "aridity",
    "temp_seasonality",
    "precip_seasonality",
)


def block_id(lat: np.ndarray, lon: np.ndarray, deg: float) -> np.ndarray:
    """Integer id of the ``deg``-degree block containing each point (unique for deg >= 0.1)."""
    row = np.floor((np.asarray(lat, dtype="float64") + 90.0) / deg).astype("int64")
    col = np.floor((np.asarray(lon, dtype="float64") + 180.0) / deg).astype("int64")
    return row * 100000 + col


def long_term_response(annual: pd.DataFrame, years: Sequence[int]) -> pd.DataFrame:
    """Mean of the valid (non-negative) annual values per band and their counts."""
    out = {}
    for band, name in BANDS.items():
        values = annual[[f"{band}_{year}" for year in years]].astype("float64")
        values = values.where(values >= 0)
        out[name] = values.mean(axis=1)
        out[f"npp_valid_years_{band}"] = values.notna().sum(axis=1).astype("int64")
    return pd.DataFrame(out, index=annual.index)


def env_clusters(frame: pd.DataFrame, n_clusters: int) -> np.ndarray:
    """k-means clusters on the standardised climate predictors (fixed seed)."""
    from sklearn.cluster import KMeans

    climate = frame[list(ENV_CLUSTER_COLUMNS)]
    climate = (climate - climate.mean()) / climate.std()
    return KMeans(n_clusters, n_init=10, random_state=0).fit_predict(climate.fillna(0))


def markdown_table(frame: pd.DataFrame) -> str:
    """Pipe table without the optional ``tabulate`` dependency."""
    frame = frame.reset_index()
    cells = frame.map(lambda v: f"{v:.4f}" if isinstance(v, float) else str(v))
    rows = [list(frame.columns), ["---"] * frame.shape[1], *cells.to_numpy().tolist()]
    return "\n".join("| " + " | ".join(map(str, row)) + " |" for row in rows)


def land_shares(path: Path = BIOME_LAND_SHARES) -> dict[int, float]:
    """Biome land shares (``BIOME_LAND_SHARES``)."""
    coverage = pd.read_csv(path)
    return dict(zip(coverage["biome"].astype(int), coverage["land_share"], strict=True))


def pool_points(
    shares: dict[int, float], total: int = POOL_TOTAL, minimum: int = POOL_MIN
) -> dict[int, int]:
    """Candidates per biome: area-proportional share of ``total``, at least ``minimum``."""
    return {biome: max(minimum, math.ceil(share * total)) for biome, share in shares.items()}


def prepare_pool(pool: pd.DataFrame, years: tuple[int, int], min_valid_years: int) -> pd.DataFrame:
    """Harmonise the pool's columns and flag eligibility (``exclusion_reason``, empty = keep)."""
    pool = pool.copy()
    if "latitude_point" in pool:
        pool = pool.drop(columns=["latitude", "longitude"], errors="ignore").rename(
            columns={"latitude_point": "latitude", "longitude_point": "longitude"}
        )
    if "sum_2014" in pool:
        span = range(years[0], years[1] + 1)
        pool = pool.join(long_term_response(pool.fillna(-9999), span))
    reasons = {
        **{
            f"worldcover_{name}": pool["worldcover"] == code
            for code, name in EXCLUDED_WORLDCOVER.items()
        },
        "worldcover_missing": ~pool["worldcover"].isin(list(WORLDCOVER_CLASSES)),
        "few_valid_years": pool[["npp_valid_years_sum", "npp_valid_years_sd"]].min(axis=1)
        < min_valid_years,
        "not_pristine": pool["pristine_mask"] != 1,
        "no_biome": ~pool["biome"].isin(list(BIOME_NAMES)),
    }
    flags = pd.DataFrame(reasons, index=pool.index)
    pool["exclusion_reason"] = flags.apply(lambda r: ";".join(flags.columns[r.to_numpy()]), axis=1)
    return pool


def select(
    pool: pd.DataFrame,
    shares: dict[int, float],
    total: int,
    min_quota: int,
    max_quota: int | None,
    min_distance_km: float,
    max_per_block: int | None,
    seed: int,
    keep_relaxed: bool = False,
    min_cluster_share: float | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Thin, allocate and draw; returns the sample and a per-biome accounting table.

    Quotas are drawn from pristine rows (``draw_design``: block cap, environmental-cluster
    minimum), weighted by the pool ``design_weight`` if present. ``env_cluster`` is fitted
    on all thinned rows.
    With ``keep_relaxed`` the sample is every thinned row that is eligible apart from
    ``not_pristine``, flagged ``reference`` (pristine) and ``selected`` (quota draw), so
    relaxed HMI sets can be derived from the same table.
    """
    allowed = ["", "not_pristine"] if keep_relaxed else [""]
    eligible = pool[pool["exclusion_reason"].isin(allowed)].copy()
    eligible["block_4deg"] = block_id(eligible["latitude"], eligible["longitude"], 4.0)
    eligible["block_10deg"] = block_id(eligible["latitude"], eligible["longitude"], 10.0)
    keep = thin_by_distance(
        eligible["latitude"].to_numpy(), eligible["longitude"].to_numpy(), min_distance_km, seed
    )
    thinned = eligible[keep].assign(env_cluster=lambda f: env_clusters(f, n_clusters=10))
    reference = thinned[thinned["exclusion_reason"] == ""]
    eligible = eligible[eligible["exclusion_reason"] == ""]
    weights = "design_weight" if "design_weight" in reference else None
    sample, quotas = draw_design(
        reference,
        shares,
        total,
        min_quota,
        max_quota,
        max_per_block,
        min_cluster_share,
        seed,
        weights,
    )

    accounting = pd.DataFrame(
        {
            "biome_name": pd.Series(BIOME_NAMES),
            "land_share": pd.Series(shares),
            "n_pool": pool["biome"].value_counts(),
            "n_eligible": eligible["biome"].value_counts(),
            "n_thinned": reference["biome"].value_counts(),
            "quota": pd.Series(quotas),
            "n_selected": sample["biome"].value_counts(),
        }
    ).fillna(0)
    counts = ["n_pool", "n_eligible", "n_thinned", "quota", "n_selected"]
    accounting[counts] = accounting[counts].astype(int)
    accounting["share"] = accounting["n_selected"] / accounting["n_selected"].sum()
    accounting["ratio"] = accounting["share"] / accounting["land_share"]
    accounting["at_floor"] = accounting["quota"] <= min_quota
    accounting.index.name = "biome"
    if keep_relaxed:
        selected = thinned.index.isin(sample.index)
        sample = thinned.assign(
            reference=thinned["exclusion_reason"] == "", selected=selected, excluded=~selected
        )
    return sample, accounting


def add_provenance(
    sample: pd.DataFrame,
    accounting: pd.DataFrame,
    years: tuple[int, int],
    min_valid_years: int,
    min_distance_km: float,
) -> pd.DataFrame:
    sample = sample.copy()
    hmi, wdpa = sample["hmi_natural_mask"] == 1, sample["wdpa_high_protection"] == 1
    sample["eligibility_source"] = np.select(
        [hmi & wdpa, hmi, wdpa], ["both", "hmi_mask", "wdpa"], "relaxed"
    )
    # Local pools (with hmi_max) compute the mask on the fly everywhere.
    on_the_fly = ("hmi_max" in sample) | (sample["latitude"] >= HMI_MASK_NORTH_EDGE)
    sample["hmi_mask_source"] = np.where(on_the_fly, "on_the_fly", "tiles")
    sample["n_valid_years"] = sample[["npp_valid_years_sum", "npp_valid_years_sd"]].min(axis=1)
    sample["sample_weight"] = sample["biome"].map(1 / accounting["ratio"])
    sample["response_source"] = RESPONSE_SOURCE
    sample["response_window"] = f"{years[0]}-{years[1]}"
    sample["response_min_valid_years"] = min_valid_years
    sample["min_distance_km"] = min_distance_km
    return sample


def report_markdown(
    title: str,
    source: str,
    pool: pd.DataFrame,
    sample: pd.DataFrame,
    accounting: pd.DataFrame,
    args: argparse.Namespace,
) -> str:
    reasons = pool["exclusion_reason"].str.split(";").explode()
    reasons = reasons[reasons != ""].value_counts().rename_axis("reason").rename("n")
    table = sample
    if "selected" in sample:
        sample = sample[sample["selected"]]
    cover = sample["worldcover"].map(WORLDCOVER_CLASSES).value_counts(normalize=True)
    contamination = sample["worldcover"].isin(list(EXCLUDED_WORLDCOVER)).mean()
    present = accounting[accounting["n_selected"] > 0]
    in_range = present[~present["at_floor"]]["ratio"].between(0.5, 2.0).all()
    response = sample["longterm_avg_npp_sum"] * NPP_STORED_SCALE
    return "\n".join(
        [
            f"# {title}",
            "",
            f"Built by `pipelines/training/sample_training_v2.py select` from `{source}`.",
            "",
            f"- Selection: total {args.total:,}, min quota {args.min_quota:,}, "
            f"max quota {args.max_quota or 'none'}, thinning {args.min_distance_km} km, "
            f"max per (biome, 4-degree block) {args.max_per_block or 'none'}, "
            f"min per (biome, env cluster) {args.min_cluster_share or 0:.0%} of quota, "
            f"seed {args.seed}.",
            f"- Response: mean of valid years {args.first_year}-{args.last_year}, "
            f">= {args.min_valid_years} valid years; min {response.min():.2f} gC m-2 yr-1, "
            f"negative responses {int((response < 0).sum())}.",
            f"- Rows: pool {len(pool):,}, eligible {int((pool['exclusion_reason'] == '').sum()):,}, "
            f"selected {len(sample):,}"
            + (f" (table {len(table):,} incl. relaxed rows)." if len(table) > len(sample) else "."),
            f"- WorldCover contamination (water, snow/ice, built, cropland) in the sample: "
            f"{contamination:.2%} (target <= 0.2 %).",
            f"- Training/area ratio within 0.5-2.0 for every biome above its floor: {in_range}. "
            "Floored biomes are over-represented by design; `sample_weight` (= 1 / ratio) "
            "restores area proportions if wanted.",
            "",
            "## Per biome",
            "",
            markdown_table(accounting),
            "",
            "## Pool exclusions (a row can have several)",
            "",
            markdown_table(reasons.to_frame()) if len(reasons) else "None.",
            "",
            "## WorldCover shares in the sample",
            "",
            markdown_table(cover.rename("share").to_frame()),
            "",
        ]
    )


def run_select(args: argparse.Namespace) -> None:
    years = (args.first_year, args.last_year)
    pool = prepare_pool(pd.read_csv(args.pool), years, args.min_valid_years)
    shares = land_shares()
    sample, accounting = select(
        pool,
        shares,
        args.total,
        args.min_quota,
        args.max_quota,
        args.min_distance_km,
        args.max_per_block or None,
        args.seed,
        args.keep_relaxed,
        args.min_cluster_share or None,
    )
    sample = add_provenance(sample, accounting, years, args.min_valid_years, args.min_distance_km)
    if args.output_csv:
        sample.to_csv(args.output_csv, index=False)
        print(f"wrote {len(sample):,} rows to {args.output_csv}")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        report_markdown(args.title, args.pool.name, pool, sample, accounting, args)
    )
    print(accounting.to_string())
    print(f"wrote {args.report}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    stages = parser.add_subparsers(dest="stage", required=True)

    sel = stages.add_parser("select", help="thinning, quotas and provenance")
    sel.add_argument("--pool", type=Path, required=True)
    sel.add_argument("--total", type=int, default=30_000)
    sel.add_argument("--min-quota", type=int, default=1_000)
    sel.add_argument("--max-quota", type=int)
    sel.add_argument("--min-distance-km", type=float, default=2.0)
    sel.add_argument("--max-per-block", type=int, default=MAX_PER_BLOCK, help="0: no cap")
    sel.add_argument(
        "--min-cluster-share", type=float, default=MIN_CLUSTER_SHARE, help="0: no minimum"
    )
    sel.add_argument("--first-year", type=int, default=SCREEN_YEARS[0])
    sel.add_argument("--last-year", type=int, default=SCREEN_YEARS[1])
    sel.add_argument("--min-valid-years", type=int, default=SCREEN_MIN_VALID_YEARS)
    sel.add_argument("--seed", type=int, default=42)
    sel.add_argument("--keep-relaxed", action="store_true", help="also keep non-pristine rows")
    sel.add_argument("--output-csv", type=Path)
    sel.add_argument("--report", type=Path, default=Path("training_selection.md"))
    sel.add_argument("--title", default="v2 potential-NPP training sample")

    args = parser.parse_args()
    run_select(args)


if __name__ == "__main__":
    main()
