"""
Training sample design v2 for the potential-NPP model.

The candidate pool is built from local rasters (``natural_areas_arrays``,
``distance_to_mask_km``, ``one_per_cell``, ``inclusion_probabilities``; driven by
``pipelines/training/build_candidate_pool_local.py``). It is then thinned against pseudo-replication at the
1 km predictor scale and drawn into per-biome quotas (``thin_by_distance``,
``allocate_quotas``, ``select_samples``, ``draw_design``).

Eligible pixels are natural (HMI mask) or highly protected (WDPA), and not water,
snow/ice, built or cropland in WorldCover 2021. Bare and sparse land stays eligible.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping

import numpy as np
import pandas as pd

# Natural-areas mask definition of the v1 hmi_masks
HMI_THRESHOLD = 0.1
HIGH_MOD_THRESHOLD = 0.4
BUFFER_PIXELS = 10
# The exported hmi_masks tiles stop at 60 N; north of it the mask is computed on the fly.
HMI_MASK_NORTH_EDGE = 60.0
EXCLUDED_WORLDCOVER = {40: "cropland", 50: "built", 70: "snow_ice", 80: "water"}
EARTH_RADIUS_KM = 6371.0088
# Candidate domain: HMI_max <= 0.3 (or high protection), banded at these edges, so relaxed
# sets come from the same pool. The reference mask itself is at HMI_THRESHOLD.
CANDIDATE_HMI_MAX = 0.3
HMI_BAND_EDGES = (0.0, 0.05, 0.1, 0.2, 0.3)
# Europe stratum for the extra quota: the box 25 W-45 E, 36-72 N (west of the Urals and
# Caspian, north of the Maghreb). It also contains western Anatolia, the Caucasus fringe,
# the northern Tunisian/Algerian coast and Iceland; it is a sampling stratum, not a border.
EUROPE_BOUNDS = (-25.0, 36.0, 45.0, 72.0)
# One candidate per cell of the 30-arc-second predictor grid (WorldClim/CHELSA, ~1 km).
CANDIDATE_CELL_DEG = 1 / 120
# Quota draw stratification: at most 40 rows per (biome,
# 4-degree block), and every environmental cluster (k-means 10 on climate) gets at least 1 %
# of the biome quota where the biome has rows in it.
MAX_PER_BLOCK = 40
MIN_CLUSTER_SHARE = 0.01
# Gap fill: strata (4-degree block x aridity class) with fewer pristine training rows than
# GAP_MIN_REFERENCE are topped up with relaxed rows (eligible apart from the pristine mask,
# HMI below GAP_HMI_MAX), fitted at GAP_WEIGHT. GAP_CHECK_SHARE of relaxed rows is never
# used, so it stays an out-of-sample check.
# Aridity edges: UNEP hyper-arid, arid (split at 0.1), semi-arid, dry sub-humid, humid.
GAP_MIN_REFERENCE = 40
GAP_HMI_MAX = 0.2
GAP_WEIGHT = 0.5
GAP_CHECK_SHARE = 0.5
GAP_BLOCK_DEG = 4.0
GAP_ARIDITY_EDGES = (0.05, 0.1, 0.2, 0.5, 0.65)


# --------------------------------------------------------------------------------------
# Offline: quotas and thinning
# --------------------------------------------------------------------------------------


def allocate_quotas(
    area_share: Mapping[int, float],
    total: int,
    min_quota: int = 1000,
    max_quota: int | None = None,
    available: Mapping[int, int] | None = None,
) -> dict[int, int]:
    """
    Area-proportional quotas with a floor and optional cap per stratum.

    Quota = ``clip(lam * area_share, min_quota, max_quota)``, with ``lam`` found by
    bisection so the quotas sum to ``total``. ``available`` caps floor and cap alike.

    Args:
        area_share: Land area (or share) per stratum.
        total: Target sample size.
        min_quota: Minimum points per stratum.
        max_quota: Optional maximum points per stratum.
        available: Optional candidate count per stratum (hard cap).

    Returns:
        Integer quota per stratum. Sums to ``total`` unless the floors exceed it (then
        every stratum gets its floor) or the caps fall short of it (every stratum capped).
    """
    keys = list(area_share)
    share = np.array([area_share[k] for k in keys], dtype="float64")
    hi = np.full(len(keys), np.inf if max_quota is None else float(max_quota))
    if available is not None:
        hi = np.minimum(hi, [available.get(k, 0) for k in keys])
    lo = np.minimum(min_quota, hi)
    target = min(total, hi.sum())
    lam_lo, lam_hi = 0.0, target / share[share > 0].min() + 1.0
    for _ in range(200):
        lam = (lam_lo + lam_hi) / 2
        if np.clip(lam * share, lo, hi).sum() < target:
            lam_lo = lam
        else:
            lam_hi = lam
    quotas = np.clip(lam_hi * share, lo, hi)
    rounded = np.floor(quotas)
    # Largest-remainder rounding, only into strata below their cap.
    for i in np.argsort(rounded - quotas):
        if rounded.sum() >= target:
            break
        if rounded[i] < hi[i]:
            rounded[i] += 1
    return {k: int(v) for k, v in zip(keys, rounded, strict=True)}


def thin_by_distance(
    lat: np.ndarray,
    lon: np.ndarray,
    min_distance_km: float,
    seed: int = 42,
) -> np.ndarray:
    """
    Random-order greedy thinning: no two kept points closer than ``min_distance_km``.

    Args:
        lat: Latitudes in degrees.
        lon: Longitudes in degrees.
        min_distance_km: Minimum great-circle distance between kept points.
        seed: Seed for the visiting order.

    Returns:
        Boolean keep mask aligned with the inputs.
    """
    from sklearn.neighbors import BallTree

    coords = np.radians(np.column_stack([lat, lon]))
    keep = np.zeros(len(coords), dtype=bool)
    if len(coords) == 0 or min_distance_km <= 0:
        return ~keep
    neighbours = BallTree(coords, metric="haversine").query_radius(
        coords, r=min_distance_km / EARTH_RADIUS_KM
    )
    blocked = np.zeros(len(coords), dtype=bool)
    for i in np.random.default_rng(seed).permutation(len(coords)):
        if not blocked[i]:
            keep[i] = True
            blocked[neighbours[i]] = True
    return keep


def select_samples(
    pool: pd.DataFrame,
    quotas: Mapping[int, int],
    stratum: str = "biome",
    max_per_block: int | None = None,
    block: str = "block_4deg",
    seed: int = 42,
    weights: str | None = None,
    min_cluster_share: float | None = None,
    cluster: str = "env_cluster",
) -> pd.DataFrame:
    """
    Draw ``quotas[s]`` random rows per stratum, optionally capped per spatial block.

    Args:
        pool: Eligible, thinned candidates.
        quotas: Rows per stratum value.
        stratum: Stratum column.
        max_per_block: Optional cap per (stratum, block); the cap can leave a quota short.
        block: Block column for the cap.
        seed: Random seed.
        weights: Optional column of draw weights, e.g. the pool design weight ``1 / p``,
            so a pool with unequal inclusion probabilities yields an area-uniform draw.
        min_cluster_share: Optional minimum per (stratum, cluster) as a share of the
            stratum quota (rounded up), met where the cluster has enough rows.
        cluster: Cluster column for the minimum.

    Returns:
        Selected rows in pool order.
    """
    if weights is None:
        shuffled = pool.sample(frac=1.0, random_state=seed)
    else:
        # Efraimidis-Spirakis order: every prefix is a weighted draw without replacement.
        keys = np.random.default_rng(seed).exponential(size=len(pool)) / pool[weights].to_numpy()
        shuffled = pool.iloc[np.argsort(keys, kind="stable")]
    if max_per_block is not None:
        rank = shuffled.groupby([stratum, block]).cumcount()
        shuffled = shuffled[rank < max_per_block]
    if min_cluster_share:
        floor = np.ceil(shuffled[stratum].map(quotas).fillna(0) * min_cluster_share)
        first = shuffled.groupby([stratum, cluster]).cumcount() < floor
        shuffled = pd.concat([shuffled[first], shuffled[~first]])
    rank = shuffled.groupby(stratum).cumcount()
    limit = shuffled[stratum].map(quotas).fillna(0)
    return shuffled[rank < limit].sort_index()


def draw_design(
    pool: pd.DataFrame,
    area_share: Mapping[int, float],
    total: int,
    min_quota: int = 1000,
    max_quota: int | None = None,
    max_per_block: int | None = MAX_PER_BLOCK,
    min_cluster_share: float | None = MIN_CLUSTER_SHARE,
    seed: int = 42,
    weights: str | None = None,
) -> tuple[pd.DataFrame, dict[int, int]]:
    """
    Quota draw: per-biome quotas (``allocate_quotas``) drawn by ``select_samples``.

    A biome's available rows count at most ``max_per_block`` per ``block_4deg``, so a
    binding block cap moves quota to other biomes and the total is kept.

    Args:
        pool: Eligible, thinned rows with ``biome``, ``block_4deg`` (and ``env_cluster``).
        area_share: Land share per biome.
        total: Target sample size.
        min_quota: Minimum rows per biome.
        max_quota: Optional maximum rows per biome.
        max_per_block: Optional cap per (biome, 4-degree block).
        min_cluster_share: Optional minimum per (biome, ``env_cluster``), share of quota.
        seed: Random seed.
        weights: Optional draw weight column (``select_samples``).

    Returns:
        Selected rows and the quota per biome.
    """
    per_block = pool.groupby(["biome", "block_4deg"]).size()
    if max_per_block is not None:
        per_block = per_block.clip(upper=max_per_block)
    available = per_block.groupby(level=0).sum().to_dict()
    quotas = allocate_quotas(area_share, total, min_quota, max_quota, available)
    sample = select_samples(
        pool,
        quotas,
        max_per_block=max_per_block,
        seed=seed,
        weights=weights,
        min_cluster_share=min_cluster_share,
    )
    return sample, quotas


def gap_fill(
    table: pd.DataFrame,
    min_reference: int = GAP_MIN_REFERENCE,
    hmi_max: float = GAP_HMI_MAX,
    check_share: float = GAP_CHECK_SHARE,
    block_deg: float = GAP_BLOCK_DEG,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Relaxed rows that top up strata short of pristine training rows, and the check rows.

    In the Sahel and on the Loess Plateau the pristine rows of a block are bare desert, while
    the semi-arid land of the same block is all in use. The potential-NPP model then
    extrapolates from the desert, is biased low, and FI scores the surplus as damage. A
    stratum is a ``block_deg`` block times a ``GAP_ARIDITY_EDGES`` class. Relaxed rows are
    ``exclusion_reason == "not_pristine"`` with ``hmi_2022`` below ``hmi_max``. A
    pseudo-random ``check_share`` of them is set aside first. In every stratum with fewer than
    ``min_reference`` reference rows, the remaining relaxed rows fill the shortfall, lowest
    HMI first.

    Args:
        table: Training table with ``reference``, ``exclusion_reason``, ``hmi_2022``,
            ``point_id``, ``aridity``, ``latitude`` and ``longitude``
            (``sample_training_v2.py select --keep-relaxed``).
        min_reference: Target training rows per stratum.
        hmi_max: Upper HMI bound (exclusive) of relaxed rows.
        check_share: Share of relaxed rows held out for checks.
        block_deg: Block size in degrees.
        seed: Seed of the check split (a hash of ``point_id``).

    Returns:
        Boolean arrays over ``table`` rows: gap-fill training rows, held-out check rows.
    """
    hmi = table["hmi_2022"].to_numpy(dtype="float64")
    relaxed = (table["exclusion_reason"].fillna("") == "not_pristine").to_numpy() & (hmi < hmi_max)
    # Hashing point_id keeps the split identical however the table was filtered or ordered.
    key = pd.util.hash_pandas_object(table["point_id"], index=False, hash_key=f"{seed:016d}")
    check = relaxed & (key.to_numpy() % 10_000 < check_share * 10_000)
    lat = np.floor(table["latitude"].to_numpy() / block_deg).astype("int64")
    lon = np.floor(table["longitude"].to_numpy() / block_deg).astype("int64")
    # NaN aridity falls in the class above the last edge, with the humid rows.
    dryness = np.digitize(table["aridity"].to_numpy(dtype="float64"), GAP_ARIDITY_EDGES)
    block = pd.Series((lat * 1000 + lon) * 10 + dryness)
    n_reference = block[table["reference"].astype(bool).to_numpy()].value_counts()
    shortfall = (min_reference - block.map(n_reference).fillna(0)).clip(lower=0).to_numpy()
    candidates = np.flatnonzero(relaxed & ~check)
    order = candidates[np.argsort(hmi[candidates], kind="stable")]
    rank = block.iloc[order].groupby(block.iloc[order].to_numpy()).cumcount().to_numpy()
    fill = np.zeros(len(table), dtype=bool)
    fill[order[rank < shortfall[order]]] = True
    return fill, check


# --------------------------------------------------------------------------------------
# Offline: candidate pool from local rasters (numpy mirrors of the Earth Engine logic)
# --------------------------------------------------------------------------------------


def _focal_max(values: np.ndarray, footprint: np.ndarray) -> np.ndarray:
    """Focal max ignoring NaN, NaN where the window has no valid value (as ee focalMax)."""
    from scipy.ndimage import maximum_filter

    filled = maximum_filter(np.where(np.isnan(values), -np.inf, values), footprint=footprint)
    return np.where(np.isneginf(filled), np.nan, filled)


def natural_areas_arrays(
    hmi_series: np.ndarray,
    hmi_threshold: float = HMI_THRESHOLD,
    high_mod_threshold: float = HIGH_MOD_THRESHOLD,
    buffer_pixels: int = BUFFER_PIXELS,
) -> tuple[np.ndarray, np.ndarray]:
    """
    ``natural_areas_mask`` on arrays of the 300 m HMI grid.

    Args:
        hmi_series: ``(years, rows, cols)`` HMI, NaN outside the data (ocean).
        hmi_threshold: Upper HMI bound for natural pixels.
        high_mod_threshold: HMI above which pixels are buffered out.
        buffer_pixels: Buffer radius in pixels (circle: centre distance <= radius).

    Returns:
        ``hmi_max`` (series max, 5x5 focal-max gap-filled, NaN where undefined) and the
        boolean natural mask (False where ``hmi_max`` is NaN). Values within
        ``buffer_pixels + 2`` of the array edge are only correct with a halo that wide.
    """
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        hmi_max = np.nanmax(hmi_series, axis=0)
    hmi_max = _focal_max(hmi_max, np.ones((5, 5), dtype=bool))
    offsets = np.arange(-buffer_pixels, buffer_pixels + 1)
    circle = offsets[:, None] ** 2 + offsets[None, :] ** 2 <= buffer_pixels**2
    high = np.where(np.isnan(hmi_max), np.nan, (hmi_max > high_mod_threshold).astype(float))
    near_high = _focal_max(high, circle) == 1
    natural = (hmi_max < hmi_threshold) & ~near_high
    return hmi_max, natural


def distance_to_mask_km(
    target: np.ndarray, lat: np.ndarray, pixel_deg: float, max_km: float
) -> np.ndarray:
    """
    Great-circle distance (km) from each pixel to the nearest ``target`` pixel.

    The nearest pixel is found with a Euclidean transform scaled at the mean latitude,
    then measured exactly; values above ``max_km`` (or with no target) become ``max_km``.

    Args:
        target: Boolean array on a regular lat/lon grid.
        lat: Latitude of each row (degrees).
        pixel_deg: Pixel size in degrees.
        max_km: Cap; set it below the halo width so edge effects stay invisible.

    Returns:
        Float32 distances in km.
    """
    from scipy.ndimage import distance_transform_edt

    if not target.any():
        return np.full(target.shape, max_km, dtype="float32")
    km = pixel_deg * np.pi / 180 * EARTH_RADIUS_KM
    sampling = (km, km * np.cos(np.radians(np.mean(lat))))
    _, (rows, cols) = distance_transform_edt(~target, sampling=sampling, return_indices=True)
    lat_r = np.radians(lat)[:, None]
    lat_t = np.radians(lat)[rows]
    dlon = np.radians((cols - np.arange(target.shape[1])[None, :]) * pixel_deg)
    hav = np.sin((lat_t - lat_r) / 2) ** 2 + np.cos(lat_r) * np.cos(lat_t) * np.sin(dlon / 2) ** 2
    dist = 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(hav, 0, 1)))
    return np.minimum(dist, max_km).astype("float32")


def hmi_band(hmi: np.ndarray, edges: tuple[float, ...] = HMI_BAND_EDGES) -> np.ndarray:
    """Band index of ``hmi``: 0 for [edges[0], edges[1]), ..., -1 above the last edge or NaN."""
    band = np.searchsorted(np.asarray(edges), hmi, side="right") - 1
    band = np.where(hmi == edges[-1], len(edges) - 2, band)
    return np.where((band >= len(edges) - 1) | (band < 0) | np.isnan(hmi), -1, band)


def in_europe(lon: np.ndarray, lat: np.ndarray, bounds: tuple = EUROPE_BOUNDS) -> np.ndarray:
    """True inside ``EUROPE_BOUNDS`` (west, south, east, north; west/south inclusive)."""
    west, south, east, north = bounds
    return (lon >= west) & (lon < east) & (lat >= south) & (lat < north)


def one_per_cell(cell_ids: np.ndarray, keys: np.ndarray) -> np.ndarray:
    """
    Indices of the row with the smallest ``key`` per cell (random pick for random keys).

    Args:
        cell_ids: Integer cell id per row.
        keys: Tie-breaking key per row, e.g. uniform random numbers.

    Returns:
        Sorted row indices, one per distinct cell.
    """
    order = np.lexsort((keys, cell_ids))
    first = np.ones(len(order), dtype=bool)
    first[1:] = cell_ids[order][1:] != cell_ids[order][:-1]
    return np.sort(order[first])


def inclusion_probabilities(
    biome: np.ndarray,
    europe: np.ndarray,
    band: np.ndarray,
    biome_targets: Mapping[int, float],
    europe_target: float = 0.0,
    band_target: float = 0.0,
) -> np.ndarray:
    """
    Poisson-sampling inclusion probability per candidate.

    Each criterion sets a rate (target / candidates in its stratum): the biome pool, the
    Europe quota and a minimum per HMI band. A candidate takes the largest rate that
    applies to it, capped at 1, so every target is met in expectation and the design
    weight is ``1 / p``.

    Args:
        biome: Biome number per candidate (strata missing from ``biome_targets`` get 0).
        europe: Europe flag per candidate.
        band: HMI band per candidate (``hmi_band``; -1 has no band minimum).
        biome_targets: Expected pool size per biome.
        europe_target: Expected pool size in Europe from the Europe criterion alone.
        band_target: Expected pool size per HMI band from the band criterion alone.

    Returns:
        Probabilities in [0, 1].
    """

    def rates(values: np.ndarray, targets: Mapping[int, float]) -> np.ndarray:
        keys, counts = np.unique(values, return_counts=True)
        lookup = {k: targets.get(k, 0.0) / n for k, n in zip(keys.tolist(), counts, strict=True)}
        return np.array([lookup[v] for v in values.tolist()], dtype="float64")

    p = rates(biome, biome_targets)
    if europe_target and europe.any():
        p = np.maximum(p, np.where(europe, europe_target / europe.sum(), 0.0))
    if band_target:
        bands = {b: band_target for b in np.unique(band).tolist() if b >= 0}
        p = np.maximum(p, rates(band, bands))
    return np.minimum(p, 1.0)
