#!/usr/bin/env python3
"""Predictors and annual CLMS NPP sampled in Earth Engine at the candidate-pool points.

Reads ``points.csv.gz`` of ``build_candidate_pool_local.py points`` and writes the candidate pool
read by ``sample_training_v2.py select``: per point the 19 predictors of the v1 model card and
annual CLMS NPP v1.1 for 2014-2025 (``sum_YEAR``, ``sd_YEAR``: stored integers, gC m-2 yr-1
and gC m-2 day-1 times 1e4; -9999 where masked or negative), sampled at the exact point. The published pool is
``npp_candidate_pool_v2_local_seed42_20261027.csv.gz``.

Both images are sampled with ``sampleRegions`` at the input point and joined on ``point_id``:
the predictor stack on its 300 m EPSG:4326 grid (``scale=300``; masked predictors are NaN),
annual NPP on the native CLMS grid. Requests of ``--chunk-size`` points run on the standard
endpoint (``getInfo``, no exports); each chunk is cached as ``<cache>/<k>.csv.gz`` and a rerun
only requests the missing chunks. The pool is written once every chunk is present.

Predictor stack (``predictor_stack``): WorldClim bio01, bio04, bio12, bio15, Global Aridity
Index, MERIT DEM elevation, Geomorpho90m slope, TPI, TRI, CTI, northness and eastness, SoilGrids
sand, clay and pH (15-30 cm, gap-filled), CHELSA NPP 1981-2010, a regional TPI (elevation minus
the 2 km focal mean of the DEM at 500 m) and the pixel longitude and latitude, declared at 90 m
and averaged to 300 m.

Sources needed to recreate (all readable by any Earth Engine user):

- Earth Engine catalog: ``WORLDCLIM/V1/BIO`` (WorldClim 1.4), ``MERIT/DEM/v1_0_3``.
- Earth Engine community catalog (awesome-gee-community-catalog):
  ``projects/sat-io/open-datasets/global_ai/global_ai_yearly`` (Global Aridity Index v3) and
  ``projects/sat-io/open-datasets/Geomorpho90m/{slope,tpi,tri,cti,northness,eastness}``.
- Landler assets, public read, as used for the published pool:
  ``projects/landler-open-data/assets/datasets/soilgrids/{sand,clay,phh2o}/
  {sand,clay,phh2o}_15-30cm_mean_gapfilled`` (SoilGrids 2.0, gap-filled),
  ``projects/landler-open-data/assets/datasets/chelsa/npp/chelsa_npp_1981_2010_v2-1``
  (CHELSA V2.1 NPP) and ``projects/landler-open-data/assets/datasets/clms_npp/v1_1/annual``
  (annual sums and SD of CLMS DMP 300 m, NPP v1.1). The NPP v1.1 collection is no longer
  published, so this step cannot be re-run as is; a rebuild starts from the published
  candidate pool.

Example:
    python pipelines/training/sample_candidate_pool_ee.py points.csv.gz pool.csv.gz \\
        --max-chunks 1   # one chunk first
    python pipelines/training/sample_candidate_pool_ee.py points.csv.gz pool.csv.gz
"""

from __future__ import annotations

import argparse
import math
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

LANDLER = "projects/landler-open-data/assets/datasets"
SOILGRIDS = LANDLER + "/soilgrids/{name}/{name}_15-30cm_mean_gapfilled"
CHELSA_NPP = LANDLER + "/chelsa/npp/chelsa_npp_1981_2010_v2-1"
CLMS_NPP_ANNUAL = LANDLER + "/clms_npp/v1_1/annual"
GEOMORPHO90M = "projects/sat-io/open-datasets/Geomorpho90m/"
PREDICTORS = [
    "mean_annual_temp",
    "temp_seasonality",
    "annual_precip",
    "precip_seasonality",
    "aridity",
    "elevation",
    "slope",
    "tpi",
    "tri",
    "cti",
    "northness",
    "eastness",
    "sand",
    "clay",
    "ph",
    "chelsa_npp",
    "tpi_regional",
    "latitude",
    "longitude",
]
YEARS = range(2014, 2026)
ANNUAL_COLUMNS = [f"{band}_{year}" for year in YEARS for band in ("sum", "sd")]
POINT_COLUMNS = [
    "longitude_point",
    "latitude_point",
    "pool_biome",
    "biome",
    "realm",
    "worldcover",
    "hmi_2022",
    "pristine_mask",
    "hmi_natural_mask",
    "wdpa_high_protection",
]
MASKED = -9999
RESOLUTION_M = 300


def predictor_stack(ee: Any) -> Any:
    """The 19 predictors in ``PREDICTORS`` order on the 300 m EPSG:4326 grid.

    Args:
        ee: The initialised ``ee`` module.

    Returns:
        ``ee.Image`` with one band per predictor.
    """
    bio = ee.Image("WORLDCLIM/V1/BIO")
    dem = ee.Image("MERIT/DEM/v1_0_3")
    elevation = dem.select("dem").rename("elevation")
    regional_mean = dem.reproject(crs="EPSG:4326", scale=500).focal_mean(2000, "circle", "meters")
    bands = [
        bio.select("bio01").rename("mean_annual_temp"),
        bio.select("bio04").rename("temp_seasonality"),
        bio.select("bio12").rename("annual_precip"),
        bio.select("bio15").rename("precip_seasonality"),
        ee.Image("projects/sat-io/open-datasets/global_ai/global_ai_yearly")
        .select("b1")
        .rename("aridity"),
        elevation,
        *[
            ee.ImageCollection(GEOMORPHO90M + name).mosaic().rename(name)
            for name in ("slope", "tpi", "tri", "cti", "northness", "eastness")
        ],
        ee.Image(CHELSA_NPP).rename("chelsa_npp"),
        ee.Image(SOILGRIDS.format(name="sand")).rename("sand"),
        ee.Image(SOILGRIDS.format(name="clay")).rename("clay"),
        ee.Image(SOILGRIDS.format(name="phh2o")).rename("ph"),
        elevation.subtract(regional_mean).rename("tpi_regional"),
        ee.Image.pixelLonLat().rename(["longitude", "latitude"]),
    ]
    return (
        ee.Image.cat(bands)
        .setDefaultProjection(crs="EPSG:4326", scale=90)
        .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=1024)
        .reproject(crs="EPSG:4326", scale=RESOLUTION_M)
        .select(PREDICTORS)
    )


def annual_npp(ee: Any) -> tuple[Any, Any]:
    """Annual NPP sum and SD per year, negative values masked and then set to ``MASKED``.

    Args:
        ee: The initialised ``ee`` module.

    Returns:
        The ``ee.Image`` and the native CLMS projection.
    """
    tiles = ee.ImageCollection(CLMS_NPP_ANNUAL)
    stack = ee.Image.cat(
        [
            tiles.filterDate(f"{year}-01-01", f"{year + 1}-01-01")
            .mosaic()
            .select(["b1", "b2"], [f"sum_{year}", f"sd_{year}"])
            for year in YEARS
        ]
    )
    return stack.updateMask(stack.gte(0)).unmask(MASKED), tiles.first().projection()


def sample_chunk(ee: Any, chunk: pd.DataFrame) -> pd.DataFrame:
    """Predictors and annual NPP at the chunk's points, indexed like ``chunk``.

    Args:
        ee: The initialised ``ee`` module.
        chunk: Points with ``longitude_point`` and ``latitude_point``, indexed by ``point_id``.

    Returns:
        ``PREDICTORS`` (NaN where masked) and ``ANNUAL_COLUMNS`` (NaN where not returned).
    """
    points = ee.FeatureCollection(
        [
            ee.Feature(ee.Geometry.Point([float(lon), float(lat)]), {"point_id": int(i)})
            for i, lon, lat in zip(
                chunk.index, chunk["longitude_point"], chunk["latitude_point"], strict=True
            )
        ]
    )
    annual, projection = annual_npp(ee)
    parts = [
        predictor_stack(ee).unmask(MASKED).sampleRegions(points, scale=RESOLUTION_M, tileScale=16),
        annual.sampleRegions(points, projection=projection, tileScale=16),
    ]
    frames = [
        pd.DataFrame([f["properties"] for f in part.getInfo()["features"]]).set_index("point_id")
        for part in parts
    ]
    out = frames[0].join(frames[1], how="outer").reindex(chunk.index)
    out = out.reindex(columns=[*PREDICTORS, *ANNUAL_COLUMNS])
    out[PREDICTORS] = out[PREDICTORS].replace(MASKED, np.nan)
    return out


def extract(
    ee: Any, points: pd.DataFrame, cache: Path, chunk_size: int, workers: int, limit: int | None
) -> None:
    """Sample the chunks of ``points`` that are not yet cached in ``cache``.

    Args:
        ee: The initialised ``ee`` module.
        points: Points indexed by ``point_id``.
        cache: Folder of the per-chunk CSVs.
        chunk_size: Points per request.
        workers: Concurrent requests.
        limit: Stop after this many chunks.
    """
    cache.mkdir(parents=True, exist_ok=True)
    n_chunks = math.ceil(len(points) / chunk_size)
    todo = [k for k in range(n_chunks) if not (cache / f"{k:05d}.csv.gz").exists()][:limit]

    def run(k: int) -> None:
        start = time.perf_counter()
        try:
            values = sample_chunk(ee, points.iloc[k * chunk_size : (k + 1) * chunk_size])
        except Exception as error:  # noqa: BLE001 - logged; the chunk stays to do
            print(f"chunk {k} failed: {error}", flush=True)
            return
        tmp = cache / f"{k:05d}.part"
        values.to_csv(tmp, compression="gzip")
        tmp.rename(cache / f"{k:05d}.csv.gz")
        print(f"chunk {k}: {len(values)} rows, {time.perf_counter() - start:.1f} s", flush=True)

    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(run, todo))


def merge(points: pd.DataFrame, cache: Path) -> pd.DataFrame | None:
    """The candidate pool from the cached chunks, or None while chunks are missing.

    Args:
        points: Points indexed by ``point_id``.
        cache: Folder of the per-chunk CSVs.

    Returns:
        Point columns, ``PREDICTORS``, ``ANNUAL_COLUMNS`` (int64, ``MASKED`` where missing) and
        the remaining point columns, indexed by ``point_id``.
    """
    chunks = sorted(cache.glob("*.csv.gz"))
    if not chunks:
        return None
    values = pd.concat([pd.read_csv(p, index_col="point_id") for p in chunks])
    if not points.index.equals(values.index.sort_values()):
        return None
    pool = points.join(values)
    pool[ANNUAL_COLUMNS] = pool[ANNUAL_COLUMNS].fillna(MASKED).astype("int64")
    head = [*POINT_COLUMNS, *PREDICTORS, *ANNUAL_COLUMNS]
    return pool[head + [c for c in points.columns if c not in head]]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("points", type=Path, help="points.csv.gz (index point_id)")
    parser.add_argument("out", type=Path, help="Output candidate pool CSV")
    parser.add_argument("--cache", type=Path, help="Chunk folder (default: next to out)")
    parser.add_argument("--chunk-size", type=int, default=2000, help="Points per request")
    parser.add_argument("--workers", type=int, default=4, help="Concurrent requests")
    parser.add_argument("--max-chunks", type=int, help="Stop after this many chunks")
    parser.add_argument("--project", help="Earth Engine Cloud project (default: from ee config)")
    args = parser.parse_args()

    import ee

    ee.Initialize(project=args.project)
    points = pd.read_csv(args.points, index_col="point_id")
    cache = args.cache or args.out.parent / f"ee_chunks_{args.chunk_size}"
    extract(ee, points, cache, args.chunk_size, args.workers, args.max_chunks)
    pool = merge(points, cache)
    if pool is None:
        print(f"chunks missing in {cache}; rerun to finish")
        return
    pool.to_csv(args.out, compression="gzip")
    print(f"wrote {len(pool):,} rows to {args.out}")


if __name__ == "__main__":
    main()
