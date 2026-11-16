#!/usr/bin/env python3
"""Equal-area global land sample of the structural-integrity calibration.

Writes the sample read by ``calibrate_structural.py calibrate`` (``lon``, ``lat``, ``biome``,
``anthrome``), published as ``si_calibration_sample.csv.gz`` (29,585 points):

1. Draw ``N_DRAWS`` points uniformly by area in ``BOX`` (longitude and sine of latitude
   uniform; ``numpy.random.default_rng(seed)``). 111,111 draws with an expected land share of
   0.27 aim at 30,000 land points.
2. Land test in Earth Engine: ESA WorldCover 10 m v200 (2021, ``ESA/WorldCover/v200``)
   sampled with ``sampleRegions`` at ``scale=1000`` in its own EPSG:4326 grid; a point is land
   where the map is valid and not permanent water (80). The 1 km value comes from the
   WorldCover pyramid (mode), which is why this step runs in Earth Engine. A land point becomes
   the centre of its 1 km sampling cell (``geometries=True``), in the order Earth Engine
   returns them; the first 30,000 are kept.
3. Complete records in Earth Engine: the sample was drawn as a table of EII v1 and its inputs,
   and ``sampleRegions`` drops a point where any sampled layer is masked. ``record_image``
   holds those layers (``RECORD_LAYERS``); points where one of them has no data at 1 km are
   dropped (about 1.4 %: ice sheets, deserts, small islands and coasts). The v1
   functional-integrity intermediates of that table are left out; their masks follow EII v1
   and natural NPP v1.
4. ``biome``: RESOLVE Ecoregions 2017 (Dinerstein et al. 2017) rasterised to ``ECO_ID`` at
   ~5 km (``datasets/ecoregions/v1_0``, as synced by ``sync_predictor_sources.py``), remapped
   by ``pipelines/preprocessing/reference/eco_id_to_biome.json``; 0 where no ECO_ID maps.
5. ``anthrome``: Anthromes v2, year 2000 (Ellis et al. 2010, ``anthro2_a2000.tif`` in
   ``anthromes_2_GeoTIFF.zip``, Harvard Dataverse doi:10.7910/DVN/6KKAWS), 5 arc-minutes,
   nearest pixel; 0 outside the anthrome land mask.

Steps 1, 4 and 5 run offline; steps 2 and 3 are interactive Earth Engine requests of up to
5,000 points (no exports, roughly 25 and 7 requests).

Sources needed to recreate:

- Earth Engine catalog: ``ESA/WorldCover/v200``, ``WORLDCLIM/V1/BIO``, ``MERIT/DEM/v1_0_3``.
- Earth Engine community catalog: ``projects/sat-io/open-datasets/GHM/HM_2022_300M`` (threat
  ``AA``), ``projects/sat-io/open-datasets/global_ai/global_ai_yearly``,
  ``projects/sat-io/open-datasets/Geomorpho90m/slope``.
- Earth Blox BII: ``projects/ebx-data/assets/earthblox/IO/BIOINTACT`` (2020 mosaic).
- Landler assets, public read: ``eii/global/eii_global_v1``,
  ``eii/intermediate/functional/predictions/natural_npp_v1``,
  ``datasets/soilgrids/{sand,clay,phh2o}/*_15-30cm_mean_gapfilled``,
  ``datasets/chelsa/npp/chelsa_npp_1981_2010_v2-1`` and
  ``datasets/ecoregions/resolve_ecoregions_rasterized`` under
  ``projects/landler-open-data/assets/``.
- Offline: the RESOLVE ECO_ID raster from ``gs://landler-opendata-gee-assets`` (requester
  pays, ``sync_predictor_sources.py``) and ``anthro2_a2000.tif`` (doi:10.7910/DVN/6KKAWS).

Example:
    python pipelines/training/build_si_calibration_sample.py si_calibration_sample.csv.gz \\
        --ecoregions predictor_stack/ecoregions/resolve_ecoregions_rasterized/<file>.tif \\
        --anthromes anthro2_a2000.tif
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import rasterio

ECO_TO_BIOME = Path(__file__).resolve().parents[1] / "preprocessing/reference/eco_id_to_biome.json"
BOX = (-180.0, -60.0, 180.0, 75.0)
N_DRAWS = 111_111
TARGET_POINTS = 30_000
SEED = 42
WORLDCOVER = "ESA/WorldCover/v200"
LANDLER = "projects/landler-open-data/assets/"
RECORD_LAYERS = {
    "eii_v1": LANDLER + "eii/global/eii_global_v1",
    "natural_npp_v1": LANDLER + "eii/intermediate/functional/predictions/natural_npp_v1",
    "hmi_2022": "projects/sat-io/open-datasets/GHM/HM_2022_300M",
    "bii_2020": "projects/ebx-data/assets/earthblox/IO/BIOINTACT",
    "worldclim": "WORLDCLIM/V1/BIO",
    "aridity": "projects/sat-io/open-datasets/global_ai/global_ai_yearly",
    "elevation": "MERIT/DEM/v1_0_3",
    "slope": "projects/sat-io/open-datasets/Geomorpho90m/slope",
    "soilgrids": LANDLER + "datasets/soilgrids/{name}/{name}_15-30cm_mean_gapfilled",
    "chelsa_npp": LANDLER + "datasets/chelsa/npp/chelsa_npp_1981_2010_v2-1",
    "ecoregions": LANDLER + "datasets/ecoregions/resolve_ecoregions_rasterized",
}
WATER = 80
SCALE_M = 1000
CHUNK_SIZE = 5000


def draw_points(
    n: int = N_DRAWS, seed: int = SEED, box: Sequence[float] = BOX
) -> tuple[np.ndarray, np.ndarray]:
    """Uniform random points by area inside ``box`` (west, south, east, north).

    Args:
        n: Number of points.
        seed: Seed of ``numpy.random.default_rng``.
        box: Bounding box in degrees.

    Returns:
        Longitudes and latitudes in degrees.
    """
    rng = np.random.default_rng(seed)
    west, south, east, north = box
    lon = rng.uniform(west, east, n)
    sin_lat = rng.uniform(np.sin(np.radians(south)), np.sin(np.radians(north)), n)
    return lon, np.degrees(np.arcsin(sin_lat))


def land_points(lon: np.ndarray, lat: np.ndarray, chunk_size: int = CHUNK_SIZE) -> pd.DataFrame:
    """Earth Engine land test at 1 km; returns the 1 km cell centres of the land points.

    Args:
        lon: Longitudes of the drawn points.
        lat: Latitudes of the drawn points.
        chunk_size: Points per request (``getInfo`` returns at most 5,000 features).

    Returns:
        Frame with ``lon`` and ``lat`` of the land points.
    """
    import ee

    cover = ee.ImageCollection(WORLDCOVER).first().select("Map")
    land = cover.mask().And(cover.neq(WATER)).rename("land")
    kept: list[list[float]] = []
    for start in range(0, len(lon), chunk_size):
        part = slice(start, start + chunk_size)
        points = ee.FeatureCollection(
            [
                ee.Feature(ee.Geometry.Point([x, y]))
                for x, y in zip(lon[part].tolist(), lat[part].tolist(), strict=True)
            ]
        )
        sampled = land.sampleRegions(
            collection=points, scale=SCALE_M, geometries=True, tileScale=4
        ).getInfo()
        kept += [
            f["geometry"]["coordinates"] for f in sampled["features"] if f["properties"].get("land")
        ]
    return pd.DataFrame(kept, columns=["lon", "lat"])


def record_image(ee: Any) -> Any:
    """The layers whose masks decide a complete record, EII v1 first (its grid is sampled).

    Args:
        ee: The initialised ``ee`` module.

    Returns:
        ``ee.Image`` with one band per layer.
    """
    layers = RECORD_LAYERS
    hmi = ee.ImageCollection(layers["hmi_2022"]).filter(ee.Filter.eq("threat_code", "AA"))
    bii = ee.ImageCollection(layers["bii_2020"]).filterDate("2020-01-01", "2020-12-31")
    return ee.Image.cat(
        [
            ee.Image(layers["eii_v1"]),
            ee.Image(layers["natural_npp_v1"]),
            ee.Image(hmi.first()).select([0]),
            bii.mosaic(),
            ee.Image(layers["worldclim"]).select(["bio01", "bio04", "bio12", "bio15"]),
            ee.Image(layers["aridity"]).select("b1"),
            ee.Image(layers["elevation"]).select("dem"),
            ee.ImageCollection(layers["slope"]).mosaic(),
            *[ee.Image(layers["soilgrids"].format(name=n)) for n in ("sand", "clay", "phh2o")],
            ee.Image(layers["chelsa_npp"]),
            ee.Image(layers["ecoregions"]).select([0]),
        ]
    )


def complete_records(points: pd.DataFrame, chunk_size: int = CHUNK_SIZE) -> pd.DataFrame:
    """The points where every ``record_image`` layer has data at 1 km, in input order.

    Args:
        points: Frame with ``lon`` and ``lat``.
        chunk_size: Points per request.

    Returns:
        The complete-record rows of ``points``.
    """
    import ee

    image = record_image(ee)
    found: set[int] = set()
    for start in range(0, len(points), chunk_size):
        part = points.iloc[start : start + chunk_size]
        features = ee.FeatureCollection(
            [
                ee.Feature(ee.Geometry.Point([x, y]), {"row": i})
                for i, x, y in zip(
                    range(start, start + len(part)), part["lon"], part["lat"], strict=True
                )
            ]
        )
        sampled = image.sampleRegions(collection=features, scale=SCALE_M, tileScale=4).getInfo()
        found.update(f["properties"]["row"] for f in sampled["features"])
    return points.iloc[sorted(found)].reset_index(drop=True)


def sample_classes(path: Path, lon: np.ndarray, lat: np.ndarray, fill: int = 0) -> np.ndarray:
    """Nearest-pixel integer class of band 1 at the points; ``fill`` for nodata or outside.

    Args:
        path: Raster in EPSG:4326.
        lon: Longitudes.
        lat: Latitudes.
        fill: Value for nodata pixels and points outside the raster.

    Returns:
        Class codes as int64.
    """
    with rasterio.open(path) as src:
        data, nodata = src.read(1), src.nodata
        rows, cols = (np.asarray(a) for a in rasterio.transform.rowcol(src.transform, lon, lat))
    inside = (rows >= 0) & (rows < data.shape[0]) & (cols >= 0) & (cols < data.shape[1])
    values = data[rows[inside], cols[inside]]
    out = np.full(len(rows), fill, dtype="int64")
    out[inside] = values if nodata is None else np.where(values == nodata, fill, values)
    return out


def attach_classes(points: pd.DataFrame, ecoregions: Path, anthromes: Path) -> pd.DataFrame:
    """Add ``biome`` (RESOLVE, via ECO_ID) and ``anthrome`` (Anthromes v2, 2000) columns.

    Args:
        points: Frame with ``lon`` and ``lat``.
        ecoregions: RESOLVE ECO_ID raster.
        anthromes: Anthromes v2 2000 raster.

    Returns:
        ``points`` with ``biome`` and ``anthrome``.
    """
    lookup = {int(k): v for k, v in json.loads(ECO_TO_BIOME.read_text()).items()}
    lon, lat = points["lon"].to_numpy(), points["lat"].to_numpy()
    eco = sample_classes(ecoregions, lon, lat, fill=-1)
    return points.assign(
        biome=np.array([lookup.get(int(e), 0) for e in eco], dtype="int64"),
        anthrome=sample_classes(anthromes, lon, lat),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("out", type=Path, help="Output CSV (lon, lat, biome, anthrome)")
    parser.add_argument("--ecoregions", type=Path, required=True, help="RESOLVE ECO_ID raster")
    parser.add_argument("--anthromes", type=Path, required=True, help="anthro2_a2000.tif")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--n-draws", type=int, default=N_DRAWS)
    parser.add_argument("--project", help="Earth Engine Cloud project (default: from ee config)")
    args = parser.parse_args()

    import ee

    ee.Initialize(project=args.project)
    points = complete_records(
        land_points(*draw_points(args.n_draws, args.seed)).head(TARGET_POINTS)
    )
    attach_classes(points, args.ecoregions, args.anthromes).to_csv(args.out, index=False)
    print(f"wrote {len(points):,} points to {args.out}")


if __name__ == "__main__":
    main()
