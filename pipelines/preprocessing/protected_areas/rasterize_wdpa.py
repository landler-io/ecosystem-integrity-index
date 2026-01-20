#!/usr/bin/env python3
"""
Rasterize WDPA high-protection areas
"""

import sys
from pathlib import Path

import ee

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.config_utils import load_config  # noqa: E402
from utils.utils_gee_assets import make_asset_dir  # noqa: E402


def rasterize_wdpa(config: dict, aoi: ee.Geometry = None) -> ee.batch.Task:
    asset_id = config["asset_id"]
    make_asset_dir("/".join(asset_id.split("/")[:-1]))

    if aoi is None:
        export_region = ee.Geometry.Rectangle([-180, -60, 180, 90], "EPSG:4326", False)
    else:
        export_region = aoi

    wdpa_full = ee.FeatureCollection(config["source_asset"])

    if aoi is not None:
        wdpa_full = wdpa_full.filterBounds(aoi)

    marine_filter = ee.Filter.Or(
        ee.Filter.eq("MARINE", "0"),
        ee.Filter.eq("MARINE", 0),
        ee.Filter.notNull(["MARINE"]).Not(),
    )

    status_filter = ee.Filter.Or(
        ee.Filter.inList("IUCN_CAT", config["iucn_categories"]),
        ee.Filter.eq("DESIG_ENG", "National Park"),
    )

    protected_areas = wdpa_full.filter(
        ee.Filter.And(
            status_filter,
            marine_filter,
            ee.Filter.gt("REP_AREA", config["min_area_km2"]),
        )
    )

    buffer_m = config.get("buffer_m", 0)
    if buffer_m != 0:

        def apply_buffer(feature):
            return feature.setGeometry(feature.geometry().buffer(buffer_m, 100))

        protected_areas = protected_areas.map(apply_buffer)
        protected_areas = protected_areas.filter(ee.Filter.notNull([".geo"]))

    # Binary Mask (1 = PA, 0 = non-PA)
    pa_mask = ee.Image(0).byte().paint(protected_areas, 1).rename("protected_class")

    pa_mask = pa_mask.clip(export_region).set(
        {
            "iucn_categories": ",".join(config["iucn_categories"]),
            "buffer_m": buffer_m,
            "description": "1=Protected, 0=Unprotected",
        }
    )

    task = ee.batch.Export.image.toAsset(
        image=pa_mask,
        description="Rasterize_WDPA_HighProt",
        assetId=asset_id,
        region=export_region,
        scale=config["scale"],
        crs="EPSG:4326",
        maxPixels=1e13,
    )
    task.start()
    return task


def main():
    ee.Initialize()
    config = load_config(Path(__file__).parent / "config.cfg")

    if isinstance(config.get("iucn_categories"), str):
        config["iucn_categories"] = [x.strip() for x in config["iucn_categories"].split(",")]

    task = rasterize_wdpa(config)
    print(f"Export task started: {task.id}")


if __name__ == "__main__":
    main()
