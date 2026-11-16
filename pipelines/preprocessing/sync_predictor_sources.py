#!/usr/bin/env python3
"""Mirror the v2 potential-NPP predictor sources from GCS to a local directory.

This is where a reproduction starts: every source below is published in the requester-pays
bucket ``gs://landler-opendata-gee-assets``.

- ``ee_export``: terrain tiles and upstream area exported by ``export_predictor_stack.py run``;
- ``chelsa_<var>``: CHELSA V2.1 1981-2010 bio1, bio4, bio12, bio15 and PET
  (``ingest_chelsa_v2_1.py``);
- ``chelsa_npp``: CHELSA V2.1 NPP 1981-2010 (``CHELSA_npp_1981-2010_V.2.1.tif`` next to the
  bioclim files of ``ingest_chelsa_v2_1.BASE_URL``; CHELSA-BIOCLIM+, Brun et al. 2022,
  doi:10.16904/envidat.332), uploaded to Earth Engine before this pipeline and exported to the
  bucket; there is no ingest script;
- ``osm_coast_distance``, ``eot20_tidal_range``, ``woa23_sss``: ``ingest_coastal_v2.py``;
- ``ecoregions``: RESOLVE Ecoregions 2017 (Dinerstein et al. 2017, Earth Engine
  ``RESOLVE/ECOREGIONS/2017``) rasterised in Earth Engine to ``ECO_ID`` (first polygon per
  pixel, 0 outside, ``scale=5000`` in EPSG:4326, -60 to 90 deg); there is no ingest script.

Re-running downloads only objects that are missing locally or differ in size. Files arrive as
``.part`` and are renamed after the client's checksum validation. Uses application default
credentials; the requester-pays reads bill ``EII_BILLING_PROJECT`` or ``--billing-project``.

Example:
    python pipelines/preprocessing/sync_predictor_sources.py predictor_stack \\
        --billing-project MY_PROJECT
"""

from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from google.cloud import storage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from local_sources import BILLING_ENV, BUCKET, billing_project  # noqa: E402

SOURCES = {
    "ee_export": "products/eii/v2_0/predictor_stack/ee_export/",
    "chelsa_npp": "datasets/chelsa_npp/v1_0/",
    **{
        f"chelsa_{var}": f"datasets/chelsa_{var}/v1_0/"
        for var in ("bio1", "bio4", "bio12", "bio15", "pet_penman_mean")
    },
    "osm_coast_distance": "datasets/osm_coast_distance/v1_0/",
    "eot20_tidal_range": "datasets/eot20_tidal_range/v1_0/",
    "woa23_sss": "datasets/woa23_sss/v1_0/",
    # RESOLVE ecoregions (biome per pixel for the interval calibration in inference)
    "ecoregions": "datasets/ecoregions/v1_0/",
}


def download(blob: storage.Blob, target: Path) -> int:
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + ".part")
    blob.download_to_filename(part, checksum="auto")
    part.rename(target)
    return blob.size


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument(
        "--billing-project",
        help="Google Cloud project billed for the reads (default: EII_BILLING_PROJECT)",
    )
    args = parser.parse_args()
    if args.billing_project:
        os.environ[BILLING_ENV] = args.billing_project

    project = billing_project()
    client = storage.Client(project=project)
    bucket = client.bucket(BUCKET, user_project=project)
    todo = []
    for name, prefix in SOURCES.items():
        blobs = [b for b in client.list_blobs(bucket, prefix=prefix) if not b.name.endswith("/")]
        missing = []
        for blob in blobs:
            target = args.root / name / blob.name[len(prefix) :]
            if not target.exists() or target.stat().st_size != blob.size:
                missing.append((blob, target))
        total, need = sum(b.size for b in blobs), sum(b.size for b, _ in missing)
        print(
            f"{name}: {len(blobs)} objects {total / 1e9:.1f} GB, "
            f"to fetch {len(missing)} ({need / 1e9:.1f} GB)",
            flush=True,
        )
        todo += missing

    done = 0
    with ThreadPoolExecutor(args.workers) as pool:
        for i, size in enumerate(pool.map(lambda job: download(*job), todo), start=1):
            done += size
            if i % 10 == 0 or i == len(todo):
                print(f"{i}/{len(todo)} files, {done / 1e9:.1f} GB", flush=True)


if __name__ == "__main__":
    main()
