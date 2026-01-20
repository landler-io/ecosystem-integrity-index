#!/bin/bash
#
# Upload GloH2O Koppen-Geiger v3 raster to Google Earth Engine.
set -euo pipefail

# Target asset path (single image)
ASSET_ID="projects/landler-open-data/assets/datasets/climatezones/gloh2o-koeppen-v3"

# GCS staging bucket (optional; required for large files or non-local uploads)
DEFAULT_GCS_BUCKET=""

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

GEOTIFF_PATH="$1"
GCS_BUCKET="$DEFAULT_GCS_BUCKET"

shift
while [[ $# -gt 0 ]]; do
    case "$1" in
        --gcs)
            GCS_BUCKET="$2"
            shift 2
            ;;
        *)
            echo -e "${RED}Unknown argument: $1${NC}"
            exit 1
            ;;
    esac
done

if [[ ! -f "$GEOTIFF_PATH" ]]; then
    echo -e "${RED}Error: File not found: $GEOTIFF_PATH${NC}"
    exit 1
fi

echo "GeoTIFF: $GEOTIFF_PATH"
echo "Asset ID: $ASSET_ID"
if [[ -n "$GCS_BUCKET" ]]; then
    echo "GCS staging: $GCS_BUCKET"
else
    echo "GCS staging: (none, upload from local path)"
fi



PARENT_PATH=$(dirname "$ASSET_ID")
earthengine create folder "$PARENT_PATH" 2>/dev/null || true


UPLOAD_SOURCE="$GEOTIFF_PATH"
if [[ -n "$GCS_BUCKET" ]]; then
    BASENAME=$(basename "$GEOTIFF_PATH")
    GCS_PATH="${GCS_BUCKET%/}/${BASENAME}"
    gsutil cp "$GEOTIFF_PATH" "$GCS_PATH"
    UPLOAD_SOURCE="$GCS_PATH"
fi


earthengine upload image \
    --asset_id="$ASSET_ID" \
    --property "source=GloH2O Koppen-Geiger v3" \
    --property "citation=Beck et al. 2023" \
    "$UPLOAD_SOURCE"
