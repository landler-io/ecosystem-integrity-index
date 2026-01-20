#!/bin/bash
#
# Upload HMI Natural Areas mask tiles to Google Earth Engine
#
# Files are staged through Google Cloud Storage before ingestion.
#
# Usage:
#   ./upload.sh /path/to/geotiffs
#   ./upload.sh /path/to/geotiffs --dry-run
#
set -euo pipefail

# =============================================================================
# CONFIGURATION - Adjust these parameters for your upload
# =============================================================================

# GEE asset path for the ImageCollection
ASSET_COLLECTION="projects/landler-open-data/assets/datasets/natural_lands/hmi_masks/v1"

# GCS bucket for staging (files are uploaded here first, then ingested to EE)
# The bucket must be accessible by your Earth Engine account
GCS_STAGING_BUCKET="gs://tlg-science-b1/hmi_masks_upload"

# HMI mask generation parameters (stored as image properties)
HMI_THRESHOLD=0.05
HIGH_MOD_THRESHOLD=0.2
BUFFER_PIXELS=10
MIN_AREA_KM2=5

# =============================================================================
# SCRIPT LOGIC
# =============================================================================

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

usage() {
    echo "Usage: $0 <source_directory> [--dry-run]"
    echo ""
    echo "Arguments:"
    echo "  source_directory   Directory containing GeoTIFF files to upload"
    echo "  --dry-run          Print commands without executing them"
    echo ""
    echo "Example:"
    echo "  $0 ~/Downloads/hmi_masks/"
    echo "  $0 ~/Downloads/hmi_masks/ --dry-run"
    echo ""
    echo "Prerequisites:"
    echo "  - gsutil configured (gcloud auth login)"
    echo "  - earthengine authenticated (earthengine authenticate)"
    echo "  - Write access to GCS_STAGING_BUCKET: $GCS_STAGING_BUCKET"
    exit 1
}

# Parse arguments
if [[ $# -lt 1 ]]; then
    usage
fi

SOURCE_DIR="$1"
DRY_RUN=false

if [[ $# -ge 2 && "$2" == "--dry-run" ]]; then
    DRY_RUN=true
    echo -e "${YELLOW}=== DRY RUN MODE ===${NC}"
fi

# Validate source directory
if [[ ! -d "$SOURCE_DIR" ]]; then
    echo -e "${RED}Error: Directory not found: $SOURCE_DIR${NC}"
    exit 1
fi

# Check for gsutil
if ! command -v gsutil &>/dev/null; then
    echo -e "${RED}Error: gsutil not found. Install Google Cloud SDK.${NC}"
    exit 1
fi

# Count GeoTIFF files
TIFF_COUNT=$(find "$SOURCE_DIR" -maxdepth 1 \( -name "*.tif" -o -name "*.tiff" \) | wc -l)
if [[ $TIFF_COUNT -eq 0 ]]; then
    echo -e "${RED}Error: No GeoTIFF files found in $SOURCE_DIR${NC}"
    exit 1
fi

echo -e "${GREEN}Found $TIFF_COUNT GeoTIFF files to upload${NC}"
echo "Target collection: $ASSET_COLLECTION"
echo "GCS staging: $GCS_STAGING_BUCKET"
echo ""

# Create the ImageCollection if it doesn't exist
echo "Checking if collection exists..."
if $DRY_RUN; then
    echo "[DRY RUN] Would create collection: $ASSET_COLLECTION"
else
    if ! earthengine asset info "$ASSET_COLLECTION" &>/dev/null; then
        echo "Creating ImageCollection: $ASSET_COLLECTION"
        # Create parent folders if needed
        PARENT_PATH=$(dirname "$ASSET_COLLECTION")
        earthengine create folder "$PARENT_PATH" 2>/dev/null || true
        earthengine create collection "$ASSET_COLLECTION"
        echo -e "${GREEN}Collection created${NC}"
    else
        echo "Collection already exists"
    fi
fi

# Build list of files to process (skip those with existing assets)
echo ""
echo "=========================================="
echo "Step 1: Checking for existing assets..."
echo "=========================================="

declare -a FILES_TO_UPLOAD
SKIPPED=0

for TIFF_FILE in "$SOURCE_DIR"/*.tif "$SOURCE_DIR"/*.tiff; do
    [[ -e "$TIFF_FILE" ]] || continue

    FILENAME=$(basename "$TIFF_FILE")

    # Extract tile ID from filename
    if [[ "$FILENAME" =~ lat_([0-9n_-]+)_lon_([0-9n_-]+) ]] || \
       [[ "$FILENAME" =~ lat__([0-9_-]+)_lon__([0-9_-]+) ]]; then
        LAT_PART="${BASH_REMATCH[1]}"
        LON_PART="${BASH_REMATCH[2]}"
        TILE_ID="lat_${LAT_PART}_lon_${LON_PART}"
    else
        TILE_ID="${FILENAME%.*}"
        TILE_ID=$(echo "$TILE_ID" | sed 's/[^a-zA-Z0-9_-]/_/g')
    fi

    ASSET_ID="${ASSET_COLLECTION}/${TILE_ID}"

    # Check if asset already exists
    if earthengine asset info "$ASSET_ID" &>/dev/null; then
        SKIPPED=$((SKIPPED + 1))
    else
        FILES_TO_UPLOAD+=("$TIFF_FILE")
    fi
done

echo "Skipping $SKIPPED files (assets already exist)"
echo "Processing ${#FILES_TO_UPLOAD[@]} new files"

if [[ ${#FILES_TO_UPLOAD[@]} -eq 0 ]]; then
    echo -e "${GREEN}All assets already exist. Nothing to do.${NC}"
    exit 0
fi

echo ""
echo "=========================================="
echo "Step 2: Uploading new files to GCS staging..."
echo "=========================================="

if $DRY_RUN; then
    echo "[DRY RUN] Would upload ${#FILES_TO_UPLOAD[@]} files to $GCS_STAGING_BUCKET/"
else
    for FILE in "${FILES_TO_UPLOAD[@]}"; do
        gsutil cp "$FILE" "$GCS_STAGING_BUCKET/"
    done
    echo -e "${GREEN}GCS upload complete${NC}"
fi

echo ""
echo "=========================================="
echo "Step 3: Ingesting from GCS to Earth Engine..."
echo "=========================================="

CURRENT=0
FAILED=0
TOTAL=${#FILES_TO_UPLOAD[@]}

for TIFF_FILE in "${FILES_TO_UPLOAD[@]}"; do
    CURRENT=$((CURRENT + 1))
    FILENAME=$(basename "$TIFF_FILE")

    # Extract tile ID from filename
    if [[ "$FILENAME" =~ lat_([0-9n_-]+)_lon_([0-9n_-]+) ]] || \
       [[ "$FILENAME" =~ lat__([0-9_-]+)_lon__([0-9_-]+) ]]; then
        LAT_PART="${BASH_REMATCH[1]}"
        LON_PART="${BASH_REMATCH[2]}"
        TILE_ID="lat_${LAT_PART}_lon_${LON_PART}"
    else
        TILE_ID="${FILENAME%.*}"
        TILE_ID=$(echo "$TILE_ID" | sed 's/[^a-zA-Z0-9_-]/_/g')
    fi

    ASSET_ID="${ASSET_COLLECTION}/${TILE_ID}"
    GCS_PATH="${GCS_STAGING_BUCKET}/${FILENAME}"

    echo -e "\n[$CURRENT/$TOTAL] ${YELLOW}$FILENAME${NC}"
    echo "  → Asset ID: $TILE_ID"

    if $DRY_RUN; then
        echo "  [DRY RUN] earthengine upload image --asset_id=$ASSET_ID ... $GCS_PATH"
    else
        if earthengine upload image \
            --asset_id="$ASSET_ID" \
            --property "hmi_threshold=$HMI_THRESHOLD" \
            --property "high_mod_threshold=$HIGH_MOD_THRESHOLD" \
            --property "buffer_pixels=$BUFFER_PIXELS" \
            --property "min_area_km2=$MIN_AREA_KM2" \
            --property "tile_id=$TILE_ID" \
            --property "source_file=$FILENAME" \
            "$GCS_PATH"; then
            echo -e "  ${GREEN}Ingestion task started${NC}"
        else
            echo -e "  ${RED}Upload failed${NC}"
            FAILED=$((FAILED + 1))
        fi
    fi
done

echo ""
echo "=========================================="
echo -e "${GREEN}All ingestion tasks submitted!${NC}"
echo "  Skipped (existing): $SKIPPED"
echo "  Processed: $TOTAL"
echo "  Failed: $FAILED"
echo ""
echo "Check ingestion status with:"
echo "  earthengine task list"
echo ""
echo "Once complete, verify the collection with:"
echo "  earthengine asset info $ASSET_COLLECTION"
echo ""
echo "To clean up GCS staging files after ingestion completes:"
echo "  gsutil -m rm ${GCS_STAGING_BUCKET}/*"
