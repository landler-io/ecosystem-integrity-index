#!/bin/bash
# This script uploads the annual NPP sums to GEE as image assets (via GCS).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config.yaml"

if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: Configuration file not found at ${CONFIG_FILE}" >&2
    exit 1
fi

. "${CONFIG_FILE}"

EE_COLLECTION_NAME="${EE_COLLECTION_NAME%/}"
LOCAL_FILE_PATH="${NPP_OUTPUT_DIR}"
EE_COLLECTION_PATH="projects/${EE_PROJECT}/assets/datasets/${EE_COLLECTION_NAME}"

make_ee_dirs() {
  local path="$1"
  local parts=(${path//\// })
  local current=""
  for i in "${!parts[@]}"; do
    current="${current}${parts[$i]}"
    if [[ "${parts[$i]}" == "assets" ]]; then
      for j in $(seq $((i+1)) $((${#parts[@]}-1))); do
        current="${current}/${parts[$j]}"
        earthengine create folder "${current}" 2>/dev/null || true
      done
      break
    fi
    current="${current}/"
  done
}

make_ee_dirs "${EE_COLLECTION_PATH}"
earthengine create collection ${EE_COLLECTION_PATH} 2>/dev/null || true

for year in $(seq ${START_YEAR} ${END_YEAR}); do
  ASSET_NAME="npp_${year}"
  LOCAL_FILENAME="annual_npp_${year}.tif"
  LOCAL_FULL_PATH="${LOCAL_FILE_PATH}/${LOCAL_FILENAME}"
  GCS_URI="gs://${GCS_BUCKET_NAME}/${EE_COLLECTION_NAME}/${LOCAL_FILENAME}"
  ASSET_ID="${EE_COLLECTION_PATH}/${ASSET_NAME}"

  MAX_RETRIES=3
  RETRY_COUNT=0
  while [ $RETRY_COUNT -lt $MAX_RETRIES ]; do
    if gsutil -m cp "${LOCAL_FULL_PATH}" "${GCS_URI}" 2>/dev/null; then
      break
    else
      RETRY_COUNT=$((RETRY_COUNT + 1))
      if [ $RETRY_COUNT -eq $MAX_RETRIES ]; then
        echo "ERROR: Failed to upload ${LOCAL_FILENAME} after ${MAX_RETRIES} attempts" >&2
        continue 2
      fi
      sleep 5
    fi
  done

  earthengine upload image --asset_id="${ASSET_ID}" \
    --pyramiding_policy=mean \
    --time_start="${year}-01-01" \
    --time_end="${year}-12-31" \
    --property "scale_factor=100" \
    --property "units=g C / m² / year (scaled)" \
    --property "source=CLMS DMP 300m 10-daily" \
    --property "processing=annual_sum_std_linear_interpolation" \
    --property "band_1=annual_sum" \
    --property "band_2=annual_sd" \
    "${GCS_URI}"
done
