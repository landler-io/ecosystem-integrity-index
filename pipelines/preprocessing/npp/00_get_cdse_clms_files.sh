#!/bin/bash

## This script downloads the raw CLMS data from the CDSE AWS endpoint to a local directory.
## Since the EODATA AWS endpoint does not recursively list deep files, we need to use a different approach.

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config.yaml"

# Check if config file exists
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: Configuration file not found at ${CONFIG_FILE}" >&2
    exit 1
fi

# Source the config file directly (no parser needed!)
# This loads all variables from the config file
. "${CONFIG_FILE}"

# Export AWS endpoint URL (not in config file as requested)
export AWS_ENDPOINT_URL=https://eodata.dataspace.copernicus.eu/

# Export variables that need to be exported
export AWS_PROFILE
export OUTPUT_DIR="${CLMS_RAW_DATA_DIR}"
export SLEEP_TIME

cd "${OUTPUT_DIR}"
for year_dir in $(aws s3 --profile "$AWS_PROFILE" ls "${S3_BASE_PATH}/" | awk '{print $2}'); do
    echo "--- Processing Year: ${year_dir} ---"
    s3_year_path="${S3_BASE_PATH}/${year_dir}"
    local_year_path="${OUTPUT_DIR}/${year_dir}"
    mkdir -p "${local_year_path}"
    for month_dir in $(aws s3 --profile "$AWS_PROFILE" ls "${s3_year_path}" | awk '{print $2}'); do
        s3_month_path="${s3_year_path}${month_dir}"
        local_month_path="${local_year_path}${month_dir}"
        mkdir -p "${local_month_path}"
        aws s3 --profile "$AWS_PROFILE" sync --exclude "*_nc*" --exclude "*RT[01234]_*" "${s3_month_path}" "${local_month_path}"
        sleep "${SLEEP_TIME}"
    done
done
