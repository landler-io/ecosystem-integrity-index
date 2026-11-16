#!/bin/bash
set -euo pipefail

## This script downloads the raw CLMS data from the CDSE AWS endpoint to a local directory.
## Since the EODATA AWS endpoint does not recursively list deep files, we need to use a different approach.
## Rerunning resumes: `aws s3 sync` skips files that are already complete.

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config.cfg"

# Check if config file exists
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: Configuration file not found at ${CONFIG_FILE}" >&2
    exit 1
fi

# Source the config file directly (no parser needed!)
# This loads all variables from the config file
. "${CONFIG_FILE}"

DOWNLOAD_START_YEAR="${DOWNLOAD_START_YEAR:-${START_YEAR}}"
DOWNLOAD_END_YEAR="${DOWNLOAD_END_YEAR:-${END_YEAR}}"

# CDSE allows 4 concurrent S3 connections per user and answers 429 above that.
# The AWS CLI default is 10, so use a generated config instead of ~/.aws/config.
# Keys still come from the shared credentials file (~/.aws/credentials).
AWS_CONFIG_FILE="$(mktemp)"
trap 'rm -f "${AWS_CONFIG_FILE}"' EXIT
cat > "${AWS_CONFIG_FILE}" <<EOF
[profile ${AWS_PROFILE}]
region = default
endpoint_url = https://eodata.dataspace.copernicus.eu
retry_mode = adaptive
max_attempts = 10
s3 =
    max_concurrent_requests = ${S3_MAX_CONCURRENT_REQUESTS}
    multipart_chunksize = 64MB
EOF

# Export variables that need to be exported
export AWS_CONFIG_FILE
export AWS_PROFILE
export OUTPUT_DIR="$(realpath -m "${CLMS_RAW_DATA_DIR}")"
export SLEEP_TIME
export PRODUCT_GLOB

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Parallel runs would exceed the connection limit and race on temporary files.
exec 9> "${OUTPUT_DIR}/.download.lock"
if ! flock -n 9; then
    echo "ERROR: another download is running in ${OUTPUT_DIR}" >&2
    exit 1
fi

# An interrupted sync leaves temporary files named like NAME.tiff.1a2B3c4D.
find "${OUTPUT_DIR}" -regextype posix-extended -type f \
    -regex '.*\.tiff?\.[0-9A-Za-z]{8}' -print -delete

# Run a command, retrying with exponential backoff when CDSE throttles or fails.
with_retries() {
    local attempt=1 delay=30
    until "$@"; do
        if [[ "$attempt" -ge "$MAX_RETRIES" ]]; then
            echo "ERROR: failed after ${attempt} attempts: $*" >&2
            return 1
        fi
        echo "Attempt ${attempt} failed; retrying in ${delay}s: $*" >&2
        sleep "$delay"
        attempt=$((attempt + 1))
        delay=$((delay * 2 > 600 ? 600 : delay * 2))
    done
}

# Capture listings before iterating so API failures stop the script.
list_prefixes() {
    local listing
    listing=$(with_retries aws s3 --profile "$AWS_PROFILE" ls "$1")
    awk '$1 == "PRE" {print $2}' <<< "$listing"
}

for year in $(seq "${DOWNLOAD_START_YEAR}" "${DOWNLOAD_END_YEAR}"); do
    year_dir="${year}/"
    echo "--- Processing Year: ${year_dir} ---"
    s3_year_path="${S3_BASE_PATH}/${year_dir}"
    local_year_path="${OUTPUT_DIR}/${year_dir}"
    mkdir -p "${local_year_path}"
    n_products=0
    month_dirs=$(list_prefixes "${s3_year_path}")
    for month_dir in $month_dirs; do
        s3_month_path="${s3_year_path}${month_dir}"
        local_month_path="${local_year_path}${month_dir}"
        mkdir -p "${local_month_path}"
        # CDSE only lists immediate children at these levels. Sync each
        # product separately, where the listing contains actual file objects.
        day_dirs=$(list_prefixes "${s3_month_path}")
        for day_dir in $day_dirs; do
            s3_day_path="${s3_month_path}${day_dir}"
            product_dirs=$(list_prefixes "${s3_day_path}")
            for product_dir in $product_dirs; do
                # shellcheck disable=SC2053  # PRODUCT_GLOB is a pattern.
                [[ "$product_dir" == $PRODUCT_GLOB ]] || continue
                with_retries aws s3 --profile "$AWS_PROFILE" sync --only-show-errors \
                    "${s3_day_path}${product_dir}" \
                    "${local_month_path}${day_dir}${product_dir}"
                n_products=$((n_products + 1))
            done
        done
        echo "${year}/${month_dir%/}: done (${n_products} products so far)"
        sleep "${SLEEP_TIME}"
    done
    echo "${year}: ${n_products} products matching ${PRODUCT_GLOB}"
    if [[ "$n_products" -eq 0 ]]; then
        echo "ERROR: no products matched for ${year}" >&2
        exit 1
    fi
done
