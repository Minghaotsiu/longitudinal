#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 RAW_DIRECTORY BATCH_DIRECTORY EXPECTED_PAIRS [SAMPLE_REGEX]" >&2
  exit 2
fi

RAW_DIRECTORY=$1
BATCH_DIRECTORY=$2
EXPECTED_PAIRS=$3
SAMPLE_REGEX=${4:-'(d[0-3]-[0-9]+-1)'}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PYTHON_BIN=${PYTHON_BIN:-python3}

: "${STAR_INDEX:?Set STAR_INDEX}"
: "${SOYBEAN_GTF:?Set SOYBEAN_GTF}"
: "${SOYBEAN_EXON_BED:?Set SOYBEAN_EXON_BED}"
: "${SOYBEAN_GENE_BED:?Set SOYBEAN_GENE_BED}"

"$PYTHON_BIN" "$SCRIPT_DIR/prepare_delivery.py" \
  "$RAW_DIRECTORY" \
  "$BATCH_DIRECTORY" \
  "$EXPECTED_PAIRS" \
  --pattern "$SAMPLE_REGEX"

cd "$BATCH_DIRECTORY"
"$PYTHON_BIN" "$SCRIPT_DIR/run_soybean_upstream.py"

BATCH_NAME=$(basename "$BATCH_DIRECTORY")
cp "fenxi/${BATCH_NAME}_raw_counts_matrix.csv" "ALL_17_raw_counts_matrix.csv"
cp "fenxi/1M/${BATCH_NAME}_1M_counts_matrix.csv" "ALL_17_1M_counts_matrix.csv"
cp "fenxi/${BATCH_NAME}_FPKM_RNA_MAP.csv" "ALL_17_RNA_MAP.csv"

