#!/usr/bin/env bash
set -euo pipefail

: "${BATCH1_RAW_DIR:?Set BATCH1_RAW_DIR}"
: "${BATCH2_RAW_DIR:?Set BATCH2_RAW_DIR}"
: "${UPSTREAM_DIR:?Set UPSTREAM_DIR}"
: "${PROCESSED_DATA_DIR:?Set PROCESSED_DATA_DIR}"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PYTHON_BIN=${PYTHON_BIN:-python3}
RSCRIPT_BIN=${RSCRIPT_BIN:-Rscript}

bash "$SCRIPT_DIR/run_delivery.sh" \
  "$BATCH1_RAW_DIR" "$UPSTREAM_DIR/batch1" 17
bash "$SCRIPT_DIR/run_delivery.sh" \
  "$BATCH2_RAW_DIR" "$UPSTREAM_DIR/batch2" 17

export BATCH1_UPSTREAM_DIR="$UPSTREAM_DIR/batch1"
export BATCH2_UPSTREAM_DIR="$UPSTREAM_DIR/batch2"
mkdir -p "$PROCESSED_DATA_DIR"

"$PYTHON_BIN" "$SCRIPT_DIR/parse_gtf_symbols.py" \
  "$SOYBEAN_GTF" \
  "$PROCESSED_DATA_DIR/soybean_gene_symbol_map_gtf.csv"
"$RSCRIPT_BIN" "$SCRIPT_DIR/combine_batches.R"
"$RSCRIPT_BIN" "$SCRIPT_DIR/collect_all_qc.R"

