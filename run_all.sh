#!/usr/bin/env bash
set -euo pipefail

if [ -f config.env ]; then
  set -a
  # shellcheck disable=SC1091
  source config.env
  set +a
fi

RSCRIPT_BIN=${RSCRIPT_BIN:-Rscript}
PYTHON_BIN=${PYTHON_BIN:-python3}

"$PYTHON_BIN" scripts/00_validate_inputs.py
"$RSCRIPT_BIN" scripts/00_check_environment.R
"$RSCRIPT_BIN" scripts/05_run_selected5_analysis.R
"$RSCRIPT_BIN" scripts/08_build_main_figures.R
"$PYTHON_BIN" scripts/09_build_report.py
"$PYTHON_BIN" tests/test_repository.py

