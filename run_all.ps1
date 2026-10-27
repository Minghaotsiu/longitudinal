$ErrorActionPreference = "Stop"

if (Test-Path -LiteralPath "config.ps1") {
    . "$PSScriptRoot/config.ps1"
}

$python = if ($env:PYTHON_BIN) { $env:PYTHON_BIN } else { "python" }
$rscript = if ($env:RSCRIPT_BIN) { $env:RSCRIPT_BIN } else { "Rscript" }

& $python "scripts/00_validate_inputs.py"
& $rscript "scripts/00_check_environment.R"
& $rscript "scripts/05_run_selected5_analysis.R"
& $rscript "scripts/08_build_main_figures.R"
& $python "scripts/09_build_report.py"
& $python "tests/test_repository.py"

