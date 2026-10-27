# Script map

| Script | Purpose |
|---|---|
| `00_check_environment.R` | Verify required R packages |
| `00_validate_inputs.py` | Verify processed matrices and metadata |
| `upstream/run_delivery.sh` | Prepare and process one delivered batch |
| `upstream/run_two_matched_batches.sh` | Process and combine the two 17-run deliveries |
| `03_run_batch_analysis.R` | Legacy single-batch analysis |
| `04_run_all34_analysis.R` | 34-run exploratory analysis |
| `05_run_selected5_analysis.R` | Main selected-5 analysis |
| `06_run_pseudotime_analysis.R` | Standalone trajectory analysis |
| `08_build_main_figures.R` | Generate main figures and source data |
| `09_build_report.py` | Generate integrated HTML report |
| `10_package_upload_bundle.py` | Assemble processed-data upload bundle |
| `11_audit_pdf_geometry.py` | Check PDF text clipping and overlap |

Run scripts from the repository root unless the script documentation says
otherwise.

