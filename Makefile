.PHONY: check validate analysis figures report test

validate:
	python scripts/00_validate_inputs.py

check:
	Rscript scripts/00_check_environment.R

analysis:
	Rscript scripts/05_run_selected5_analysis.R

figures:
	Rscript scripts/08_build_main_figures.R

report:
	python scripts/09_build_report.py

test:
	python tests/test_repository.py

