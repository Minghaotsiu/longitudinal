# Contributing

This repository accompanies a defined dataset and manuscript. Reproducibility
fixes are welcome through issues or pull requests.

Before submitting a change:

```bash
python tests/test_repository.py
Rscript scripts/00_check_environment.R
```

Do not commit raw sequencing data, credentials, private server paths or large
intermediate files. Any change to sample identifiers must be checked against
`metadata/scientific_data_sample_crosswalk.csv`.

