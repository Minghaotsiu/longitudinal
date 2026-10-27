# Reproducibility notes

## Starting points

There are two valid starting points:

1. Raw FASTQ from GSA BioProject `PRJCA078283`.
2. The processed count matrices derived from those reads.

The upstream section starts from FASTQ. The downstream section starts from:

```text
data/doudou_all34/
  ALL_34_raw_counts_matrix.csv
  ALL_34_1M_counts_matrix.csv
  ALL_34_sample_metadata.csv
  ALL_34_QC_metrics.csv
  soybean_gene_symbol_map.csv
```

## Environment

Use a clean R installation and set `R_LIBS_USER` if packages are installed in a
non-default library. The analysis does not depend on a hard-coded user account
or server path.

For Linux upstream processing, activate an environment containing Python,
R, FastQC, cutadapt, BBTools, STAR and samtools. Then export the variables in
`config.example.env`.

## Randomness

R scripts set `set.seed(123)`. Random downsampling in the fixed-depth counting
scripts uses Python random sampling. The analysis should therefore be repeated
with the same software versions when exact identity of downsampled counts is
required.

## Pseudotime

Calendar Day is measured metadata. Slingshot, TSCAN and SCORPIUS infer an
ordering from expression and are evaluated against, but never substituted for,
the recorded Day.

## Figure source data

`results/source_data/` contains the small tabular data used by the main figures.
The figures themselves can be regenerated from the processed matrices and the
analysis tables.

## Known non-reproducible local details

The original working environment contained local and server-specific paths.
These have been replaced by environment variables. No credentials or private
keys are included.

