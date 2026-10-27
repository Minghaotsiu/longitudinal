# Longitudinal soybean seed-germination transcriptomics

Analysis code and reusable metadata for the study:

**Longitudinal soybean seed-germination transcriptomes acquired by active microneedle sampling**

The dataset combines 20 active-microneedle samples collected from five
individual soybean seeds at Day0-Day3 with 14 independent destructive bulk
samples. Raw reads are deposited under Genome Sequence Archive BioProject
[`PRJCA078283`](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA078283).

This repository contains the code used for read processing, expression
quantification, quality control, differential-expression analysis, temporal
clustering, functional enrichment, trajectory inference and manuscript figures.
It does not redistribute raw FASTQ files or large intermediate BAM/SAM files.

## What is included

- `scripts/upstream/`: FASTQ preparation, trimming, STAR alignment, gene
  counting, batch merging and QC collection.
- `scripts/03_run_batch_analysis.R`: downstream analysis of one 17-run batch.
- `scripts/04_run_all34_analysis.R`: analysis of the complete 34-run matrix.
- `scripts/05_run_selected5_analysis.R`: final selected-5 time-course analysis,
  including DESeq2, candidate genes, Mfuzz, enrichment and pseudotime.
- `scripts/06_run_pseudotime_analysis.R`: standalone Slingshot, TSCAN and
  SCORPIUS trajectory analysis.
- `scripts/08_build_main_figures.R`: generation and audit of the six main
  figures.
- `scripts/09_build_report.py`: integrated HTML report.
- `scripts/10_package_upload_bundle.py`: processed-data and figure archive
  bundle for a public repository.
- `metadata/`: sample crosswalk, analysis metadata and complete QC table.
- `results/source_data/`: numerical source data for the main figures.
- `env/`: recorded R session, package requirements and installation helper.
- `tests/`: lightweight repository and metadata checks.

## What is not included

- Raw FASTQ: obtain from GSA BioProject `PRJCA078283`.
- Processed expression matrices: obtain from the linked GSA/OMIX record after
  its accession is assigned.
- The soybean GTF, STAR index, BED files and mitochondrial/ribosomal gene lists:
  these are external reference resources and must be identified through
  environment variables.
- Credentials, server connection files, temporary logs and manuscript DOCX
  files.

## Software

The final R analysis was run with R 4.3.1 and DESeq2 1.40.2. The complete
recorded package set is in `env/R-sessionInfo.txt`.

Command-line tools used by the upstream workflow:

- FastQC
- cutadapt
- BBTools `repair.sh`
- STAR 2.7.11b
- samtools
- Python 3 with pandas, numpy, openpyxl, pysam and tqdm

R package installation:

```bash
Rscript env/install_R_packages.R
```

Python package installation:

```bash
python -m pip install -r env/requirements.txt
```

See `docs/software.md` for the exact role of each tool.

## Quick start

Clone the repository and enter it:

```bash
git clone <YOUR-REPOSITORY-URL>
cd soybean-microneedle-transcriptome
```

Create a local environment file and edit every path:

```bash
cp config.example.env config.env
set -a
source config.env
set +a
```

On Windows PowerShell, copy `config.example.ps1` to `config.ps1` and dot-source
that file instead.

Check the R environment:

```bash
Rscript scripts/00_check_environment.R
```

Place the processed matrices in the directory named by
`PROCESSED_DATA_DIR`. Required files include:

- `ALL_34_raw_counts_matrix.csv`
- `ALL_34_1M_counts_matrix.csv`
- `ALL_34_sample_metadata.csv`
- `ALL_34_QC_metrics.csv`
- `soybean_gene_symbol_map.csv`

Validate the inputs:

```bash
python scripts/00_validate_inputs.py
```

Run the final selected-5 analysis:

```bash
Rscript scripts/05_run_selected5_analysis.R
```

Build the main figures after the selected-5 analysis:

```bash
Rscript scripts/08_build_main_figures.R
```

Build the HTML report:

```bash
python scripts/09_build_report.py
```

## Upstream processing

The upstream workflow expects Linux and the external reference paths listed in
`config.example.env`. A single delivered batch can be prepared and processed
with:

```bash
bash scripts/upstream/run_delivery.sh \
  /path/to/raw_delivery \
  work/upstream/batch1 \
  17
```

For the two-delivery design used in this project:

```bash
export BATCH1_RAW_DIR=/path/to/first_delivery
export BATCH2_RAW_DIR=/path/to/second_delivery
export UPSTREAM_DIR=work/upstream
export PROCESSED_DATA_DIR=data/doudou_all34

bash scripts/upstream/run_two_matched_batches.sh
```

The upstream code writes only derived files into the batch directory. The
preparation step uses hard links when possible so original FASTQ files are not
modified.

## Important interpretation boundaries

- Calendar Day is experimental metadata.
- Slingshot, TSCAN and SCORPIUS produce data-derived orderings and are not
  measurements of chronological time.
- Pseudotime is not used as proof of sample identity.
- The 14 destructive samples provide a technical comparator; they are not
  additional longitudinal observations of the five microneedle-sampled seeds.
- Differential-expression and enrichment results are descriptive analyses of
  reuse value and do not establish causal regulation.

## Data availability

Raw RNA-seq reads:

`GSA BioProject PRJCA078283`

Processed matrices, source data and repository metadata will be linked from the
manuscript and the BioProject record. The exact OMIX accession is not invented
in this repository.

## Citation

Use `CITATION.cff`. The repository is versioned as `v1.0.0` for the initial
Scientific Data submission.

## Licence

Code is released under the MIT License. Data records are distributed separately
under the licence specified in their repository record.

