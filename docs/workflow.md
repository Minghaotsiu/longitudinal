# Workflow

## 1. Delivered FASTQ

Each delivery is organized into a batch directory. R1 and R2 files are paired
and linked under a unique sample name; original delivery files are not modified.

## 2. Read quality control

Raw reads are assessed with FastQC. cutadapt removes sequencing adapters and
low-quality tails. BBTools `repair.sh` restores synchronized paired reads.

## 3. Alignment and quantification

Clean paired reads are aligned with STAR. Gene-level and fixed-depth counts are
generated against the *Glycine max* v2.1 annotation.

## 4. Batch assembly

The two 17-run deliveries are combined into a 34-run matrix. Sample identity and
Day remain in the metadata table. The code also supports one-batch analysis.

## 5. Downstream analysis

The analysis includes:

- expression filtering and VST transformation;
- PCA, UMAP and t-SNE;
- six pairwise differential-expression comparisons;
- candidate-gene ANOVA and Tukey comparisons;
- Mfuzz temporal clustering;
- g:Profiler enrichment;
- Slingshot, TSCAN and SCORPIUS pseudotime;
- tradeSeq association testing.

## 6. Figures and archive

The figure builder exports SVG, PDF, TIFF and PNG files. Source data are saved as
CSV. The report builder embeds the figures and major tables in a portable HTML
file.

