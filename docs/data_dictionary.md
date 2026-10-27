# Data dictionary

## Sample identifiers

- `CorrectedSampleID`: analysis label such as `MN_Seed2_Day0`.
- `SamplingMethod`: `MN` or `Bulk`.
- `SeedID`: `Seed1` to `Seed5` for longitudinal microneedle samples.
- `Day`: `Day0`, `Day1`, `Day2` or `Day3`, counted from the start of
  imbibition.
- `RunAlias`: unique sequencing-run identifier.

## Matrices

- Raw counts: integer gene-level counts.
- 1M counts: counts recalculated after fixed-depth downsampling.
- Normalized counts: DESeq2 size-factor-normalized counts.
- VST: variance-stabilizing transformation from DESeq2.
- FPKM: fragments per kilobase of transcript per million mapped reads.

## Analysis fields

- `padj`: Benjamini-Hochberg-adjusted P value.
- `log2FoldChange`: DESeq2 effect size.
- `cluster`: Mfuzz soft-cluster assignment by maximum membership.
- `Slingshot`, `TSCAN`, `SCORPIUS`: data-derived pseudotime values.
- `SpearmanRho`: rank correlation used to compare orderings.

