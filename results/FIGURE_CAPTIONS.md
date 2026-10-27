# Scientific Data figure captions

## Figure 1 | Data-acquisition workflow for longitudinal microneedle and destructive bulk transcriptomics

Five individual soybean seeds were sampled repeatedly at Day0, Day1, Day2 and Day3 using active microneedles. A separate cross-sectional comparator comprised 14 independent destructive bulk samples collected across the same time course. Daily bulk sample counts are intentionally omitted from the figure. The lower workflow shows material release, RNA extraction and library preparation for both sampling routes, followed by sequencing, quality control, quantification and construction of the reusable expression and Seed-by-Day metadata tables. The diagram defines the sample identities used throughout the dataset; original sequencing-run labels are not used as biological groupings. Source data are provided as a Source Data file.

## Figure 2 | Technical metrics without significant MN–Bulk differences

Technical characteristics of the 34 sequencing runs for metrics with unadjusted one-way ANOVA P ≥ 0.05 between MN and Bulk: 20 microneedle samples and 14 independent destructive bulk samples. Each point is one run, and horizontal marks indicate group medians. P values are shown in each panel. Metrics with P < 0.05 are intentionally not included in this display; the complete statistical table, including significant metrics, is provided in Source Data. Q20/Q30 are read-level mean-quality summaries, not base-level aggregates. Mapping and too-short fractions are retained to document the technical limits of the dataset. Source data are provided as a Source Data file.

## Figure 3 | Global expression structure of microneedle and bulk samples

**a**, Principal component analysis (PCA), **b**, uniform manifold approximation and projection (UMAP), and **c**, t-distributed stochastic neighbour embedding (t-SNE) of the top 1,000 variable genes. All samples are shown as circles, and colour encodes sampling method only. Batch identifiers are not used as plotting variables or colours. These panels describe overall expression structure and do not constitute a significance test. Source data are provided as a Source Data file.

## Figure 4 | Concordance between microneedle and destructive bulk expression profiles

**a**, Correlation between microneedle and bulk group-mean VST expression across 2,372 filtered genes. The line and shaded region show a linear fit and its 95% confidence interval. **b**, Mean within-method Pearson correlation calculated separately for each run; one value per run is used to avoid treating pairwise correlations as independent observations. The P value is from one-way ANOVA comparing MN and Bulk. **c**, Row-scaled expression of the 80 most variable genes, with samples ordered by sampling method only. These analyses describe concordance and heterogeneity; they do not establish equivalence. Source data are provided as a Source Data file.

## Figure 5 | Data-derived ordering compared with recorded sampling Day

**a**, Microneedle samples in UMAP space, coloured by recorded Day. Arrows connect Day centroids ordered by median Slingshot pseudotime and therefore summarize a data-derived ordering rather than a fitted single-cell trajectory. **b**, Within-seed scaled pseudotime for Slingshot, TSCAN and SCORPIUS across Day0–Day3; lines connect samples assigned to the same seed. **c**, Spearman correlation between each method and recorded Day. **d**, Pairwise Spearman agreement among the three methods. Pseudotime is inferred from expression and is not equivalent to experimental time. These panels assess method behaviour and do not identify a unique true trajectory. Source data are provided as a Source Data file.

## Figure 6 | Developmental and reserve-mobilization candidate genes across the germination time course

Expression of ICL1, NAC2, AS1 and TOP6B is shown for the microneedle time series. Points are individual samples (n = 5 per Day). One-way ANOVA P values are shown in each panel, and Tukey comparisons are displayed for the selected Day pairs. The figure demonstrates reuse of the processed expression matrix and is not intended to establish causal regulation. Source data are provided as a Source Data file.
