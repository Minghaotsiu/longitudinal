#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "data.table",
  "dplyr",
  "ggh4x",
  "gghalves",
  "ggalluvial",
  "ggbeeswarm",
  "ggplot2",
  "ggrepel",
  "httr",
  "jsonlite",
  "MASS",
  "openxlsx",
  "patchwork",
  "pheatmap",
  "png",
  "ragg",
  "Rtsne",
  "readr",
  "scales",
  "svglite",
  "tibble",
  "tidyr",
  "uwot"
)

bioc_packages <- c(
  "Biobase",
  "DESeq2",
  "DynDoc",
  "Mfuzz",
  "SCORPIUS",
  "SingleCellExperiment",
  "slingshot",
  "tradeSeq",
  "TSCAN"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran)) {
  install.packages(missing_cran)
}

missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

cat("R dependencies checked.\n")

