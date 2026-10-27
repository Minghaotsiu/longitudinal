#!/usr/bin/env Rscript

packages <- c(
  "data.table",
  "DESeq2",
  "dplyr",
  "ggh4x",
  "gghalves",
  "ggalluvial",
  "ggbeeswarm",
  "ggplot2",
  "ggrepel",
  "jsonlite",
  "MASS",
  "Mfuzz",
  "patchwork",
  "pheatmap",
  "png",
  "ragg",
  "readr",
  "Rtsne",
  "scales",
  "SCORPIUS",
  "SingleCellExperiment",
  "slingshot",
  "svglite",
  "tibble",
  "tidyr",
  "tradeSeq",
  "TSCAN",
  "uwot"
)

available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
print(data.frame(package = packages, available = available), row.names = FALSE)

if (any(!available)) {
  stop(
    "Missing R packages: ",
    paste(packages[!available], collapse = ", "),
    call. = FALSE
  )
}

cat("\nR version:", R.version.string, "\n")
cat("Environment check passed.\n")

