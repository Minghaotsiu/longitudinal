options(error = function() { traceback(20); quit(status = 1) })
suppressPackageStartupMessages({
  library(TSCAN)
  library(SCORPIUS)
  library(dplyr)
})

base <- file.path(
  normalizePath(
    Sys.getenv("PROCESSED_DATA_DIR", "data/doudou_all34"),
    winslash = "/",
    mustWork = TRUE
  ),
  "selected5_pseudotime",
  "tables"
)
counts_df <- read.csv(file.path(base, "selected_n5_raw_counts.csv"), check.names = FALSE)
vst_df <- read.csv(file.path(base, "selected_n5_VST_expression.csv"), check.names = FALSE)
meta <- read.csv(file.path(base, "selected_microseed_n5.csv"), check.names = FALSE)
counts <- as.matrix(counts_df[, -c(1, 2), drop = FALSE])
rownames(counts) <- counts_df$Gene_ID
vst <- as.matrix(vst_df[, -c(1, 2), drop = FALSE])
rownames(vst) <- vst_df$Gene_ID

keep <- rowSums(counts >= 10) >= 5
counts <- round(counts[keep, , drop = FALSE])
vst <- vst[rownames(vst) %in% rownames(counts), , drop = FALSE]
var <- apply(vst, 1, var)
top <- names(sort(var, decreasing = TRUE))[seq_len(min(1000, length(var)))]

counts_top <- counts[rownames(vst[top[seq_len(min(500, length(top)))], , drop = FALSE]), , drop = FALSE]
proc_tscan <- TSCAN::preprocess(counts_top)
clust_tscan <- TSCAN::exprmclust(
  proc_tscan,
  clusternum = 2:4,
  modelNames = "EII"
)
ord_tscan <- TSCAN::TSCANorder(clust_tscan, orderonly = FALSE, flip = TRUE)
cat("TSCAN_ORDER\n")
print(head(ord_tscan))
print(str(ord_tscan))

expr_scorpius <- t(vst[top, , drop = FALSE])
space_scorpius <- SCORPIUS::reduce_dimensionality(expr_scorpius, "spearman")
traj_scorpius <- SCORPIUS::infer_trajectory(space_scorpius)
cat("SCORPIUS_TIME\n")
print(head(traj_scorpius$time))
print(length(traj_scorpius$time))
