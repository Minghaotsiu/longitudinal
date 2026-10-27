#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

old_dir <- Sys.getenv("BATCH1_UPSTREAM_DIR")
new_dir <- Sys.getenv("BATCH2_UPSTREAM_DIR")
out_dir <- Sys.getenv("PROCESSED_DATA_DIR", "data/doudou_all34")
if (!nzchar(old_dir) || !nzchar(new_dir)) {
  stop("Set BATCH1_UPSTREAM_DIR and BATCH2_UPSTREAM_DIR.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

old_counts <- read.csv(
  file.path(old_dir, "ALL_17_raw_counts_matrix.csv"),
  check.names = FALSE
)
new_counts <- read.csv(
  file.path(new_dir, "ALL_17_raw_counts_matrix.csv"),
  check.names = FALSE
)
old_1m <- read.csv(
  file.path(old_dir, "ALL_17_1M_counts_matrix.csv"),
  check.names = FALSE
)
new_1m <- read.csv(
  file.path(new_dir, "ALL_17_1M_counts_matrix.csv"),
  check.names = FALSE
)

stopifnot(identical(old_counts$Gene_ID, new_counts$Gene_ID))
stopifnot(identical(old_counts$Gene_ID, old_1m$Gene_ID))
stopifnot(identical(old_counts$Gene_ID, new_1m$Gene_ID))

old_samples <- colnames(old_counts)[-1]
new_samples <- colnames(new_counts)[-1]
stopifnot(setequal(old_samples, new_samples))
new_counts <- new_counts[, c("Gene_ID", old_samples), drop = FALSE]
new_1m <- new_1m[, c("Gene_ID", old_samples), drop = FALSE]

colnames(old_counts)[-1] <- paste0(old_samples, "_B1")
colnames(new_counts)[-1] <- paste0(new_samples, "_B2")
colnames(old_1m)[-1] <- paste0(old_samples, "_B1")
colnames(new_1m)[-1] <- paste0(new_samples, "_B2")

combined <- cbind(old_counts, new_counts[, -1, drop = FALSE])
write.csv(
  combined,
  file.path(out_dir, "ALL_34_raw_counts_matrix.csv"),
  row.names = FALSE
)
combined_1m <- cbind(old_1m, new_1m[, -1, drop = FALSE])
write.csv(
  combined_1m,
  file.path(out_dir, "ALL_34_1M_counts_matrix.csv"),
  row.names = FALSE
)

metadata <- data.frame(
  RunSample = c(paste0(old_samples, "_B1"), paste0(old_samples, "_B2")),
  Sample = c(old_samples, new_samples),
  Batch = factor(rep(c("B1", "B2"), each = length(old_samples)), levels = c("B1", "B2")),
  Timepoint = factor(
    sub(".*-d([0-3])-.*", "d\\1", c(old_samples, new_samples)),
    levels = c("d0", "d1", "d2", "d3")
  ),
  stringsAsFactors = FALSE
)
write.csv(
  metadata,
  file.path(out_dir, "ALL_34_sample_metadata.csv"),
  row.names = FALSE
)

summed <- data.frame(Gene_ID = old_counts$Gene_ID, check.names = FALSE)
for (sample_id in old_samples) {
  summed[[sample_id]] <-
    old_counts[[paste0(sample_id, "_B1")]] +
    new_counts[[paste0(sample_id, "_B2")]]
}
write.csv(
  summed,
  file.path(out_dir, "ALL_17_summed_counts_matrix.csv"),
  row.names = FALSE
)

symbol_file <- file.path(out_dir, "soybean_gene_symbol_map_gtf.csv")
if (!file.exists(symbol_file)) {
  stop("Run parse_gtf_symbols.py before combining batches.")
}
symbol_map <- read.csv(symbol_file, check.names = FALSE)
symbol_map <- symbol_map[
  match(old_counts$Gene_ID, symbol_map$Gene_ID),
  c("Gene_ID", "Symbol", "Gene_has_GTF_name", "Gene_biotype"),
  drop = FALSE
]
if (anyNA(symbol_map$Gene_ID)) {
  stop("The GTF-derived symbol map does not cover every soybean gene.")
}
write.csv(
  symbol_map,
  file.path(out_dir, "soybean_gene_symbol_map.csv"),
  row.names = FALSE
)

cat("combined_dimensions=", paste(dim(combined), collapse = "x"), "\n", sep = "")
cat("combined_1M_dimensions=", paste(dim(combined_1m), collapse = "x"), "\n", sep = "")
cat("summed_dimensions=", paste(dim(summed), collapse = "x"), "\n", sep = "")
cat("samples=", length(old_samples), " base; ", nrow(metadata), " run-level\n", sep = "")
