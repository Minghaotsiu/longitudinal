root <- normalizePath(
  Sys.getenv("PROCESSED_DATA_DIR", "data/doudou_all34"),
  winslash = "/",
  mustWork = TRUE
)
vst <- read.csv(
  file.path(root, "downstream_all34/tables/VST_expression_34.csv"),
  check.names = FALSE
)
manifest <- read.csv(
  file.path(root, "selected5_pseudotime/tables/sample_selection_manifest.csv"),
  check.names = FALSE
)
symbol_map <- read.csv(
  file.path(root, "soybean_gene_symbol_map.csv"),
  check.names = FALSE
)
selected <- manifest$RunSample[manifest$Selection == "Microseed_analysis_n5"]
excluded <- manifest$RunSample[
  manifest$Selection == "Excluded_bulk_simulation_pool"
]
gene_ids <- symbol_map$Gene_ID[symbol_map$Gene_biotype == "protein_coding"]
vst <- vst[vst$Gene_ID %in% gene_ids, , drop = FALSE]
matrix_vst <- as.matrix(vst[, -c(1, 2), drop = FALSE])
x <- rowMeans(matrix_vst[, selected, drop = FALSE])
y <- rowMeans(matrix_vst[, excluded, drop = FALSE])
fit <- lm(y ~ x)
test <- cor.test(x, y, method = "pearson")
out <- data.frame(
  comparison = "Excluded bulk pool vs microseed selected",
  genes = length(x),
  slope = coefficients(fit)[2],
  intercept = coefficients(fit)[1],
  R2 = summary(fit)$r.squared,
  model_p = coefficients(summary(fit))[2, 4],
  pearson_r = unname(test$estimate),
  pearson_p = test$p.value
)
write.csv(
  out,
  file.path(root, "selected5_pseudotime/tables/bulk_vs_micro_expression_profile_regression.csv"),
  row.names = FALSE
)
print(out)
