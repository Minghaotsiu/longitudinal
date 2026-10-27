#!/usr/bin/env Rscript

# Soybean germination RNA-seq downstream analysis.
# Plot structure and colors are adapted from the user's fixed F4 scripts.

options(stringsAsFactors = FALSE)
set.seed(123)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(gghalves)
  library(ggh4x)
  library(patchwork)
  library(pheatmap)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(Mfuzz)
  library(Biobase)
  library(uwot)
  library(Rtsne)
  library(MASS)
  library(httr)
  library(jsonlite)
  library(ggalluvial)
})

root_dir <- normalizePath(
  Sys.getenv("PROCESSED_DATA_DIR", "data/doudou_all34"),
  winslash = "/",
  mustWork = TRUE
)
output_dir <- normalizePath(
  Sys.getenv("ALL34_OUTPUT_DIR", "outputs/doudou_all34_downstream"),
  winslash = "/",
  mustWork = FALSE
)
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

timeline_levels <- c("d0", "d1", "d2", "d3")
day_levels <- c("Day0", "Day1", "Day2", "Day3")
comparisons <- c(
  Day1_vs_Day0 = "Day1",
  Day2_vs_Day0 = "Day2",
  Day3_vs_Day0 = "Day3"
)

# Colors reused from F4UMAP分群图.R and 7.mfuzz/mfuzz.R.
time_colors <- c(
  Day0 = "#7895C1",
  Day1 = "#7BC0CD",
  Day2 = "#ECB66C",
  Day3 = "#ED8D5A"
)
traditional_colors <- c(
  "#BFDFD2", "#51999F", "#4198AC", "#7BC0CD",
  "#DBCB92", "#ECB66C", "#EA9E58", "#ED8D5A"
)

theme_f4 <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = base_size),
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold")
    )
}

save_plot <- function(plot, file_name, width = 8, height = 6, dpi = 800) {
  ggsave(
    file.path(figure_dir, file_name),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

symbol_from_id <- function(gene_id) {
  ifelse(
    grepl("^GLYMA_", gene_id),
    sub("^GLYMA_", "Glyma.", gene_id),
    gene_id
  )
}

message("[1/10] Reading 34-run count matrix")
count_file <- file.path(root_dir, "ALL_34_raw_counts_matrix.csv")
metadata_file <- file.path(root_dir, "ALL_34_sample_metadata.csv")
symbol_file <- file.path(root_dir, "soybean_gene_symbol_map.csv")
qc_file <- file.path(root_dir, "ALL_34_QC_metrics.csv")
summed_file <- file.path(root_dir, "ALL_17_summed_counts_matrix.csv")
if (!file.exists(summed_file)) {
  legacy_summed_file <- file.path(
    root_dir,
    "ALL_17_B1_B2_summed_counts_matrix.csv"
  )
  if (file.exists(legacy_summed_file)) {
    summed_file <- legacy_summed_file
  }
}
stopifnot(
  file.exists(count_file),
  file.exists(metadata_file),
  file.exists(symbol_file),
  file.exists(qc_file)
)

raw_counts <- read.csv(count_file, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)
symbol_map <- read.csv(symbol_file, check.names = FALSE)

if (nrow(metadata) != 34 || ncol(raw_counts) != 35) {
  stop("Expected 34 run-level samples and one Gene_ID column.")
}
if (!identical(raw_counts$Gene_ID, symbol_map$Gene_ID)) {
  stop("Gene_ID order differs between the count matrix and symbol map.")
}

metadata <- metadata[
  match(colnames(raw_counts)[-1], metadata$RunSample),
  ,
  drop = FALSE
]
if (anyNA(metadata$RunSample)) {
  stop("Some count-matrix columns are absent from the run metadata.")
}
rownames(metadata) <- metadata$RunSample
metadata$Batch <- factor(metadata$Batch, levels = c("B1", "B2"))
metadata$Timepoint <- factor(metadata$Timepoint, levels = timeline_levels)
metadata$Day <- factor(
  sub("d", "Day", as.character(metadata$Timepoint)),
  levels = day_levels
)
metadata$DisplaySample <- sub(
  "^XM-QMH-d([0-3])-",
  "Day\\1-",
  metadata$RunSample
)
write.csv(
  metadata,
  file.path(table_dir, "sample_metadata_34.csv"),
  row.names = FALSE
)

count_mat <- as.matrix(raw_counts[, -1, drop = FALSE])
suppressWarnings(storage.mode(count_mat) <- "double")
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)
rownames(count_mat) <- raw_counts$Gene_ID
symbol_lookup <- setNames(symbol_map$Symbol, symbol_map$Gene_ID)
gtf_name_lookup <- setNames(symbol_map$Gene_has_GTF_name, symbol_map$Gene_ID)

message("[2/10] Filtering and DESeq2 with batch correction")
# At least 10 counts in six samples, the smallest run-level timepoint group.
keep <- rowSums(count_mat >= 10) >= min(table(metadata$Timepoint))
filtered_counts <- count_mat[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = metadata,
  design = ~ Batch + Day
)
dds <- DESeq(dds, parallel = FALSE)
normalized_counts <- counts(dds, normalized = TRUE)
vst_mat <- assay(vst(dds, blind = FALSE))

write.csv(
  data.frame(
    Gene_ID = rownames(normalized_counts),
    Symbol = unname(symbol_lookup[rownames(normalized_counts)]),
    normalized_counts,
    check.names = FALSE
  ),
  file.path(table_dir, "DESeq2_normalized_counts_34.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    Gene_ID = rownames(vst_mat),
    Symbol = unname(symbol_lookup[rownames(vst_mat)]),
    vst_mat,
    check.names = FALSE
  ),
  file.path(table_dir, "VST_expression_34.csv"),
  row.names = FALSE
)

message("[3/10] Six-panel QC in F4 half-violin style")
qc_source <- read.csv(qc_file, check.names = FALSE)
qc_df <- qc_source[match(colnames(count_mat), qc_source$RunSample), , drop = FALSE]
qc_df$Day <- factor(qc_df$Day, levels = day_levels)
qc_df$InputReadsM <- qc_df$InputReads / 1e6
qc_df$CleanReadsM <- qc_df$CleanReads / 1e6
qc_df$TotalCountsM <- qc_df$TotalCounts / 1e6
qc_df$TotalCounts1MM <- qc_df$TotalCounts1M / 1e6
qc_df$AllRuns <- "All runs"
write.csv(
  qc_df,
  file.path(table_dir, "QC_summary_34_six_metrics.csv"),
  row.names = FALSE
)

half_panel <- function(data, y_col, y_label, fill_value, x_col = "AllRuns") {
  ggplot(
    data,
    aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[x_col]])
  ) +
    geom_half_violin(side = "r", color = NA, alpha = 0.48) +
    geom_half_boxplot(
      side = "r",
      errorbar.draw = FALSE,
      width = 0.13,
      linewidth = 0.65,
      outlier.shape = NA
    ) +
    geom_half_point_panel(
      side = "l",
      shape = 22,
      size = 1.9,
      color = "white",
      alpha = 0.92
    ) +
    scale_fill_manual(values = setNames(rep(fill_value, length(unique(data[[x_col]]))),
                                        unique(data[[x_col]])), guide = "none") +
    labs(x = NULL, y = y_label) +
    theme_f4(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank()
    )
}

reads_long <- qc_df %>%
  dplyr::select(RunSample, InputReadsM, CleanReadsM) %>%
  pivot_longer(-RunSample, names_to = "ReadStage", values_to = "ReadsM") %>%
  mutate(
    ReadStage = factor(
      ReadStage,
      levels = c("InputReadsM", "CleanReadsM"),
      labels = c("Input reads", "Clean reads")
    )
  )

p_reads <- ggplot(
  reads_long,
  aes(x = ReadStage, y = ReadsM, fill = ReadStage)
) +
  geom_half_violin(side = "r", color = NA, alpha = 0.48) +
  geom_half_boxplot(
    side = "r",
    errorbar.draw = FALSE,
    width = 0.13,
    linewidth = 0.65,
    outlier.shape = NA
  ) +
  geom_half_point_panel(
    side = "l",
    shape = 22,
    size = 1.9,
    color = "white",
    alpha = 0.92
  ) +
  scale_fill_manual(
    values = c(`Input reads` = "#992224", `Clean reads` = "#EF8B67"),
    guide = "none"
  ) +
  labs(x = NULL, y = "Reads (million)") +
  theme_f4(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

plot_list_qc <- list(
  p_reads,
  half_panel(qc_df, "Q30ReadPct", "Q30 reads (%)", "#F0C284"),
  half_panel(qc_df, "TotalCountsM", "Assigned counts (million)", "#D6EFF4"),
  half_panel(qc_df, "nGeneRaw", "Detected genes", "#7895C1"),
  half_panel(qc_df, "TotalCounts1MM", "Counts at 1M depth (million)", "#A8CBDF"),
  half_panel(qc_df, "nGene1M", "Detected genes at 1M", "#8074C8")
)

p_qc <- wrap_plots(plot_list_qc, ncol = 3) +
  plot_annotation(
    title = "Sequencing quality and transcript capture across 34 runs",
    subtitle = "Six-panel half-violin layout follows the F4 QC script; no grouping by Day",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30")
    )
  )
save_plot(p_qc, "01_QC_six_panel.png", width = 10, height = 7)

message("[4/10] PCA, UMAP, t-SNE and LDA in the fixed F4 style")
gene_variance <- apply(vst_mat, 1, var)
top_variable_n <- min(1000, sum(is.finite(gene_variance) & gene_variance > 0))
top_variable <- names(sort(gene_variance, decreasing = TRUE))[seq_len(top_variable_n)]
pca_fit <- prcomp(t(vst_mat[top_variable, , drop = FALSE]), scale. = FALSE)
percent_var <- 100 * pca_fit$sdev^2 / sum(pca_fit$sdev^2)
pca_df <- data.frame(
  RunSample = rownames(pca_fit$x),
  DisplaySample = metadata[rownames(pca_fit$x), "DisplaySample"],
  PC1 = pca_fit$x[, 1],
  PC2 = pca_fit$x[, 2],
  Day = metadata[rownames(pca_fit$x), "Day"],
  Batch = metadata[rownames(pca_fit$x), "Batch"]
)

axis_truncated <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

dimension_plot <- function(data, x_col, y_col, title, x_label, y_label) {
  ggplot(data, aes(
    x = .data[[x_col]],
    y = .data[[y_col]],
    color = Day,
    shape = Batch
  )) +
    geom_point(size = 3.5, stroke = 0.8) +
    geom_text_repel(
      aes(label = DisplaySample),
      color = "black",
      size = 2.55,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.35,
      segment.size = 0.3,
      show.legend = FALSE
    ) +
    theme_f4(base_size = 16) +
    theme(
      aspect.ratio = 1,
      axis.line = element_line(
        arrow = arrow(type = "closed", length = unit(0.18, "inches")),
        linewidth = 0.9
      ),
      axis.title = element_text(hjust = 0.5, face = "italic", size = 16),
      axis.text = element_blank(),
      legend.text = element_text(size = 13),
      plot.title = element_text(size = 17, hjust = 0.5)
    ) +
    scale_color_manual(values = time_colors) +
    scale_shape_manual(values = c(B1 = 16, B2 = 17)) +
    guides(x = axis_truncated, y = axis_truncated) +
    labs(
      x = x_label,
      y = y_label,
      color = "Day",
      shape = "Batch",
      title = title
    )
}

p_pca <- dimension_plot(
  pca_df,
  "PC1",
  "PC2",
  "PCA",
  sprintf("PC1 (%.1f%%)", percent_var[1]),
  sprintf("PC2 (%.1f%%)", percent_var[2])
)
save_plot(p_pca, "02a_PCA_original_style.png", width = 9, height = 8)

set.seed(123)
umap_fit <- uwot::umap(
  pca_fit$x[, seq_len(min(20, ncol(pca_fit$x))), drop = FALSE],
  n_neighbors = 8,
  min_dist = 0.25,
  metric = "cosine",
  n_components = 2,
  verbose = FALSE
)
colnames(umap_fit) <- c("UMAP1", "UMAP2")
umap_df <- data.frame(
  RunSample = rownames(pca_fit$x),
  DisplaySample = metadata[rownames(pca_fit$x), "DisplaySample"],
  umap_fit,
  Day = metadata[rownames(pca_fit$x), "Day"],
  Batch = metadata[rownames(pca_fit$x), "Batch"]
)
p_umap <- dimension_plot(
  umap_df,
  "UMAP1",
  "UMAP2",
  "UMAP",
  "UMAP1",
  "UMAP2"
)
save_plot(p_umap, "02b_UMAP_original_style.png", width = 9, height = 8)

set.seed(123)
tsne_fit <- Rtsne::Rtsne(
  t(vst_mat[top_variable, , drop = FALSE]),
  dims = 2,
  perplexity = 8,
  theta = 0.0,
  check_duplicates = FALSE,
  pca = TRUE,
  max_iter = 1500
)$Y
colnames(tsne_fit) <- c("tSNE1", "tSNE2")
tsne_df <- data.frame(
  RunSample = colnames(vst_mat),
  DisplaySample = metadata[colnames(vst_mat), "DisplaySample"],
  tsne_fit,
  Day = metadata[colnames(vst_mat), "Day"],
  Batch = metadata[colnames(vst_mat), "Batch"]
)
p_tsne <- dimension_plot(
  tsne_df,
  "tSNE1",
  "tSNE2",
  "t-SNE",
  "t-SNE1",
  "t-SNE2"
)
save_plot(p_tsne, "02c_tSNE_original_style.png", width = 9, height = 8)

lda_fit <- MASS::lda(
  x = pca_fit$x[, seq_len(min(20, ncol(pca_fit$x))), drop = FALSE],
  grouping = metadata[rownames(pca_fit$x), "Day"]
)
lda_coords <- predict(lda_fit)$x
if (ncol(lda_coords) < 2) {
  lda_coords <- cbind(lda_coords, LD2 = 0)
}
colnames(lda_coords)[1:2] <- c("LD1", "LD2")
lda_df <- data.frame(
  RunSample = rownames(pca_fit$x),
  DisplaySample = metadata[rownames(pca_fit$x), "DisplaySample"],
  lda_coords[, 1:2, drop = FALSE],
  Day = metadata[rownames(pca_fit$x), "Day"],
  Batch = metadata[rownames(pca_fit$x), "Batch"]
)
p_lda <- dimension_plot(
  lda_df,
  "LD1",
  "LD2",
  "LDA",
  "LD1",
  "LD2"
)
save_plot(p_lda, "02d_LDA_original_style.png", width = 9, height = 8)

p_dimensions <- (p_pca + p_umap) / (p_tsne + p_lda)
save_plot(
  p_dimensions,
  "02_dimensional_reduction_4panel.png",
  width = 15,
  height = 13
)

write.csv(pca_df, file.path(table_dir, "PCA_coordinates_34.csv"), row.names = FALSE)
write.csv(umap_df, file.path(table_dir, "UMAP_coordinates_34.csv"), row.names = FALSE)
write.csv(tsne_df, file.path(table_dir, "tSNE_coordinates_34.csv"), row.names = FALSE)
write.csv(lda_df, file.path(table_dir, "LDA_coordinates_34.csv"), row.names = FALSE)

message("[5/10] Differential expression with batch correction")
de_results <- list()
de_summary <- list()
significant_union <- character(0)

for (contrast_name in names(comparisons)) {
  target <- comparisons[[contrast_name]]
  res <- results(
    dds,
    contrast = c("Day", target, "Day0"),
    alpha = 0.05,
    independentFiltering = TRUE
  )
  res_df <- as.data.frame(res)
  res_df$Gene_ID <- rownames(res_df)
  res_df$Symbol <- unname(symbol_lookup[res_df$Gene_ID])
  res_df$direction <- "ns"
  res_df$direction[
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      res_df$log2FoldChange >= 0.3
  ] <- "up"
  res_df$direction[
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      res_df$log2FoldChange <= -0.3
  ] <- "down"
  res_df$comparison <- contrast_name
  res_df <- res_df[, c(
    "Gene_ID", "Symbol", "baseMean", "log2FoldChange", "lfcSE", "stat",
    "pvalue", "padj", "direction", "comparison"
  )]
  write.csv(
    res_df,
    file.path(table_dir, paste0("DE_", contrast_name, "_all.csv")),
    row.names = FALSE
  )
  sig <- res_df[res_df$direction != "ns", , drop = FALSE]
  write.csv(
    sig,
    file.path(table_dir, paste0("DE_", contrast_name, "_significant.csv")),
    row.names = FALSE
  )
  de_results[[contrast_name]] <- res_df
  significant_union <- union(significant_union, sig$Gene_ID)
  de_summary[[contrast_name]] <- data.frame(
    comparison = contrast_name,
    up = sum(res_df$direction == "up"),
    down = sum(res_df$direction == "down"),
    tested = sum(!is.na(res_df$padj)),
    stringsAsFactors = FALSE
  )

  plot_df <- res_df %>%
    filter(!is.na(pvalue), is.finite(pvalue)) %>%
    mutate(neglog10p = -log10(pmax(pvalue, 1e-50)))
  label_df <- plot_df %>%
    filter(direction != "ns") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 10)

  p_volcano <- ggplot(
    plot_df,
    aes(x = log2FoldChange, y = neglog10p, color = log2FoldChange, size = neglog10p)
  ) +
    geom_point(alpha = 0.8) +
    geom_text_repel(
      data = label_df,
      aes(label = Symbol),
      color = "black",
      size = 4.2,
      box.padding = 0.8,
      max.overlaps = 50,
      min.segment.length = 0
    ) +
    scale_color_gradientn(colours = c("#7BC0CD", "#ED8D5A")) +
    scale_size(range = c(1, 7)) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = 2, linewidth = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = 4, linewidth = 0.7) +
    labs(
      x = expression(log[2]~fold~change),
      y = expression(-log[10](italic(p))),
      color = "log2FC",
      size = "-log10(p)",
      title = gsub("_", " ", contrast_name),
      subtitle = "DESeq2 Wald test, BH padj<0.05, |log2FC|>=0.3; batch-adjusted"
    ) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1),
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13)
    )
  save_plot(
    p_volcano,
    paste0("03_volcano_", contrast_name, ".png"),
    width = 8,
    height = 6.5
  )
}

de_summary_df <- bind_rows(de_summary)
write.csv(
  de_summary_df,
  file.path(table_dir, "DE_summary_34.csv"),
  row.names = FALSE
)

p_de_counts <- de_summary_df %>%
  pivot_longer(c(up, down), names_to = "direction", values_to = "count") %>%
  mutate(
    comparison = factor(comparison, levels = names(comparisons)),
    direction = factor(direction, levels = c("down", "up"))
  ) %>%
  ggplot(aes(comparison, count, fill = direction)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.65,
    color = "black"
  ) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    size = 4
  ) +
  scale_fill_manual(values = c(down = "#7BC0CD", up = "#ED8D5A")) +
  labs(
    x = NULL,
    y = "Number of genes",
    fill = "Direction",
    title = "Differentially expressed genes",
    subtitle = "DESeq2 Wald test, BH padj<0.05, |log2FC|>=0.3; n=6/10/10/8"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    plot.title = element_text(hjust = 0.5, size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30")
  )
save_plot(p_de_counts, "04_DE_counts_original_style.png", width = 8, height = 5)

message("[6/10] Bar plots for ranked genes")
union_table <- bind_rows(de_results)
gene_rank <- union_table %>%
  filter(Gene_ID %in% significant_union) %>%
  group_by(Gene_ID, Symbol) %>%
  summarise(
    min_padj = min(padj, na.rm = TRUE),
    max_abs_lfc = max(abs(log2FoldChange), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(min_padj, desc(max_abs_lfc))
gene_rank$Gene_has_GTF_name <- unname(gtf_name_lookup[gene_rank$Gene_ID])

# Curated, development- or hormone-related families with conserved plant
# evidence. Statistical support is reported separately so nominal candidates
# are not presented as FDR-significant findings.
curated_candidates <- data.frame(
  Gene_ID = c(
    "GLYMA_13G069900",
    "GLYMA_19G006400",
    "GLYMA_06G303200",
    "GLYMA_06G114000",
    "GLYMA_18G061100",
    "GLYMA_13G342900",
    "GLYMA_07G132600",
    "GLYMA_17G027400"
  ),
  Candidate_category = c(
    "GA-responsive cell expansion",
    "Auxin signaling",
    "Glyoxylate cycle / lipid reserve mobilization",
    "NAC developmental transcription factor",
    "AS1-like developmental regulator",
    "TOP6B organ growth",
    "Expansin / cell-wall loosening",
    "ABA/desiccation-associated LEA protein"
  ),
  Evidence_note = c(
    "GASA/GAST1-family member; GA-regulated growth and cell elongation",
    "Auxin response factor family; auxin-dependent growth regulation",
    "Isocitrate lyase; glyoxylate cycle supports lipid-to-carbohydrate conversion during germination",
    "NAC transcription-factor family; developmental and stress-responsive regulator",
    "ASYMMETRIC LEAVES1-like developmental regulator",
    "DNA topoisomerase VI subunit B; organ growth and development",
    "Expansin family; cell-wall loosening during growth",
    "Late-embryogenesis-abundant family; ABA/desiccation response"
  ),
  Evidence_reference = c(
    "Aubert et al., Plant Mol Biol 1998; doi:10.1023/A:1005938624418",
    "Guilfoyle & Hagen, Curr Opin Plant Biol 2007; doi:10.1016/j.pbi.2007.08.014",
    "Eastmond et al., PNAS 2000; doi:10.1073/pnas.97.10.5669",
    "Olsen et al., Trends Plant Sci 2005; doi:10.1016/j.tplants.2004.12.010",
    "Byrne et al., Nature 2000; doi:10.1038/35050091",
    "Hartung et al., Current Biology 2002; doi:10.1016/S0960-9822(02)01218-6",
    "Cosgrove, Nature 2000; doi:10.1038/35030000",
    "Battaglia et al., Plant Physiol 2008; doi:10.1104/pp.108.120725"
  ),
  stringsAsFactors = FALSE
)
curated_candidates$Symbol <- unname(symbol_lookup[curated_candidates$Gene_ID])
all_de <- bind_rows(de_results)
curated_stats <- all_de %>%
  filter(Gene_ID %in% curated_candidates$Gene_ID) %>%
  group_by(Gene_ID) %>%
  summarise(
    min_pvalue = min(pvalue, na.rm = TRUE),
    min_padj = ifelse(
      all(is.na(padj)),
      NA_real_,
      min(padj, na.rm = TRUE)
    ),
    max_abs_lfc = max(abs(log2FoldChange), na.rm = TRUE),
    significant_any = any(direction != "ns"),
    significant_contrasts = paste(
      unique(comparison[direction != "ns"]),
      collapse = "; "
    ),
    .groups = "drop"
  )
candidate_genes <- curated_candidates %>%
  left_join(curated_stats, by = "Gene_ID") %>%
  mutate(
    Statistical_support = ifelse(
      significant_any,
      "FDR-significant in at least one comparison",
      ifelse(min_pvalue < 0.05,
             "Nominal p<0.05 only; not FDR-significant",
             "Not significant")
    )
  )
write.csv(
  candidate_genes,
  file.path(table_dir, "development_ABA_candidate_genes_with_evidence.csv"),
  row.names = FALSE
)
write.csv(
  gene_rank,
  file.path(table_dir, "DE_ranked_gene_union_with_symbols.csv"),
  row.names = FALSE
)

top_genes <- head(candidate_genes$Gene_ID, 8)
bar_long <- data.frame(
  Gene_ID = rep(top_genes, each = ncol(vst_mat)),
  RunSample = rep(colnames(vst_mat), times = length(top_genes)),
  Day = rep(metadata[colnames(vst_mat), "Day"], times = length(top_genes)),
  Expression = as.vector(t(vst_mat[top_genes, , drop = FALSE]))
)
bar_long$Symbol <- unname(symbol_lookup[bar_long$Gene_ID])
bar_summary <- bar_long %>%
  group_by(Gene_ID, Symbol, Day) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression),
    n = n(),
    .groups = "drop"
  )

p_bars <- ggplot(
  bar_summary,
  aes(Day, mean_expression, fill = Day)
) +
  geom_col(width = 0.65, color = "black", linewidth = 0.45) +
  geom_errorbar(
    aes(ymin = mean_expression, ymax = mean_expression + sd_expression),
    width = 0.18,
    linewidth = 0.55
  ) +
  geom_jitter(
    data = bar_long,
    aes(Day, Expression),
    width = 0.13,
    size = 1,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.35,
    inherit.aes = FALSE
  ) +
  facet_wrap(~Symbol, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = time_colors, guide = "none") +
  labs(
    x = NULL,
    y = "VST-normalized expression",
    title = "Development- and hormone-related candidate genes",
    subtitle = "Bars: mean; error bars: +SD; points: individual runs"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30")
  )
save_plot(p_bars, "05_candidate_gene_bars.png", width = 11, height = 8)
write.csv(bar_long, file.path(table_dir, "ranked_gene_sample_expression.csv"), row.names = FALSE)

message("[7/10] Mfuzz clustering in the fixed style")
time_mean_vst <- sapply(
  day_levels,
  function(day_label) rowMeans(vst_mat[, metadata$Day == day_label, drop = FALSE])
)
colnames(time_mean_vst) <- day_levels
time_var <- apply(time_mean_vst, 1, var)
mfuzz_n <- min(2000, sum(is.finite(time_var) & time_var > 0))
mfuzz_genes <- names(sort(time_var, decreasing = TRUE))[seq_len(mfuzz_n)]
mfuzz_gene_mat <- time_mean_vst[mfuzz_genes, , drop = FALSE]

mfuzz_class <- new("ExpressionSet", exprs = mfuzz_gene_mat)
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = "mean")
mfuzz_class <- filter.std(mfuzz_class, min.std = 0.1, visu = FALSE)
mfuzz_class <- standardise(mfuzz_class)
m_value <- mestimate(mfuzz_class)
cluster_num <- 6
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = m_value)

cluster_size <- mfuzz_cluster$size
names(cluster_size) <- seq_along(cluster_size)
write.csv(
  data.frame(cluster = names(cluster_size), genes = as.integer(cluster_size)),
  file.path(table_dir, "Mfuzz_cluster_sizes.csv"),
  row.names = FALSE
)

memship <- mfuzz_cluster$membership
colnames(memship) <- paste0("membership", seq_len(ncol(memship)))
exp <- data.frame(
  exprs(mfuzz_class),
  Identifier = rownames(mfuzz_class),
  clusterindex = mfuzz_cluster$cluster,
  memship,
  check.names = FALSE
)
exp <- exp %>%
  pivot_longer(
    cols = all_of(day_levels),
    names_to = "sample",
    values_to = "expression"
  ) %>%
  mutate(
    sample = factor(sample, levels = day_levels),
    clusterindex = factor(clusterindex, levels = seq_len(cluster_num))
  )
exp$maxMembership <- apply(
  exp[, grep("^membership", colnames(exp)), drop = FALSE],
  1,
  max
)

gene_cluster <- data.frame(
  Gene_ID = names(mfuzz_cluster$cluster),
  Symbol = unname(symbol_lookup[names(mfuzz_cluster$cluster)]),
  cluster = as.integer(mfuzz_cluster$cluster),
  maxMembership = apply(memship, 1, max),
  check.names = FALSE
)
write.csv(
  gene_cluster,
  file.path(table_dir, "Mfuzz_gene_clusters.csv"),
  row.names = FALSE
)

png(
  file.path(figure_dir, "06a_Mfuzz_plot2.png"),
  width = 3000,
  height = 2200,
  res = 300
)
mfuzz.plot2(
  mfuzz_class,
  cl = mfuzz_cluster,
  mfrow = c(2, 3),
  time.labels = day_levels,
  xlab = "",
  x11 = FALSE
)
dev.off()

set.seed(123)
line_ids <- exp %>%
  distinct(Identifier) %>%
  slice_sample(n = min(800, nrow(.))) %>%
  pull(Identifier)
exp_lines <- exp %>% filter(Identifier %in% line_ids)

p_mfuzz <- ggplot(
  exp_lines,
  aes(x = sample, y = expression)
) +
  geom_line(
    aes(group = Identifier, color = maxMembership),
    linewidth = 0.35,
    alpha = 0.18
  ) +
  facet_grid(clusterindex ~ ., switch = "y") +
  scale_color_gradientn(
    "Membership",
    colours = traditional_colors
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Standardized expression",
    title = "Mfuzz temporal expression clusters",
    subtitle = paste0("Fuzzy c-means, c=6, m=", round(m_value, 2),
                      "; top ", mfuzz_n, " dynamic genes; no p-value")
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    plot.title = element_text(hjust = 0.5, size = 17),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30")
  )
save_plot(p_mfuzz, "06b_Mfuzz_trajectory.png", width = 8, height = 10, dpi = 400)

heat_genes <- gene_cluster %>%
  group_by(cluster) %>%
  arrange(desc(maxMembership)) %>%
  slice_head(n = 12) %>%
  ungroup()
heat_matrix <- exprs(mfuzz_class)[heat_genes$Gene_ID, , drop = FALSE]
heat_df <- as.data.frame(heat_matrix) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(all_of(day_levels), names_to = "Day", values_to = "Z") %>%
  left_join(heat_genes[, c("Gene_ID", "Symbol", "cluster")], by = "Gene_ID") %>%
  mutate(
    Day = factor(Day, levels = day_levels),
    cluster = factor(cluster, levels = seq_len(cluster_num))
  )

p_mfuzz_heat <- ggplot(heat_df, aes(Day, Symbol, fill = Z)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradientn(
    colours = traditional_colors,
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Z-score",
    title = "Mfuzz cluster heatmap",
    subtitle = "Top 12 membership genes per cluster"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 5),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16)
  )
save_plot(p_mfuzz_heat, "06c_Mfuzz_heatmap.png", width = 8, height = 11, dpi = 400)

message("[8/10] g:Profiler enrichment for Glycine max")
gprofiler_enrich <- function(genes, contrast_name, direction) {
  if (length(genes) == 0) {
    return(NULL)
  }
  body <- list(
    organism = "gmax",
    query = unname(genes),
    sources = c("GO:BP", "GO:CC", "GO:MF", "KEGG"),
    user_threshold = 0.05,
    domain_scope = "annotated"
  )
  response <- tryCatch(
    httr::POST(
      "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
      body = jsonlite::toJSON(body, auto_unbox = TRUE),
      httr::content_type_json(),
      httr::timeout(120)
    ),
    error = function(error) error
  )
  if (inherits(response, "error")) {
    warning(contrast_name, " ", direction, ": ", response$message)
    return(NULL)
  }
  if (httr::status_code(response) != 200) {
    warning(contrast_name, " ", direction, ": HTTP ", httr::status_code(response))
    return(NULL)
  }
  parsed <- jsonlite::fromJSON(
    httr::content(response, as = "text", encoding = "UTF-8"),
    simplifyDataFrame = TRUE
  )
  result <- parsed$result
  if (is.null(result) || !is.data.frame(result) || nrow(result) == 0) {
    return(NULL)
  }
  result$contrast <- contrast_name
  result$direction <- direction
  result
}

enrichment_list <- list()
for (contrast_name in names(de_results)) {
  for (direction_label in c("up", "down")) {
    genes <- de_results[[contrast_name]]$Gene_ID[
      de_results[[contrast_name]]$direction == direction_label
    ]
    message("  ", contrast_name, " ", direction_label, ": ", length(genes), " genes")
    result <- gprofiler_enrich(genes, contrast_name, direction_label)
    if (!is.null(result)) {
      enrichment_list[[paste(contrast_name, direction_label, sep = "_")]] <- result
    }
  }
}

enrichment_df <- if (length(enrichment_list)) bind_rows(enrichment_list) else data.frame()
if (nrow(enrichment_df) > 0) {
  enrichment_df[] <- lapply(
    enrichment_df,
    function(column) {
      if (is.list(column)) {
        vapply(column, function(value) paste(value, collapse = ";"), character(1))
      } else {
        column
      }
    }
  )
}
write.csv(
  enrichment_df,
  file.path(table_dir, "GO_KEGG_enrichment_gmax.csv"),
  row.names = FALSE
)

message("[9/10] Enrichment plots")
if (nrow(enrichment_df) > 0 && "p_value" %in% colnames(enrichment_df)) {
  enrichment_plot <- enrichment_df %>%
    mutate(
      comparison_direction = paste(contrast, direction, sep = "\n"),
      gene_ratio = intersection_size / query_size,
      minus_log10_p = -log10(pmax(p_value, 1e-50))
    ) %>%
    group_by(comparison_direction) %>%
    arrange(p_value, desc(intersection_size)) %>%
    slice_head(n = 8) %>%
    ungroup()

  p_bubble <- ggplot(
    enrichment_plot,
    aes(gene_ratio, reorder(name, gene_ratio))
  ) +
    geom_point(
      aes(size = intersection_size, color = minus_log10_p)
    ) +
    facet_wrap(~comparison_direction, scales = "free_y", ncol = 2) +
    scale_color_gradientn(colours = c("#7BC0CD", "#ED8D5A")) +
    labs(
      x = "Gene ratio",
      y = NULL,
      size = "Count",
      color = "-log10(p)",
      title = "GO and KEGG enrichment",
      subtitle = "g:Profiler gmax; hypergeometric test; g:SCS correction; p<0.05"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(size = 9, color = "black"),
      plot.title = element_text(hjust = 0.5, size = 17),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30")
    )
  save_plot(p_bubble, "07_GO_KEGG_bubble_original_style.png", width = 11, height = 9)

  # Fixed term-to-gene fan plot in the visual language of the F4 enrichment script.
  term_genes <- enrichment_df %>%
    filter(contrast == "Day2_vs_Day0") %>%
    arrange(p_value) %>%
    slice_head(n = 5)
  if (nrow(term_genes) > 0) {
    links <- term_genes %>%
      mutate(genes = strsplit(intersections, ";", fixed = TRUE)) %>%
      dplyr::select(source, name, direction, genes) %>%
      unnest(genes) %>%
      mutate(
        Symbol = unname(symbol_lookup[genes]),
        Symbol = ifelse(is.na(Symbol) | Symbol == "", genes, Symbol)
      ) %>%
      distinct() %>%
      group_by(name) %>%
      slice_head(n = 5) %>%
      ungroup() %>%
      mutate(
        term_rank = dense_rank(name),
        gene_rank = dense_rank(genes)
      )
    max_rank <- max(links$gene_rank)
    term_nodes <- links %>%
      distinct(name, term_rank) %>%
      transmute(label = name, x = 1, y = max_rank - term_rank + 1)
    gene_nodes <- links %>%
      distinct(genes, gene_rank) %>%
      transmute(label = genes, x = 0, y = max_rank - gene_rank + 1)

    p_fan <- ggplot(links) +
      geom_segment(
        aes(
          x = 1,
          xend = 0,
          y = max_rank - term_rank + 1,
          yend = max_rank - gene_rank + 1
        ),
        color = "#ED8D5A",
        alpha = 0.45,
        linewidth = 0.8
      ) +
      geom_point(
        data = term_nodes,
        aes(x, y),
        shape = 15,
        size = 4,
        color = "#7BC0CD"
      ) +
      geom_point(
        data = gene_nodes,
        aes(x, y),
        shape = 16,
        size = 2.5,
        color = "#ED8D5A"
      ) +
      geom_text(
        data = term_nodes,
        aes(x = 1.01, y = y, label = label),
        hjust = 0,
        size = 3.2,
        color = "black"
      ) +
      geom_text(
        data = gene_nodes,
        aes(x = -0.01, y = y, label = label),
        hjust = 1,
        size = 2.8,
        color = "black"
      ) +
      coord_cartesian(xlim = c(-0.75, 1.75), clip = "off") +
      theme_void() +
      theme(
        plot.margin = margin(18, 18, 18, 18),
        plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30")
      ) +
      labs(
        title = "Top enriched terms and their genes",
        subtitle = "Day2 vs Day0; g:Profiler gmax, g:SCS corrected p<0.05"
      )
    save_plot(p_fan, "08a_GO_KEGG_term_gene_links.png", width = 12, height = 8)

    sankey_long <- data.frame(
      Comparison = "Day2 vs Day0",
      Direction = factor(
        as.character(links[["direction"]]),
        levels = c("up", "down")
      ),
      Term = ifelse(
        nchar(as.character(links[["name"]])) > 48,
        paste0(substr(as.character(links[["name"]]), 1, 45), "..."),
        as.character(links[["name"]])
      ),
      Gene = as.character(links[["Symbol"]]),
      stringsAsFactors = FALSE
    ) %>%
      group_by(Comparison, Direction, Term, Gene) %>%
      summarise(Weight = n(), .groups = "drop")

    p_sankey <- ggplot(
      sankey_long,
      aes(
        axis1 = Comparison,
        axis2 = Term,
        axis3 = Gene,
        y = Weight
      )
    ) +
      geom_alluvium(aes(fill = Direction), width = 0.20, alpha = 0.72) +
      geom_stratum(width = 0.20, fill = "grey96", color = "white") +
      geom_text(
        stat = "stratum",
        aes(label = after_stat(stratum)),
        size = 2.3,
        color = "black"
      ) +
      scale_x_discrete(
        limits = c("Comparison", "Enriched term", "Gene"),
        expand = c(0.04, 0.04)
      ) +
      scale_fill_manual(
        values = c(up = "#ED8D5A", down = "#7BC0CD"),
        drop = FALSE
      ) +
      labs(
        x = NULL,
        y = "Gene-term links",
        fill = "Direction",
        title = "Supplementary enrichment Sankey",
        subtitle = "Top five Day2/Day0 terms; up to five genes per term"
      ) +
      theme_void(base_size = 10) +
      theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30")
      )
    save_plot(p_sankey, "08b_GO_KEGG_sankey.png", width = 13, height = 9)
  }
}

message("[9b/10] Mfuzz cluster functional enrichment")
cluster_enrichment_list <- list()
for (cluster_id in seq_len(cluster_num)) {
  cluster_genes <- gene_cluster$Gene_ID[gene_cluster$cluster == cluster_id]
  message("  Cluster", cluster_id, ": ", length(cluster_genes), " genes")
  cluster_result <- gprofiler_enrich(
    cluster_genes,
    paste0("Cluster", cluster_id),
    "cluster"
  )
  if (!is.null(cluster_result)) {
    cluster_enrichment_list[[paste0("Cluster", cluster_id)]] <- cluster_result
  }
}
cluster_enrichment_df <- if (length(cluster_enrichment_list) > 0) {
  bind_rows(cluster_enrichment_list)
} else {
  data.frame()
}
if (nrow(cluster_enrichment_df) > 0) {
  cluster_enrichment_df[] <- lapply(
    cluster_enrichment_df,
    function(column) {
      if (is.list(column)) {
        vapply(column, function(value) paste(value, collapse = ";"), character(1))
      } else {
        column
      }
    }
  )
}
write.csv(
  cluster_enrichment_df,
  file.path(table_dir, "Mfuzz_cluster_GO_KEGG.csv"),
  row.names = FALSE
)
if (nrow(cluster_enrichment_df) > 0 && "p_value" %in% colnames(cluster_enrichment_df)) {
  cluster_enrichment_plot <- cluster_enrichment_df %>%
    mutate(minus_log10_p = -log10(pmax(p_value, 1e-50))) %>%
    group_by(contrast) %>%
    arrange(p_value) %>%
    slice_head(n = 5) %>%
    ungroup()
  p_cluster_enrichment <- ggplot(
    cluster_enrichment_plot,
    aes(minus_log10_p, reorder(name, minus_log10_p))
  ) +
    geom_point(aes(size = intersection_size, color = source)) +
    facet_wrap(~contrast, scales = "free_y", ncol = 3) +
    scale_color_manual(
      values = c(
        `GO:BP` = "#7895C1",
        `GO:CC` = "#7BC0CD",
        `GO:MF` = "#F0C284",
        KEGG = "#ED8D5A"
      )
    ) +
    labs(
      x = "-log10(p)",
      y = NULL,
      size = "Genes",
      color = "Source",
      title = "Functional enrichment of Mfuzz clusters",
      subtitle = "g:Profiler gmax; g:SCS corrected p<0.05"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 6.8),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30")
    )
  save_plot(p_cluster_enrichment, "09_Mfuzz_cluster_enrichment.png", width = 12, height = 8)
}

message("[10/10] Aggregated 17-sample sensitivity analysis")
summed <- read.csv(
  summed_file,
  check.names = FALSE
)
summed_mat <- as.matrix(summed[, -1, drop = FALSE])
suppressWarnings(storage.mode(summed_mat) <- "double")
summed_mat <- round(summed_mat)
rownames(summed_mat) <- summed$Gene_ID
base_samples <- colnames(summed_mat)
base_meta <- data.frame(
  Sample = base_samples,
  Day = factor(
    paste0("Day", sub(".*-d([0-3])-.*", "\\1", base_samples)),
    levels = day_levels
  ),
  row.names = base_samples
)

summed_filter <- rowSums(summed_mat >= 10) >= min(table(base_meta$Day))
dds_summed <- DESeqDataSetFromMatrix(
  countData = summed_mat[summed_filter, , drop = FALSE],
  colData = base_meta,
  design = ~ Day
)
dds_summed <- DESeq(dds_summed, parallel = FALSE)

sensitivity_summary <- list()
for (contrast_name in names(comparisons)) {
  target <- comparisons[[contrast_name]]
  res <- results(
    dds_summed,
    contrast = c("Day", target, "Day0"),
    alpha = 0.05
  )
  res_df <- as.data.frame(res)
  res_df$Gene_ID <- rownames(res_df)
  res_df$direction <- "ns"
  res_df$direction[
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      res_df$log2FoldChange >= 0.3
  ] <- "up"
  res_df$direction[
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      res_df$log2FoldChange <= -0.3
  ] <- "down"
  write.csv(
    res_df,
    file.path(table_dir, paste0("SENSITIVITY_17_DE_", contrast_name, ".csv")),
    row.names = FALSE
  )

  common <- intersect(
    de_results[[contrast_name]]$Gene_ID,
    res_df$Gene_ID
  )
  run_level <- de_results[[contrast_name]]$log2FoldChange[
    match(common, de_results[[contrast_name]]$Gene_ID)
  ]
  sensitivity_summary[[contrast_name]] <- data.frame(
    comparison = contrast_name,
    run34_up = sum(de_results[[contrast_name]]$direction == "up"),
    run34_down = sum(de_results[[contrast_name]]$direction == "down"),
    summed17_up = sum(res_df$direction == "up"),
    summed17_down = sum(res_df$direction == "down"),
    common_genes = length(common),
    spearman_log2FC = cor(run_level, res_df$log2FoldChange[match(common, res_df$Gene_ID)],
                           method = "spearman", use = "complete.obs"),
    stringsAsFactors = FALSE
  )
}
write.csv(
  bind_rows(sensitivity_summary),
  file.path(table_dir, "sensitivity_34_vs_17.csv"),
  row.names = FALSE
)

figure_stats <- bind_rows(
  data.frame(
    figure = "01_QC_six_panel.png",
    n = "34 runs; Day0=6, Day1=10, Day2=10, Day3=8",
    test = "Descriptive; no group test",
    threshold_or_result = "Input/Clean reads, Q30, counts, nGene, 1M counts and nGene",
    stringsAsFactors = FALSE
  ),
  data.frame(
    figure = c(
      "02_dimensional_reduction_4panel.png",
      "02a_PCA_original_style.png",
      "02b_UMAP_original_style.png",
      "02c_tSNE_original_style.png",
      "02d_LDA_original_style.png"
    ),
    n = "34 runs",
    test = "Unsupervised dimension reduction; no significance test",
    threshold_or_result = c(
      "PCA/UMAP/t-SNE/LDA with sample labels",
      "Top 1,000 variable genes",
      "8 nearest neighbours",
      "perplexity=12",
      "20 PCs as LDA predictors"
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    figure = paste0("03_volcano_", names(comparisons), ".png"),
    n = c("Day0 n=6 vs Day1 n=10", "Day0 n=6 vs Day2 n=10", "Day0 n=6 vs Day3 n=8"),
    test = "DESeq2 Wald test; BH correction; batch-adjusted",
    threshold_or_result = "padj<0.05; |log2FC|>=0.3",
    stringsAsFactors = FALSE
  ),
  data.frame(
    figure = c(
      "04_DE_counts_original_style.png",
      "05_candidate_gene_bars.png",
      "06a_Mfuzz_plot2.png",
      "06b_Mfuzz_trajectory.png",
      "06c_Mfuzz_heatmap.png",
      "07_GO_KEGG_bubble_original_style.png",
      "08a_GO_KEGG_term_gene_links.png",
      "08b_GO_KEGG_sankey.png"
    ),
    n = c(
      "34 runs",
      "34 runs; 8 curated candidate genes",
      "4 Day means; top 2,000 genes",
      "4 Day means; top 2,000 genes",
      "Top 12 genes per cluster",
      "DEGs by comparison and direction",
      "Day2 DEGs",
      "Day2 DEGs"
    ),
    test = c(
      "DESeq2 Wald test; BH correction",
      "Mean and SD; descriptive",
      "Fuzzy c-means; no p-value",
      "Fuzzy c-means; no p-value",
      "Fuzzy c-means; descriptive",
      "g:Profiler hypergeometric; g:SCS correction",
      "g:Profiler hypergeometric; g:SCS correction",
      "g:Profiler hypergeometric; g:SCS correction"
    ),
    threshold_or_result = c(
      "padj<0.05; |log2FC|>=0.3",
      "VST expression; points are runs",
      "c=6",
      "c=6; membership coloring",
      "Z-score",
      "p<0.05; top 8 per comparison/direction",
      "p<0.05; top 5 Day2 terms",
      "p<0.05; top 5 Day2 terms"
    ),
    stringsAsFactors = FALSE
  )
)
if (exists("p_cluster_enrichment")) {
  figure_stats <- bind_rows(
    figure_stats,
    data.frame(
      figure = "09_Mfuzz_cluster_enrichment.png",
      n = "Mfuzz clusters 1-6",
      test = "g:Profiler hypergeometric; g:SCS correction",
      threshold_or_result = "p<0.05; top 5 terms per cluster",
      stringsAsFactors = FALSE
    )
  )
}
write.csv(
  figure_stats,
  file.path(table_dir, "figure_statistics.csv"),
  row.names = FALSE
)

summary_lines <- c(
  "Soybean germination RNA-seq, two sequencing batches",
  "Run-level samples: 34 (17 base samples x 2 batches)",
  paste0("Raw genes: ", nrow(count_mat)),
  paste0("Filtered genes: ", nrow(filtered_counts)),
  "Filter: counts >= 10 in at least six runs",
  "Primary model: DESeq2 ~ Batch + Day",
  "Comparisons: Day1 vs Day0, Day2 vs Day0, Day3 vs Day0",
  "DE threshold: BH padj < 0.05 and |log2FoldChange| >= 0.3",
  "Gene labels: GTF gene_name where available; Glyma.xxGxxxxxx fallback",
  "Enrichment: g:Profiler organism gmax, GO BP/CC/MF and KEGG",
  "Sensitivity model: same-timepoint B1+B2 counts summed to 17 biological samples",
  "Interpretation: exploratory; independent validation still required"
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
writeLines(as.character(packageVersion("DESeq2")), file.path(output_dir, "DESeq2_version.txt"))

message("Completed: ", output_dir)
