#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(123)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(gghalves)
  library(ggh4x)
  library(patchwork)
  library(DESeq2)
  library(Mfuzz)
  library(Biobase)
  library(uwot)
  library(Rtsne)
  library(slingshot)
  library(TSCAN)
  library(SCORPIUS)
  library(SingleCellExperiment)
  library(tradeSeq)
  library(ggbeeswarm)
  library(rtracklayer)
  library(httr)
  library(jsonlite)
  library(ggalluvial)
  library(scales)
})

project_dir <- normalizePath(
  Sys.getenv("PROCESSED_DATA_DIR", "data/doudou_all34"),
  winslash = "/",
  mustWork = TRUE
)
output_dir <- normalizePath(
  Sys.getenv(
    "PSEUDOTIME_OUTPUT_DIR",
    "outputs/doudou_selected5_pseudotime"
  ),
  winslash = "/",
  mustWork = FALSE
)
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

day_levels <- c("Day0", "Day1", "Day2", "Day3")
comparisons <- c(
  Day1_vs_Day0 = "Day1",
  Day2_vs_Day0 = "Day2",
  Day3_vs_Day0 = "Day3"
)
day_colors <- c(
  Day0 = "#7895C1",
  Day1 = "#7BC0CD",
  Day2 = "#ECB66C",
  Day3 = "#ED8D5A"
)
traditional_colors <- c(
  "#BFDFD2", "#51999F", "#4198AC", "#7BC0CD",
  "#DBCB92", "#ECB66C", "#EA9E58", "#ED8D5A"
)

theme_f4 <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black")
    )
}

save_plot <- function(plot, file_name, width = 8, height = 6, dpi = 600) {
  ggsave(
    file.path(figure_dir, file_name),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

message("[1/10] Building selected and excluded sample manifests")
raw_counts <- read.csv(
  file.path(project_dir, "ALL_34_raw_counts_matrix.csv"),
  check.names = FALSE
)
counts_1m <- read.csv(
  file.path(project_dir, "ALL_34_1M_counts_matrix.csv"),
  check.names = FALSE
)
metadata <- read.csv(
  file.path(project_dir, "ALL_34_sample_metadata.csv"),
  check.names = FALSE
)
qc_table <- read.csv(
  file.path(project_dir, "ALL_34_QC_metrics.csv"),
  check.names = FALSE
)
symbol_map <- read.csv(
  file.path(project_dir, "soybean_gene_symbol_map.csv"),
  check.names = FALSE
)

annotation_audit <- data.frame(
  Category = c(
    "All GTF gene models",
    "Protein-coding gene models",
    "Non-coding RNA gene models"
  ),
  Count = c(
    nrow(symbol_map),
    sum(symbol_map$Gene_biotype == "protein_coding"),
    sum(symbol_map$Gene_biotype != "protein_coding")
  )
)
write.csv(
  annotation_audit,
  file.path(table_dir, "soybean_gene_model_audit.csv"),
  row.names = FALSE
)

protein_coding_ids <- symbol_map$Gene_ID[
  symbol_map$Gene_biotype == "protein_coding"
]
raw_counts <- raw_counts[raw_counts$Gene_ID %in% protein_coding_ids, , drop = FALSE]
counts_1m <- counts_1m[counts_1m$Gene_ID %in% protein_coding_ids, , drop = FALSE]
symbol_map <- symbol_map[symbol_map$Gene_ID %in% protein_coding_ids, , drop = FALSE]

metadata$Day <- factor(
  sub("d", "Day", metadata$Timepoint),
  levels = day_levels
)
metadata$Batch <- factor(metadata$Batch, levels = c("B1", "B2"))
metadata$DisplaySample <- sub(
  "^XM-QMH-d([0-3])-",
  "Day\\1-",
  metadata$RunSample
)
metadata$Delivery <- ifelse(
  metadata$Batch == "B1",
  "RJYY-20260923-L-02",
  "RJYY-20261015-L-01"
)

selected_runs <- c(
  "XM-QMH-d0-1_B1", "XM-QMH-d0-1_B2",
  "XM-QMH-d0-4_B2",
  "XM-QMH-d0-5_B1", "XM-QMH-d0-5_B2",
  "XM-QMH-d1-2_B1", "XM-QMH-d1-2_B2",
  "XM-QMH-d1-4_B1", "XM-QMH-d1-4_B2",
  "XM-QMH-d1-5_B1",
  "XM-QMH-d2-2_B1", "XM-QMH-d2-2_B2",
  "XM-QMH-d2-3_B1", "XM-QMH-d2-3_B2",
  "XM-QMH-d2-5_B2",
  "XM-QMH-d3-2_B2",
  "XM-QMH-d3-3_B1", "XM-QMH-d3-3_B2",
  "XM-QMH-d3-5_B1", "XM-QMH-d3-5_B2"
)
if (any(!selected_runs %in% metadata$RunSample)) {
  stop("Selected sample list contains unknown run IDs.")
}
if (any(table(metadata$Day[metadata$RunSample %in% selected_runs]) != 5)) {
  stop("Selected micro-needle group must contain exactly five runs per Day.")
}

manifest <- metadata %>%
  left_join(
    qc_table %>%
      dplyr::select(
        RunSample,
        RawFastQC_R1,
        RawFastQC_R2,
        CleanFastQC_R1,
        CleanFastQC_R2
      ),
    by = "RunSample"
  ) %>%
  mutate(
    Selection = ifelse(
      RunSample %in% selected_runs,
      "Microseed_analysis_n5",
      "Excluded_bulk_simulation_pool"
    )
  )
excluded_runs <- manifest$RunSample[
  manifest$Selection == "Excluded_bulk_simulation_pool"
]
write.csv(
  manifest,
  file.path(table_dir, "sample_selection_manifest.csv"),
  row.names = FALSE
)
write.csv(
  manifest %>% filter(Selection == "Excluded_bulk_simulation_pool"),
  file.path(table_dir, "excluded_sample_pool.csv"),
  row.names = FALSE
)
write.csv(
  manifest %>% filter(Selection == "Microseed_analysis_n5"),
  file.path(table_dir, "selected_microseed_n5.csv"),
  row.names = FALSE
)
cat("Selected:", length(selected_runs), " Excluded:", length(excluded_runs), "\n")

symbol_lookup <- setNames(symbol_map$Symbol, symbol_map$Gene_ID)
gtf_name_lookup <- setNames(symbol_map$Gene_has_GTF_name, symbol_map$Gene_ID)

message("[2/10] Subsetting selected samples and running DESeq2")
raw_mat <- as.matrix(raw_counts[, -1, drop = FALSE])
storage.mode(raw_mat) <- "double"
raw_mat[is.na(raw_mat)] <- 0
raw_mat <- round(raw_mat)
rownames(raw_mat) <- raw_counts$Gene_ID

one_m_mat <- as.matrix(counts_1m[, -1, drop = FALSE])
storage.mode(one_m_mat) <- "double"
rownames(one_m_mat) <- counts_1m$Gene_ID
one_m_mat[is.na(one_m_mat)] <- 0

qc_table$TotalCounts <- colSums(raw_mat[, qc_table$RunSample, drop = FALSE])
qc_table$nGeneRaw <- colSums(raw_mat[, qc_table$RunSample, drop = FALSE] > 0)
qc_table$nGeneGe10 <- colSums(raw_mat[, qc_table$RunSample, drop = FALSE] >= 10)
qc_table$TotalCounts1M <- colSums(one_m_mat[, qc_table$RunSample, drop = FALSE])
qc_table$nGene1M <- colSums(one_m_mat[, qc_table$RunSample, drop = FALSE] > 0)

selected_meta <- metadata[match(selected_runs, metadata$RunSample), , drop = FALSE]
rownames(selected_meta) <- selected_meta$RunSample
selected_counts <- raw_mat[, selected_runs, drop = FALSE]
selected_1m <- one_m_mat[, selected_runs, drop = FALSE]
qc_selected <- qc_table[match(selected_runs, qc_table$RunSample), , drop = FALSE]

keep <- rowSums(selected_counts >= 10) >= 5
selected_counts_f <- selected_counts[keep, , drop = FALSE]
dds <- DESeqDataSetFromMatrix(
  countData = selected_counts_f,
  colData = selected_meta,
  design = ~ Batch + Day
)
dds <- DESeq(dds, parallel = FALSE)
selected_vst <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
selected_norm <- counts(dds, normalized = TRUE)

write.csv(
  data.frame(
    Gene_ID = rownames(selected_counts),
    Symbol = unname(symbol_lookup[rownames(selected_counts)]),
    selected_counts,
    check.names = FALSE
  ),
  file.path(table_dir, "selected_n5_raw_counts.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    Gene_ID = rownames(selected_vst),
    Symbol = unname(symbol_lookup[rownames(selected_vst)]),
    selected_vst,
    check.names = FALSE
  ),
  file.path(table_dir, "selected_n5_VST_expression.csv"),
  row.names = FALSE
)

message("[3/10] Selected-sample QC")
write.csv(
  qc_selected,
  file.path(table_dir, "selected_n5_QC_metrics.csv"),
  row.names = FALSE
)
qc_selected$InputReadsM <- qc_selected$InputReads / 1e6
qc_selected$CleanReadsM <- qc_selected$CleanReads / 1e6
qc_selected$TotalCountsM <- qc_selected$TotalCounts / 1e6
qc_selected$TotalCounts1MM <- qc_selected$TotalCounts1M / 1e6

qc_metric_panel <- function(data, y, label, color) {
  data$.all_runs <- "All runs"
  complete_n <- sum(is.finite(data[[y]]))
  ggplot(data, aes(x = .all_runs, y = .data[[y]])) +
    geom_violin(width = 0.48, fill = color, color = NA, alpha = 0.28) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
    geom_quasirandom(width = 0.19, size = 1.8, color = "black", alpha = 0.8) +
    labs(x = NULL, y = label) +
    theme_f4(11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, size = 8)
    ) +
    labs(subtitle = paste0("visible runs = ", complete_n))
}

p_qc_selected <- wrap_plots(
  list(
    qc_metric_panel(qc_selected, "InputReadsM", "Input reads (million)", "#992224"),
    qc_metric_panel(qc_selected, "CleanReadsM", "Clean reads (million)", "#EF8B67"),
    qc_metric_panel(qc_selected, "CleanReadRetentionPct", "Clean-read retention (%)", "#F0C284"),
    qc_metric_panel(qc_selected, "Q30ReadPct", "Q30 reads (%)", "#F5EBAE"),
    qc_metric_panel(qc_selected, "TotalCountsM", "Assigned counts (million)", "#D6EFF4"),
    qc_metric_panel(qc_selected, "nGeneRaw", "Detected genes", "#7895C1"),
    qc_metric_panel(qc_selected, "TotalCounts1MM", "Counts at 1M (million)", "#A8CBDF"),
    qc_metric_panel(qc_selected, "nGene1M", "Detected genes at 1M", "#8074C8")
  ),
  ncol = 4
) +
  plot_annotation(
    title = "Micro-needle selected samples: eight QC panels",
    subtitle = "Single-metric panels only; 20 selected runs; no grouping by Day"
  )
save_plot(p_qc_selected, "01_selected_QC_eight_panel.png", 13, 7)

message("[4/10] PCA, UMAP, t-SNE and LDA for selected samples")
gene_var <- apply(selected_vst, 1, var)
top_n <- min(1000, sum(is.finite(gene_var) & gene_var > 0))
top_genes_var <- names(sort(gene_var, decreasing = TRUE))[seq_len(top_n)]

pca <- prcomp(t(selected_vst[top_genes_var, , drop = FALSE]), scale. = FALSE)
pca_df <- data.frame(
  RunSample = rownames(pca$x),
  DisplaySample = selected_meta[rownames(pca$x), "DisplaySample"],
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Day = selected_meta[rownames(pca$x), "Day"],
  Batch = selected_meta[rownames(pca$x), "Batch"]
)

set.seed(123)
umap_fit <- uwot::umap(
  selected_vst[top_genes_var, , drop = FALSE] %>% t(),
  n_neighbors = 5,
  min_dist = 0.25,
  metric = "cosine",
  n_components = 2,
  verbose = FALSE
)
colnames(umap_fit) <- c("UMAP1", "UMAP2")
umap_df <- data.frame(
  RunSample = rownames(umap_fit),
  DisplaySample = selected_meta[rownames(umap_fit), "DisplaySample"],
  umap_fit,
  Day = selected_meta[rownames(umap_fit), "Day"],
  Batch = selected_meta[rownames(umap_fit), "Batch"]
)

set.seed(123)
tsne_fit <- Rtsne::Rtsne(
  t(selected_vst[top_genes_var, , drop = FALSE]),
  dims = 2,
  perplexity = 5,
  theta = 0,
  check_duplicates = FALSE,
  pca = TRUE,
  max_iter = 1500
)$Y
colnames(tsne_fit) <- c("tSNE1", "tSNE2")
rownames(tsne_fit) <- colnames(selected_vst)
tsne_df <- data.frame(
  RunSample = rownames(tsne_fit),
  DisplaySample = selected_meta[rownames(tsne_fit), "DisplaySample"],
  tsne_fit,
  Day = selected_meta[rownames(tsne_fit), "Day"],
  Batch = selected_meta[rownames(tsne_fit), "Batch"]
)

lda_fit <- MASS::lda(
  x = umap_fit,
  grouping = selected_meta[rownames(umap_fit), "Day"]
)
lda_coords <- predict(lda_fit)$x[, 1:2, drop = FALSE]
colnames(lda_coords) <- c("LD1", "LD2")
lda_df <- data.frame(
  RunSample = rownames(lda_coords),
  DisplaySample = selected_meta[rownames(lda_coords), "DisplaySample"],
  lda_coords,
  Day = selected_meta[rownames(lda_coords), "Day"],
  Batch = selected_meta[rownames(lda_coords), "Batch"]
)

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(2.4, "cm")
)
dim_plot <- function(data, x, y, title) {
  ggplot(data, aes(.data[[x]], .data[[y]], color = Day, shape = Batch)) +
    geom_point(size = 3.2) +
    geom_text_repel(
      aes(label = DisplaySample),
      size = 2.4,
      max.overlaps = Inf,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    scale_color_manual(values = day_colors) +
    scale_shape_manual(values = c(B1 = 16, B2 = 17)) +
    guides(x = axis, y = axis) +
    labs(x = x, y = y, title = title) +
    theme_f4(15) +
    theme(
      aspect.ratio = 1,
      axis.line = element_line(
        arrow = arrow(type = "closed", length = unit(0.18, "inches")),
        linewidth = 0.8
      ),
      axis.text = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5)
    )
}

p_pca <- dim_plot(pca_df, "PC1", "PC2", "PCA (visualization only)")
p_umap <- dim_plot(umap_df, "UMAP1", "UMAP2", "UMAP")
p_tsne <- dim_plot(tsne_df, "tSNE1", "tSNE2", "t-SNE")
p_lda <- dim_plot(lda_df, "LD1", "LD2", "LDA from UMAP")
p_dims <- (p_pca | p_umap) / (p_tsne | p_lda)
save_plot(p_dims, "02_selected_dimension_reduction_4panel.png", 14, 12)
write.csv(pca_df, file.path(table_dir, "selected_PCA_coordinates.csv"), row.names = FALSE)
write.csv(umap_df, file.path(table_dir, "selected_UMAP_coordinates.csv"), row.names = FALSE)
write.csv(tsne_df, file.path(table_dir, "selected_tSNE_coordinates.csv"), row.names = FALSE)
write.csv(lda_df, file.path(table_dir, "selected_LDA_coordinates.csv"), row.names = FALSE)

message("[5/10] Differential expression in selected samples")
de_results <- list()
de_summary <- list()
for (contrast_name in names(comparisons)) {
  target <- comparisons[[contrast_name]]
  res <- results(dds, contrast = c("Day", target, "Day0"), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$Gene_ID <- rownames(res_df)
  res_df$Symbol <- unname(symbol_lookup[res_df$Gene_ID])
  res_df$SymbolLabel <- ifelse(
    unname(gtf_name_lookup[res_df$Gene_ID]) == "yes",
    res_df$Symbol,
    NA_character_
  )
  res_df$direction <- "ns"
  res_df$direction[
    !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 0.3
  ] <- "up"
  res_df$direction[
    !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -0.3
  ] <- "down"
  res_df$comparison <- contrast_name
  res_df$neglog10p <- ifelse(
    res_df$pvalue > 0,
    -log10(res_df$pvalue),
    -log10(2) - pnorm(-abs(res_df$stat), log.p = TRUE) / log(10)
  )
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
  de_summary[[contrast_name]] <- data.frame(
    comparison = contrast_name,
    up = sum(res_df$direction == "up"),
    down = sum(res_df$direction == "down"),
    tested = sum(!is.na(res_df$padj))
  )

  plot_df <- res_df %>%
    filter(
      !is.na(log2FoldChange),
      !is.na(neglog10p),
      is.finite(log2FoldChange),
      is.finite(neglog10p)
    )
  labels <- plot_df %>%
    filter(!is.na(SymbolLabel), direction != "ns") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 12)
  p <- ggplot(plot_df, aes(log2FoldChange, neglog10p, color = log2FoldChange, size = neglog10p)) +
    geom_point(alpha = 0.75) +
    geom_text_repel(
      data = labels,
      aes(label = SymbolLabel),
      color = "black",
      size = 3.5,
      max.overlaps = Inf,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 4) +
    scale_color_gradientn(colours = c("#7BC0CD", "#ED8D5A")) +
    scale_size(range = c(1, 7)) +
    labs(
      x = expression(log[2]~fold~change),
      y = expression(-log[10](italic(p))),
      color = "log2FC",
      size = "-log10(p)",
      title = gsub("_", " ", contrast_name),
      subtitle = "DESeq2 Wald test; BH padj<0.05; |log2FC|>=0.3; n=5/group"
    ) +
    theme_bw(15) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey30")
    )
  save_plot(p, paste0("03_volcano_", contrast_name, ".png"), 8, 6.5)
}
de_summary_df <- bind_rows(de_summary)
write.csv(de_summary_df, file.path(table_dir, "DE_summary_selected_n5.csv"), row.names = FALSE)

significant_union <- unique(unlist(lapply(
  de_results,
  function(x) x$Gene_ID[x$direction != "ns"]
)))

message("[6/10] Candidate genes and ANOVA")
curated <- data.frame(
  Gene_ID = c(
    "GLYMA_13G069900", "GLYMA_19G006400", "GLYMA_06G303200",
    "GLYMA_06G114000", "GLYMA_18G061100", "GLYMA_13G342900",
    "GLYMA_07G132600", "GLYMA_17G027400"
  ),
  Category = c(
    "GA-responsive cell expansion", "Auxin signaling",
    "Glyoxylate cycle / lipid mobilization",
    "NAC developmental transcription factor",
    "AS1-like developmental regulator", "TOP6B organ growth",
    "Expansin / cell-wall loosening", "ABA/desiccation-associated LEA"
  )
)
curated$Symbol <- unname(symbol_lookup[curated$Gene_ID])
curated$Label <- ifelse(
  unname(gtf_name_lookup[curated$Gene_ID]) == "yes",
  curated$Symbol,
  NA_character_
)
curated <- curated[!is.na(curated$Label) & curated$Gene_ID %in% rownames(selected_vst), , drop = FALSE]

anova_rows <- list()
candidate_plot_rows <- list()
for (i in seq_len(nrow(curated))) {
  gene_id <- curated$Gene_ID[i]
  values <- as.numeric(selected_vst[gene_id, ])
  plot_df <- data.frame(
    Expression = values,
    Day = selected_meta[colnames(selected_vst), "Day"],
    RunSample = colnames(selected_vst)
  )
  fit <- aov(Expression ~ Day, data = plot_df)
  p_anova <- summary(fit)[[1]][["Pr(>F)"]][1]
  anova_rows[[i]] <- data.frame(
    Gene_ID = gene_id,
    Symbol = curated$Symbol[i],
    Category = curated$Category[i],
    ANOVA_F = summary(fit)[[1]][["F value"]][1],
    ANOVA_p = p_anova,
    n_per_group = 5
  )
  label <- if (is.na(p_anova)) {
    "ANOVA unavailable"
  } else {
    paste0("ANOVA p = ", formatC(p_anova, format = "e", digits = 3))
  }
  y_pos <- max(plot_df$Expression, na.rm = TRUE)
  candidate_plot_rows[[i]] <- ggplot(plot_df, aes(Day, Expression, fill = Day)) +
    stat_summary(fun = mean, geom = "col", width = 0.62, color = "black") +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = "errorbar", width = 0.15) +
    geom_jitter(width = 0.12, size = 1.4, shape = 21, fill = "white", color = "black") +
    annotate("text", x = 2.5, y = y_pos, label = label, size = 2.8, vjust = -0.5) +
    scale_fill_manual(values = day_colors, guide = "none") +
    labs(x = NULL, y = "VST expression", title = curated$Symbol[i]) +
    theme_bw(10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6),
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}
anova_df <- bind_rows(anova_rows)
write.csv(anova_df, file.path(table_dir, "candidate_gene_one_way_ANOVA.csv"), row.names = FALSE)
save_plot(wrap_plots(candidate_plot_rows, ncol = 4), "04_candidate_gene_ANOVA.png", 12, 7)

pairwise_rows <- list()
for (i in seq_len(nrow(curated))) {
  gene_id <- curated$Gene_ID[i]
  plot_df <- data.frame(
    Expression = as.numeric(selected_vst[gene_id, ]),
    Day = selected_meta[colnames(selected_vst), "Day"]
  )
  fit <- aov(Expression ~ Day, data = plot_df)
  tukey <- TukeyHSD(fit, "Day")$Day
  pairwise_rows[[i]] <- data.frame(
    Gene_ID = gene_id,
    Symbol = curated$Symbol[i],
    Category = curated$Category[i],
    Pair = rownames(tukey),
    Difference = tukey[, "diff"],
    Lower95 = tukey[, "lwr"],
    Upper95 = tukey[, "upr"],
    Tukey_p = tukey[, "p adj"],
    n_per_group = 5,
    check.names = FALSE
  )
}
pairwise_df <- bind_rows(pairwise_rows)
pairwise_df$Significance <- case_when(
  pairwise_df$Tukey_p < 0.001 ~ "***",
  pairwise_df$Tukey_p < 0.01 ~ "**",
  pairwise_df$Tukey_p < 0.05 ~ "*",
  TRUE ~ "n.s."
)
write.csv(
  pairwise_df,
  file.path(table_dir, "candidate_gene_pairwise_Tukey.csv"),
  row.names = FALSE
)
p_tukey <- ggplot(
  pairwise_df,
  aes(Pair, Symbol, fill = -log10(pmax(Tukey_p, .Machine$double.xmin)))
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(Tukey_p < 0.05, Significance, "")),
    size = 3
  ) +
  scale_fill_gradientn(colours = c("#D6EFF4", "#7BC0CD", "#ED8D5A")) +
  labs(
    x = "Pairwise Day comparison",
    y = NULL,
    fill = "-log10(Tukey p)",
    title = "Pairwise Tukey tests for candidate genes",
    subtitle = "n=5 per Day; BH-adjusted within the six pairwise comparisons"
  ) +
  theme_bw(11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
save_plot(p_tukey, "04b_candidate_pairwise_Tukey.png", 9, 6)

message("[7a/10] Mfuzz on selected samples with symbol labels")
gprofiler_enrich <- function(genes, label) {
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
  if (inherits(response, "error") || httr::status_code(response) != 200) {
    warning(label, ": g:Profiler unavailable or returned an error")
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
  result$cluster <- label
  result
}

time_mean_selected <- sapply(
  day_levels,
  function(day_label) {
    rowMeans(selected_vst[, selected_meta$Day == day_label, drop = FALSE])
  }
)
colnames(time_mean_selected) <- day_levels
mfuzz_variance <- apply(time_mean_selected, 1, var)
mfuzz_gene_ids <- names(sort(mfuzz_variance, decreasing = TRUE))[
  seq_len(min(2000, sum(is.finite(mfuzz_variance) & mfuzz_variance > 0)))
]
mfuzz_eset <- new("ExpressionSet", exprs = time_mean_selected[mfuzz_gene_ids, , drop = FALSE])
mfuzz_eset <- filter.NA(mfuzz_eset, thres = 0.25)
mfuzz_eset <- fill.NA(mfuzz_eset, mode = "mean")
mfuzz_eset <- filter.std(mfuzz_eset, min.std = 0.1, visu = FALSE)
mfuzz_eset <- standardise(mfuzz_eset)
mfuzz_m <- mestimate(mfuzz_eset)
set.seed(123)
mfuzz_fit <- mfuzz(mfuzz_eset, c = 6, m = mfuzz_m)

png(file.path(figure_dir, "04a_Mfuzz_plot2.png"), width = 3000, height = 2200, res = 300)
mfuzz.plot2(
  mfuzz_eset,
  cl = mfuzz_fit,
  mfrow = c(2, 3),
  time.labels = day_levels,
  xlab = "",
  x11 = FALSE
)
dev.off()

mfuzz_membership <- mfuzz_fit$membership
colnames(mfuzz_membership) <- paste0("membership", seq_len(ncol(mfuzz_membership)))
mfuzz_exp <- data.frame(
  exprs(mfuzz_eset),
  Gene_ID = rownames(mfuzz_eset),
  cluster = mfuzz_fit$cluster,
  mfuzz_membership,
  check.names = FALSE
) %>%
  pivot_longer(
    cols = all_of(day_levels),
    names_to = "Day",
    values_to = "Z"
  ) %>%
  mutate(
    Day = factor(Day, levels = day_levels),
    cluster = factor(cluster, levels = 1:6)
  )
mfuzz_exp$Membership <- apply(
  mfuzz_exp[, grepl("^membership", colnames(mfuzz_exp)), drop = FALSE],
  1,
  max
)

set.seed(123)
draw_ids <- mfuzz_exp %>%
  distinct(Gene_ID) %>%
  slice_sample(n = min(800, nrow(.))) %>%
  pull(Gene_ID)
p_mfuzz <- ggplot(
  mfuzz_exp %>% filter(Gene_ID %in% draw_ids),
  aes(Day, Z)
) +
  geom_line(aes(group = Gene_ID, color = Membership), alpha = 0.22, linewidth = 0.35) +
  facet_grid(cluster ~ ., switch = "y") +
  scale_color_gradientn(colours = traditional_colors) +
  labs(
    x = NULL,
    y = "Standardized expression",
    title = "Mfuzz trajectories: selected micro-needle samples",
    subtitle = paste0("c=6; m=", round(mfuzz_m, 2), "; no p-value")
  ) +
  theme_bw(12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
save_plot(p_mfuzz, "04b_Mfuzz_trajectory.png", 8, 10, dpi = 400)

mfuzz_gene_table <- data.frame(
  Gene_ID = names(mfuzz_fit$cluster),
  Symbol = unname(symbol_lookup[names(mfuzz_fit$cluster)]),
  Gene_has_GTF_name = unname(gtf_name_lookup[names(mfuzz_fit$cluster)]),
  cluster = as.integer(mfuzz_fit$cluster),
  maxMembership = apply(mfuzz_fit$membership, 1, max),
  significant_DE = names(mfuzz_fit$cluster) %in% significant_union,
  curated_development_or_hormone = names(mfuzz_fit$cluster) %in% curated$Gene_ID
)
mfuzz_gene_table$Label <- ifelse(
  mfuzz_gene_table$Gene_has_GTF_name == "yes" &
    mfuzz_gene_table$significant_DE &
    mfuzz_gene_table$curated_development_or_hormone,
  mfuzz_gene_table$Symbol,
  ""
)
write.csv(
  mfuzz_gene_table,
  file.path(table_dir, "Mfuzz_selected_gene_clusters.csv"),
  row.names = FALSE
)

heat_top <- mfuzz_gene_table %>%
  group_by(cluster) %>%
  arrange(desc(maxMembership)) %>%
  slice_head(n = 12) %>%
  ungroup()
heat_curated <- mfuzz_gene_table %>%
  filter(significant_DE, curated_development_or_hormone)
heat_genes <- bind_rows(heat_top, heat_curated) %>%
  distinct(Gene_ID, .keep_all = TRUE)
heat_matrix <- exprs(mfuzz_eset)[heat_genes$Gene_ID, , drop = FALSE]
heat_df <- as.data.frame(heat_matrix) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(all_of(day_levels), names_to = "Day", values_to = "Z") %>%
  left_join(heat_genes[, c("Gene_ID", "Label", "cluster")], by = "Gene_ID") %>%
  mutate(
    Day = factor(Day, levels = day_levels),
    cluster = factor(cluster, levels = 1:6),
    Label = ifelse(is.na(Label) | Label == "", "", Label)
  )
p_heat <- ggplot(heat_df, aes(Day, reorder(Gene_ID, Z), fill = Z)) +
  geom_tile() +
  facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradientn(colours = traditional_colors, limits = c(-2.5, 2.5), oob = scales::squish) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Z-score",
    title = "Mfuzz heatmap",
    subtitle = "Only significant curated development/hormone genes are labeled with symbols"
  ) +
  theme_bw(10) +
  theme(
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
save_plot(p_heat, "04c_Mfuzz_heatmap.png", 8, 11, dpi = 400)

cluster_enrichment_list <- list()
for (cluster_id in 1:6) {
  genes <- mfuzz_gene_table$Gene_ID[mfuzz_gene_table$cluster == cluster_id]
  cluster_result <- gprofiler_enrich(genes, paste0("Cluster", cluster_id))
  if (!is.null(cluster_result)) {
    cluster_enrichment_list[[paste0("Cluster", cluster_id)]] <- cluster_result
  }
}
cluster_enrichment <- if (length(cluster_enrichment_list)) bind_rows(cluster_enrichment_list) else data.frame()
if (nrow(cluster_enrichment) > 0) {
  cluster_enrichment[] <- lapply(
    cluster_enrichment,
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
  cluster_enrichment,
  file.path(table_dir, "Mfuzz_selected_cluster_enrichment.csv"),
  row.names = FALSE
)
if (nrow(cluster_enrichment) > 0 && "p_value" %in% colnames(cluster_enrichment)) {
  p_cluster_enrich <- cluster_enrichment %>%
    mutate(minus_log10_p = -log10(pmax(as.numeric(p_value), .Machine$double.xmin))) %>%
    group_by(cluster) %>%
    arrange(p_value) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    ggplot(aes(minus_log10_p, reorder(name, minus_log10_p))) +
    geom_point(aes(size = intersection_size, color = source)) +
    facet_wrap(~cluster, scales = "free_y", ncol = 3) +
    scale_color_manual(
      values = c(`GO:BP` = "#7895C1", `GO:CC` = "#7BC0CD", `GO:MF` = "#F0C284", KEGG = "#ED8D5A")
    ) +
    labs(
      x = "-log10(p)",
      y = NULL,
      size = "Genes",
      color = "Source",
      title = "Mfuzz cluster enrichment",
      subtitle = "g:Profiler gmax; g:SCS correction; p<0.05"
    ) +
    theme_bw(10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(hjust = 0.5)
    )
  save_plot(p_cluster_enrich, "04d_Mfuzz_cluster_enrichment.png", 12, 8)
  p_mfuzz_complete <- (p_mfuzz | p_heat) / p_cluster_enrich +
    plot_annotation(
      title = "Mfuzz complete figure",
      subtitle = "Trajectory + cluster heatmap + cluster functional enrichment",
      tag_levels = "A"
    )
  save_plot(p_mfuzz_complete, "04e_Mfuzz_complete.png", 16, 14, dpi = 300)
}

message("[7/10] Excluded-pool bulk simulation QC/correlation/clustering")
excluded_meta <- metadata[match(excluded_runs, metadata$RunSample), , drop = FALSE]
rownames(excluded_meta) <- excluded_meta$RunSample
excluded_counts <- raw_mat[, excluded_runs, drop = FALSE]
excluded_1m <- one_m_mat[, excluded_runs, drop = FALSE]
qc_excluded <- qc_table[match(excluded_runs, qc_table$RunSample), , drop = FALSE]
write.csv(
  data.frame(
    Gene_ID = rownames(excluded_counts),
    Symbol = unname(symbol_lookup[rownames(excluded_counts)]),
    excluded_counts,
    check.names = FALSE
  ),
  file.path(table_dir, "excluded_pool_raw_counts.csv"),
  row.names = FALSE
)

excluded_keep <- rowSums(excluded_counts >= 10) >= 2
excluded_counts_f <- excluded_counts[excluded_keep, , drop = FALSE]
dds_excluded <- DESeqDataSetFromMatrix(
  countData = excluded_counts_f,
  colData = excluded_meta,
  design = ~ Batch + Day
)
dds_excluded <- DESeq(dds_excluded, parallel = FALSE)
excluded_vst <- assay(varianceStabilizingTransformation(dds_excluded, blind = FALSE))

qc_excluded$InputReadsM <- qc_excluded$InputReads / 1e6
qc_excluded$CleanReadsM <- qc_excluded$CleanReads / 1e6
qc_excluded$TotalCountsM <- qc_excluded$TotalCounts / 1e6
qc_excluded$TotalCounts1MM <- qc_excluded$TotalCounts1M / 1e6
qc_excluded$AllRuns <- "Excluded pool"
p_qc_excluded <- wrap_plots(
  list(
    qc_metric_panel(qc_excluded, "InputReadsM", "Input reads (million)", "#992224"),
    qc_metric_panel(qc_excluded, "CleanReadsM", "Clean reads (million)", "#EF8B67"),
    qc_metric_panel(qc_excluded, "CleanReadRetentionPct", "Clean-read retention (%)", "#F0C284"),
    qc_metric_panel(qc_excluded, "Q30ReadPct", "Q30 reads (%)", "#F5EBAE"),
    qc_metric_panel(qc_excluded, "TotalCountsM", "Assigned counts (million)", "#D6EFF4"),
    qc_metric_panel(qc_excluded, "nGeneRaw", "Detected genes", "#7895C1"),
    qc_metric_panel(qc_excluded, "TotalCounts1MM", "Counts at 1M (million)", "#A8CBDF"),
    qc_metric_panel(qc_excluded, "nGene1M", "Detected genes at 1M", "#8074C8")
  ),
  ncol = 4
) +
  plot_annotation(
    title = "Excluded-sample pool: bulk-like QC",
    subtitle = "14 excluded runs; descriptive only; no Day grouping"
  )
save_plot(p_qc_excluded, "05_excluded_pool_QC.png", 13, 7)

# Compare selected micro-needle runs with the excluded bulk-like pool.
qc_compare <- qc_table %>%
  mutate(
    Selection = ifelse(
      RunSample %in% selected_runs,
      "Microseed selected",
      "Excluded bulk pool"
    ),
    InputReadsM = InputReads / 1e6,
    CleanReadsM = CleanReads / 1e6,
    TotalCountsM = TotalCounts / 1e6,
    TotalCounts1MM = TotalCounts1M / 1e6
  )
qc_compare_metrics <- c(
  "InputReadsM", "CleanReadsM", "CleanReadRetentionPct", "Q30ReadPct",
  "TotalCountsM", "nGeneRaw", "TotalCounts1MM", "nGene1M"
)
qc_metric_labels <- c(
  "Input reads (million)", "Clean reads (million)",
  "Clean-read retention (%)", "Q30 reads (%)",
  "Assigned counts (million)", "Detected genes",
  "Counts at 1M (million)", "Detected genes at 1M"
)
point_audit <- data.frame(
  Metric = qc_compare_metrics,
  Selected_total_rows = vapply(
    qc_compare_metrics,
    function(metric) sum(qc_compare$Selection == "Microseed selected"),
    integer(1)
  ),
  Selected_nonmissing = vapply(
    qc_compare_metrics,
    function(metric) {
      sum(
        qc_compare$Selection == "Microseed selected" &
          !is.na(qc_compare[[metric]])
      )
    },
    integer(1)
  ),
  Selected_finite = vapply(
    qc_compare_metrics,
    function(metric) {
      sum(
        qc_compare$Selection == "Microseed selected" &
          is.finite(qc_compare[[metric]])
      )
    },
    integer(1)
  ),
  Selected_unique_values = vapply(
    qc_compare_metrics,
    function(metric) {
      length(unique(
        qc_compare[[metric]][
          qc_compare$Selection == "Microseed selected"
        ]
      ))
    },
    integer(1)
  ),
  Excluded_nonmissing = vapply(
    qc_compare_metrics,
    function(metric) {
      sum(
        qc_compare$Selection == "Excluded bulk pool" &
          !is.na(qc_compare[[metric]])
      )
    },
    integer(1)
  )
)
write.csv(
  point_audit,
  file.path(table_dir, "QC_point_count_audit.csv"),
  row.names = FALSE
)

comparison_rows <- lapply(seq_along(qc_compare_metrics), function(i) {
  metric <- qc_compare_metrics[i]
  selected_values <- qc_compare[[metric]][
    qc_compare$Selection == "Microseed selected"
  ]
  excluded_values <- qc_compare[[metric]][
    qc_compare$Selection == "Excluded bulk pool"
  ]
  welch <- t.test(selected_values, excluded_values)
  wilcox <- wilcox.test(selected_values, excluded_values, exact = FALSE)
  pooled_sd <- sqrt(
    ((length(selected_values) - 1) * var(selected_values) +
       (length(excluded_values) - 1) * var(excluded_values)) /
      (length(selected_values) + length(excluded_values) - 2)
  )
  data.frame(
    Metric = metric,
    Label = qc_metric_labels[i],
    Selected_n = length(selected_values),
    Excluded_n = length(excluded_values),
    Selected_mean = mean(selected_values),
    Excluded_mean = mean(excluded_values),
    Welch_p = welch$p.value,
    Wilcoxon_p = wilcox$p.value,
    Cohen_d = (mean(selected_values) - mean(excluded_values)) / pooled_sd,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
write.csv(
  comparison_rows,
  file.path(table_dir, "QC_selected_vs_excluded_tests.csv"),
  row.names = FALSE
)

qc_compare_panel <- function(metric, label) {
  ggplot(
    qc_compare,
    aes(Selection, .data[[metric]], fill = Selection)
  ) +
    geom_violin(width = 0.55, alpha = 0.28, color = NA) +
    geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA) +
    geom_quasirandom(
      aes(shape = Batch),
      width = 0.15,
      size = 1.6,
      alpha = 0.85
    ) +
    scale_fill_manual(
      values = c(
        `Microseed selected` = "#7BC0CD",
        `Excluded bulk pool` = "#F0C284"
      ),
      guide = "none"
    ) +
    scale_shape_manual(values = c(B1 = 16, B2 = 17)) +
    labs(x = NULL, y = label, subtitle = "n = 20 vs 14") +
    theme_f4(10) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      axis.ticks.x = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, size = 8)
    )
}
qc_compare_plots <- Map(qc_compare_panel, qc_compare_metrics, qc_metric_labels)
p_qc_compare <- wrap_plots(qc_compare_plots, ncol = 4) +
  plot_annotation(
    title = "QC comparison: selected micro-needle runs vs excluded bulk pool",
    subtitle = "Welch and Wilcoxon tests are reported in QC_selected_vs_excluded_tests.csv"
  )
save_plot(p_qc_compare, "05a_selected_vs_excluded_QC.png", 14, 8)

full_vst_file <- file.path(
  project_dir,
  "downstream_all34/tables/VST_expression_34.csv"
)
full_vst <- read.csv(full_vst_file, check.names = FALSE)
full_vst <- full_vst[full_vst$Gene_ID %in% protein_coding_ids, , drop = FALSE]
full_vst_mat <- as.matrix(full_vst[, -c(1, 2), drop = FALSE])
rownames(full_vst_mat) <- full_vst$Gene_ID
all_selection <- ifelse(
  colnames(full_vst_mat) %in% selected_runs,
  "Microseed selected",
  "Excluded bulk pool"
)
names(all_selection) <- colnames(full_vst_mat)
all_day <- metadata$Day[match(colnames(full_vst_mat), metadata$RunSample)]

all_var <- apply(full_vst_mat, 1, var)
all_top <- names(sort(all_var, decreasing = TRUE))[
  seq_len(min(1000, sum(is.finite(all_var) & all_var > 0)))
]
all_cor <- cor(full_vst_mat[all_top, , drop = FALSE], method = "pearson")
all_cor_df <- as.data.frame(as.table(all_cor))
colnames(all_cor_df) <- c("Sample1", "Sample2", "Correlation")
all_cor_df$Sample1 <- factor(all_cor_df$Sample1, levels = colnames(full_vst_mat))
all_cor_df$Sample2 <- factor(all_cor_df$Sample2, levels = colnames(full_vst_mat))
p_all_cor <- ggplot(all_cor_df, aes(Sample1, Sample2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradientn(colours = traditional_colors, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, title = "All-run expression correlation") +
  theme_minimal(8) +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 5.5),
    axis.text.y = element_text(size = 5.5),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

selected_profile <- rowMeans(full_vst_mat[, names(all_selection)[all_selection == "Microseed selected"], drop = FALSE])
excluded_profile <- rowMeans(full_vst_mat[, names(all_selection)[all_selection == "Excluded bulk pool"], drop = FALSE])
profile_fit <- lm(excluded_profile ~ selected_profile)
profile_r2 <- summary(profile_fit)$r.squared
profile_p <- coefficients(summary(profile_fit))[2, 4]
profile_cor <- cor.test(selected_profile, excluded_profile, method = "pearson")
write.csv(
  data.frame(
    comparison = "Excluded bulk pool vs microseed selected",
    genes = length(selected_profile),
    slope = coefficients(profile_fit)[2],
    intercept = coefficients(profile_fit)[1],
    R2 = profile_r2,
    model_p = profile_p,
    pearson_r = unname(profile_cor$estimate),
    pearson_p = profile_cor$p.value
  ),
  file.path(table_dir, "bulk_vs_micro_expression_profile_regression.csv"),
  row.names = FALSE
)
p_profile_fit <- ggplot(
  data.frame(Selected = selected_profile, Excluded = excluded_profile),
  aes(Selected, Excluded)
) +
  geom_point(size = 1.1, alpha = 0.5, color = "#4198AC") +
  geom_smooth(method = "lm", formula = y ~ x, color = "#ED8D5A", se = TRUE) +
  labs(
    x = "Mean VST, selected micro-needle group",
    y = "Mean VST, excluded bulk pool",
    title = "Group-mean expression linear fit",
    subtitle = sprintf("R² = %.3f; p = %.3e", profile_r2, profile_p)
  ) +
  theme_bw(11) +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

set.seed(123)
all_umap <- uwot::umap(
  t(full_vst_mat[all_top, , drop = FALSE]),
  n_neighbors = 8,
  min_dist = 0.25,
  metric = "cosine",
  verbose = FALSE
)
colnames(all_umap) <- c("UMAP1", "UMAP2")
all_umap_df <- data.frame(
  RunSample = rownames(all_umap),
  UMAP1 = all_umap[, 1],
  UMAP2 = all_umap[, 2],
  Selection = all_selection[rownames(all_umap)],
  Day = all_day[match(rownames(all_umap), colnames(full_vst_mat))],
  Batch = metadata$Batch[match(rownames(all_umap), metadata$RunSample)]
)
p_all_umap <- ggplot(
  all_umap_df,
  aes(UMAP1, UMAP2, color = Selection, shape = Day)
) +
  geom_point(size = 3.2) +
  scale_color_manual(values = c(
    `Microseed selected` = "#7BC0CD",
    `Excluded bulk pool` = "#F0C284"
  )) +
  scale_shape_manual(values = c(Day0 = 16, Day1 = 15, Day2 = 17, Day3 = 18)) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    title = "Combined sample clustering",
    subtitle = "Selected vs excluded; all 34 runs"
  ) +
  theme_f4(11) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))

p_bulk_comparison <- (p_all_cor | p_profile_fit | p_all_umap) +
  plot_annotation(
    title = "Selected micro-needle samples versus excluded bulk pool",
    subtitle = "Bulk-like comparison uses all 34 runs; descriptive and model-based"
  )
save_plot(p_bulk_comparison, "05b_bulk_vs_micro_expression_structure.png", 16, 7)

excluded_var <- apply(excluded_vst, 1, var)
excluded_top <- names(sort(excluded_var, decreasing = TRUE))[seq_len(min(1000, sum(excluded_var > 0)))]
corr_mat <- cor(excluded_vst[excluded_top, , drop = FALSE], method = "pearson")
corr_df <- as.data.frame(as.table(corr_mat))
colnames(corr_df) <- c("Sample1", "Sample2", "Correlation")
corr_df$Sample1 <- factor(corr_df$Sample1, levels = excluded_runs)
corr_df$Sample2 <- factor(corr_df$Sample2, levels = excluded_runs)
p_corr <- ggplot(corr_df, aes(Sample1, Sample2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradientn(colours = traditional_colors, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, title = "Excluded-pool expression correlation") +
  theme_minimal(9) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

paired_bases <- excluded_meta %>%
  group_by(Sample) %>%
  filter(n() == 2, all(c("B1", "B2") %in% Batch)) %>%
  summarise(
    B1 = RunSample[Batch == "B1"][1],
    B2 = RunSample[Batch == "B2"][1],
    .groups = "drop"
  )
pair_metrics <- lapply(seq_len(nrow(paired_bases)), function(i) {
  x <- excluded_vst[excluded_top, paired_bases$B1[i]]
  y <- excluded_vst[excluded_top, paired_bases$B2[i]]
  fit <- lm(y ~ x)
  data.frame(
    Sample = paired_bases$Sample[i],
    B1 = paired_bases$B1[i],
    B2 = paired_bases$B2[i],
    R2 = summary(fit)$r.squared,
    p = coefficients(summary(fit))[2, 4]
  )
}) %>% bind_rows()
write.csv(pair_metrics, file.path(table_dir, "excluded_pool_pair_regression.csv"), row.names = FALSE)
chosen_pair <- pair_metrics[which.min(abs(pair_metrics$R2 - median(pair_metrics$R2))), ]
reg_df <- data.frame(
  B1 = excluded_vst[excluded_top, chosen_pair$B1],
  B2 = excluded_vst[excluded_top, chosen_pair$B2]
)
p_reg <- ggplot(reg_df, aes(B1, B2)) +
  geom_point(size = 1.2, alpha = 0.55, color = "#4198AC") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#ED8D5A") +
  labs(
    x = paste0(chosen_pair$Sample, " B1 (VST)"),
    y = paste0(chosen_pair$Sample, " B2 (VST)"),
    title = "Within-sample B1 vs B2 regression",
    subtitle = sprintf("R2 = %.3f; linear-model p = %.3e", chosen_pair$R2, chosen_pair$p)
  ) +
  theme_bw(11) +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

set.seed(123)
excluded_umap <- uwot::umap(
  t(excluded_vst[excluded_top, , drop = FALSE]),
  n_neighbors = 5,
  min_dist = 0.25,
  metric = "cosine",
  verbose = FALSE
)
colnames(excluded_umap) <- c("UMAP1", "UMAP2")
excluded_umap_df <- data.frame(
  RunSample = rownames(excluded_umap),
  DisplaySample = excluded_meta[rownames(excluded_umap), "DisplaySample"],
  excluded_umap,
  Day = excluded_meta[rownames(excluded_umap), "Day"],
  Batch = excluded_meta[rownames(excluded_umap), "Batch"]
)
p_excluded_umap <- dim_plot(
  excluded_umap_df,
  "UMAP1",
  "UMAP2",
  "Excluded-pool clustering"
)
write.csv(excluded_umap_df, file.path(table_dir, "excluded_pool_UMAP.csv"), row.names = FALSE)

p_excluded_all <- p_qc_excluded /
  ((p_corr | p_reg) / p_excluded_umap + plot_layout(heights = c(1.1, 1))) +
  plot_annotation(tag_levels = "A")
save_plot(p_excluded_all, "05_bulk_pool_QC_regression_clustering.png", 14, 16)

message("[8/10] Pseudotime on UMAP, without PCA coordinates")
umap_pseudo <- uwot::umap(
  t(selected_vst[top_genes_var, , drop = FALSE]),
  n_neighbors = 5,
  min_dist = 0.20,
  metric = "cosine",
  n_components = 2,
  verbose = FALSE
)
colnames(umap_pseudo) <- c("UMAP1", "UMAP2")
set.seed(123)
km <- kmeans(umap_pseudo, centers = 4, nstart = 100)$cluster
names(km) <- rownames(umap_pseudo)

pseudotime_sce <- SingleCellExperiment(
  assays = list(counts = selected_counts_f),
  colData = selected_meta,
  reducedDims = list(UMAP = umap_pseudo)
)
colData(pseudotime_sce)$Kmeans <- factor(km[colnames(selected_counts_f)])
pseudotime_sce <- slingshot(
  pseudotime_sce,
  clusterLabels = "Day",
  reducedDim = "UMAP",
  start.clus = "Day0",
  reweight = FALSE
)

pt_mat <- slingPseudotime(pseudotime_sce, na = TRUE)
day_numeric <- as.numeric(factor(selected_meta$Day, levels = day_levels))
lineage_cors <- apply(pt_mat, 2, function(values) {
  if (sum(!is.na(values)) < 3) {
    return(NA_real_)
  }
  cor(values, day_numeric, method = "spearman", use = "complete.obs")
})

# Branch-aware consensus pseudotime: orient every lineage from Day0 toward later
# Days, scale each lineage to [0, 1], and average across lineages covering a run.
pt_scaled <- pt_mat
for (lineage_id in seq_len(ncol(pt_mat))) {
  values <- pt_mat[, lineage_id]
  if (!is.na(lineage_cors[lineage_id]) && lineage_cors[lineage_id] < 0) {
    values <- max(values, na.rm = TRUE) - values
  }
  valid <- !is.na(values)
  if (any(valid)) {
    value_range <- range(values[valid])
    if (diff(value_range) > 0) {
      pt_scaled[valid, lineage_id] <-
        (values[valid] - value_range[1]) / diff(value_range)
    } else {
      pt_scaled[valid, lineage_id] <- 0
    }
  }
}
pseudotime <- rowMeans(pt_scaled, na.rm = TRUE)
names(pseudotime) <- rownames(pt_mat)
lineage_coverage <- rowSums(!is.na(pt_mat))
coverage_correlation <- suppressWarnings(cor(
  pseudotime,
  day_numeric,
  method = "spearman",
  use = "complete.obs"
))
pt_df <- data.frame(
  RunSample = rownames(pt_mat),
  DisplaySample = selected_meta[rownames(pt_mat), "DisplaySample"],
  Day = selected_meta[rownames(pt_mat), "Day"],
  Batch = selected_meta[rownames(pt_mat), "Batch"],
  Pseudotime = pseudotime,
  UMAP1 = umap_pseudo[, 1],
  UMAP2 = umap_pseudo[, 2],
  Kmeans = km,
  LineageCoverage = lineage_coverage
)
write.csv(pt_df, file.path(table_dir, "selected_n5_pseudotime.csv"), row.names = FALSE)

p_pt_umap <- ggplot(pt_df, aes(UMAP1, UMAP2, color = Pseudotime)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = DisplaySample), color = "black", size = 2.5,
                  max.overlaps = Inf, show.legend = FALSE) +
  scale_color_gradientn(colours = traditional_colors) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "Pseudotime",
    title = "Pseudotime trajectory on UMAP",
    subtitle = sprintf(
      "Consensus across %d lineages; Spearman with Day = %.3f",
      ncol(pt_mat),
      coverage_correlation
    )
  ) +
  theme_f4(13) +
  theme(aspect.ratio = 1, plot.subtitle = element_text(hjust = 0.5, size = 9))
for (curve in slingCurves(pseudotime_sce)$curves) {
  p_pt_umap <- p_pt_umap +
    annotate(
      "path",
      x = curve$s[, 1],
      y = curve$s[, 2],
      color = "#B84D28",
      linewidth = 1.1
    )
}
save_plot(p_pt_umap, "06_pseudotime_UMAP.png", 9, 7.5)

pt_complete <- pt_df %>% filter(!is.na(Pseudotime))
p_pt_day <- ggplot(pt_complete, aes(Day, Pseudotime, color = Day)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.25) +
  geom_quasirandom(width = 0.15, size = 2.5) +
  scale_color_manual(values = day_colors, guide = "none") +
  labs(
    x = NULL,
    y = "Pseudotime",
    title = "Pseudotime progression by Day",
    subtitle = sprintf(
      "n = %s; Spearman rho=%.3f",
      paste(table(pt_complete$Day), collapse = " / "),
      coverage_correlation
    )
  ) +
  theme_f4(13) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 9))
save_plot(p_pt_day, "07_pseudotime_by_Day.png", 6.5, 5.5)

message("[9/10] Pseudotime-associated genes with tradeSeq")
model_genes <- names(sort(apply(selected_counts_f, 1, var), decreasing = TRUE))[1:min(500, nrow(selected_counts_f))]
model_samples <- names(pseudotime)[!is.na(pseudotime)]
model_counts <- selected_counts_f[model_genes, model_samples, drop = FALSE]
model_pt <- matrix(pseudotime[model_samples], ncol = 1)
rownames(model_pt) <- model_samples
colnames(model_pt) <- "Lineage1"
model_weights <- matrix(1, nrow = length(model_samples), ncol = 1)
rownames(model_weights) <- colnames(model_counts)
colnames(model_weights) <- "Lineage1"

association_df <- tryCatch({
  set.seed(123)
  gam_sce <- fitGAM(
    counts = model_counts,
    pseudotime = model_pt,
    cellWeights = model_weights,
    nknots = 4,
    verbose = FALSE
  )
  assoc <- associationTest(gam_sce)
  assoc$Gene_ID <- rownames(assoc)
  assoc$Symbol <- unname(symbol_lookup[assoc$Gene_ID])
  assoc$Has_GTF_name <- unname(gtf_name_lookup[assoc$Gene_ID])
  assoc
}, error = function(error) {
  warning("tradeSeq fitGAM failed: ", error$message)
  data.frame()
})
write.csv(
  association_df,
  file.path(table_dir, "pseudotime_association_tradeSeq.csv"),
  row.names = FALSE
)

if (nrow(association_df) > 0) {
  assoc_plot <- association_df %>%
    filter(Has_GTF_name == "yes") %>%
    arrange(pvalue) %>%
    slice_head(n = 6)
  trend_rows <- list()
  for (i in seq_len(nrow(assoc_plot))) {
    gene_id <- assoc_plot$Gene_ID[i]
    trend_rows[[i]] <- data.frame(
      Gene_ID = gene_id,
      Symbol = assoc_plot$Symbol[i],
      Pseudotime = pseudotime,
      Expression = as.numeric(selected_vst[gene_id, names(pseudotime)])
    )
  }
  trend_df <- bind_rows(trend_rows)
  p_trend <- ggplot(trend_df, aes(Pseudotime, Expression)) +
    geom_point(size = 1.8, color = "#4198AC", alpha = 0.7) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE,
                color = "#ED8D5A") +
    facet_wrap(~Symbol, scales = "free_y", ncol = 3) +
    labs(
      x = "Pseudotime",
      y = "VST expression",
      title = "Top named genes associated with pseudotime",
      subtitle = "tradeSeq association test; named GTF symbols only"
    ) +
    theme_bw(11) +
    theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
  save_plot(p_trend, "08_pseudotime_top_gene_trends.png", 10, 7)
}

message("[9b/10] Cross-checking trajectories with TSCAN and SCORPIUS")
tscan_genes <- top_genes_var[seq_len(min(500, length(top_genes_var)))]
tscan_counts <- selected_counts_f[tscan_genes, , drop = FALSE]
tscan_proc <- TSCAN::preprocess(tscan_counts)
tscan_clusters <- TSCAN::exprmclust(
  tscan_proc,
  clusternum = 2:4,
  modelNames = "EII"
)
png(file.path(figure_dir, "09a_TSCAN_state_tree.png"),
    width = 2400, height = 1800, res = 300)
TSCAN::plotmclust(tscan_clusters)
dev.off()
tscan_order <- TSCAN::TSCANorder(
  tscan_clusters,
  orderonly = FALSE,
  flip = TRUE
)
tscan_pt <- setNames(tscan_order$Pseudotime, tscan_order$sample_name)
tscan_pt <- tscan_pt[selected_runs]
tscan_day_cor <- cor(
  tscan_pt,
  as.numeric(selected_meta[names(tscan_pt), "Day"]),
  method = "spearman",
  use = "complete.obs"
)
if (tscan_day_cor < 0) {
  tscan_pt <- max(tscan_pt, na.rm = TRUE) + min(tscan_pt, na.rm = TRUE) - tscan_pt
}
write.csv(
  data.frame(
    RunSample = names(tscan_pt),
    Day = selected_meta[names(tscan_pt), "Day"],
    TSCAN_pseudotime = as.numeric(tscan_pt),
    TSCAN_state = tscan_order$State[match(names(tscan_pt), tscan_order$sample_name)]
  ),
  file.path(table_dir, "TSCAN_pseudotime.csv"),
  row.names = FALSE
)

scorpius_expression <- t(selected_vst[top_genes_var, , drop = FALSE])
scorpius_space <- SCORPIUS::reduce_dimensionality(scorpius_expression, "spearman")
scorpius_traj <- SCORPIUS::infer_trajectory(scorpius_space)
scorpius_pt <- scorpius_traj$time
names(scorpius_pt) <- rownames(scorpius_space)
scorpius_day_cor <- cor(
  scorpius_pt,
  as.numeric(selected_meta[names(scorpius_pt), "Day"]),
  method = "spearman",
  use = "complete.obs"
)
if (scorpius_day_cor < 0) {
  scorpius_pt <- max(scorpius_pt, na.rm = TRUE) - scorpius_pt
}
png(file.path(figure_dir, "09b_SCORPIUS_trajectory.png"),
    width = 2400, height = 1800, res = 300)
print(
  SCORPIUS::draw_trajectory_plot(
    scorpius_space,
    progression_group = selected_meta[rownames(scorpius_space), "Day"],
    progression_group_palette = day_colors,
    path = scorpius_traj$path,
    contour = TRUE
  )
)
dev.off()
write.csv(
  data.frame(
    RunSample = names(scorpius_pt),
    Day = selected_meta[names(scorpius_pt), "Day"],
    SCORPIUS_pseudotime = as.numeric(scorpius_pt)
  ),
  file.path(table_dir, "SCORPIUS_pseudotime.csv"),
  row.names = FALSE
)

method_df <- data.frame(
  RunSample = selected_runs,
  Day = selected_meta[selected_runs, "Day"],
  Slingshot = pseudotime[selected_runs],
  TSCAN = as.numeric(tscan_pt[selected_runs]),
  SCORPIUS = as.numeric(scorpius_pt[selected_runs]),
  UMAP1 = umap_pseudo[selected_runs, 1],
  UMAP2 = umap_pseudo[selected_runs, 2]
)
write.csv(method_df, file.path(table_dir, "pseudotime_method_comparison.csv"), row.names = FALSE)
method_cor <- cor(
  method_df[, c("Slingshot", "TSCAN", "SCORPIUS")],
  method = "spearman",
  use = "pairwise.complete.obs"
)
write.csv(
  data.frame(method = rownames(method_cor), method_cor, check.names = FALSE),
  file.path(table_dir, "pseudotime_method_correlation.csv"),
  row.names = FALSE
)

method_long <- method_df %>%
  pivot_longer(
    cols = c(Slingshot, TSCAN, SCORPIUS),
    names_to = "Method",
    values_to = "Pseudotime"
  ) %>%
  group_by(Method) %>%
  mutate(
    Pseudotime01 = ifelse(
      diff(range(Pseudotime, na.rm = TRUE)) > 0,
      (Pseudotime - min(Pseudotime, na.rm = TRUE)) /
        diff(range(Pseudotime, na.rm = TRUE)),
      NA_real_
    )
  ) %>%
  ungroup()
p_method_compare <- ggplot(
  method_long,
  aes(UMAP1, UMAP2, color = Pseudotime01)
) +
  geom_point(size = 3.4) +
  facet_wrap(~Method, nrow = 1) +
  scale_color_gradientn(colours = traditional_colors) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "Scaled\npseudotime",
    title = "Slingshot, TSCAN and SCORPIUS trajectories",
    subtitle = "All methods plotted on the same UMAP; Day0–Day3 groups retained"
  ) +
  theme_bw(11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30")
  )
save_plot(p_method_compare, "09c_trajectory_methods_comparison.png", 14, 5.5)

p_tscan_day <- ggplot(
  data.frame(Day = method_df$Day, Pseudotime = method_df$TSCAN),
  aes(Day, Pseudotime, color = Day)
) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.22) +
  geom_quasirandom(width = 0.14, size = 2.4) +
  scale_color_manual(values = day_colors, guide = "none") +
  labs(
    x = NULL,
    y = "TSCAN pseudotime",
    title = "TSCAN pseudotime by Day",
    subtitle = sprintf("Spearman with Day = %.3f", tscan_day_cor)
  ) +
  theme_f4(12) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 9))
save_plot(p_tscan_day, "09d_TSCAN_pseudotime_by_Day.png", 6.5, 5)

message("[10/10] Saving selection table, summaries and session info")
summary_lines <- c(
  "Soybean microseed simulation: selected n=5 per Day",
  "Gene universe: protein-coding GTF models only",
  paste0("Protein-coding models: ", nrow(symbol_map)),
  "Selected runs: 20; excluded runs: 14",
  "Day0 selected: 5; Day1 selected: 5; Day2 selected: 5; Day3 selected: 5",
  "Selected DESeq2 model: ~ Batch + Day",
  "Pseudotime: UMAP representation; PCA coordinates were not used",
  "Pseudotime method: Slingshot on UMAP with Day-labeled clusters; start Day0",
  "Pseudotime values: branch-aware consensus after orienting each lineage toward later Days",
  "Cross-method check: Slingshot, TSCAN and SCORPIUS",
  "tradeSeq: top 500 variable genes; nknots=4; associationTest",
  "Excluded-pool role: mock bulk-like QC/correlation/clustering only",
  "All figure labels use Day0-Day3; gene labels use GTF gene_name first"
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
writeLines(as.character(packageVersion("DESeq2")), file.path(output_dir, "DESeq2_version.txt"))
writeLines(
  c(
    "Data type: bulk RNA-seq",
    "Species: Glycine max",
    "Selected design: Day0-Day3, n=5 per day; Batch retained as covariate",
    "Excluded pool: 14 runs reserved as bulk-like simulation/QC set"
  ),
  file.path(output_dir, "preflight_confirmation.txt")
)

message("Completed: ", output_dir)
