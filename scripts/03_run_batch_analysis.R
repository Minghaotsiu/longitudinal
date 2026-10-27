#!/usr/bin/env Rscript

# Soybean germination time-course bulk RNA-seq downstream analysis.
# The DESeq2, volcano and Mfuzz structure follows the user's F4 scripts,
# with organism, sample grouping, thresholds and output paths adapted.

options(stringsAsFactors = FALSE)
set.seed(123)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  library(patchwork)
  library(Mfuzz)
  library(Biobase)
  library(uwot)
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

input_dir <- normalizePath(
  Sys.getenv("BATCH_UPSTREAM_DIR", "data/doudou_v2/upstream_run"),
  winslash = "/",
  mustWork = TRUE
)
output_dir <- normalizePath(
  Sys.getenv("BATCH_OUTPUT_DIR", "outputs/doudou_v2_downstream"),
  winslash = "/",
  mustWork = FALSE
)
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

palette_time <- c(
  d0 = "#7895C1",
  d1 = "#7BC0CD",
  d2 = "#F0C284",
  d3 = "#ED8D5A"
)

timeline_levels <- c("d0", "d1", "d2", "d3")
comparisons <- c(
  d1_vs_d0 = "d1",
  d2_vs_d0 = "d2",
  d3_vs_d0 = "d3"
)

theme_manuscript <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(color = "black", face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      plot.subtitle = element_text(color = "grey25"),
      legend.title = element_text(face = "bold")
    )
}

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(
    filename = file.path(figure_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

message("[1/10] Reading count matrix and metadata")
count_file <- file.path(input_dir, "ALL_17_raw_counts_matrix.csv")
map_file <- file.path(input_dir, "ALL_17_RNA_MAP.csv")
stopifnot(file.exists(count_file), file.exists(map_file))

raw_counts <- read.csv(count_file, check.names = FALSE)
if (!identical(colnames(raw_counts)[1], "Gene_ID")) {
  stop("The first column of the count matrix is not Gene_ID.")
}

gene_ids <- trimws(raw_counts$Gene_ID)
if (anyDuplicated(gene_ids)) {
  stop("Duplicated Gene_ID values are present.")
}

count_mat <- as.matrix(raw_counts[, -1, drop = FALSE])
suppressWarnings(storage.mode(count_mat) <- "double")
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)
rownames(count_mat) <- gene_ids

sample_annotation <- read.csv(map_file, check.names = FALSE)
sample_ids <- colnames(count_mat)
if (!all(sample_ids %in% sample_annotation$Sample)) {
  stop("Sample metadata does not cover all matrix columns.")
}

meta <- sample_annotation[
  match(sample_ids, sample_annotation$Sample),
  c("Sample", "Filename"),
  drop = FALSE
]
meta$timepoint <- sub(".*-d([0-3])-.*", "d\\1", meta$Sample)
if (any(!meta$timepoint %in% timeline_levels)) {
  stop("Unexpected sample naming; timepoint could not be parsed.")
}
meta$timepoint <- factor(meta$timepoint, levels = timeline_levels)
rownames(meta) <- meta$Sample
write.csv(meta, file.path(table_dir, "sample_metadata.csv"), row.names = FALSE)

message("[2/10] Filtering low-count genes and running DESeq2")
# Retain genes with >=10 reads in at least the smallest biological group (n=3).
keep <- rowSums(count_mat >= 10) >= min(table(meta$timepoint))
filtered_counts <- count_mat[keep, , drop = FALSE]

if (nrow(filtered_counts) < 100) {
  stop("Too few genes remain after filtering.")
}

dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = meta,
  design = ~ timepoint
)
dds <- DESeq(dds, parallel = FALSE)

normalized_counts <- counts(dds, normalized = TRUE)
vst_mat <- assay(vst(dds, blind = FALSE))

write.csv(
  data.frame(Gene_ID = rownames(normalized_counts), normalized_counts, check.names = FALSE),
  file.path(table_dir, "DESeq2_normalized_counts.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Gene_ID = rownames(vst_mat), vst_mat, check.names = FALSE),
  file.path(table_dir, "VST_expression.csv"),
  row.names = FALSE
)

message("[3/10] QC summaries and figures")
qc_all <- data.frame(
  Sample = colnames(count_mat),
  timepoint = meta$timepoint,
  library_size = colSums(count_mat),
  genes_detected_gt0 = colSums(count_mat > 0),
  genes_detected_ge10 = colSums(count_mat >= 10),
  stringsAsFactors = FALSE
)
write.csv(
  qc_all,
  file.path(table_dir, "QC_summary_counts.csv"),
  row.names = FALSE
)

qc_long <- qc_all %>%
  select(Sample, timepoint, library_size, genes_detected_ge10) %>%
  pivot_longer(
    cols = c(library_size, genes_detected_ge10),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("library_size", "genes_detected_ge10"),
      labels = c("Library size (reads)", "Genes with counts >= 10")
    )
  )

p_qc <- ggplot(qc_long, aes(x = timepoint, y = value, fill = timepoint)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.12, size = 2.2, shape = 21, fill = "white", stroke = 0.8) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = palette_time, guide = "none") +
  labs(
    x = "Germination timepoint",
    y = NULL,
    title = "Sample-level RNA-seq yield and gene detection",
    subtitle = "Descriptive statistics; box = median and IQR; points = individual samples"
  ) +
  theme_manuscript()
save_plot(p_qc, "01_QC_counts.png", width = 9, height = 4.8)

message("[4/10] PCA, UMAP and sample correlation")
gene_variances <- apply(vst_mat, 1, var)
top_variable_n <- min(1000, sum(is.finite(gene_variances) & gene_variances > 0))
top_variable <- names(sort(gene_variances, decreasing = TRUE))[seq_len(top_variable_n)]

pca_fit <- prcomp(t(vst_mat[top_variable, , drop = FALSE]), scale. = FALSE)
percent_var <- 100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2))
pca_df <- data.frame(
  Sample = rownames(pca_fit$x),
  timepoint = meta[rownames(pca_fit$x), "timepoint"],
  PC1 = pca_fit$x[, 1],
  PC2 = pca_fit$x[, 2]
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = timepoint, label = Sample)) +
  geom_point(size = 4.2, alpha = 0.95) +
  geom_text_repel(size = 3.1, max.overlaps = 30, show.legend = FALSE) +
  scale_color_manual(values = palette_time) +
  labs(
    x = sprintf("PC1 (%.1f%%)", percent_var[1]),
    y = sprintf("PC2 (%.1f%%)", percent_var[2]),
    color = "Timepoint",
    title = "PCA of DESeq2 VST-normalized expression",
    subtitle = sprintf("Top %d variable genes; unsupervised, descriptive only", top_variable_n)
  ) +
  theme_manuscript()
save_plot(p_pca, "02_PCA.png", width = 7.2, height = 5.6)

set.seed(123)
pc_input <- pca_fit$x[, seq_len(min(10, ncol(pca_fit$x))), drop = FALSE]
umap_fit <- uwot::umap(
  pc_input,
  n_neighbors = 5,
  min_dist = 0.25,
  metric = "cosine",
  n_components = 2,
  verbose = FALSE
)
colnames(umap_fit) <- c("UMAP1", "UMAP2")
umap_df <- data.frame(
  Sample = rownames(pca_fit$x),
  timepoint = meta[rownames(pca_fit$x), "timepoint"],
  umap_fit
)

p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = timepoint)) +
  geom_point(size = 4.5, alpha = 0.95) +
  scale_color_manual(values = palette_time) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "Timepoint",
    title = "Exploratory UMAP of sample expression profiles",
    subtitle = "17 samples; 5 nearest neighbours; unsupervised, no significance test"
  ) +
  theme_manuscript()
save_plot(p_umap, "03_UMAP.png", width = 6.2, height = 5.2)

cor_mat <- cor(vst_mat[top_variable, , drop = FALSE], method = "pearson")
annotation_col <- data.frame(
  Timepoint = meta[colnames(cor_mat), "timepoint"],
  row.names = colnames(cor_mat)
)
annotation_colors <- list(Timepoint = palette_time)

png(file.path(figure_dir, "04_sample_correlation.png"),
    width = 2200, height = 1900, res = 300)
pheatmap(
  cor_mat,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("#7895C1", "white", "#ED8D5A"))(100),
  border_color = "white",
  fontsize = 8,
  main = "Pearson correlation of VST-normalized expression"
)
dev.off()

message("[5/10] Differential expression for day1/day2/day3 versus day0")
de_results <- list()
de_summary <- list()
significant_union <- character(0)

for (contrast_name in names(comparisons)) {
  target <- comparisons[[contrast_name]]
  res <- results(
    dds,
    contrast = c("timepoint", target, "d0"),
    alpha = 0.05,
    independentFiltering = TRUE
  )
  res_df <- as.data.frame(res)
  res_df$Gene_ID <- rownames(res_df)
  res_df <- res_df[, c(
    "Gene_ID", "baseMean", "log2FoldChange", "lfcSE", "stat",
    "pvalue", "padj"
  )]
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

  write.csv(
    res_df,
    file.path(table_dir, paste0("DE_", contrast_name, "_all.csv")),
    row.names = FALSE
  )
  sig_df <- res_df[res_df$direction != "ns", ]
  write.csv(
    sig_df,
    file.path(table_dir, paste0("DE_", contrast_name, "_significant.csv")),
    row.names = FALSE
  )
  de_results[[contrast_name]] <- res_df
  significant_union <- union(significant_union, sig_df$Gene_ID)

  de_summary[[contrast_name]] <- data.frame(
    comparison = contrast_name,
    up = sum(res_df$direction == "up"),
    down = sum(res_df$direction == "down"),
    tested = sum(!is.na(res_df$padj)),
    stringsAsFactors = FALSE
  )

  plot_df <- res_df %>%
    filter(!is.na(pvalue), is.finite(pvalue))
  plot_df$neglog10p <- -log10(pmax(plot_df$pvalue, 1e-50))
  label_df <- plot_df %>%
    filter(direction != "ns") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 10)

  p_volcano <- ggplot(plot_df, aes(log2FoldChange, neglog10p, color = direction)) +
    geom_point(size = 1.35, alpha = 0.65) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = 2, linewidth = 0.45) +
    geom_hline(yintercept = -log10(0.05), linetype = 4, linewidth = 0.45) +
    geom_text_repel(
      data = label_df,
      aes(label = Gene_ID),
      size = 2.7,
      max.overlaps = 25,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(down = "#7BC0CD", ns = "grey72", up = "#ED8D5A"),
      labels = c(down = "Down", ns = "Not significant", up = "Up")
    ) +
    labs(
      x = expression(log[2]~fold~change),
      y = expression(-log[10](italic(p)~value)),
      color = NULL,
      title = paste0(target, " versus d0"),
      subtitle = "DESeq2 Wald test; BH padj < 0.05 and |log2FC| >= 0.3; p capped at 1e-50 for display"
    ) +
    theme_manuscript()
  save_plot(
    p_volcano,
    paste0("05_volcano_", contrast_name, ".png"),
    width = 7.2,
    height = 5.8
  )
}

de_summary_df <- bind_rows(de_summary)
write.csv(
  de_summary_df,
  file.path(table_dir, "DE_summary.csv"),
  row.names = FALSE
)

p_de_counts <- de_summary_df %>%
  pivot_longer(c(up, down), names_to = "direction", values_to = "count") %>%
  mutate(
    comparison = factor(comparison, levels = names(comparisons)),
    direction = factor(direction, levels = c("down", "up"))
  ) %>%
  ggplot(aes(x = comparison, y = count, fill = direction)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.72),
    vjust = -0.3,
    size = 3.1
  ) +
  scale_fill_manual(values = c(down = "#7BC0CD", up = "#ED8D5A")) +
  labs(
    x = NULL,
    y = "Number of differentially expressed genes",
    fill = "Direction",
    title = "Differential expression across germination timepoints",
    subtitle = "DESeq2 Wald test, BH adjustment, padj < 0.05 and |log2FC| >= 0.3"
  ) +
  theme_manuscript()
save_plot(p_de_counts, "06_DE_counts.png", width = 7.2, height = 4.8)

message("[6/10] Heatmap and key dynamic genes")
if (length(significant_union) > 0) {
  union_table <- bind_rows(de_results)
  gene_score <- union_table %>%
    filter(Gene_ID %in% significant_union) %>%
    group_by(Gene_ID) %>%
    summarise(
      min_padj = min(padj, na.rm = TRUE),
      max_abs_lfc = max(abs(log2FoldChange), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(min_padj, desc(max_abs_lfc))
  heat_genes <- head(gene_score$Gene_ID, 50)
} else {
  gene_score <- data.frame(Gene_ID = character())
  heat_genes <- rownames(vst_mat)[order(gene_variances, decreasing = TRUE)][1:40]
}

write.csv(
  gene_score,
  file.path(table_dir, "DE_ranked_gene_union.csv"),
  row.names = FALSE
)

heat_mat <- vst_mat[heat_genes, , drop = FALSE]
heat_mat <- t(scale(t(heat_mat)))
heat_mat[!is.finite(heat_mat)] <- 0
heat_annot <- data.frame(
  Timepoint = meta[colnames(heat_mat), "timepoint"],
  row.names = colnames(heat_mat)
)

png(file.path(figure_dir, "07_top50_heatmap.png"),
    width = 2400, height = 2500, res = 300)
pheatmap(
  heat_mat,
  annotation_col = heat_annot,
  annotation_colors = annotation_colors,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  fontsize_row = 6,
  fontsize_col = 8,
  color = colorRampPalette(c("#7895C1", "white", "#ED8D5A"))(100),
  border_color = NA,
  main = "Top 50 differential genes across germination timepoints",
  angle_col = 45
)
dev.off()

key_genes <- head(gene_score$Gene_ID, 12)
if (length(key_genes) == 0) {
  key_genes <- names(sort(gene_variances, decreasing = TRUE))[1:12]
}

trajectory_long <- data.frame(
  Gene_ID = rep(key_genes, each = ncol(vst_mat)),
  Sample = rep(colnames(vst_mat), times = length(key_genes)),
  timepoint = rep(meta[colnames(vst_mat), "timepoint"], times = length(key_genes)),
  expression = as.vector(t(vst_mat[key_genes, , drop = FALSE]))
)
write.csv(
  trajectory_long,
  file.path(table_dir, "key_gene_sample_expression.csv"),
  row.names = FALSE
)

trajectory_summary <- trajectory_long %>%
  group_by(Gene_ID, timepoint) %>%
  summarise(
    mean_expression = mean(expression),
    sd_expression = sd(expression),
    n = n(),
    .groups = "drop"
  )

p_key <- ggplot(trajectory_long, aes(timepoint, expression)) +
  geom_point(
    aes(color = timepoint),
    position = position_jitter(width = 0.09, height = 0),
    size = 1.8,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  geom_line(
    data = trajectory_summary,
    aes(timepoint, mean_expression, group = 1),
    color = "black",
    linewidth = 0.8
  ) +
  geom_point(
    data = trajectory_summary,
    aes(timepoint, mean_expression),
    color = "black",
    size = 2.7
  ) +
  facet_wrap(~Gene_ID, scales = "free_y", ncol = 4) +
  scale_color_manual(values = palette_time) +
  labs(
    x = "Germination timepoint",
    y = "VST-normalized expression",
    title = "Expression trajectories of leading differential genes",
    subtitle = "Black line = group mean; points = individual samples; descriptive only"
  ) +
  theme_manuscript(base_size = 9)
save_plot(p_key, "08_key_gene_trajectories.png", width = 11, height = 7.5)

message("[7/10] Mfuzz temporal clustering")
time_mean_vst <- sapply(
  timeline_levels,
  function(tp) rowMeans(vst_mat[, meta$timepoint == tp, drop = FALSE])
)
colnames(time_mean_vst) <- timeline_levels
mean_variance <- apply(time_mean_vst, 1, var)
mfuzz_n <- min(2000, sum(is.finite(mean_variance) & mean_variance > 0))
mfuzz_genes <- names(sort(mean_variance, decreasing = TRUE))[seq_len(mfuzz_n)]
mfuzz_mat <- time_mean_vst[mfuzz_genes, , drop = FALSE]

mfuzz_object <- new("ExpressionSet", exprs = mfuzz_mat)
mfuzz_object <- filter.NA(mfuzz_object, thres = 0.25)
mfuzz_object <- fill.NA(mfuzz_object, mode = "mean")
mfuzz_object <- filter.std(mfuzz_object, min.std = 0.1, visu = FALSE)
mfuzz_object <- standardise(mfuzz_object)

cluster_number <- 6
m_value <- mestimate(mfuzz_object)
mfuzz_result <- mfuzz(
  mfuzz_object,
  c = cluster_number,
  m = m_value
)

mfuzz_assignments <- data.frame(
  Gene_ID = names(mfuzz_result$cluster),
  cluster = as.integer(mfuzz_result$cluster),
  max_membership = apply(mfuzz_result$membership, 1, max),
  check.names = FALSE
)
write.csv(
  mfuzz_assignments,
  file.path(table_dir, "Mfuzz_gene_clusters.csv"),
  row.names = FALSE
)

cluster_sizes <- as.data.frame(table(cluster = mfuzz_assignments$cluster))
write.csv(
  cluster_sizes,
  file.path(table_dir, "Mfuzz_cluster_sizes.csv"),
  row.names = FALSE
)

standardized_expr <- exprs(mfuzz_object)
mfuzz_long <- as.data.frame(standardized_expr) %>%
  rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, names_to = "timepoint", values_to = "z_expression") %>%
  left_join(mfuzz_assignments, by = "Gene_ID") %>%
  mutate(
    timepoint = factor(timepoint, levels = timeline_levels),
    cluster = factor(cluster, levels = seq_len(cluster_number))
  )

set.seed(123)
line_genes <- mfuzz_long %>%
  distinct(Gene_ID) %>%
  slice_sample(n = min(500, nrow(.)))
mfuzz_lines <- mfuzz_long %>% filter(Gene_ID %in% line_genes$Gene_ID)

centers <- as.data.frame(mfuzz_result$centers) %>%
  rownames_to_column("cluster") %>%
  pivot_longer(-cluster, names_to = "timepoint", values_to = "center") %>%
  mutate(
    cluster = factor(as.integer(cluster), levels = seq_len(cluster_number)),
    timepoint = factor(timepoint, levels = timeline_levels)
  )

p_mfuzz <- ggplot(mfuzz_lines, aes(timepoint, z_expression)) +
  geom_line(aes(group = Gene_ID), color = "#7BC0CD", alpha = 0.08, linewidth = 0.35) +
  geom_line(
    data = centers,
    aes(timepoint, center, group = 1),
    color = "#ED8D5A",
    linewidth = 1.15,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = centers,
    aes(timepoint, center),
    color = "#B84D28",
    size = 2.2,
    inherit.aes = FALSE
  ) +
  facet_wrap(~cluster, scales = "fixed", ncol = 3) +
  labs(
    x = "Germination timepoint",
    y = "Standardized expression",
    title = "Mfuzz clustering of dynamic expression trajectories",
    subtitle = paste0(
      "Fuzzy c-means, c=6, m=", round(m_value, 2),
      "; up to ", mfuzz_n, " variable genes; no p-value"
    )
  ) +
  theme_manuscript(base_size = 10)
save_plot(p_mfuzz, "09_Mfuzz_trajectories.png", width = 9.2, height = 6.5)

message("[8/10] Functional enrichment using g:Profiler (Glycine max, gmax)")
all_background_genes <- rownames(filtered_counts)
enrichment_list <- list()

gprofiler_enrich <- function(genes, label, direction) {
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
  resp <- tryCatch(
    httr::POST(
      "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
      body = jsonlite::toJSON(body, auto_unbox = TRUE),
      httr::content_type_json(),
      httr::timeout(120)
    ),
    error = function(e) e
  )
  if (inherits(resp, "error")) {
    warning(label, ": g:Profiler request failed: ", resp$message)
    return(NULL)
  }
  if (httr::status_code(resp) != 200) {
    warning(label, ": g:Profiler returned HTTP ", httr::status_code(resp))
    return(NULL)
  }
  parsed <- jsonlite::fromJSON(
    httr::content(resp, as = "text", encoding = "UTF-8"),
    simplifyDataFrame = TRUE
  )
  result <- parsed$result
  if (is.null(result) || length(result) == 0 || !is.data.frame(result) || nrow(result) == 0) {
    return(NULL)
  }
  result$contrast <- label
  result$direction <- direction
  result
}

for (contrast_name in names(de_results)) {
  for (direction_label in c("up", "down")) {
    genes <- de_results[[contrast_name]] %>%
      filter(direction == direction_label) %>%
      pull(Gene_ID)
    message(
      "  g:Profiler ", contrast_name, " ", direction_label, ": ", length(genes), " genes"
    )
    result <- gprofiler_enrich(
      genes,
      contrast_name,
      direction_label
    )
    if (!is.null(result)) {
      enrichment_list[[paste(contrast_name, direction_label, sep = "_")]] <- result
    }
  }
}

enrichment_df <- if (length(enrichment_list) > 0) {
  bind_rows(enrichment_list)
} else {
  data.frame()
}

# g:Profiler returns nested fields such as intersections and parents. Keep all
# information, but serialize nested list columns so the table can be written to CSV.
if (nrow(enrichment_df) > 0) {
  enrichment_df[] <- lapply(
    enrichment_df,
    function(x) {
      if (is.list(x)) {
        vapply(
          x,
          function(value) paste(value, collapse = ";"),
          character(1)
        )
      } else {
        x
      }
    }
  )
}

write.csv(
  enrichment_df,
  file.path(table_dir, "GO_KEGG_enrichment_gprofiler.csv"),
  row.names = FALSE
)

if (nrow(enrichment_df) > 0 && "p_value" %in% colnames(enrichment_df)) {
  enrich_plot_df <- enrichment_df %>%
    mutate(
      comparison_direction = paste(contrast, direction, sep = " "),
      neglog10p = -log10(pmax(p_value, 1e-50))
    ) %>%
    group_by(comparison_direction) %>%
    arrange(p_value, desc(intersection_size)) %>%
    slice_head(n = 10) %>%
    ungroup()

  p_enrich <- ggplot(
    enrich_plot_df,
    aes(
      x = neglog10p,
      y = reorder(paste(source, name, sep = " | "), neglog10p),
      size = intersection_size,
      color = source
    )
  ) +
    geom_point(alpha = 0.85) +
    facet_wrap(~comparison_direction, scales = "free_y", ncol = 2) +
    scale_color_manual(
      values = c(
        `GO:BP` = "#7895C1",
        `GO:CC` = "#7BC0CD",
        `GO:MF` = "#F0C284",
        KEGG = "#ED8D5A"
      )
    ) +
    labs(
      x = expression(-log[10](italic(p))),
      y = NULL,
      size = "Genes",
      color = "Source",
      title = "Functional enrichment of differential genes",
      subtitle = "g:Profiler hypergeometric test; g:SCS correction; p < 0.05"
    ) +
    theme_manuscript(base_size = 9)
  save_plot(p_enrich, "10_GO_KEGG_enrichment.png", width = 12, height = 8)
}

message("[9/10] Writing figure statistics table")
figure_stats <- data.frame(
  figure = c(
    "01_QC_counts.png",
    "02_PCA.png",
    "03_UMAP.png",
    "04_sample_correlation.png",
    paste0("05_volcano_", names(comparisons), ".png"),
    "06_DE_counts.png",
    "07_top50_heatmap.png",
    "08_key_gene_trajectories.png",
    "09_Mfuzz_trajectories.png",
    if (exists("p_enrich")) "10_GO_KEGG_enrichment.png" else character(0)
  ),
  n = c(
    "17 samples; n=3/5/5/4 by d0-d3",
    "17 samples; n=3/5/5/4",
    "17 samples; n=3/5/5/4",
    "17 samples; n=3/5/5/4",
    rep("d0 n=3 vs target n=5/5/4", 3),
    "17 samples",
    "17 samples",
    "17 samples",
    paste0("4 timepoint means; up to ", mfuzz_n, " genes"),
    if (exists("p_enrich")) "DEGs per contrast/direction" else character(0)
  ),
  test = c(
    "Descriptive; no group test",
    "Unsupervised PCA; no significance test",
    "Unsupervised UMAP; no significance test",
    "Pearson correlation; descriptive",
    rep("DESeq2 Wald test; BH correction", 3),
    "DESeq2 Wald test; BH correction",
    "Unsupervised clustering; descriptive",
    "Group mean and SD; descriptive",
    "Fuzzy c-means; no p-value",
    if (exists("p_enrich")) "g:Profiler hypergeometric; g:SCS correction" else character(0)
  ),
  result_or_threshold = c(
    "Counts >=10; library sizes and detected genes shown",
    paste0("Top ", top_variable_n, " variable genes"),
    "5 nearest neighbours; exploratory",
    paste0("Top ", top_variable_n, " variable genes"),
    rep("padj<0.05 and |log2FC|>=0.3", 3),
    "padj<0.05 and |log2FC|>=0.3",
    "Top 50 differential or variable genes",
    "Top 12 ranked genes",
    "c=6; membership-based assignment",
    if (exists("p_enrich")) "p<0.05; top 10 terms per contrast/direction" else character(0)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  figure_stats,
  file.path(table_dir, "figure_statistics.csv"),
  row.names = FALSE
)

message("[10/10] Saving analysis summary and session information")
summary_lines <- c(
  "Soybean germination time-course bulk RNA-seq",
  paste0("Samples: ", ncol(count_mat)),
  paste0("Genes in raw matrix: ", nrow(count_mat)),
  paste0("Genes after filter: ", nrow(filtered_counts)),
  paste0("Filter: rowSums(counts >= 10) >= ", min(table(meta$timepoint))),
  paste0("Group sizes: ", paste(names(table(meta$timepoint)), table(meta$timepoint), sep = "=", collapse = ", ")),
  "Comparisons: d1 vs d0, d2 vs d0, d3 vs d0",
  "DE model: DESeq2, design ~ timepoint, Wald test, BH correction",
  "DE threshold: padj < 0.05 and |log2FoldChange| >= 0.3",
  "Enrichment: g:Profiler gmax, GO BP/CC/MF and KEGG, g:SCS correction",
  "Interpretation: exploratory due to uneven sequencing depth and incomplete genome alignment"
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))
writeLines(as.character(packageVersion("DESeq2")), file.path(output_dir, "DESeq2_version.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

message("Completed: ", output_dir)
