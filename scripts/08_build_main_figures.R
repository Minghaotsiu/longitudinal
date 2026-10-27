#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(123)

workspace_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", "."),
  winslash = "/",
  mustWork = TRUE
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(slingshot)
  library(SingleCellExperiment)
  library(uwot)
  library(Rtsne)
  library(scales)
})

root <- workspace_root
data_dir <- normalizePath(
  Sys.getenv("PROCESSED_DATA_DIR", file.path(root, "data", "doudou_all34")),
  winslash = "/",
  mustWork = TRUE
)
all34_downstream_dir <- normalizePath(
  Sys.getenv("ALL34_OUTPUT_DIR", file.path(root, "outputs", "doudou_all34_downstream")),
  winslash = "/",
  mustWork = TRUE
)
selected5_dir <- normalizePath(
  Sys.getenv("SELECTED5_OUTPUT_DIR", file.path(root, "outputs", "doudou_selected5_v2")),
  winslash = "/",
  mustWork = TRUE
)
pseudotime_dir <- normalizePath(
  Sys.getenv(
    "PSEUDOTIME_OUTPUT_DIR",
    file.path(root, "outputs", "doudou_selected5_pseudotime")
  ),
  winslash = "/",
  mustWork = TRUE
)
out_dir <- normalizePath(
  Sys.getenv("SCIENTIFIC_DATA_OUTPUT_DIR", file.path(root, "outputs", "doudou_scientific_data")),
  winslash = "/",
  mustWork = FALSE
)
fig_dir <- file.path(out_dir, "figures")
source_dir <- file.path(out_dir, "source_data")
qa_dir <- file.path(out_dir, "qa")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

alignment_script <- Sys.getenv("PANEL_ALIGN_AUDIT", "")
alignment_helper <- file.path(root, "scripts", "lib", "panel_alignment_compat.R")
audit_python <- Sys.getenv("PYTHON_BIN", Sys.which("python3"))
source(alignment_helper)

analysis_colors <- c(MN = "#356C82", Bulk = "#B46A4F")
day_colors <- c(
  Day0 = "#66758A",
  Day1 = "#5F8FA6",
  Day2 = "#78A58F",
  Day3 = "#C57C5A"
)
seed_colors <- c(
  Seed1 = "#2F657A",
  Seed2 = "#6B8E6E",
  Seed3 = "#C08A48",
  Seed4 = "#A76358",
  Seed5 = "#746A9A"
)

theme_descriptor <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#272727"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#272727"),
      axis.title = element_text(size = base_size, colour = "#202020"),
      axis.text = element_text(size = base_size - 0.4, colour = "#303030"),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.6),
      legend.key.size = unit(3.4, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.1, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, colour = "#555555", hjust = 0),
      plot.tag = element_text(size = 8, face = "bold"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(3, 5, 3, 4)
    )
}
theme_set(theme_descriptor())

save_plot <- function(
    plot,
    stem,
    width_mm,
    height_mm,
    panel_ids = NULL,
    row_groups = NULL,
    column_groups = NULL,
    exemptions = list(),
    dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  if (inherits(plot, "patchwork")) {
    require_patchwork_panel_alignment(
      plot = plot,
      manifest_path = file.path(qa_dir, paste0(stem, ".alignment-layout.json")),
      report_path = file.path(qa_dir, paste0(stem, ".alignment.json")),
      overlay_svg = file.path(qa_dir, paste0(stem, ".alignment.svg")),
      audit_script = alignment_script,
      python = audit_python,
      panel_ids = panel_ids,
      row_groups = row_groups,
      column_groups = column_groups,
      exemptions = exemptions,
      width_in = width_in,
      height_in = height_in,
      tolerance_pt = 1.5,
      gutter_tolerance_pt = 1.5,
      strict = TRUE
    )
  }

  svglite::svglite(
    file.path(fig_dir, paste0(stem, ".svg")),
    width = width_in,
    height = height_in,
    bg = "white"
  )
  print(plot)
  dev.off()

  grDevices::cairo_pdf(
    file.path(fig_dir, paste0(stem, ".pdf")),
    width = width_in,
    height = height_in,
    family = "sans",
    bg = "white"
  )
  print(plot)
  dev.off()

  ragg::agg_tiff(
    file.path(fig_dir, paste0(stem, ".tiff")),
    width = width_in,
    height = height_in,
    units = "in",
    res = dpi,
    background = "white"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    file.path(fig_dir, paste0(stem, ".png")),
    width = width_in,
    height = height_in,
    units = "in",
    res = 300,
    background = "white"
  )
  print(plot)
  dev.off()
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
           formatC(p, format = "f", digits = 3))
  )
}

crosswalk <- read.csv(
  file.path(root, "metadata", "scientific_data_sample_crosswalk.csv"),
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
) %>%
  mutate(
    Day = factor(Day, levels = c("Day0", "Day1", "Day2", "Day3")),
    SamplingMethod = factor(SamplingMethod, levels = c("MN", "Bulk")),
    AnalysisRunID = sprintf("Run%02d", row_number())
  )

qc_raw <- read.csv(
  file.path(data_dir, "ALL_34_QC_metrics.csv"),
  check.names = FALSE
) %>%
  transmute(
    OriginalRunSample = RunSample,
    InputReads,
    CleanReads,
    CleanReadRetentionPct,
    MeanReadQuality,
    TotalCounts,
    nGeneGe10,
    nGene1M
  )
q20 <- read.delim(
  file.path(data_dir, "ALL_34_Q20_Q30_readlevel.tsv"),
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
) %>%
  transmute(
    OriginalRunSample = RunSample,
    ReadQ20Pct = Q20ReadPct,
    ReadQ30Pct = Q30ReadPct
  )
star <- read.csv(
  file.path(data_dir, "ALL_34_STAR_metrics_with_source.csv"),
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
) %>%
  dplyr::select(
    OriginalRunSample,
    UniquePct,
    MultiPct,
    TooShortPct,
    OtherUnmappedPct,
    MismatchPct,
    AvgReadLength,
    AvgMappedLength
  )

raw_counts <- read.csv(
  file.path(data_dir, "ALL_34_raw_counts_matrix.csv"),
  check.names = FALSE
)
raw_matrix <- as.matrix(raw_counts[, -1, drop = FALSE])
storage.mode(raw_matrix) <- "double"
rownames(raw_matrix) <- raw_counts$Gene_ID

top100_metrics <- data.frame(
  OriginalRunSample = colnames(raw_matrix),
  Top100Pct = apply(raw_matrix, 2, function(values) {
    ordered <- sort(values, decreasing = TRUE)
    total <- sum(ordered)
    if (total == 0) NA_real_ else 100 * sum(head(ordered, 100)) / total
  })
)

qc_all <- crosswalk %>%
  left_join(qc_raw, by = "OriginalRunSample") %>%
  left_join(q20, by = "OriginalRunSample") %>%
  left_join(star, by = "OriginalRunSample") %>%
  left_join(top100_metrics, by = "OriginalRunSample") %>%
  mutate(
    InputReadsM = InputReads / 1e6,
    CleanReadsM = CleanReads / 1e6,
    TotalCountsM = TotalCounts / 1e6
  )

source_columns <- c(
  "AnalysisRunID", "CorrectedSampleID", "SamplingMethod", "SeedID", "Day",
  "InputReadsM", "CleanReadsM", "CleanReadRetentionPct",
  "ReadQ20Pct", "ReadQ30Pct", "MeanReadQuality",
  "UniquePct", "MultiPct", "TooShortPct",
  "TotalCountsM", "nGeneGe10", "nGene1M", "Top100Pct"
)
write.csv(
  qc_all[, source_columns],
  file.path(source_dir, "Figure2_technical_QC_source_data.csv"),
  row.names = FALSE
)

# Figure 1: data-acquisition workflow
day_x <- c(Day0 = 9, Day1 = 31, Day2 = 53, Day3 = 75)
mn_design <- crosswalk %>%
  filter(SamplingMethod == "MN") %>%
  mutate(
    x = unname(day_x[as.character(Day)]),
    y = 82 - 5.2 * (as.integer(sub("Seed", "", SeedID)) - 1L),
    Label = "MN"
  )
if (nrow(mn_design) != 20 || any(table(mn_design$SeedID, mn_design$Day) != 1)) {
  stop("MN crosswalk must contain one sample for every Seed x Day cell.")
}

bulk_day_counts <- crosswalk %>%
  filter(SamplingMethod == "Bulk") %>%
  dplyr::count(Day, name = "n") %>%
  complete(
    Day = factor(c("Day0", "Day1", "Day2", "Day3"),
                 levels = c("Day0", "Day1", "Day2", "Day3")),
    fill = list(n = 0)
  ) %>%
  mutate(
    x = unname(day_x[as.character(Day)]),
    y = 45,
    Label = paste0(Day, "\nBulk")
  )

workflow_nodes <- data.frame(
  x = c(12, 35, 58, 12, 35, 58, 82, 82),
  y = c(31, 31, 31, 19, 19, 19, 25, 12),
  label = c(
    "Microneedle\nsample",
    "Dissolve needle\nrelease material",
    "RNA extraction\nlibrary preparation",
    "Bulk seed tissue\n(total n=14)",
    "Freeze, grind\nTRIzol RNA",
    "RNA extraction\nlibrary preparation",
    "QC, alignment\ncount matrix",
    "Seed-by-Day metadata\nlongitudinal matrix"
  ),
  group = c(rep("MN processing", 3), rep("Bulk processing", 3),
            rep("Shared data products", 2))
)

workflow_edges <- data.frame(
  x = c(17.5, 40.5, 17.5, 40.5, 63.5, 63.5, 82),
  y = c(31, 31, 19, 19, 31, 19, 21.5),
  xend = c(29.5, 52.5, 29.5, 52.5, 75.5, 75.5, 82),
  yend = c(31, 31, 19, 19, 25, 25, 15.5)
)

figure1 <- ggplot() +
  annotate(
    "text", x = 2, y = 97, label = "Longitudinal same-seed sampling",
    hjust = 0, size = 2.8, fontface = "bold"
  ) +
  annotate(
    "segment", x = 6, xend = 81, y = 91, yend = 91,
    linewidth = 0.8, colour = "#356C82",
    arrow = arrow(type = "closed", length = unit(1.7, "mm"))
  ) +
  annotate(
    "text", x = unname(day_x), y = 93.3, label = names(day_x),
    size = 2.35, fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2,
    y = 82 - 5.2 * (0:4),
    label = paste0("Seed", 1:5),
    hjust = 0, size = 2.1
  ) +
  geom_segment(
    data = mn_design,
    aes(x = x - 4.0, xend = x + 4.0, y = y, yend = y),
    colour = "#B9CDD6",
    linewidth = 0.55
  ) +
  geom_point(
    data = mn_design,
    aes(x, y),
    size = 3.1,
    shape = 21,
    stroke = 0.35,
    fill = "#DCEAF0",
    colour = "#356C82"
  ) +
  geom_point(
    data = mn_design,
    aes(x, y + 1.7),
    size = 1.25,
    shape = 25,
    fill = "#356C82",
    colour = "#356C82"
  ) +
  annotate(
    "text", x = 84, y = 76,
    label = "same five seeds\nsampled at every Day",
    hjust = 0, size = 2.1, colour = "#356C82"
  ) +
  annotate(
    "text", x = 2, y = 51, label = "Destructive cross-sectional comparator",
    hjust = 0, size = 2.7, fontface = "bold"
  ) +
  geom_segment(
    data = bulk_day_counts,
    aes(x = x - 3.8, xend = x + 3.8, y = y, yend = y),
    colour = "#D8C2B8",
    linewidth = 0.55
  ) +
  geom_label(
    data = bulk_day_counts,
    aes(x, y, label = Label),
    size = 1.95,
    label.size = 0.22,
    fill = "#F5EDE9",
    colour = "#704637",
    lineheight = 0.92
  ) +
  annotate(
    "text", x = 84, y = 45,
    label = "independent seeds\nfor each time point",
    hjust = 0, size = 2.1, colour = "#9B614B"
  ) +
  geom_segment(
    data = workflow_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "#5B5B5B",
    linewidth = 0.45,
    arrow = arrow(type = "closed", length = unit(1.25, "mm"))
  ) +
  geom_label(
    data = workflow_nodes,
    aes(x, y, label = label, fill = group),
    size = 1.95,
    label.size = 0.22,
    colour = "#252525",
    lineheight = 0.95
  ) +
  scale_fill_manual(
    values = c(
      "MN processing" = "#DCEAF0",
      "Bulk processing" = "#F5EDE9",
      "Shared data products" = "#E8ECE5"
    ),
    guide = "none"
  ) +
  annotate(
    "text", x = 2, y = 36.5,
    label = "Sample processing",
    hjust = 0, size = 2.55, fontface = "bold"
  ) +
  annotate(
    "text", x = 2, y = 4.2,
    label = "Data products: FASTQ, QC, count matrices, sample metadata and reusable longitudinal analyses",
    hjust = 0, size = 1.95, colour = "#4D4D4D"
  ) +
  coord_cartesian(xlim = c(0, 102), ylim = c(1, 100), clip = "off") +
  labs(
    title = "Data-acquisition workflow for longitudinal microneedle and destructive bulk transcriptomics",
    subtitle = "Active sampling, material release, library preparation and reusable data products",
    caption = "MN samples follow five individual seeds from Day0 to Day3. Bulk samples are independent destructive collections."
  ) +
  theme_void(base_size = 7, base_family = "sans") +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 6.2, colour = "#555555", hjust = 0),
    plot.caption = element_text(size = 5.6, colour = "#666666", hjust = 0),
    plot.margin = margin(5, 7, 5, 7)
  )

save_plot(figure1, "Figure1_sampling_workflow", 183, 120)
write.csv(
  mn_design %>%
    dplyr::select(CorrectedSampleID, SeedID, Day, x, y),
  file.path(source_dir, "Figure1_MN_design_source_data.csv"),
  row.names = FALSE
)
write.csv(
  bulk_day_counts %>% dplyr::select(Day, n, x, y),
  file.path(source_dir, "Figure1_Bulk_counts_source_data.csv"),
  row.names = FALSE
)
write.csv(
  workflow_nodes,
  file.path(source_dir, "Figure1_workflow_nodes_source_data.csv"),
  row.names = FALSE
)
write.csv(
  workflow_edges,
  file.path(source_dir, "Figure1_workflow_edges_source_data.csv"),
  row.names = FALSE
)

# Replace the schematic drawing with the author-supplied sampling image as the
# upper panel. The lower panel retains only the data-generating workflow and
# omits per-day bulk sample counts.
sampling_image_path <- file.path(root, "inputs", "sampling_flowchart.png")
sampling_image <- png::readPNG(sampling_image_path)
image_ratio <- dim(sampling_image)[1] / dim(sampling_image)[2]

p_sampling_source <- ggplot() +
  annotation_raster(
    sampling_image,
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = image_ratio,
    interpolate = TRUE
  ) +
  coord_fixed(
    xlim = c(0, 1),
    ylim = c(0, image_ratio),
    expand = FALSE
  ) +
  theme_void()

p_workflow <- ggplot() +
  annotate(
    "text",
    x = 2,
    y = 38,
    label = "Sample processing and reusable data products",
    hjust = 0,
    size = 2.8,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 12,
    y = 36,
    label = "Microneedle route",
    size = 2.15,
    colour = "#356C82",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 12,
    y = 24,
    label = "Destructive bulk route",
    size = 2.15,
    colour = "#9B614B",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 82,
    y = 31,
    label = "Shared data products",
    size = 2.15,
    colour = "#536554",
    fontface = "bold"
  ) +
  geom_segment(
    data = workflow_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "#5B5B5B",
    linewidth = 0.45,
    arrow = arrow(type = "closed", length = unit(1.25, "mm"))
  ) +
  geom_label(
    data = workflow_nodes,
    aes(x, y, label = label, fill = group),
    size = 1.95,
    label.size = 0.22,
    colour = "#252525",
    lineheight = 0.95
  ) +
  scale_fill_manual(
    values = c(
      "MN processing" = "#DCEAF0",
      "Bulk processing" = "#F5EDE9",
      "Shared data products" = "#E8ECE5"
    ),
    guide = "none"
  ) +
  annotate(
    "text",
    x = 2,
    y = 1.5,
    label = "RNA-seq outputs include FASTQ, quality-control summaries, expression matrices, sample metadata and reusable analysis tables.",
    hjust = 0,
    size = 1.95,
    colour = "#4D4D4D"
  ) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 40), clip = "off") +
  theme_void(base_size = 7, base_family = "sans") +
  theme(plot.margin = margin(0, 7, 2, 7))

figure1_composite <- wrap_plots(
  p_sampling_source,
  p_workflow,
  ncol = 1,
  heights = c(1.20, 0.80)
) +
  plot_annotation(
    title = "Data-acquisition workflow for longitudinal microneedle and destructive bulk transcriptomics",
    subtitle = "The upper panel is the supplied sampling schematic; the lower panel summarizes material processing and data products",
    caption = "Daily bulk sample counts are intentionally omitted from the figure.",
    tag_levels = "a"
  )
save_plot(
  figure1_composite,
  "Figure1_sampling_workflow",
  183,
  190,
  panel_ids = c("a", "b")
)
write.csv(
  data.frame(
    SourceImage = "inputs/sampling_flowchart.png",
    WidthPixels = dim(sampling_image)[2],
    HeightPixels = dim(sampling_image)[1],
    Note = "Author-supplied sampling schematic used as the upper panel; daily bulk n values are not displayed."
  ),
  file.path(source_dir, "Figure1_source_image_provenance.csv"),
  row.names = FALSE
)

# Figure 2: technical quality and limitations
metric_specs <- tribble(
  ~Variable, ~Label,
  "InputReadsM", "Input reads (M)",
  "CleanReadsM", "Clean reads (M)",
  "ReadQ20Pct", "Read Q20 (%)",
  "ReadQ30Pct", "Read Q30 (%)",
  "CleanReadRetentionPct", "Clean retention (%)",
  "UniquePct", "Unique mapping (%)",
  "MultiPct", "Multi-mapping (%)",
  "TooShortPct", "Unmapped: too short (%)",
  "TotalCountsM", "Assigned counts (M)",
  "nGeneGe10", "Genes with counts >= 10",
  "nGene1M", "Genes at 1M counts",
  "Top100Pct", "Top-100 gene share (%)"
)

anova_p <- function(data, variable) {
  model <- aov(reformulate("SamplingMethod", response = variable), data = data)
  summary(model)[[1]][["Pr(>F)"]][1]
}

metric_specs <- metric_specs %>%
  mutate(PanelLabel = factor(
    Label,
    levels = Label
  ))

qc_long_all <- qc_all %>%
  dplyr::select(
    CorrectedSampleID,
    SamplingMethod,
    all_of(metric_specs$Variable)
  ) %>%
  pivot_longer(
    cols = all_of(metric_specs$Variable),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  left_join(
    metric_specs %>% dplyr::select(Variable, PanelLabel),
    by = "Variable"
  ) %>%
  filter(is.finite(Value))
write.csv(
  qc_long_all,
  file.path(source_dir, "Figure2_technical_QC_long_source_data.csv"),
  row.names = FALSE
)

qc_p_all <- qc_long_all %>%
  group_by(Variable, PanelLabel) %>%
  summarise(PValue = anova_p(pick(everything()), "Value"), .groups = "drop") %>%
  mutate(
    LabelText = paste0("ANOVA P = ", format_p(PValue)),
    Displayed = PValue >= 0.05
  )
write.csv(
  qc_p_all,
  file.path(source_dir, "Figure2_ANOVA_statistics.csv"),
  row.names = FALSE
)
qc_p <- qc_p_all %>%
  filter(Displayed)
qc_long <- qc_long_all %>%
  filter(Variable %in% qc_p$Variable) %>%
  mutate(PanelLabel = factor(PanelLabel, levels = qc_p$PanelLabel))

figure2 <- ggplot(
  qc_long,
  aes(SamplingMethod, Value, colour = SamplingMethod)
) +
  geom_jitter(
    width = 0.085,
    height = 0,
    size = 1.12,
    alpha = 0.72,
    shape = 16
  ) +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "errorbar",
    width = 0.30,
    linewidth = 0.55,
    colour = "#202020"
  ) +
  geom_text(
    data = qc_p,
    aes(x = 1.5, y = Inf, label = LabelText),
    inherit.aes = FALSE,
    vjust = 1.45,
    hjust = 0.5,
    size = 1.85,
    colour = "#303030"
  ) +
  facet_wrap(~PanelLabel, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = analysis_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0.06, 0.18)),
    labels = label_number_auto()
  ) +
  labs(
    x = NULL,
    y = NULL,
    colour = NULL,
    title = "Technical metrics without significant MN-Bulk differences",
    subtitle = "All displayed metrics have one-way ANOVA P >= 0.05; horizontal marks are medians",
    caption = "Metrics with ANOVA P < 0.05 are not shown in this panel. Full results, including significant metrics, are provided in Source Data."
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 5.7)
  )
save_plot(
  figure2,
  "Figure2_technical_QC",
  183,
  190
)

# Figure 3: expression structure without batch encoding
vst <- read.csv(
  file.path(all34_downstream_dir, "tables",
            "VST_expression_34.csv"),
  check.names = FALSE
)
vst_matrix <- as.matrix(vst[, -(1:2), drop = FALSE])
rownames(vst_matrix) <- vst$Gene_ID

gene_var <- apply(vst_matrix, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(1000, sum(gene_var > 0)))]
expr_top <- t(vst_matrix[top_genes, , drop = FALSE])

pca_fit <- prcomp(expr_top, scale. = FALSE)
pca_var <- 100 * pca_fit$sdev^2 / sum(pca_fit$sdev^2)
set.seed(123)
umap_fit <- uwot::umap(
  expr_top,
  n_neighbors = 8,
  min_dist = 0.25,
  metric = "cosine",
  n_components = 2,
  verbose = FALSE
)
set.seed(123)
tsne_fit <- Rtsne::Rtsne(
  expr_top,
  dims = 2,
  perplexity = 7,
  theta = 0,
  check_duplicates = FALSE,
  pca = TRUE,
  max_iter = 1500
)$Y
structure_df <- bind_rows(
  data.frame(
    OriginalRunSample = rownames(pca_fit$x),
    Method = "PCA",
    Axis1 = pca_fit$x[, 1],
    Axis2 = pca_fit$x[, 2]
  ),
  data.frame(
    OriginalRunSample = rownames(expr_top),
    Method = "UMAP",
    Axis1 = umap_fit[, 1],
    Axis2 = umap_fit[, 2]
  ),
  data.frame(
    OriginalRunSample = rownames(expr_top),
    Method = "t-SNE",
    Axis1 = tsne_fit[, 1],
    Axis2 = tsne_fit[, 2]
  ),
) %>%
  left_join(
    crosswalk %>% dplyr::select(OriginalRunSample, CorrectedSampleID, SamplingMethod),
    by = "OriginalRunSample"
  ) %>%
  mutate(Method = factor(Method, levels = c("PCA", "UMAP", "t-SNE")))

figure3 <- ggplot(
  structure_df,
  aes(Axis1, Axis2, colour = SamplingMethod)
) +
  geom_point(size = 2.15, shape = 16, alpha = 0.86) +
  facet_wrap(~Method, scales = "free", ncol = 3) +
  scale_colour_manual(values = analysis_colors) +
  labs(
    x = "Dimension 1",
    y = "Dimension 2",
    colour = NULL,
    title = "Global expression structure of MN and bulk samples",
    subtitle = "PCA, UMAP and t-SNE; all samples are circles",
    caption = "Colour encodes sampling method only. These ordinations are descriptive and are not significance tests."
  ) +
  theme(
    aspect.ratio = 1,
    legend.position = "top"
  )
save_plot(
  figure3,
  "Figure3_expression_structure",
  183,
  90
)

structure_export <- structure_df %>%
  dplyr::select(CorrectedSampleID, SamplingMethod, Method, Axis1, Axis2)
write.csv(
  structure_export,
  file.path(source_dir, "Figure3_expression_structure_source_data.csv"),
  row.names = FALSE
)

# Figure 4: comparison between MN and destructive bulk profiles
mn_samples <- crosswalk$OriginalRunSample[crosswalk$SamplingMethod == "MN"]
bulk_samples <- crosswalk$OriginalRunSample[crosswalk$SamplingMethod == "Bulk"]
profile_df <- data.frame(
  Gene_ID = rownames(vst_matrix),
  Mean_MN = rowMeans(vst_matrix[, mn_samples, drop = FALSE]),
  Mean_Bulk = rowMeans(vst_matrix[, bulk_samples, drop = FALSE])
)
profile_r <- cor(profile_df$Mean_MN, profile_df$Mean_Bulk)

p_profile <- ggplot(profile_df, aes(Mean_MN, Mean_Bulk)) +
  geom_point(size = 0.70, alpha = 0.30, colour = analysis_colors[["MN"]]) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = analysis_colors[["Bulk"]],
    fill = "#DDB7A7",
    alpha = 0.24,
    linewidth = 0.65
  ) +
  labs(
    x = "Group-mean VST expression, MN",
    y = "Group-mean VST expression, Bulk",
    title = "Group-level expression concordance",
    subtitle = sprintf("Pearson r = %.3f; %s genes", profile_r, comma(nrow(profile_df)))
  ) +
  theme(aspect.ratio = 1)

same_method_corr <- vst_matrix[top_genes, , drop = FALSE]
corr_all <- cor(same_method_corr)
sample_similarity <- crosswalk %>%
  dplyr::select(OriginalRunSample, CorrectedSampleID, SamplingMethod) %>%
  rowwise() %>%
  mutate(
    MeanWithinMethodR = mean(
      corr_all[
        OriginalRunSample,
        crosswalk$OriginalRunSample[
          crosswalk$SamplingMethod == SamplingMethod &
            crosswalk$OriginalRunSample != OriginalRunSample
        ]
      ]
    )
  ) %>%
  ungroup()
similarity_p <- anova_p(sample_similarity, "MeanWithinMethodR")

p_similarity <- ggplot(
  sample_similarity,
  aes(SamplingMethod, MeanWithinMethodR, colour = SamplingMethod)
) +
  geom_jitter(width = 0.08, size = 1.4, alpha = 0.78, shape = 16) +
  stat_summary(
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "errorbar",
    width = 0.30,
    linewidth = 0.65,
    colour = "#202020"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    label = paste0("ANOVA P = ", format_p(similarity_p)),
    vjust = 1.55,
    size = 2.1
  ) +
  scale_colour_manual(values = analysis_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.20))) +
  labs(
    x = NULL,
    y = "Mean within-method Pearson r",
    title = "Sample concordance within each method",
    subtitle = "One value per run"
  ) +
  theme(aspect.ratio = 1)

heat_genes <- names(sort(apply(vst_matrix, 1, var), decreasing = TRUE))[1:80]
heat_matrix <- t(scale(t(vst_matrix[heat_genes, , drop = FALSE])))
heat_order <- c(mn_samples, bulk_samples)
heat_df <- as.data.frame(heat_matrix[, heat_order, drop = FALSE]) %>%
  tibble::rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, names_to = "OriginalRunSample", values_to = "Z") %>%
  left_join(
    crosswalk %>% dplyr::select(OriginalRunSample, SamplingMethod),
    by = "OriginalRunSample"
  ) %>%
  mutate(OriginalRunSample = factor(OriginalRunSample, levels = heat_order))

p_heat <- ggplot(heat_df, aes(OriginalRunSample, reorder(Gene_ID, Z), fill = Z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#4F7896",
    mid = "white",
    high = "#B65E49",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = squish
  ) +
  annotate(
    "segment",
    x = 0.5,
    xend = length(mn_samples) + 0.5,
    y = 81.2,
    yend = 81.2,
    colour = analysis_colors[["MN"]],
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = (1 + length(mn_samples)) / 2,
    y = 83.2,
    label = "MN",
    colour = analysis_colors[["MN"]],
    size = 2.3,
    fontface = "bold"
  ) +
  annotate(
    "segment",
    x = length(mn_samples) + 0.5,
    xend = ncol(vst_matrix) + 0.5,
    y = 81.2,
    yend = 81.2,
    colour = analysis_colors[["Bulk"]],
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = length(mn_samples) + (length(bulk_samples) + 1) / 2,
    y = 83.2,
    label = "Bulk",
    colour = analysis_colors[["Bulk"]],
    size = 2.3,
    fontface = "bold"
  ) +
  coord_cartesian(ylim = c(0.5, 84.5), clip = "off") +
  labs(
    x = NULL,
    y = NULL,
    fill = "Z-score",
    title = "Top 80 variable genes"
  ) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

figure4 <- (wrap_plots(p_profile, p_similarity, p_heat) +
  plot_layout(design = "AB\nCC", guides = "keep") +
  plot_annotation(
    title = "Concordance between micro-needle and destructive bulk expression profiles",
    subtitle = "Descriptive comparison of group-level structure and within-method sample similarity",
    caption = "Pairwise-sample dependence is avoided by calculating one mean within-method correlation per run before ANOVA.",
    tag_levels = "a"
  ))
save_plot(
  figure4,
  "Figure4_MN_Bulk_comparison",
  183,
  190,
  panel_ids = letters[1:3]
)

write.csv(
  profile_df,
  file.path(source_dir, "Figure4_group_mean_expression_source_data.csv"),
  row.names = FALSE
)
write.csv(
  sample_similarity %>%
    dplyr::select(CorrectedSampleID, SamplingMethod, MeanWithinMethodR),
  file.path(source_dir, "Figure4_within_method_similarity_source_data.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    Comparison = c(
      "Group-mean expression concordance",
      "Within-method sample concordance"
    ),
    Statistic = c("Pearson r", "ANOVA P"),
    Value = c(profile_r, similarity_p),
    Definition = c(
      "Correlation of MN and Bulk group means across filtered genes",
      "One-way ANOVA of one mean within-method correlation per run"
    )
  ),
  file.path(source_dir, "Figure4_statistics.csv"),
  row.names = FALSE
)

# Figure 5: data-derived ordering relative to recorded sampling Day
pt <- read.csv(
  file.path(pseudotime_dir, "tables",
            "pseudotime_method_comparison.csv"),
  check.names = FALSE
) %>%
  left_join(
    crosswalk %>%
      dplyr::select(OriginalRunSample, CorrectedSampleID, SeedID, SamplingMethod, DayCrosswalk = Day),
    by = c("RunSample" = "OriginalRunSample")
  ) %>%
  filter(SamplingMethod == "MN") %>%
  mutate(Day = factor(Day, levels = c("Day0", "Day1", "Day2", "Day3")))

methods <- c("Slingshot", "TSCAN", "SCORPIUS")
pt_long <- pt %>%
  dplyr::select(CorrectedSampleID, SeedID, Day, all_of(methods)) %>%
  pivot_longer(all_of(methods), names_to = "Method", values_to = "Pseudotime") %>%
  group_by(Method) %>%
  mutate(
    ScaledPseudotime = (Pseudotime - min(Pseudotime, na.rm = TRUE)) /
      diff(range(Pseudotime, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(Method = factor(Method, levels = methods))

day_rank <- sapply(methods, function(method) {
  cor(pt[[method]], as.integer(pt$Day), method = "spearman")
})
method_rank <- cor(pt[, methods], method = "spearman")

p_seed <- ggplot(
  pt_long,
  aes(as.integer(Day) - 1L, ScaledPseudotime, colour = SeedID, group = SeedID)
) +
  geom_line(linewidth = 0.50, alpha = 0.82) +
  geom_point(size = 1.65, shape = 16) +
  facet_wrap(~Method, nrow = 3) +
  scale_colour_manual(values = seed_colors) +
  scale_x_continuous(breaks = 0:3, labels = paste0("D", 0:3)) +
  labs(
    x = "Recorded sampling Day",
    y = "Scaled data-derived pseudotime",
    colour = NULL,
    title = "Pseudotime trajectories by seed"
  ) +
  theme(legend.position = "top")

rank_df <- data.frame(
  Method = factor(methods, levels = methods),
  SpearmanRho = as.numeric(day_rank[methods])
)
p_day_rank <- ggplot(rank_df, aes(Method, SpearmanRho, fill = Method)) +
  geom_col(width = 0.58, colour = "white", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#555555") +
  geom_text(
    aes(label = sprintf("%.2f", SpearmanRho)),
    hjust = ifelse(rank_df$SpearmanRho >= 0, -0.18, 1.18),
    size = 2.2
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c(
    Slingshot = "#356C82",
    TSCAN = "#78A58F",
    SCORPIUS = "#C57C5A"
  ), guide = "none") +
  scale_y_continuous(limits = c(-0.10, 0.60), expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = NULL,
    y = "Spearman rho with recorded Day",
    title = "Agreement with calendar time"
  )

method_heat <- as.data.frame(as.table(method_rank)) %>%
  setNames(c("Method1", "Method2", "SpearmanRho")) %>%
  mutate(
    Method1 = factor(Method1, levels = methods),
    Method2 = factor(Method2, levels = methods)
  )
p_method_heat <- ggplot(method_heat, aes(Method1, Method2, fill = SpearmanRho)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.2f", SpearmanRho)), size = 2.2) +
  scale_fill_gradient2(
    low = "#4F7896",
    mid = "white",
    high = "#B65E49",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Spearman rho",
    title = "Agreement among pseudotime methods"
  ) +
  theme(aspect.ratio = 1)

selected_vst <- read.csv(
  file.path(pseudotime_dir, "tables",
            "selected_n5_VST_expression.csv"),
  check.names = FALSE
)
selected_matrix <- as.matrix(selected_vst[, -(1:2), drop = FALSE])
rownames(selected_matrix) <- selected_vst$Gene_ID
selected_umap <- read.csv(
  file.path(pseudotime_dir, "tables",
            "selected_UMAP_coordinates.csv"),
  check.names = FALSE
)
umap_mat <- as.matrix(selected_umap[, c("UMAP1", "UMAP2")])
rownames(umap_mat) <- selected_umap$RunSample

umap_points <- selected_umap %>%
  left_join(
    crosswalk %>% dplyr::select(OriginalRunSample, CorrectedSampleID, SeedID),
    by = c("RunSample" = "OriginalRunSample")
  ) %>%
  left_join(
    pt %>% dplyr::select(RunSample, Slingshot),
    by = "RunSample"
  )
day_centroids <- umap_points %>%
  group_by(Day) %>%
  summarise(
    UMAP1 = median(UMAP1, na.rm = TRUE),
    UMAP2 = median(UMAP2, na.rm = TRUE),
    MedianSlingshot = median(Slingshot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(MedianSlingshot) %>%
  mutate(
    UMAP1End = lead(UMAP1),
    UMAP2End = lead(UMAP2)
  ) %>%
  filter(!is.na(UMAP1End), !is.na(UMAP2End))
p_umap_path <- ggplot() +
  geom_point(
    data = umap_points,
    aes(UMAP1, UMAP2, colour = Day),
    size = 2.15,
    shape = 16,
    alpha = 0.90
  ) +
  geom_segment(
    data = day_centroids,
    aes(
      x = UMAP1,
      y = UMAP2,
      xend = UMAP1End,
      yend = UMAP2End
    ),
    inherit.aes = FALSE,
    colour = "#B46A4F",
    linewidth = 0.70,
    arrow = arrow(type = "closed", length = unit(1.7, "mm"))
  ) +
  scale_colour_manual(values = day_colors) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    colour = "Day",
    title = "Day centroids ordered by median Slingshot pseudotime"
  ) +
  theme(aspect.ratio = 1)

figure5 <- (wrap_plots(
  p_umap_path, p_seed, p_day_rank, p_method_heat,
  ncol = 2,
  guides = "keep"
) +
  plot_annotation(
    title = "Data-derived ordering compared with recorded sampling Day",
    subtitle = "Pseudotime is inferred from expression and is not equivalent to experimental time",
    caption = "Spearman correlations describe ordering only. The analysis does not establish a unique true trajectory.",
    tag_levels = "a"
  )
)
save_plot(
  figure5,
  "Figure5_pseudotime_comparison",
  160,
  160,
  panel_ids = letters[1:4]
)

write.csv(
  pt_long %>%
    dplyr::select(CorrectedSampleID, SeedID, Day, Method, Pseudotime, ScaledPseudotime),
  file.path(source_dir, "Figure5_pseudotime_by_seed_source_data.csv"),
  row.names = FALSE
)
write.csv(
  rank_df,
  file.path(source_dir, "Figure5_pseudotime_vs_day_source_data.csv"),
  row.names = FALSE
)
write.csv(
  method_heat,
  file.path(source_dir, "Figure5_method_agreement_source_data.csv"),
  row.names = FALSE
)

# Figure 6: examples of analyses supported by the dataset
de_summary <- read.csv(
  file.path(selected5_dir, "tables",
            "DE_summary_selected_n5.csv"),
  check.names = FALSE
) %>%
  mutate(
    Comparison = factor(
      comparison,
      levels = c(
        "Day1_vs_Day0", "Day2_vs_Day0", "Day3_vs_Day0",
        "Day2_vs_Day1", "Day3_vs_Day1", "Day3_vs_Day2"
      )
    )
  ) %>%
  select(Comparison, Up = up, Down = down) %>%
  pivot_longer(c(Up, Down), names_to = "Direction", values_to = "Genes")

p_de <- ggplot(de_summary, aes(Comparison, Genes, fill = Direction)) +
  geom_col(width = 0.68, position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_fill_manual(values = c(Up = "#B46A4F", Down = "#6E91A6")) +
  labs(
    x = NULL,
    y = "Significant genes",
    fill = NULL,
    title = "Pairwise differential expression"
  ) +
  theme(legend.position = "top")

candidate_anova <- read.csv(
  file.path(selected5_dir, "tables",
            "candidate_gene_one_way_ANOVA.csv"),
  check.names = FALSE
) %>%
  mutate(
    PText = paste0("P = ", formatC(ANOVA_p, format = "e", digits = 2))
  )

candidate_expr <- selected_vst %>%
  filter(Symbol %in% candidate_anova$Symbol) %>%
  select(Gene_ID, Symbol, all_of(mn_samples)) %>%
  pivot_longer(all_of(mn_samples), names_to = "OriginalRunSample", values_to = "VST") %>%
  left_join(
    crosswalk %>% dplyr::select(OriginalRunSample, CorrectedSampleID, Day),
    by = "OriginalRunSample"
  ) %>%
  left_join(
    candidate_anova %>% dplyr::select(Symbol),
    by = "Symbol"
  )

p_candidate <- ggplot(
  candidate_expr,
  aes(Day, VST, colour = Day)
) +
  geom_boxplot(
    width = 0.52,
    fill = "white",
    outlier.shape = NA,
    linewidth = 0.35
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 1.15,
    alpha = 0.76
  ) +
  geom_text(
    data = candidate_anova,
    aes(x = 3.68, y = Inf, label = PText),
    inherit.aes = FALSE,
    hjust = 0.98,
    vjust = 1.45,
    size = 1.75,
    colour = "#303030"
  ) +
  facet_wrap(~Symbol, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = day_colors, guide = "none") +
  labs(
    x = NULL,
    y = "VST expression",
    title = "Candidate-gene expression"
  )

mfuzz_assign <- read.csv(
  file.path(selected5_dir, "tables",
            "Mfuzz_selected_gene_clusters.csv"),
  check.names = FALSE
) %>%
  select(Gene_ID, cluster)

mfuzz_int <- read.csv(
  file.path(selected5_dir, "tables",
            "Mfuzz_cluster_interpretation.csv"),
  check.names = FALSE
)

selected_z <- t(scale(t(selected_matrix)))
selected_z[!is.finite(selected_z)] <- 0
mfuzz_long <- as.data.frame(selected_z[, mn_samples, drop = FALSE]) %>%
  tibble::rownames_to_column("Gene_ID") %>%
  pivot_longer(-Gene_ID, names_to = "OriginalRunSample", values_to = "Z") %>%
  left_join(mfuzz_assign, by = "Gene_ID") %>%
  left_join(
    crosswalk %>% dplyr::select(OriginalRunSample, Day),
    by = "OriginalRunSample"
  ) %>%
  group_by(cluster, Day) %>%
  summarise(MeanZ = mean(Z, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    mfuzz_int %>% dplyr::select(cluster = Cluster, Peak_day),
    by = "cluster"
  )

p_mfuzz <- ggplot(
  mfuzz_long,
  aes(as.integer(Day) - 1L, MeanZ, colour = factor(cluster), group = cluster)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#999999") +
  geom_line(linewidth = 0.58) +
  geom_point(size = 1.45) +
  facet_wrap(
    ~paste0("C", cluster, " | peak ", Peak_day),
    ncol = 2,
    scales = "fixed"
  ) +
  scale_x_continuous(breaks = 0:3, labels = paste0("D", 0:3)) +
  scale_colour_manual(
    values = c("#356C82", "#78A58F", "#C57C5A", "#8B98A8", "#B08B57", "#736A99"),
    guide = "none"
  ) +
  labs(
    x = "Day index",
    y = "Mean within-gene Z-score",
    title = "Mfuzz temporal modules"
  )

mfuzz_enrichment <- read.csv(
  file.path(selected5_dir, "tables",
            "Mfuzz_selected_cluster_enrichment.csv"),
  check.names = FALSE
) %>%
  filter(
    significant %in% c(TRUE, "TRUE"),
    is.finite(p_value),
    source %in% c("GO:BP", "GO:CC", "GO:MF", "KEGG")
  ) %>%
  group_by(cluster, source) %>%
  slice_min(order_by = p_value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    NegLog10P = -log10(pmax(p_value, .Machine$double.xmin)),
    Label = paste0(substr(name, 1, 30), " [", cluster, "]")
  ) %>%
  arrange(desc(NegLog10P)) %>%
  slice_head(n = 6)

p_enrichment <- ggplot(
  mfuzz_enrichment,
  aes(NegLog10P, reorder(Label, NegLog10P), fill = source)
) +
  geom_col(width = 0.60) +
  scale_fill_manual(
    values = c(
      "GO:BP" = "#356C82",
      "GO:CC" = "#78A58F",
      "GO:MF" = "#B08B57",
      "KEGG" = "#C57C5A"
    )
  ) +
  labs(
    x = "-log10(P)",
    y = NULL,
    fill = "Source",
    title = "Leading enriched terms"
  )

figure6 <- wrap_plots(
  p_de, p_candidate, p_mfuzz, p_enrichment,
  ncol = 2,
  guides = "keep"
) +
  plot_annotation(
    title = "Examples of transcriptomic analyses supported by the dataset",
    subtitle = "Differential expression, candidate-gene profiles, temporal modules and functional enrichment",
    caption = "These panels show reuse of the processed matrices and are not used to make causal claims.",
    tag_levels = "a"
  )
save_plot(
  figure6,
  "Figure6_dataset_supported_analyses",
  183,
  180,
  panel_ids = letters[1:4],
  exemptions = lapply(letters[1:4], function(id) {
    list(
      panels = id,
      checks = c("panel-width"),
      reason = "Mixed plot types have intentionally different y-axis label widths."
    )
  })
)

# The user-selected main-data figure retains four polished panels from the
# earlier analysis report. The generated quantitative summary above remains
# available as source data but is not used as the final Figure 6.
legacy_figure_paths <- file.path(root, "assets", "04_candidate_gene_ANOVA.png")
if (!file.exists(legacy_figure_paths)) {
  stop("The finalized candidate-gene figure is missing.")
}

legacy_image_panel <- function(path) {
  image <- png::readPNG(path)
  ratio <- dim(image)[1] / dim(image)[2]
  ggplot() +
    annotation_raster(
      image,
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = ratio,
      interpolate = TRUE
    ) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, ratio), expand = FALSE) +
    theme_void()
}

figure6_legacy <- legacy_image_panel(legacy_figure_paths)
save_plot(
  figure6_legacy,
  "Figure6_dataset_supported_analyses",
  183,
  153
)

write.csv(
  data.frame(
    Panel = "Figure 6",
    SourceImage = legacy_figure_paths,
    SourceType = "Finalized candidate-gene ANOVA and Tukey figure from the earlier selected-5 report"
  ),
  file.path(source_dir, "Figure6_legacy_source_manifest.csv"),
  row.names = FALSE
)

write.csv(de_summary, file.path(source_dir, "Figure6_DE_summary_source_data.csv"), row.names = FALSE)
write.csv(candidate_expr, file.path(source_dir, "Figure6_candidate_genes_source_data.csv"), row.names = FALSE)
write.csv(mfuzz_long, file.path(source_dir, "Figure6_Mfuzz_trajectories_source_data.csv"), row.names = FALSE)
write.csv(mfuzz_enrichment, file.path(source_dir, "Figure6_enrichment_source_data.csv"), row.names = FALSE)

writeLines(
  c(
    "Figure generation completed.",
    sprintf("R version: %s", R.version.string),
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "MN/Bulk colours are the only sampling-method encoding in Figures 2-4.",
    "Original B1/B2 labels are not used as plot mappings or in exported figure source data.",
    "Figure 5 shows data-derived ordering as an exploratory comparison, not as sample-identity evidence."
  ),
  file.path(qa_dir, "figure_generation_notes.txt")
)
