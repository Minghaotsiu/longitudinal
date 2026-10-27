#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root_old <- Sys.getenv("BATCH1_UPSTREAM_DIR")
root_new <- Sys.getenv("BATCH2_UPSTREAM_DIR")
project_dir <- Sys.getenv("PROCESSED_DATA_DIR", "data/doudou_all34")
if (!nzchar(root_old) || !nzchar(root_new)) {
  stop("Set BATCH1_UPSTREAM_DIR and BATCH2_UPSTREAM_DIR.")
}
metadata <- read.csv(
  file.path(project_dir, "ALL_34_sample_metadata.csv"),
  check.names = FALSE
)

parse_fastqc <- function(zip_path) {
  if (length(zip_path) == 0 || is.na(zip_path[1])) {
    return(c(total_sequences = NA_real_, q30_read_pct = NA_real_, mean_quality = NA_real_))
  }
  members <- unzip(zip_path[1], list = TRUE)$Name
  data_member <- grep("fastqc_data\\.txt$", members, value = TRUE)[1]
  if (is.na(data_member)) {
    return(c(total_sequences = NA_real_, q30_read_pct = NA_real_, mean_quality = NA_real_))
  }
  lines <- readLines(unz(zip_path[1], data_member), warn = FALSE)

  total_line <- grep("^Total Sequences\\t", lines, value = TRUE)[1]
  total_sequences <- suppressWarnings(as.numeric(sub(".*\\t", "", total_line)))

  start <- grep("^>>Per sequence quality scores", lines)[1]
  end <- grep("^>>END_MODULE", lines)
  end <- end[end > start][1]
  if (is.na(start) || is.na(end)) {
    return(c(
      total_sequences = total_sequences,
      q30_read_pct = NA_real_,
      mean_quality = NA_real_
    ))
  }
  block <- lines[(start + 1):(end - 1)]
  block <- block[
    !grepl("^#", block) &
      !grepl("^>>", block) &
      nzchar(block)
  ]
  values <- do.call(
    rbind,
    strsplit(block, "\t", fixed = TRUE)
  )
  quality <- as.numeric(values[, 1])
  counts <- as.numeric(values[, 2])
  q30_read_pct <- 100 * sum(counts[quality >= 30], na.rm = TRUE) /
    sum(counts, na.rm = TRUE)
  mean_quality <- weighted.mean(quality, counts, na.rm = TRUE)
  c(
    total_sequences = total_sequences,
    q30_read_pct = q30_read_pct,
    mean_quality = mean_quality
  )
}

find_one <- function(root, pattern) {
  files <- list.files(
    root,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(files) == 0) {
    NA_character_
  } else {
    files[1]
  }
}

input_old <- read.csv(
  file.path(root_old, "ALL_17_RNA_MAP.csv"),
  check.names = FALSE
)
input_new <- read.csv(
  file.path(root_new, "ALL_17_RNA_MAP.csv"),
  check.names = FALSE
)
input_all <- rbind(
  data.frame(
    Sample = input_old$Sample,
    Batch = "B1",
    InputReads = input_old$`Number of input reads`,
    check.names = FALSE
  ),
  data.frame(
    Sample = input_new$Sample,
    Batch = "B2",
    InputReads = input_new$`Number of input reads`,
    check.names = FALSE
  )
)

raw_counts <- read.csv(
  file.path(project_dir, "ALL_34_raw_counts_matrix.csv"),
  check.names = FALSE
)
counts_1m <- read.csv(
  file.path(project_dir, "ALL_34_1M_counts_matrix.csv"),
  check.names = FALSE
)
stopifnot(identical(raw_counts$Gene_ID, counts_1m$Gene_ID))

raw_mat <- as.matrix(raw_counts[, -1, drop = FALSE])
one_m_mat <- as.matrix(counts_1m[, -1, drop = FALSE])
rownames(raw_mat) <- raw_counts$Gene_ID
rownames(one_m_mat) <- counts_1m$Gene_ID

qc_rows <- vector("list", nrow(metadata))
for (i in seq_len(nrow(metadata))) {
  run_sample <- metadata$RunSample[i]
  base_sample <- metadata$Sample[i]
  roots <- if (metadata$Batch[i] == "B1") root_old else root_new

  raw_r1 <- find_one(
    roots,
    paste0("^", base_sample, "_R1_fastqc\\.zip$")
  )
  raw_r2 <- find_one(
    roots,
    paste0("^", base_sample, "_R2_fastqc\\.zip$")
  )
  clean_r1 <- find_one(
    roots,
    paste0("^", base_sample, "_repair_R1_fastqc\\.zip$")
  )
  clean_r2 <- find_one(
    roots,
    paste0("^", base_sample, "_repair_R2_fastqc\\.zip$")
  )

  raw_r1_stats <- parse_fastqc(raw_r1)
  raw_r2_stats <- parse_fastqc(raw_r2)
  clean_r1_stats <- parse_fastqc(clean_r1)
  clean_r2_stats <- parse_fastqc(clean_r2)

  star_input_reads <- input_all$InputReads[
    input_all$Sample == base_sample & input_all$Batch == metadata$Batch[i]
  ][1]
  input_reads <- mean(
    c(raw_r1_stats["total_sequences"], raw_r2_stats["total_sequences"]),
    na.rm = TRUE
  )
  clean_reads <- mean(
    c(clean_r1_stats["total_sequences"], clean_r2_stats["total_sequences"]),
    na.rm = TRUE
  )
  q30_read_pct <- mean(
    c(raw_r1_stats["q30_read_pct"], raw_r2_stats["q30_read_pct"]),
    na.rm = TRUE
  )
  mean_quality <- mean(
    c(raw_r1_stats["mean_quality"], raw_r2_stats["mean_quality"]),
    na.rm = TRUE
  )

  qc_rows[[i]] <- data.frame(
    RunSample = run_sample,
    Sample = base_sample,
    Batch = metadata$Batch[i],
    Day = sub("d", "Day", metadata$Timepoint[i]),
    InputReads = input_reads,
    CleanReads = clean_reads,
    STARInputReads = star_input_reads,
    CleanReadRetentionPct = 100 * clean_reads / input_reads,
    Q30ReadPct = q30_read_pct,
    MeanReadQuality = mean_quality,
    TotalCounts = sum(raw_mat[, run_sample]),
    nGeneRaw = sum(raw_mat[, run_sample] > 0),
    nGeneGe10 = sum(raw_mat[, run_sample] >= 10),
    TotalCounts1M = sum(one_m_mat[, run_sample]),
    nGene1M = sum(one_m_mat[, run_sample] > 0),
    RawFastQC_R1 = raw_r1,
    RawFastQC_R2 = raw_r2,
    CleanFastQC_R1 = clean_r1,
    CleanFastQC_R2 = clean_r2,
    stringsAsFactors = FALSE
  )
}

qc_table <- do.call(rbind, qc_rows)
write.csv(
  qc_table,
  file.path(project_dir, "ALL_34_QC_metrics.csv"),
  row.names = FALSE
)
cat("rows=", nrow(qc_table), "\n", sep = "")
cat("missing_Q30=", sum(is.na(qc_table$Q30ReadPct)), "\n", sep = "")
cat("missing_clean=", sum(is.na(qc_table$CleanReads)), "\n", sep = "")
print(summary(qc_table[, c(
  "InputReads", "CleanReads", "Q30ReadPct", "TotalCounts",
  "nGeneRaw", "nGene1M", "TotalCounts1M"
)]))
