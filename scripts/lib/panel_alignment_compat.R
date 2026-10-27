# Compatibility wrapper for patchwork layouts whose null units are not
# resolved by older patchwork/grid combinations.

alignment_library <- Sys.getenv("PANEL_ALIGNMENT_R", "")
if (nzchar(alignment_library) && file.exists(alignment_library)) {
  source(alignment_library)
}

if (!exists("require_patchwork_panel_alignment", mode = "function")) {
  require_patchwork_panel_alignment <- function(...) {
    message("Panel-alignment audit skipped; set PANEL_ALIGNMENT_R to enable it.")
    invisible(FALSE)
  }
}

.resolve_null_units_pt <- function(units, total_pt) {
  labels <- as.character(units)
  values <- grid::convertUnit(units, "pt", valueOnly = TRUE)
  null_index <- grepl("null$", labels)
  if (!any(null_index)) {
    return(values)
  }
  shares <- suppressWarnings(as.numeric(sub("null$", "", labels[null_index])))
  if (any(!is.finite(shares)) || any(shares <= 0)) {
    stop("Patchwork layout contains an invalid null unit.", call. = FALSE)
  }
  remaining <- total_pt - sum(values[!null_index])
  if (!is.finite(remaining) || remaining <= 0) {
    stop("Patchwork layout leaves no finite space for null units.", call. = FALSE)
  }
  values[null_index] <- remaining * shares / sum(shares)
  values
}

write_patchwork_panel_layout <- function(
    plot,
    manifest_path,
    width_in,
    height_in,
    panel_ids = NULL,
    row_groups = NULL,
    column_groups = NULL,
    exemptions = list()) {
  for (package in c("patchwork", "grid", "jsonlite")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for patchwork alignment QA", package), call. = FALSE)
    }
  }
  if (!is.numeric(width_in) || length(width_in) != 1 || !is.finite(width_in) || width_in <= 0 ||
      !is.numeric(height_in) || length(height_in) != 1 || !is.finite(height_in) || height_in <= 0) {
    stop("width_in and height_in must be positive finite numbers", call. = FALSE)
  }

  grob <- patchwork::patchworkGrob(plot)
  panel_rows <- grob$layout[
    grepl(
      "^panel(?:;|$)|^panel-(?:[0-9]+|nested-patchwork-[0-9]+)$",
      grob$layout$name,
      perl = TRUE
    ),
    ,
    drop = FALSE
  ]
  if (nrow(panel_rows) < 2) {
    stop("At least two patchwork panel cells are required for alignment QA", call. = FALSE)
  }
  panel_rows <- panel_rows[order(panel_rows$t, panel_rows$l), , drop = FALSE]

  if (is.null(panel_ids)) {
    panel_ids <- vapply(
      seq_len(nrow(panel_rows)),
      function(index) if (index <= length(letters)) letters[index] else paste0("panel-", index),
      character(1)
    )
  }
  panel_ids <- as.character(panel_ids)
  if (length(panel_ids) != nrow(panel_rows) || anyDuplicated(panel_ids) || any(!nzchar(panel_ids))) {
    stop("panel_ids must be unique, non-empty, and match the measured patchwork panels", call. = FALSE)
  }

  probe_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(
    probe_path,
    width = width_in,
    height = height_in,
    useDingbats = FALSE,
    family = "sans"
  )
  on.exit({
    if (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
    unlink(probe_path)
  }, add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(grob)
  grid::grid.force()

  widths_pt <- .resolve_null_units_pt(grob$widths, width_in * 72)
  heights_pt <- .resolve_null_units_pt(grob$heights, height_in * 72)
  x_edges <- c(0, cumsum(widths_pt))
  top_edges <- c(0, cumsum(heights_pt))
  total_height_pt <- sum(heights_pt)

  panels <- lapply(seq_len(nrow(panel_rows)), function(index) {
    row <- panel_rows[index, ]
    left <- x_edges[row$l]
    right <- x_edges[row$r + 1]
    top <- total_height_pt - top_edges[row$t]
    bottom <- total_height_pt - top_edges[row$b + 1]
    list(
      id = panel_ids[index],
      bbox_pt = unname(c(left, bottom, right, top)),
      grid_id = "patchwork-grid-1",
      row_start = as.integer(row$t - 1),
      row_stop = as.integer(row$b),
      col_start = as.integer(row$l - 1),
      col_stop = as.integer(row$r)
    )
  })

  manifest <- list(
    schema_version = 1L,
    backend = "r-patchwork",
    figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
    panels = panels,
    exemptions = exemptions
  )
  normalized_rows <- .nature_alignment_groups(row_groups)
  normalized_columns <- .nature_alignment_groups(column_groups)
  if (!is.null(normalized_rows)) {
    manifest$row_groups <- normalized_rows
  }
  if (!is.null(normalized_columns)) {
    manifest$column_groups <- normalized_columns
  }

  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    pattern = ".panel-layout-",
    tmpdir = dirname(manifest_path),
    fileext = ".json"
  )
  jsonlite::write_json(
    manifest,
    temporary_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    digits = NA
  )
  if (!file.rename(temporary_path, manifest_path)) {
    unlink(temporary_path)
    stop("Could not activate the patchwork alignment manifest atomically", call. = FALSE)
  }
  grDevices::dev.off()
  unlink(probe_path)
  invisible(manifest)
}
