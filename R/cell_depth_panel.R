# Per-sample sequencing-depth panel for the sample summaries (github #43).
#
# What this is for
# ----------------
# Nothing in a sample summary currently says how deeply the sample was
# sequenced, so a noisy heatmap region or a thin clone call cannot be told apart
# from a shallow one. Depth across this cohort spans 23-fold (median
# nCount_gene 831 to 18,980), and the single lowest-depth sample (SRX10031191,
# median 831, 60% of cells under 1000 counts) is also the only one whose numbat
# run left no phylogeny. That is exactly the kind of thing this panel is meant
# to make obvious at a glance.
#
# Two things worth knowing about what is being drawn:
#
#  * The pipeline's own filter is `nCount_gene > 1000` (filter_sample_qc(), and
#    the annotation filter in generate_filtering_cell_counts()). It is drawn as
#    a reference line so a sample sitting mostly below it is visible.
#  * numbat does NOT apply that filter. run_numbat.R takes its cell set from the
#    UNFILTERED *_seu.rds and subsets on cell type only, and numbat's own
#    min_depth defaults to 0. So the cells drawn here are the cells numbat
#    actually saw -- which is the point of plotting the unfiltered stage.
#
# Data comes from the cell_qc_values table in batch_hashes.sqlite, which already
# holds per-cell nCount_gene / nFeature_gene / percent_mt for the whole cohort.
# Re-reading the Seurat objects to get the same numbers would cost minutes per
# sample for no gain.

#' Per-cell sequencing depth for one sample, from the QC database
#'
#' @param sample_id Sample accession.
#' @param sqlite_path Path to the metadata SQLite database.
#' @param stage_pattern Regex picking the processing stage to report. The
#'   default matches the unfiltered `<sample>_seu.rds`, which is the cell set
#'   numbat is run on.
#' @return A data frame with `cell`, `nCount_gene`, `nFeature_gene`,
#'   `percent_mt`, or a zero-row frame when the sample is not in the database.
#' @export
get_sample_cell_depth <- function(sample_id,
                                  sqlite_path   = "batch_hashes.sqlite",
                                  stage_pattern = "/[A-Z]{3}[0-9]+_seu\\.rds$") {
  empty <- data.frame(cell = character(), nCount_gene = numeric(),
                      nFeature_gene = numeric(), percent_mt = numeric(),
                      stringsAsFactors = FALSE)
  if (!file.exists(sqlite_path)) return(empty)

  con <- tryCatch(connect_hash_db(sqlite_path), error = function(e) NULL)
  if (is.null(con)) return(empty)
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

  d <- tryCatch(
    db_retry(DBI::dbGetQuery(con, "
      SELECT filepath, cell, nCount_gene, nFeature_gene, percent_mt
      FROM   cell_qc_values
      WHERE  sample_id = ?", params = list(sample_id))),
    error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(empty)

  # A sample can have rows for several processing stages. Prefer the unfiltered
  # object; fall back to whichever stage has the most cells rather than
  # returning nothing, so a renamed path does not blank the panel.
  hit <- d[grepl(stage_pattern, d$filepath), , drop = FALSE]
  if (nrow(hit) == 0) {
    biggest <- names(sort(table(d$filepath), decreasing = TRUE))[1]
    hit <- d[d$filepath == biggest, , drop = FALSE]
  } else if (length(unique(hit$filepath)) > 1) {
    biggest <- names(sort(table(hit$filepath), decreasing = TRUE))[1]
    hit <- hit[hit$filepath == biggest, , drop = FALSE]
  }
  hit[, c("cell", "nCount_gene", "nFeature_gene", "percent_mt")]
}

#' Draw the per-sample sequencing-depth panel
#'
#' A log-scaled distribution of `nCount_gene` with the pipeline's own
#' `> 1000` filter marked, annotated with the sample's median depth and the
#' fraction of cells below the line.
#'
#' @param sample_id Sample accession.
#' @param out_dir Directory for the PDF.
#' @param sqlite_path Path to the metadata SQLite database.
#' @param threshold The count threshold to mark. Defaults to the pipeline's own.
#' @param width,height PDF dimensions in inches.
#' @return Path to the PDF, or `NA_character_` when the sample has no QC rows.
#' @export
plot_cell_depth_panel <- function(sample_id,
                                  out_dir     = "results/cell_depth",
                                  sqlite_path = "batch_hashes.sqlite",
                                  threshold   = 1000,
                                  width  = 6,
                                  height = 4) {

  d <- get_sample_cell_depth(sample_id, sqlite_path = sqlite_path)
  if (nrow(d) == 0) return(NA_character_)

  v <- d$nCount_gene[!is.na(d$nCount_gene) & d$nCount_gene > 0]
  if (length(v) == 0) return(NA_character_)

  med    <- stats::median(v)
  n_below <- sum(v < threshold)
  pct_below <- 100 * n_below / length(v)

  # The floor is worth showing explicitly: this cohort was not filtered
  # uniformly upstream -- 13 samples retain cells down to ~500 counts while 22
  # never go below ~1250 -- so the minimum is a provenance signal, not noise.
  sub <- sprintf(
    "n = %s cells   median = %s   min = %s   below %s: %s (%.1f%%)",
    format(length(v), big.mark = ","), format(round(med), big.mark = ","),
    format(round(min(v)), big.mark = ","), format(threshold, big.mark = ","),
    format(n_below, big.mark = ","), pct_below)

  p <- ggplot2::ggplot(data.frame(nCount = v), ggplot2::aes(x = nCount)) +
    ggplot2::geom_histogram(bins = 60, fill = "grey35", colour = NA) +
    ggplot2::geom_vline(xintercept = threshold, colour = "firebrick",
                        linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = med, colour = "#2166AC",
                        linetype = "solid", linewidth = 0.5) +
    ggplot2::scale_x_log10(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    ggplot2::annotation_logticks(sides = "b") +
    ggplot2::labs(
      title    = paste0(sample_id, " - sequencing depth (unfiltered cells)"),
      subtitle = sub,
      caption  = paste0("dashed = pipeline filter (nCount_gene > ", threshold,
                        "); solid = median. numbat is run on these cells,",
                        " before that filter."),
      x = "nCount_gene (log scale)", y = "cells") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 7))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(out_dir, paste0(sample_id, "_cell_depth.pdf"))

  grDevices::pdf(pdf_path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  err <- tryCatch({ print(p); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(err)) {
    grid::grid.newpage()
    grid::grid.text(paste0(sample_id, "\n\ndepth panel could not be rendered\n\n",
                           paste(strwrap(err, width = 70), collapse = "\n")),
                    gp = grid::gpar(fontsize = 11, col = "firebrick"))
  }
  pdf_path
}
