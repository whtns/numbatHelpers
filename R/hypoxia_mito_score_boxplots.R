# Per-cluster mean HYPOXIA and MITOCHONDRIAL module scores, side by side.
#
# Companion to the composite-score boxplot (src/plot_hypoxia_mean_score_boxplots.R,
# results/hypoxia_cluster_split/hypoxia_mean_score_boxplots.pdf), which plots the
# single statistic the split actually thresholds on: mean_score, the per-cluster
# mean of `hypoxia_score`.
#
# `hypoxia_score` is a COMPOSITE. add_hypoxia_score() builds it as
#     hypoxia_score = rescale(rowMeans(hypoxia, MT), 0..1),   with MT stored NEGATED
# so a cluster can reach the median + 3*MAD fence two different ways: genuinely
# high HALLMARK_HYPOXIA expression, or merely low mitochondrial content. The
# composite cannot be decomposed after the fact, so this plot shows the two
# ingredients separately on a shared axis. A cluster whose hypoxia box sits at
# the bulk median while its mito box sits far below it was flagged by the
# mitochondrial half of the score, not by hypoxia.
#
# DATA SOURCE -- two paths, preferred first:
#   1. hypoxia_split_log_all.csv, if it carries mean_hypoxia / mean_mito.
#      identify_hypoxia_clusters() records them, so every split run from that
#      change onward logs them directly.
#   2. Otherwise replay the split's round-1 sweep from the *_seu_hypoxia.rds
#      objects (replay_cluster_score_means()) and reconcile against the log.
#      Logs written before the column was added have no other way to be read,
#      and hypoxia_partition_paths is pinned with cue(command = FALSE), so
#      editing the split does not by itself refresh the log.

#' Replay the hypoxia split's round-1 sweep to recover per-cluster score means
#'
#' Reproduces the clustering that [split_hypoxia_by_clusters()] performs in its
#' round-1 resolution sweep and returns the per-cluster mean of the raw hypoxia
#' and mitochondrial module scores — the two ingredients of `hypoxia_score` that
#' the older split logs do not record.
#'
#' The clustering is reproduced, not re-derived: same assay, same
#' `FindNeighbors(dims = 1:30, k.param = min(20, n - 1))` on the stored `pca`,
#' same `FindClusters` per resolution. Seurat's Louvain is seeded, so given the
#' same object this returns the same cluster assignments. Only the expensive
#' `presto::wilcoxauc()` marker scan is skipped — it feeds the marker gate, not
#' the score means. Verify rather than trust: the caller reconciles `n_cells` and
#' `mean_score` against the log (see [plot_hypoxia_mito_score_boxplots()]).
#'
#' @param hypoxia_seu_paths Character vector of `*_seu_hypoxia.rds` paths (the
#'   `hypoxia_seus` target). `NA`/missing entries are skipped.
#' @param resolutions Resolutions swept, matching the pipeline's call to
#'   [split_hypoxia_by_clusters()] (default `seq(0.2, 1.2, by = 0.2)`).
#' @param split_assay Assay the split clusters on (default `"gene"`).
#' @return data.frame with `sample_id`, `resolution`, `cluster`, `n_cells`,
#'   `mean_hypoxia`, `mean_mito`, `mean_score`. Zero rows if nothing could be
#'   replayed. Never errors on a single bad sample — it is skipped with a message.
#' @export
replay_cluster_score_means <- function(hypoxia_seu_paths,
                                       resolutions = seq(0.2, 1.2, by = 0.2),
                                       split_assay = "gene") {
  paths <- unlist(hypoxia_seu_paths, use.names = FALSE)
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  paths <- paths[file.exists(paths)]
  resolutions <- sort(unique(round(as.numeric(resolutions), 3)))

  empty <- data.frame(sample_id = character(0), resolution = numeric(0),
                      cluster = character(0), n_cells = integer(0),
                      mean_hypoxia = numeric(0), mean_mito = numeric(0),
                      mean_score = numeric(0), stringsAsFactors = FALSE)
  if (length(paths) == 0) return(empty)

  out <- lapply(paths, function(p) {
    sample_id <- stringr::str_extract(p, "SR[RX][0-9]+")
    tryCatch({
      seu <- readRDS(p)

      # Same assay resolution as split_hypoxia_by_clusters().
      da <- if (split_assay %in% names(seu@assays)) split_assay else
        Seurat::DefaultAssay(seu)
      Seurat::DefaultAssay(seu) <- da

      # Same guards. NormalizeData is deliberately NOT repeated: clustering runs
      # off the stored pca, so the data layer cannot change the assignments.
      if (!"pca" %in% names(seu@reductions)) {
        message(sample_id, ": no pca reduction -> cannot replay")
        return(NULL)
      }
      if (ncol(seu) < 20) {
        message(sample_id, ": fewer than 20 cells -> split skipped it too")
        return(NULL)
      }
      if (!all(c("hypoxia", "MT", "hypoxia_score") %in% colnames(seu@meta.data))) {
        message(sample_id, ": object lacks the hypoxia/MT score columns -> skip")
        return(NULL)
      }

      n_cells <- ncol(seu)
      sub <- Seurat::FindNeighbors(
        seu, dims = 1:30, reduction = "pca",
        graph.name = paste0(da, c("_nn", "_snn")),
        k.param = min(20L, n_cells - 1L), verbose = FALSE)
      snn_name <- paste0(da, "_snn")

      per_res <- lapply(resolutions, function(res) {
        sub <- Seurat::FindClusters(sub, graph.name = snn_name,
                                    resolution = res, verbose = FALSE)
        grpv     <- as.character(sub$seurat_clusters)
        clusters <- sort(unique(grpv))
        f        <- factor(grpv, levels = clusters)
        cl_mean  <- function(x) unname(as.numeric(
          tapply(x, f, mean, na.rm = TRUE)[clusters]))

        data.frame(
          sample_id    = sample_id,
          resolution   = res,
          cluster      = clusters,
          n_cells      = as.integer(table(f)[clusters]),
          mean_hypoxia = cl_mean(sub$hypoxia),
          # MT is stored negated so it can be averaged into the composite; flip
          # it back so mean_mito reads in its natural direction (high = more mito).
          mean_mito    = cl_mean(-sub$MT),
          mean_score   = cl_mean(sub$hypoxia_score),
          stringsAsFactors = FALSE)
      })
      message("replayed ", sample_id, ": ", length(resolutions), " resolutions, ",
              n_cells, " cells")
      dplyr::bind_rows(per_res)
    }, error = function(e) {
      message("!! replay failed for ", sample_id, ": ", conditionMessage(e))
      NULL
    })
  })

  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(empty)
  dplyr::bind_rows(out)
}


#' Boxplot of per-cluster hypoxia and mitochondrial score by resolution
#'
#' Same layout as the composite-score boxplot — one facet per sample, clustering
#' resolution on x, one point per cluster coloured by its split decision — but
#' with two dodged boxes per resolution: the raw HALLMARK_HYPOXIA module score
#' and the raw mitochondrial module score. Sharing one y axis is the point: it
#' shows which of the two halves pushed a flagged cluster over the fence.
#'
#' Reads `mean_hypoxia` / `mean_mito` straight from the split log when present,
#' otherwise recovers them with [replay_cluster_score_means()] and keeps only
#' the clusters that reconcile with the log on both `n_cells` and `mean_score`,
#' so a cluster mapping that drifted cannot be plotted as if it were the split's.
#'
#' @param split_log_csv Path to the collated split log
#'   (`hypoxia_split_log_collated`).
#' @param hypoxia_seu_paths `hypoxia_seus` — only read if the log lacks the two
#'   columns. `NULL` disables the replay path.
#' @param out_pdf Output PDF path.
#' @param resolutions,split_assay Passed to [replay_cluster_score_means()].
#' @param tol Absolute tolerance when reconciling replayed `mean_score` against
#'   the log.
#' @return `out_pdf` on success, `NA_character_` if there was nothing to plot.
#' @export
plot_hypoxia_mito_score_boxplots <- function(
    split_log_csv,
    hypoxia_seu_paths = NULL,
    out_pdf = "results/hypoxia_cluster_split/hypoxia_mito_score_boxplots.pdf",
    resolutions = seq(0.2, 1.2, by = 0.2),
    split_assay = "gene",
    tol = 1e-6) {

  if (is.null(split_log_csv) || length(split_log_csv) != 1 ||
      is.na(split_log_csv) || !file.exists(split_log_csv)) {
    message("no split log at '", split_log_csv, "' -> nothing to plot")
    return(NA_character_)
  }

  d <- readr::read_csv(split_log_csv, show_col_types = FALSE)
  d$resolution <- round(as.numeric(d$resolution), 3)

  # Round 1 only: the confirmatory rounds re-cluster a shrinking pool at
  # resolutions computed as r_flag + k*recluster_step, so their cluster ids do
  # not share a grid with the sweep. Restricting to round 1 also makes the log
  # and replay paths cover exactly the same rows.
  d <- d[!is.na(d$round) & d$round == 1L, , drop = FALSE]
  if (nrow(d) == 0) {
    message("split log has no round-1 rows -> nothing to plot")
    return(NA_character_)
  }

  has_cols <- all(c("mean_hypoxia", "mean_mito") %in% colnames(d)) &&
    any(!is.na(d$mean_hypoxia)) && any(!is.na(d$mean_mito))

  if (has_cols) {
    message("using mean_hypoxia / mean_mito straight from the split log")
  } else {
    if (is.null(hypoxia_seu_paths)) {
      message("split log lacks mean_hypoxia / mean_mito and no hypoxia_seus ",
              "given -> nothing to plot")
      return(NA_character_)
    }
    message("split log predates mean_hypoxia / mean_mito -> replaying the ",
            "round-1 sweep from the hypoxia objects")
    rep_df <- replay_cluster_score_means(hypoxia_seu_paths,
                                         resolutions = resolutions,
                                         split_assay = split_assay)
    if (nrow(rep_df) == 0) {
      message("replay produced no rows -> nothing to plot")
      return(NA_character_)
    }

    # A log that HAS the columns but leaves them all NA would otherwise survive
    # the join and blank the plot -- the replayed values must win outright.
    d <- d[, setdiff(colnames(d), c("mean_hypoxia", "mean_mito")), drop = FALSE]
    d$cluster      <- as.character(d$cluster)
    rep_df$cluster <- as.character(rep_df$cluster)
    j <- dplyr::inner_join(
      d, rep_df[, c("sample_id", "resolution", "cluster", "n_cells",
                    "mean_hypoxia", "mean_mito", "mean_score")],
      by = c("sample_id", "resolution", "cluster"),
      suffix = c("", "_replay"))

    # Verify rather than trust. A cluster only survives if the replay reproduced
    # both its size and its composite mean -- together those pin the assignment.
    ok <- !is.na(j$mean_score_replay) &
      j$n_cells == j$n_cells_replay &
      abs(j$mean_score - j$mean_score_replay) < tol
    frac <- if (nrow(j) > 0) mean(ok) else 0
    message("reconciled ", sum(ok), "/", nrow(d), " logged round-1 clusters (",
            round(100 * frac), "% of the ", nrow(j), " joined)")
    if (sum(ok) == 0) {
      warning("replay reproduced none of the logged clusters -- refusing to ",
              "plot means that may belong to a different clustering")
      return(NA_character_)
    }
    if (frac < 0.95) {
      warning("only ", round(100 * frac), "% of joined clusters reconciled; ",
              "the unmatched ones are dropped from the plot")
    }
    d <- j[ok, , drop = FALSE]
  }

  d <- d[!is.na(d$mean_hypoxia) | !is.na(d$mean_mito), , drop = FALSE]
  if (nrow(d) == 0) {
    message("no clusters with both score means -> nothing to plot")
    return(NA_character_)
  }

  # Same three-state status encoding as the composite-score boxplot: being a
  # score outlier is necessary but not sufficient to be moved to the high subset.
  d$is_outlier <- as.logical(d$is_outlier)
  d$flagged    <- as.logical(d$flagged)
  d$status <- factor(
    dplyr::case_when(d$flagged ~ "flagged", d$is_outlier ~ "spared",
                     TRUE ~ "cluster"),
    levels = c("cluster", "spared", "flagged"))

  res_levels <- sort(unique(d$resolution))
  d$res_lab  <- factor(paste0("res ", res_levels[match(d$resolution, res_levels)]),
                       levels = paste0("res ", res_levels))
  d$lab <- ifelse(d$is_outlier,
                  paste0("c", d$cluster, " (n=", d$n_cells, ")"), NA_character_)

  n_by <- d |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::summarise(n_out  = sum(.data$is_outlier, na.rm = TRUE),
                     n_flag = sum(.data$flagged, na.rm = TRUE),
                     .groups = "drop")
  d <- dplyr::left_join(d, n_by, by = "sample_id")
  d$strip <- paste0(d$sample_id, " (", d$n_flag, "/", d$n_out, " flagged)")

  n_flag_tot <- sum(d$flagged, na.rm = TRUE)
  n_out_tot  <- sum(d$is_outlier, na.rm = TRUE)
  samples    <- sort(unique(d$sample_id))

  # ONE PAGE PER SCORE, each with its own y range, and a FREE y per sample.
  #
  # Unlike the composite-score boxplot, the axis is not shared across samples.
  # That plot can share one because hypoxia_score is rescaled to [0,1] within
  # each sample, so every facet already lives on the same axis. These are raw
  # module scores: their absolute level tracks library depth and composition, so
  # a shared axis is neither comparable across samples nor readable -- one
  # high-range sample flattens the other 32 into the bottom sliver of their
  # panels. The comparison that matters is within a sample anyway: where do the
  # flagged clusters sit relative to their OWN sample's bulk?
  #
  # The two scores are not on a comparable scale and cannot share an axis: the
  # mitochondrial module is 5 very highly expressed genes and lands around
  # 0.5-2.5, while HALLMARK_HYPOXIA is ~200 genes and lands near 0.05. Dodging
  # them into one panel flattens the hypoxia boxes to a line. Faceting sample x
  # score with free scales would fix that but costs the cross-sample comparison
  # the composite-score boxplot is built around, so instead each score keeps the
  # familiar shared-y layout on its own page. Facets, order, colours and labels
  # are identical between the pages, so flipping between them compares like
  # with like -- only the y range differs.
  page <- function(col, score_label, box_fill) {
    pd <- d
    pd$value <- pd[[col]]
    pd <- pd[!is.na(pd$value), , drop = FALSE]
    if (nrow(pd) == 0) return(NULL)

    ggplot2::ggplot(pd, ggplot2::aes(x = .data$res_lab, y = .data$value)) +
      ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6,
                            fill = box_fill, alpha = 0.35, colour = "grey45") +
      ggplot2::geom_jitter(ggplot2::aes(colour = .data$status),
                           width = 0.12, height = 0, size = 1.1, alpha = 0.85) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = .data$lab, colour = .data$status), size = 2.0,
        min.segment.length = 0, box.padding = 0.25, max.overlaps = Inf,
        na.rm = TRUE, seed = 1, show.legend = FALSE) +
      ggplot2::scale_colour_manual(
        values = c(cluster = "grey35", spared = "#e08214", flagged = "firebrick"),
        labels = c(cluster = "not an outlier",
                   spared  = "outlier, spared by gates (stays low)",
                   flagged = "flagged \u2192 high subset"),
        name = NULL, drop = FALSE) +
      ggplot2::facet_wrap(~ strip, ncol = 6, scales = "free_y") +
      ggplot2::labs(
        title = paste0("Cluster mean ", score_label,
                       " score by resolution \u2014 all samples (y free per sample)"),
        subtitle = paste0(
          n_flag_tot, " of ", n_out_tot, " score-outlier clusters flagged across ",
          length(samples), " samples (round-1 sweep). Clusters are flagged on the ",
          "COMPOSITE hypoxia_score = rescale(mean(hypoxia, -mito)), so a cluster ",
          "can clear the median + 3*MAD fence by being\nhigh on hypoxia genes OR ",
          "low on mitochondrial content. This page shows the ", score_label,
          " half alone; the companion page shows the other. Compare WITHIN a ",
          "panel, not across them:\nthese are raw module scores, so the y axis ",
          "is free per sample and the two pages do not share a range. ",
          "Mitochondrial is un-negated here, so higher = more mitochondrial."),
        x = "clustering resolution",
        y = paste0("cluster mean ", score_label, " score")) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        legend.position = "top",
        plot.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_text(size = 8))
  }

  # Same hues as the collage column annotations for these two scores, so a
  # reader moving between the figures keeps the colour association.
  pages <- list(page("mean_hypoxia", "HALLMARK_HYPOXIA", "#762A83"),
                page("mean_mito",    "mitochondrial",    "#8C510A"))
  pages <- pages[!vapply(pages, is.null, logical(1))]
  if (length(pages) == 0) {
    message("neither score had plottable values -> nothing to plot")
    return(NA_character_)
  }

  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(out_pdf, width = 20, height = 24)
  for (pg in pages) print(pg)
  invisible(grDevices::dev.off())
  message("wrote ", out_pdf)
  out_pdf
}
