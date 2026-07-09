# Cluster-based hypoxia split (single round, outlier mean-score rule).
#
# Alternative to the per-cell hypoxia_score threshold used by
# subset_seu_by_expression(): cluster the cells (at several resolutions in ONE
# round) and assign WHOLE clusters to the high-hypoxia subset when their mean
# hypoxia_score is a high-side statistical outlier among that resolution's
# per-cluster means. Hits are unioned across resolutions.
#
# The core rule flags a cluster "high" when its mean hypoxia_score is a
# high-side outlier by a robust-z rule, median(means) + 3*MAD. (Robust z is used
# instead of a Tukey fence because, when 2+ of a sample's few clusters are
# hypoxic, those clusters inflate Q3 and the fence rises above them; median/MAD
# track the bulk.) The score outlier is then gated by two filters (both on by
# default in split_hypoxia_by_clusters):
#   1. Marker gate (min_hyp_markers = 2): the cluster must carry >= 2 hypoxia
#      markers in its top-N. Drops pure cell-cycle-phase clusters (histone/PLK2
#      S-phase, CENPF/MKI67 G2M) whose hypoxia_score is inflated by the
#      snoRNA-host genes in the marker set but carry no real hypoxia markers.
#   2. Direct phase gate (max_dom_frac = 0.70): the cluster must be phase-MIXED
#      (dominant cell-cycle-phase fraction < 0.70). This catches what the marker
#      gate cannot — clusters that DO carry hypoxia markers (snoRNA-host
#      ride-along) yet are still single-phase, so filtering them would delete a
#      fine-grained phase state. Phase is scored once via
#      annotate_cell_cycle_without_1q(). See doc/hypoxia_split_outlier_rule.md
#      and the offline diagnostic src/eval_hypoxia_cluster_phase.R.
# Set min_hyp_markers = 0, max_dom_frac = NULL for the original score-only rule.
# See doc/hypoxia_split_outlier_rule.md for the rule and alternatives.

# snoRNA-host genes (the ZFAS1/GAS5 family that drives these clusters) plus
# classic canonical-hypoxia genes.
.hypoxia_marker_set <- c(
  "ZFAS1", "GAS5", "SNHG5", "SNHG6", "SNHG7", "SNHG8", "SNHG15", "SNHG29",
  "VEGFA", "NDRG1", "BNIP3", "BNIP3L", "SLC2A1", "CA9", "PGK1", "LDHA",
  "ENO1", "ANKRD37", "ADM", "P4HA1", "PLOD2"
)

#' Identify outlier-hypoxia clusters at one grouping
#'
#' A cluster is flagged hypoxia if its mean `hypoxia_score` is a high-side
#' outlier among the per-cluster means, using a robust z rule:
#' `mean_score > median(means) + mad_k * mad(means)` (MAD scaled to be
#' SD-consistent). This is robust to several elevated clusters because the
#' median and MAD are dominated by the bulk majority — unlike a Tukey fence,
#' whose Q3 is inflated by the very clusters it should catch. Flagging is
#' independent of hypoxia marker-gene presence; the marker columns
#' (`top_markers`, `n_hyp_markers`, `matched_genes`) are reported for
#' interpretation only.
#'
#' @param seu Seurat object with the grouping column `grp` and a `hypoxia_score`
#'   metadata column.
#' @param grp Name of the cluster grouping column in `seu@meta.data`.
#' @param hypoxia_genes Character vector of marker genes (default
#'   [.hypoxia_marker_set]); used only to populate the informational
#'   `n_hyp_markers` / `matched_genes` columns.
#' @param top_n Window of top markers (by logFC, padj < 0.5) scanned for the
#'   hypoxia-gene overlap columns.
#' @param top_markers_n Number of top markers per cluster reported in the
#'   `top_markers` column (default 5).
#' @param mad_k Robust-z multiplier for the outlier threshold (default 3).
#' @param min_hyp_markers Minimum number of hypoxia marker genes a cluster must
#'   carry in its top-`top_n` markers to be flagged, in addition to being a
#'   score outlier (default 0 = score-only, marker-independent). Set > 0 to
#'   suppress pure cell-cycle-phase clusters that score high on the
#'   snoRNA-host-gene-laden hypoxia module without genuine hypoxia markers.
#' @param max_dom_frac Direct phase-purity gate. When non-NULL and the object
#'   carries a `Phase` metadata column, a cluster is flagged only if it is
#'   phase-MIXED — i.e. its dominant cell-cycle-phase fraction is `< max_dom_frac`
#'   (default `NULL` = no phase gate). Set e.g. 0.70 so phase-restricted clusters
#'   (a single dominant G1/S/G2M phase) are NOT swept into the high-hypoxia
#'   subset even when they carry hypoxia markers (the snoRNA-host ride-along).
#'   Clusters with no usable Phase info (`dom_frac` NA) are treated permissively
#'   (the gate does not block them; the marker gate still governs).
#' @param return_detail If `TRUE`, return the full per-cluster decision table
#'   (one row per cluster) instead of just the flagged cluster ids.
#' @return By default a character vector of outlier cluster ids (possibly
#'   empty). With `return_detail = TRUE`, a data.frame with columns `cluster`,
#'   `n_cells`, `top_markers`, `n_hyp_markers`, `matched_genes`, `mean_score`,
#'   `dom_frac`, `outlier_fence`, `is_outlier`, `phase_mixed`, `flagged`. Never
#'   errors — returns an empty result of the requested shape on failure.
#' @export
identify_hypoxia_clusters <- function(seu, grp,
                                      hypoxia_genes = .hypoxia_marker_set,
                                      top_n = 20, top_markers_n = 5,
                                      mad_k = 3, min_hyp_markers = 0,
                                      max_dom_frac = NULL,
                                      return_detail = FALSE) {
  empty <- if (return_detail) {
    data.frame(cluster = character(0), n_cells = integer(0),
               top_markers = character(0), n_hyp_markers = integer(0),
               matched_genes = character(0), mean_score = numeric(0),
               dom_frac = numeric(0),
               outlier_fence = numeric(0), is_outlier = logical(0),
               phase_mixed = logical(0), flagged = logical(0),
               stringsAsFactors = FALSE)
  } else character(0)

  tryCatch({
    da       <- Seurat::DefaultAssay(seu)
    grpv     <- as.character(seu@meta.data[[grp]])
    clusters <- sort(unique(grpv))
    cl_n     <- as.integer(table(factor(grpv, levels = clusters)))

    # marker overlap + top-N markers per cluster (informational only)
    n_hyp   <- stats::setNames(rep(0L, length(clusters)), clusters)
    matched <- stats::setNames(rep("", length(clusters)), clusters)
    topmk   <- stats::setNames(rep("", length(clusters)), clusters)
    if (length(clusters) >= 2) {
      top <- presto::wilcoxauc(seu, grp, seurat_assay = da) |>
        dplyr::filter(.data$padj < 0.5) |>
        dplyr::group_by(.data$group) |>
        dplyr::arrange(dplyr::desc(.data$logFC), .by_group = TRUE) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::summarise(
          tm = paste(utils::head(.data$feature, top_markers_n), collapse = ","),
          n  = sum(.data$feature %in% hypoxia_genes),
          g  = paste(.data$feature[.data$feature %in% hypoxia_genes],
                     collapse = ","),
          .groups = "drop")
      idx <- as.character(top$group)
      topmk[idx]   <- top$tm
      n_hyp[idx]   <- as.integer(top$n)
      matched[idx] <- top$g
    }

    # per-cluster mean hypoxia_score
    if ("hypoxia_score" %in% colnames(seu@meta.data)) {
      sc      <- seu$hypoxia_score
      cl_mean <- tapply(sc, factor(grpv, levels = clusters), mean, na.rm = TRUE)
    } else {
      cl_mean <- stats::setNames(rep(NA_real_, length(clusters)), clusters)
    }

    # per-cluster dominant cell-cycle-phase fraction (for the direct phase gate).
    # NA when no Phase column or no non-NA phase calls in the cluster.
    if ("Phase" %in% colnames(seu@meta.data)) {
      phv      <- as.character(seu@meta.data[["Phase"]])
      dom_frac <- vapply(clusters, function(cl) {
        ph <- phv[grpv == cl]; ph <- ph[!is.na(ph)]
        if (length(ph) == 0) return(NA_real_)
        max(table(ph)) / length(ph)
      }, numeric(1))
    } else {
      dom_frac <- stats::setNames(rep(NA_real_, length(clusters)), clusters)
    }

    detail <- data.frame(
      cluster       = clusters,
      n_cells       = cl_n,
      top_markers   = unname(topmk[clusters]),
      n_hyp_markers = unname(n_hyp[clusters]),
      matched_genes = unname(matched[clusters]),
      mean_score    = unname(as.numeric(cl_mean[clusters])),
      dom_frac      = unname(as.numeric(dom_frac[clusters])),
      stringsAsFactors = FALSE)

    # Robust-z (median + mad_k * MAD) over the per-cluster means: flag high-side
    # outliers. Robust to several elevated clusters because median/MAD track the
    # bulk majority. mad() scales by 1.4826 (SD-consistent). A zero/NA MAD means
    # no usable spread (e.g. <= 2 clusters or many identical means) -> flag none.
    med   <- stats::median(detail$mean_score, na.rm = TRUE)
    md    <- stats::mad(detail$mean_score, na.rm = TRUE)
    fence <- as.numeric(med + mad_k * md)
    detail$outlier_fence <- fence
    detail$is_outlier    <- !is.na(detail$mean_score) & is.finite(fence) &
                            md > 0 & detail$mean_score > fence
    # Optional marker gate: with min_hyp_markers > 0 a cluster must ALSO carry at
    # least that many hypoxia marker genes in its top-N to be flagged. This
    # suppresses pure cell-cycle-phase clusters (e.g. histone/PLK2 S-phase or
    # CENPF/MKI67 G2M) whose mean hypoxia_score is inflated by the snoRNA-host
    # genes in the marker set but that are not genuinely hypoxic. Default 0
    # preserves the score-only rule (used by the phase-evaluation diagnostic).
    #
    # Direct phase gate (max_dom_frac): the marker gate cannot catch clusters
    # that carry hypoxia markers yet are still single-phase (the snoRNA-host
    # ride-along in proliferating cells). When max_dom_frac is set, a cluster is
    # flagged only if it is phase-MIXED (dominant-phase fraction < max_dom_frac).
    # Clusters with no usable Phase info (dom_frac NA) are treated permissively
    # so the gate never blocks on absent phase annotation.
    detail$phase_mixed <- if (is.null(max_dom_frac)) {
      rep(TRUE, nrow(detail))
    } else {
      is.na(detail$dom_frac) | detail$dom_frac < max_dom_frac
    }
    detail$flagged       <- detail$is_outlier &
                            detail$n_hyp_markers >= min_hyp_markers &
                            detail$phase_mixed

    if (return_detail) detail else as.character(detail$cluster[detail$flagged])
  }, error = function(e) {
    warning("identify_hypoxia_clusters failed for ", grp, ": ",
            conditionMessage(e))
    empty
  })
}

# QC PDF for one sample: UMAP by partition, UMAP by round, and a marker dotplot
# of the final partition. Best-effort — never aborts the split.
.write_hypoxia_split_qc <- function(seu, sample_id, out_dir = "results/hypoxia_cluster_split") {
  tryCatch({
    fs::dir_create(out_dir)
    pdf_path <- file.path(out_dir, glue::glue("{sample_id}_hypoxia_split.pdf"))
    has_umap <- "umap" %in% names(seu@reductions)
    panels <- list()
    if (has_umap) {
      panels[[length(panels) + 1]] <- Seurat::DimPlot(
        seu, group.by = "hypoxia_partition", reduction = "umap",
        cols = c(low = "grey70", high = "#d73027")) +
        ggplot2::labs(title = glue::glue("{sample_id} · partition"))
      seu$hypoxia_round_f <- factor(ifelse(is.na(seu$hypoxia_round), "low",
                                           paste0("round_", seu$hypoxia_round)))
      panels[[length(panels) + 1]] <- Seurat::DimPlot(
        seu, group.by = "hypoxia_round_f", reduction = "umap") +
        ggplot2::labs(title = glue::glue("{sample_id} · round peeled"))
    }
    if (exists("hypoxia_grid_marker_dotplot", mode = "function")) {
      panels[[length(panels) + 1]] <- hypoxia_grid_marker_dotplot(
        seu, "hypoxia_partition",
        glue::glue("{sample_id} · partition markers"))
    }
    panels <- purrr::compact(panels)
    if (length(panels) > 0) {
      grid <- patchwork::wrap_plots(panels, ncol = 1)
      ggplot2::ggsave(pdf_path, grid, height = 4.5 * length(panels), width = 8,
                      units = "in", limitsize = FALSE)
    }
    pdf_path
  }, error = function(e) {
    warning("hypoxia split QC failed for ", sample_id, ": ", conditionMessage(e))
    NA_character_
  })
}

#' Split a hypoxia-scored Seurat object into low/high subsets by clustering
#'
#' Clusters cells at several resolutions in a single round, flags clusters whose
#' mean `hypoxia_score` is a high-side outlier (see [identify_hypoxia_clusters()]),
#' and assigns their cells to the high subset. A cell is high if it lands in a
#' flagged cluster at ANY resolution in `resolutions` (hits unioned across
#' resolutions). After labelling, the established [subset_seu_by_expression()]
#' path produces the `*_hypoxia_low_seu.rds` / `*_hypoxia_high_seu.rds` outputs
#' (re-clustering, markers, UMAP, hash metadata, barcode DB) so every
#' downstream target is unchanged.
#'
#' @param hypoxia_seu_path Path to a `*_seu_hypoxia.rds` (has `hypoxia_score`,
#'   a `pca` reduction, and ideally a `umap`).
#' @param resolutions Louvain resolutions to cluster at (hits unioned across
#'   them; default `c(0.2, 0.6, 1.0)`).
#' @param n_iter Number of peeling rounds (default 1 — single round).
#' @param hypoxia_genes,top_n,top_markers_n,mad_k,min_hyp_markers,max_dom_frac
#'   Passed to [identify_hypoxia_clusters()]. `min_hyp_markers` defaults to 2
#'   here (marker gate on) and `max_dom_frac` to 0.70 (direct phase gate on),
#'   unlike the diagnostic defaults of 0 / NULL. Together they flag a cluster
#'   only when it is a score outlier, carries >= 2 hypoxia markers, AND is
#'   phase-mixed (dominant-phase fraction < 0.70). Phase is computed once up
#'   front via [annotate_cell_cycle_without_1q()] if not already present.
#' @param split_assay Assay used for the partition clustering (default "gene").
#' @param low_assay,high_assay Assays passed to [subset_seu_by_expression()]
#'   for the low/high subset re-clustering (defaults match the prior pipeline:
#'   gene for low, SCT for high).
#' @param write_qc Whether to emit a per-sample QC PDF.
#' @return Named character vector `c(low = ..., high = ...)` of output RDS
#'   paths, or `NA_character_` for a missing input.
#' @export
split_hypoxia_by_clusters <- function(hypoxia_seu_path,
                                      resolutions = c(0.2, 0.6, 1.0),
                                      n_iter = 1,
                                      hypoxia_genes = .hypoxia_marker_set,
                                      top_n = 20, top_markers_n = 5, mad_k = 3,
                                      min_hyp_markers = 2,
                                      max_dom_frac = 0.70,
                                      split_assay = "gene",
                                      low_assay = "gene", high_assay = "SCT",
                                      write_qc = TRUE) {
  if (is.null(hypoxia_seu_path) || length(hypoxia_seu_path) == 0 ||
      is.na(hypoxia_seu_path)) {
    return(NA_character_)
  }

  sample_id <- stringr::str_extract(hypoxia_seu_path, "SR[RX][0-9]+")
  seu <- readRDS(hypoxia_seu_path)

  da <- if (split_assay %in% names(seu@assays)) split_assay else Seurat::DefaultAssay(seu)
  Seurat::DefaultAssay(seu) <- da
  if (!"data" %in% SeuratObject::Layers(seu[[da]])) {
    seu <- Seurat::NormalizeData(seu, assay = da, verbose = FALSE)
  }

  # Cell-cycle Phase for the direct phase-purity gate. Compute once up front
  # (Phase is a per-cell property, invariant to clustering resolution). The
  # hypoxia_seu inputs do not reliably carry Phase, so annotate here when the
  # gate is on. Best-effort: if scoring fails the gate degrades to permissive
  # (identify_hypoxia_clusters treats NA dom_frac as not-blocking).
  if (!is.null(max_dom_frac) && !"Phase" %in% colnames(seu@meta.data)) {
    seu <- tryCatch(
      annotate_cell_cycle_without_1q(seu),
      error = function(e) {
        warning(sample_id, " cell-cycle scoring failed (phase gate off): ",
                conditionMessage(e))
        seu
      })
  }

  seu$hypoxia_round <- NA_integer_
  pool <- colnames(seu)  # current low pool
  log_rows <- list()     # per round x resolution x cluster decision trace

  if (!"pca" %in% names(seu@reductions)) {
    warning("No pca reduction in ", sample_id, "; cannot cluster — all cells -> low.")
  } else {
    for (round in seq_len(n_iter)) {
      if (length(pool) < 20) break
      sub <- seu[, pool]
      n_cells <- ncol(sub)
      k_param <- min(20L, n_cells - 1L)
      graph_names <- paste0(da, c("_nn", "_snn"))
      sub <- tryCatch(
        Seurat::FindNeighbors(sub, dims = 1:30, reduction = "pca",
                              graph.name = graph_names, k.param = k_param,
                              verbose = FALSE),
        error = function(e) { warning("FindNeighbors failed round ", round,
                                      " for ", sample_id, ": ",
                                      conditionMessage(e)); NULL })
      if (is.null(sub)) break

      hyp_cells <- character(0)
      for (res in resolutions) {
        grp <- glue::glue("hypsplit_r{round}_res.{res}")
        sub <- Seurat::FindClusters(sub, graph.name = paste0(da, "_snn"),
                                    resolution = res, verbose = FALSE)
        sub@meta.data[[grp]] <- sub$seurat_clusters
        det <- identify_hypoxia_clusters(
          sub, grp, hypoxia_genes = hypoxia_genes, top_n = top_n,
          top_markers_n = top_markers_n, mad_k = mad_k,
          min_hyp_markers = min_hyp_markers, max_dom_frac = max_dom_frac,
          return_detail = TRUE)
        if (is.data.frame(det) && nrow(det) > 0) {
          det <- cbind(sample_id = sample_id, round = round, resolution = res,
                       pool_n = n_cells, det, stringsAsFactors = FALSE)
          log_rows[[length(log_rows) + 1]] <- det
          hyp_cl <- as.character(det$cluster[det$flagged])
          if (length(hyp_cl) > 0) {
            cells_res <- colnames(sub)[as.character(sub@meta.data[[grp]]) %in% hyp_cl]
            hyp_cells <- union(hyp_cells, cells_res)
          }
        }
      }
      if (length(hyp_cells) == 0) break
      seu$hypoxia_round[colnames(seu) %in% hyp_cells] <- round
      pool <- setdiff(pool, hyp_cells)
    }
  }

  # Per-sample decision log: round x resolution x cluster, with the top-5 marker
  # genes, mean hypoxia_score, the robust-z (median+3*MAD) outlier_fence, and the
  # is_outlier flag.
  if (length(log_rows) > 0) {
    tryCatch({
      fs::dir_create("results/hypoxia_cluster_split")
      log_df <- dplyr::bind_rows(log_rows)
      readr::write_csv(
        log_df,
        glue::glue("results/hypoxia_cluster_split/{sample_id}_hypoxia_split_log.csv"))
    }, error = function(e) warning("Failed to write split log for ", sample_id,
                                   ": ", conditionMessage(e)))
  }

  seu$hypoxia_partition <- ifelse(is.na(seu$hypoxia_round), "low", "high")
  message(sample_id, " hypoxia partition: ",
          paste(names(table(seu$hypoxia_partition)),
                table(seu$hypoxia_partition), sep = "=", collapse = ", "))

  if (isTRUE(write_qc)) .write_hypoxia_split_qc(seu, sample_id)

  # Persist the labelled object, then reuse the established subset/recluster
  # path. str_replace("_seu.*.rds", ...) collapses "_seu_hypoxia_labeled.rds"
  # to the expected "_hypoxia_low_seu.rds" / "_hypoxia_high_seu.rds" outputs.
  labeled_path <- sub("_seu_hypoxia\\.rds$", "_seu_hypoxia_labeled.rds",
                      hypoxia_seu_path)
  add_hash_metadata(labeled_path, seu = seu)

  low_path <- tryCatch(
    subset_seu_by_expression(
      labeled_path, run_hypoxia_clustering = TRUE,
      hypoxia_expr = "hypoxia_partition == 'low'",
      slug = "hypoxia_low", assay = low_assay),
    error = function(e) {
      warning(sample_id, " low partition failed: ", conditionMessage(e))
      NA_character_
    })
  high_path <- tryCatch(
    subset_seu_by_expression(
      labeled_path, run_hypoxia_clustering = TRUE,
      hypoxia_expr = "hypoxia_partition == 'high'",
      slug = "hypoxia_high", assay = high_assay),
    error = function(e) {
      warning(sample_id, " high partition failed: ", conditionMessage(e))
      NA_character_
    })

  c(low = low_path, high = high_path)
}
