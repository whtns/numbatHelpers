# Cluster-based hypoxia split (lowest-resolution-first drop + confirmatory
# clean-up).
#
# Alternative to the per-cell hypoxia_score threshold used by
# subset_seu_by_expression(): cluster the cells and assign WHOLE clusters to the
# high-hypoxia subset when their mean hypoxia_score is a high-side statistical
# outlier among that resolution's per-cluster means.
#
# The sweep runs over `resolutions` from coarse to fine, logging every
# resolution, but the DROP is anchored at r_flag — the LOWEST resolution at which
# any cluster flags, i.e. the coarsest scale at which the hypoxic population
# first becomes separable. All cells of every cluster flagged at r_flag go high.
# The survivors are then re-clustered at r_flag + recluster_step (default +0.4)
# and any cluster still flagging there is dropped too — and this confirmatory
# clean-up ITERATES (r_flag + k*recluster_step, k = 1, 2, ...) until a pass flags
# nothing. One pass is not sufficient: residual hypoxia can stay merged until well
# above r_flag + recluster_step, so the loop runs to quiescence to keep the low
# subset (the SCNA/clone input) clean.
#
# This replaces the earlier rule, which unioned hits across a fixed resolution
# grid and was hard to reason about and prone to over-dropping.
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

#' Pick the anchor (lowest stable) resolution from a sweep's flagged-cell counts
#'
#' Shared by [split_hypoxia_by_clusters()] (which records the choice in the split
#' log) and by readers of that log, so the rule lives in exactly one place and
#' cannot drift between producer and consumer.
#'
#' Returns the index of the lowest resolution whose flagged-cell count agrees,
#' within `stability_tol`, with the next `stability_window` flagging resolutions.
#' A single successor is not enough: a monotonically rising count can contain a
#' spurious one-step flat spot (SRX22868102 goes 385, 480, 437, 434, 538, 566 —
#' 0.6/0.8 agree to 0.7% while the population is still growing), and anchoring
#' there under-drops by ~20%. If nothing stabilises, falls back to the
#' most-cells index (with attribute `fallback = TRUE`) so the split never
#' knowingly under-drops.
#'
#' @param counts Integer vector of flagged-cell counts, one per resolution,
#'   ordered coarse -> fine. Zeros mean that resolution flagged nothing.
#' @param stability_tol Relative tolerance (default 0.05).
#' @param stability_window Number of successive resolutions the anchor must agree
#'   with (default 2).
#' @return Integer index into `counts`, or `NA_integer_` if nothing flagged.
#'   Carries attribute `fallback = TRUE` when no resolution stabilised.
#' @export
pick_hypoxia_anchor <- function(counts, stability_tol = 0.05,
                                stability_window = 2) {
  counts   <- as.integer(counts)
  idx_flag <- which(counts > 0)
  if (length(idx_flag) == 0) return(NA_integer_)

  for (i in idx_flag) {
    js <- (i + 1L):(i + stability_window)
    js <- js[js <= length(counts)]
    if (length(js) == 0 || any(counts[js] == 0)) next
    ok <- vapply(js, function(j) {
      abs(counts[j] - counts[i]) / max(counts[i], counts[j]) <= stability_tol
    }, logical(1))
    if (all(ok)) return(as.integer(i))
  }

  pick <- as.integer(idx_flag[which.max(counts[idx_flag])])
  attr(pick, "fallback") <- TRUE
  pick
}

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

# One marker-heatmap/phase-scatter collage for an intermediate state of the split
# (a clustering that only exists inside split_hypoxia_by_clusters). Best-effort —
# never aborts the split.
#
# plot_seu_marker_heatmap() takes a path, and both it and dummy_cluster_order()
# key on the "SCT_snn_res.0.6" column, so the split clustering is copied into
# that column and the object is written to a temp RDS named `slug`. Because
# `slug` does not contain "_filtered_seu", the plotter keeps the whole name and
# emits results/{slug}__filtered_heatmap_phase_scatter_patchwork.pdf.
.write_hypoxia_stage_heatmap <- function(seu, group_col, slug) {
  out <- tryCatch({
    seu@meta.data[["SCT_snn_res.0.6"]] <- factor(
      as.character(seu@meta.data[[group_col]]))
    # plot_seu_marker_heatmap() computes markers on SCT but filters them by
    # VariableFeatures(seu), which reads the DEFAULT assay. The split object's
    # default is "gene", so the two disagree and the marker set can come out empty
    # ("No marker heatmap features" -> collage silently skipped, as on
    # SRX22868102). Point the default assay at SCT so they match.
    if ("SCT" %in% names(seu@assays)) {
      Seurat::DefaultAssay(seu) <- "SCT"
      if (length(Seurat::VariableFeatures(seu)) == 0) {
        seu <- Seurat::FindVariableFeatures(seu, assay = "SCT", verbose = FALSE)
      }
    }
    tmp_path <- file.path(tempdir(), slug)
    saveRDS(seu, tmp_path)
    plot_seu_marker_heatmap(tmp_path, cluster_order = NULL, nb_paths = NULL,
                            label = "_filtered_")
  }, error = function(e) {
    # message(), not warning(): Rscript defers warnings to the end of the script,
    # which hides a failure here until long after (or never, if the run is killed)
    message("!! hypoxia stage heatmap FAILED for ", slug, ": ",
            conditionMessage(e))
    NA_character_
  })
  expected <- glue::glue("results/{slug}__filtered_heatmap_phase_scatter_patchwork.pdf")
  if (!file.exists(expected)) {
    message("!! hypoxia stage heatmap produced NO pdf at ", expected)
  } else {
    message("wrote stage heatmap ", expected)
  }
  out
}

#' Split a hypoxia-scored Seurat object into low/high subsets by clustering
#'
#' Sweeps `resolutions` coarse-to-fine, flagging clusters whose mean
#' `hypoxia_score` is a high-side outlier (see [identify_hypoxia_clusters()]).
#' Every resolution is logged, but the drop is anchored at `r_flag` — the LOWEST
#' resolution at which anything flags. All cells of every cluster flagged at
#' `r_flag` go high. The survivors are then re-clustered at
#' `r_flag + recluster_step` and any cluster still flagging there is dropped too;
#' this confirmatory clean-up **iterates** at `r_flag + k * recluster_step`,
#' k = 1, 2, …, until a pass flags nothing, so residual hypoxia that only becomes
#' separable at a much finer scale is still caught. `hypoxia_round` records which
#' pass took a cell (1 = primary drop at `r_flag`; 2, 3, … = the confirmatory
#' passes, so the log shows which resolution finally separated each pocket).
#'
#' After labelling, the established [subset_seu_by_expression()] path produces
#' the `*_hypoxia_low_seu.rds` / `*_hypoxia_high_seu.rds` outputs (re-clustering
#' over its own full resolution grid, markers, UMAP, hash metadata, barcode DB)
#' so every downstream target is unchanged. The `r_flag + recluster_step`
#' clustering is diagnostic only and is not what the persisted low object
#' carries.
#'
#' @param hypoxia_seu_path Path to a `*_seu_hypoxia.rds` (has `hypoxia_score`,
#'   a `pca` reduction, and ideally a `umap`).
#' @param resolutions Louvain resolutions to sweep, coarse to fine (default
#'   `seq(0.2, 1.2, by = 0.2)`). All are logged; the drop happens at the lowest
#'   *stable* one (`r_flag`, see `stability_tol`).
#' @param stability_tol Relative tolerance (default 0.05) used to pick `r_flag`:
#'   the anchor is the lowest resolution whose flagged-cell count is within this
#'   fraction of the next `stability_window` flagging resolutions'. This rejects a
#'   too-coarse resolution that *detects* the hypoxic population but under-resolves
#'   it — capturing only part of it and merging the rest into neighbouring clusters
#'   (SRX10264518: 507 cells at res 0.2 vs a stable ~660 from 0.4 up). If no
#'   resolution stabilises, falls back to the one flagging the most cells so the
#'   split never knowingly under-drops.
#' @param stability_window How many successive resolutions the anchor's count must
#'   agree with (default 2). One is not enough: a monotonically rising count can
#'   contain a spurious one-step flat spot (SRX22868102 goes 385, 480, 437, 434,
#'   538, 566 — 0.6 and 0.8 agree to 0.7% but the population is still growing), and
#'   anchoring there under-drops by ~20%. Requiring two successive agreements
#'   rejects that flat spot and advances to the real plateau.
#' @param recluster_step Resolution increment per confirmatory pass; pass `k`
#'   re-clusters the survivors at `r_flag + k * recluster_step` (default 0.4).
#' @param confirmatory_drop Whether clusters still flagging in a confirmatory pass
#'   are also assigned to the high subset (default `TRUE`, and the passes then
#'   iterate to convergence). Set `FALSE` to run a single pass that records the
#'   residual in the log without dropping it.
#' @param max_confirmatory_rounds Safety cap on the number of confirmatory drop
#'   passes (default 5). Hitting it warns — the low subset may still hold hypoxia.
#' @param max_recluster_res Ceiling for the confirmatory ladder: stop once
#'   `r_flag + k * recluster_step` would exceed it. Defaults to `max(resolutions)`
#'   so the confirmatory passes never explore beyond the swept range — above it,
#'   clusters fragment and the MAD fence gets easier to clear, so an unbounded
#'   ladder manufactures drops. Note this means a sample whose `r_flag` is already
#'   the top resolution gets the primary drop only.
#' @param write_diagnostics Whether to emit the three per-sample stage heatmap
#'   PDFs (at `r_flag` before the drop, at `r_flag` after it, and at
#'   `r_flag + recluster_step`).
#' @param n_iter Deprecated, retained for call-site back-compatibility; the split
#'   is now always a single pass plus the confirmatory clean-up.
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
                                      resolutions = seq(0.2, 1.2, by = 0.2),
                                      stability_tol = 0.05,
                                      stability_window = 2,
                                      recluster_step = 0.4,
                                      confirmatory_drop = TRUE,
                                      max_confirmatory_rounds = 5,
                                      max_recluster_res = NULL,
                                      write_diagnostics = TRUE,
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
  log_rows <- list()      # per resolution x cluster decision trace
  resolutions   <- sort(unique(as.numeric(resolutions)))
  # The confirmatory ladder must not explore beyond the swept range: at high
  # resolution clusters fragment and the MAD fence gets easier to clear, so an
  # unbounded ladder manufactures drops at resolutions never asked for
  # (SRX10264520 flags 74 cells at 1.2, then another 105 at 1.6 — outside the
  # sweep). Default the ceiling to the top of `resolutions`.
  max_recluster_res <- max_recluster_res %||% max(resolutions)
  r_flag        <- NA_real_       # lowest resolution that flags anything
  r_recluster   <- NA_real_       # r_flag + recluster_step
  flagged_grp   <- NULL           # cluster column at r_flag
  recluster_col <- NULL           # cluster column at r_recluster
  hyp_cells_1   <- character(0)   # primary drop (clusters flagged at r_flag)
  hyp_cells_2   <- character(0)   # all confirmatory drops, unioned
  confirm_rounds <- list()        # round index -> cells taken by that pass
  sub           <- NULL
  retained_sub  <- NULL

  graph_names <- paste0(da, c("_nn", "_snn"))
  snn_name    <- paste0(da, "_snn")

  if (!"pca" %in% names(seu@reductions)) {
    warning("No pca reduction in ", sample_id, "; cannot cluster — all cells -> low.")
  } else if (ncol(seu) < 20) {
    warning("Fewer than 20 cells in ", sample_id, "; skipping split — all cells -> low.")
  } else {
    n_cells <- ncol(seu)
    sub <- tryCatch(
      Seurat::FindNeighbors(seu, dims = 1:30, reduction = "pca",
                            graph.name = graph_names,
                            k.param = min(20L, n_cells - 1L), verbose = FALSE),
      error = function(e) { warning("FindNeighbors failed for ", sample_id, ": ",
                                    conditionMessage(e)); NULL })

    # Sweep coarse -> fine, recording the flagged cell set at EVERY resolution
    # (all are logged, so the full picture stays inspectable).
    sweep <- list()
    if (!is.null(sub)) {
      for (res in resolutions) {
        grp <- glue::glue("hypsplit_res.{res}")
        sub <- Seurat::FindClusters(sub, graph.name = snn_name,
                                    resolution = res, verbose = FALSE)
        sub@meta.data[[grp]] <- sub$seurat_clusters
        det <- identify_hypoxia_clusters(
          sub, grp, hypoxia_genes = hypoxia_genes, top_n = top_n,
          top_markers_n = top_markers_n, mad_k = mad_k,
          min_hyp_markers = min_hyp_markers, max_dom_frac = max_dom_frac,
          return_detail = TRUE)
        if (!is.data.frame(det) || nrow(det) == 0) next
        log_rows[[length(log_rows) + 1]] <- cbind(
          sample_id = sample_id, round = 1L, resolution = res,
          pool_n = n_cells, det, stringsAsFactors = FALSE)

        hyp_cl <- as.character(det$cluster[det$flagged])
        cells  <- if (length(hyp_cl) > 0) {
          colnames(sub)[as.character(sub@meta.data[[grp]]) %in% hyp_cl]
        } else character(0)
        sweep[[length(sweep) + 1]] <- list(res = res, grp = grp, cells = cells)
      }
    }

    # Anchor: the lowest STABLE resolution — the coarsest scale at which the
    # hypoxic population is actually RESOLVED, not merely first detected.
    #
    # "First resolution that flags anything" is not enough. A too-coarse
    # clustering can detect the hypoxic population while under-resolving it,
    # capturing only part of it and merging the rest into neighbouring clusters.
    # SRX10264518 is the worked example: res 0.2 flags a 507-cell cluster, but
    # from 0.4 up the same population reads a stable ~660 cells — so anchoring at
    # 0.2 would leave ~24% of it in the low subset. The confirmatory pass cannot
    # rescue those cells either: once the 507-cell core is removed, the remainder
    # is diluted and no longer forms an outlier cluster (verified — the pass
    # flagged nothing).
    #
    # So require the anchor to AGREE with the next resolution up: pick the lowest
    # resolution whose flagged-cell count is within stability_tol of the next
    # flagging resolution's. For SRX10264518 that selects 0.4 (663 vs 667 at 0.6,
    # 0.6% apart). Where the count is already stable at the first flagging
    # resolution (15 of 18 flagging samples in this cohort) the anchor is
    # unchanged. If nothing ever stabilises (e.g. only the finest resolution
    # flags), fall back to the resolution flagging the MOST cells, so we never
    # knowingly under-drop.
    counts <- vapply(sweep, function(x) length(x$cells), integer(1))
    pick   <- pick_hypoxia_anchor(counts, stability_tol = stability_tol,
                                  stability_window = stability_window)
    if (!is.na(pick) && isTRUE(attr(pick, "fallback"))) {
      warning(sample_id, " no resolution stabilised within stability_tol (",
              stability_tol, "); falling back to the max-extent resolution (",
              sweep[[pick]]$res, ", ", counts[pick], " cells).")
    }
    if (!is.na(pick)) {
      r_flag      <- sweep[[pick]]$res
      flagged_grp <- sweep[[pick]]$grp
      hyp_cells_1 <- sweep[[pick]]$cells
      message(sample_id, " anchor resolution ", r_flag, " (", counts[pick],
              " cells); sweep counts: ",
              paste(vapply(sweep, function(x) x$res, numeric(1)), counts,
                    sep = "=", collapse = ", "))
    }

    # Confirmatory clean-up, ITERATED TO CONVERGENCE: re-cluster the survivors one
    # step finer (r_flag + k * recluster_step), drop anything that still flags, and
    # repeat at the next step up until a pass flags NOTHING. A single pass is not
    # enough — residual hypoxia may only become separable well above
    # r_flag + recluster_step (e.g. SRX10264518 flags at 0.2, is clean at 0.6, but
    # carries a further hypoxic pocket that only separates at 1.0). Iterating until
    # quiescence is what keeps the low subset — the SCNA/clone input — clean.
    #
    # Terminates on the first no-flag pass, or when the pool drops below 20 cells,
    # the resolution passes max_recluster_res, or max_confirmatory_rounds is hit.
    if (!is.na(r_flag)) {
      if (isTRUE(write_diagnostics)) {
        .write_hypoxia_stage_heatmap(
          sub, flagged_grp,
          glue::glue("{sample_id}_hypoxia_rflag{r_flag}_before_seu.rds"))
      }
      retained <- setdiff(colnames(sub), hyp_cells_1)
      if (isTRUE(write_diagnostics) && length(retained) >= 20) {
        .write_hypoxia_stage_heatmap(
          sub[, retained], flagged_grp,
          glue::glue("{sample_id}_hypoxia_rflag{r_flag}_after_seu.rds"))
      }

      k <- 1L
      repeat {
        if (length(retained) < 20) {
          message(sample_id, " confirmatory stop: fewer than 20 cells left.")
          break
        }
        if (k > max_confirmatory_rounds) {
          warning(sample_id, " confirmatory pass hit max_confirmatory_rounds (",
                  max_confirmatory_rounds, ") without converging; residual ",
                  "hypoxia may remain in the low subset.")
          break
        }
        r_recluster <- round(r_flag + k * recluster_step, 1)  # dodge float noise
        if (r_recluster > max_recluster_res) {
          message(sample_id, " confirmatory stop: next resolution (", r_recluster,
                  ") exceeds max_recluster_res (", max_recluster_res, ").")
          break
        }
        recluster_col <- glue::glue("hypsplit_recluster_res.{r_recluster}")

        retained_sub <- tryCatch({
          rs <- sub[, retained]
          # the graph is stale once cells are removed — recompute each pass
          rs <- Seurat::FindNeighbors(rs, dims = 1:30, reduction = "pca",
                                      graph.name = graph_names,
                                      k.param = min(20L, ncol(rs) - 1L),
                                      verbose = FALSE)
          rs <- Seurat::FindClusters(rs, graph.name = snn_name,
                                     resolution = r_recluster, verbose = FALSE)
          rs@meta.data[[recluster_col]] <- rs$seurat_clusters
          rs
        }, error = function(e) {
          warning("Confirmatory recluster failed for ", sample_id, " at res ",
                  r_recluster, ": ", conditionMessage(e))
          NULL
        })
        if (is.null(retained_sub)) break

        det2 <- identify_hypoxia_clusters(
          retained_sub, recluster_col, hypoxia_genes = hypoxia_genes,
          top_n = top_n, top_markers_n = top_markers_n, mad_k = mad_k,
          min_hyp_markers = min_hyp_markers, max_dom_frac = max_dom_frac,
          return_detail = TRUE)
        if (!is.data.frame(det2) || nrow(det2) == 0) break

        log_rows[[length(log_rows) + 1]] <- cbind(
          sample_id = sample_id, round = k + 1L, resolution = r_recluster,
          pool_n = ncol(retained_sub), det2, stringsAsFactors = FALSE)

        # Collage of this pass's clustering, with any flagged cluster still
        # visible — i.e. the state the plot shows is BEFORE this round's drop.
        if (isTRUE(write_diagnostics)) {
          .write_hypoxia_stage_heatmap(
            retained_sub, recluster_col,
            glue::glue("{sample_id}_hypoxia_recluster{r_recluster}_seu.rds"))
        }

        hyp_cl2 <- as.character(det2$cluster[det2$flagged])
        if (length(hyp_cl2) == 0) {
          message(sample_id, " confirmatory converged at res ", r_recluster,
                  " (nothing flagged) after ", k - 1L, " drop round(s).")
          break
        }

        cells_k <- colnames(retained_sub)[
          as.character(retained_sub@meta.data[[recluster_col]]) %in% hyp_cl2]

        if (!isTRUE(confirmatory_drop)) {
          # log-only mode: record the residual, take nothing, stop.
          hyp_cells_2 <- cells_k
          message(sample_id, " residual hypoxia at res ", r_recluster, ": ",
                  length(cells_k), " cells flagged but NOT dropped ",
                  "(confirmatory_drop = FALSE); see round = 2 rows in the log.")
          break
        }

        confirm_rounds[[as.character(k + 1L)]] <- cells_k
        hyp_cells_2 <- union(hyp_cells_2, cells_k)
        retained    <- setdiff(retained, cells_k)
        k <- k + 1L
      }
    }
  }

  seu$hypoxia_round[colnames(seu) %in% hyp_cells_1] <- 1L
  if (isTRUE(confirmatory_drop)) {
    # each confirmatory pass gets its own round index (2, 3, ...) so the log shows
    # which resolution finally separated each pocket
    for (rn in names(confirm_rounds)) {
      seu$hypoxia_round[colnames(seu) %in% confirm_rounds[[rn]]] <- as.integer(rn)
    }
  }

  # Per-sample decision log: round x resolution x cluster, with the top-5 marker
  # genes, mean hypoxia_score, the robust-z (median+3*MAD) outlier_fence, and the
  # is_outlier flag. round = 1 rows are the full coarse-to-fine sweep; round = 2
  # rows are the confirmatory pass on the survivors at r_flag + recluster_step.
  if (length(log_rows) > 0) {
    tryCatch({
      fs::dir_create("results/hypoxia_cluster_split")
      log_df <- dplyr::bind_rows(log_rows)
      # Record the chosen anchor on every row so downstream consumers (e.g.
      # plot_hypoxia_low_res_collages()) can read it straight off the log instead
      # of re-deriving the stability rule and risking drifting out of sync with it.
      # NA = nothing flagged for this sample.
      log_df$anchor_resolution <- r_flag
      readr::write_csv(
        log_df,
        glue::glue("results/hypoxia_cluster_split/{sample_id}_hypoxia_split_log.csv"))
    }, error = function(e) warning("Failed to write split log for ", sample_id,
                                   ": ", conditionMessage(e)))
  }

  seu$hypoxia_partition <- ifelse(is.na(seu$hypoxia_round), "low", "high")
  message(sample_id, " hypoxia partition (r_flag=",
          if (is.na(r_flag)) "none" else r_flag, ", primary=",
          length(hyp_cells_1), " cells, confirmatory=", length(hyp_cells_2),
          " cells over ", length(confirm_rounds), " pass(es)): ",
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
      slug = "hypoxia_low", assay = low_assay,
      # recompute the PCA on the low-hypoxia cells so the persisted clusters
      # (SCT_snn_res.* read by the collages, and gene_snn_res.*) are defined on
      # the retained population, not the inherited pre-split PCA
      recompute_pca = TRUE),
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
