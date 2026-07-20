# Marker-heatmap collages of the LOW-hypoxia object at two clustering resolutions.
#
# The standard low-hypoxia collage (heatmap_collages_low_hypoxia ->
# plot_seu_marker_heatmap) is pinned to SCT_snn_res.0.6 for every sample, because
# dummy_cluster_order() hardcodes that column. That is fine as a fixed reference
# but it does not show the clustering the hypoxia split actually acted on.
#
# This emits TWO collages per sample, on the persisted *_hypoxia_low_seu.rds:
#   1. the anchor resolution r_flag -- the resolution whose clusters were dropped
#   2. r_flag + step (default +0.4) -- the same cells one step finer
# When the sample flagged nothing (no clusters removed, r_flag = NA) it falls back
# to a fixed 0.6 / 1.0 pair so every sample still gets a comparable pair.

#' Two-resolution marker-heatmap collages for a low-hypoxia object
#'
#' Emits the `*_heatmap_phase_scatter_patchwork.pdf` collage for a
#' `*_hypoxia_low_seu.rds` at two clustering resolutions: the hypoxia split's
#' anchor resolution (`r_flag`, the one whose clusters were dropped) and
#' `r_flag + step`. If the sample flagged nothing, defaults to `fallback` (0.6 and
#' 1.0), so every sample yields a comparable pair.
#'
#' The anchor is read from the split log's `anchor_resolution` column (written by
#' [split_hypoxia_by_clusters()]) rather than re-derived, so this cannot drift out
#' of sync with the split's own stability rule.
#'
#' Output paths follow [plot_seu_marker_heatmap()]'s naming with the resolution in
#' the label segment, e.g.
#' `results/SRX11133592_hypoxia_low_seu.rds__filtered_res0.2_heatmap_phase_scatter_patchwork.pdf`.
#'
#' @param low_seu_path Path to a `*_hypoxia_low_seu.rds`.
#' @param nb_paths Numbat RDS paths (for the clone tree panel); pass `NULL` to
#'   omit the clone tree.
#' @param clone_simplifications Passed through to [plot_seu_marker_heatmap()].
#' @param split_log_dir Directory holding `<SRX>_hypoxia_split_log.csv`.
#' @param step Resolution increment above the anchor (default 0.4).
#' @param fallback Resolutions used when the sample flagged nothing
#'   (default `c(0.6, 1.0)`).
#' @param resolutions Optional explicit resolution vector. When supplied, the
#'   anchor logic is bypassed entirely and a collage is emitted at each of these
#'   resolutions (clustering on demand when the `<assay>_snn_res.<r>` column is
#'   absent). Use e.g. `seq(0.2, 1.2, 0.2)` for a full sweep. Default `NULL`
#'   keeps the anchor + `step` (or `fallback`) behaviour.
#' @param assay Assay whose `<assay>_snn_res.<r>` columns are used (default "SCT").
#' @param recompute_pca When `TRUE`, recompute a fresh PCA on the (already
#'   hypoxia-subsetted) low object and re-cluster EVERY requested resolution on
#'   it, instead of using the persisted `<assay>_snn_res.<r>` columns. Those
#'   persisted columns were clustered in the PCA inherited from the full pre-split
#'   object (PCA is not recomputed by `subset_seu_by_expression()`), so the
#'   surviving cells are clustered on axes partly defined by the removed
#'   high-hypoxia cells; a fresh PCA re-orients the axes to the retained
#'   population and yields cleaner clusters at higher resolution. Falls back to
#'   the persisted-column path if the PCA recompute fails. Default `FALSE`.
#' @return Character vector of the PDF paths written (`NA_character_` per failed
#'   resolution). Never errors — collages are best-effort.
#' @export
plot_hypoxia_low_res_collages <- function(low_seu_path,
                                          nb_paths = NULL,
                                          clone_simplifications = NULL,
                                          split_log_dir = "results/hypoxia_cluster_split",
                                          step = 0.4,
                                          fallback = c(0.6, 1.0),
                                          resolutions = NULL,
                                          assay = "SCT",
                                          recompute_pca = FALSE) {
  if (is.null(low_seu_path) || length(low_seu_path) == 0 ||
      is.na(low_seu_path) || !file.exists(low_seu_path)) {
    return(NA_character_)
  }
  sample_id <- stringr::str_extract(low_seu_path, "SR[RX][0-9]+")

  # Explicit resolution override: skip the split-log anchor lookup entirely.
  resolutions_override <- if (!is.null(resolutions)) {
    as.numeric(resolutions)
  } else NULL

  # Anchor from the split log. Prefer the recorded anchor_resolution column; for
  # logs written before that column existed, recompute it from the round-1 sweep
  # with pick_hypoxia_anchor() — the SAME function split_hypoxia_by_clusters()
  # uses, so the two cannot disagree. NA = the sample flagged nothing.
  anchor <- tryCatch({
    lg_path <- file.path(split_log_dir,
                         glue::glue("{sample_id}_hypoxia_split_log.csv"))
    if (!file.exists(lg_path)) NA_real_ else {
      lg <- readr::read_csv(lg_path, show_col_types = FALSE)
      if ("anchor_resolution" %in% names(lg)) {
        suppressWarnings(as.numeric(lg$anchor_resolution[1]))
      } else {
        sweep <- lg |>
          dplyr::filter(.data$round == 1) |>
          dplyr::group_by(.data$resolution) |>
          dplyr::summarise(cells = sum(.data$n_cells[.data$flagged]),
                           .groups = "drop") |>
          dplyr::arrange(.data$resolution)
        pick <- pick_hypoxia_anchor(sweep$cells)
        if (is.na(pick)) NA_real_ else as.numeric(sweep$resolution[pick])
      }
    }
  }, error = function(e) NA_real_)

  resolutions <- if (!is.null(resolutions_override)) {
    resolutions_override
  } else if (is.na(anchor)) {
    as.numeric(fallback)
  } else {
    c(anchor, round(anchor + step, 1))  # round: 0.2 + 0.4 = 0.6000000000000001
  }
  message(sample_id, " low-hypoxia collages at res ",
          paste(resolutions, collapse = " and "),
          if (!is.null(resolutions_override)) " (explicit resolutions)" else
          if (is.na(anchor)) " (no clusters dropped -> fallback pair)" else
            glue::glue(" (anchor {anchor})"))

  seu <- readRDS(low_seu_path)
  if (assay %in% names(seu@assays)) Seurat::DefaultAssay(seu) <- assay

  # Recompute the PCA on the surviving (low-hypoxia) cells and rebuild the SNN
  # graph on it, so every resolution below is clustered in an embedding defined by
  # the retained population rather than the inherited pre-split PCA. Guarded: on
  # any failure fall back to the persisted-column / inherited-PCA path.
  snn_name <- glue::glue("{assay}_snn")
  if (isTRUE(recompute_pca)) {
    recompute_pca <- tryCatch({
      seu <- Seurat::RunPCA(seu, assay = assay, npcs = 30, verbose = FALSE)
      seu <- Seurat::FindNeighbors(seu, dims = 1:30, reduction = "pca",
                                   graph.name = paste0(assay, c("_nn", "_snn")),
                                   verbose = FALSE)
      message(sample_id, ": recomputed PCA on ", ncol(seu),
              " low-hypoxia cells; re-clustering all resolutions on it")
      TRUE
    }, error = function(e) {
      message("!! ", sample_id, ": PCA recompute failed (", conditionMessage(e),
              "); falling back to persisted clusters")
      FALSE
    })
  }

  purrr::map_chr(resolutions, function(res) {
    tryCatch({
      col <- glue::glue("{assay}_snn_res.{res}")
      s   <- seu
      if (isTRUE(recompute_pca)) {
        # Fresh PCA + graph already on `seu`: cluster this resolution on it,
        # ignoring any inherited persisted column (which used the old PCA).
        s <- Seurat::FindClusters(s, graph.name = snn_name, resolution = res,
                                  verbose = FALSE)
        s@meta.data[[col]] <- s$seurat_clusters
      } else if (!col %in% colnames(s@meta.data)) {
        # subset_seu_by_expression() clusters the low object over seq(0.2, 1, 0.2),
        # so a resolution above 1.0 has no column -- cluster on demand rather than
        # silently falling back to a resolution we did not ask for.
        if (!snn_name %in% names(s@graphs)) {
          s <- Seurat::FindNeighbors(s, dims = 1:30, reduction = "pca",
                                     graph.name = paste0(assay, c("_nn", "_snn")),
                                     verbose = FALSE)
        }
        s <- Seurat::FindClusters(s, graph.name = snn_name, resolution = res,
                                  verbose = FALSE)
        s@meta.data[[col]] <- s$seurat_clusters
      }
      # dummy_cluster_order() and group.by both key on SCT_snn_res.0.6, so the
      # requested resolution has to be copied into that column (same trick as
      # .write_hypoxia_stage_heatmap).
      s@meta.data[["SCT_snn_res.0.6"]] <- factor(as.character(s@meta.data[[col]]))

      # Temp copy keeps the *_hypoxia_low_seu.rds basename: plot_seu_marker_heatmap
      # derives the output name from it, and the resolution rides in `label`.
      tmp <- file.path(tempdir(), basename(low_seu_path))
      saveRDS(s, tmp)

      out <- plot_seu_marker_heatmap(
        tmp, cluster_order = NULL, nb_paths = nb_paths,
        clone_simplifications = clone_simplifications,
        label = glue::glue("_filtered_res{res}_"))

      expected <- glue::glue(
        "results/{basename(low_seu_path)}__filtered_res{res}_heatmap_phase_scatter_patchwork.pdf")
      if (!file.exists(expected)) {
        message("!! low-hypoxia collage produced NO pdf at ", expected)
        return(NA_character_)
      }
      message("wrote low-hypoxia collage ", expected)
      as.character(expected)
    }, error = function(e) {
      message("!! low-hypoxia collage FAILED for ", sample_id, " res ", res, ": ",
              conditionMessage(e))
      NA_character_
    })
  })
}
