# All-clone SCNA marker-heatmap collages across a clustering-resolution sweep,
# with the cells removed during hypoxia splitting marked on the clustree (node
# colour = removed fraction) and on a dedicated UMAP panel.
#
# Sibling of plot_scna_two_clone_res_collages(): same per-resolution collage via
# plot_seu_marker_heatmap(), same guarded PCA-recompute + full sweep so the
# clustree spans every requested resolution -- but it keeps ALL clones (no
# two-clone subset, bar_var = "clone") and carries a per-cell "removed during
# hypoxia splitting" flag so the marking lands on the clustree and UMAP.
#
# The removed set is defined per CELL, taken from the sample's
# *_hypoxia_high_seu.rds barcodes (the high-hypoxia partition = the cells dropped
# by split_hypoxia_by_clusters). It is deliberately NOT re-derived from the
# split's own clusters: those live on a different object at a different resolution
# grid and do not map 1:1 to this object's clusters. Carrying the flag per barcode
# and aggregating to a node-level MEAN in the clustree is the only sound bridge.

#' All-clone SCNA collages across a resolution sweep, hypoxia-removed cells marked
#'
#' Builds one marker-heatmap + phase + UMAP + tree + clustree collage per
#' requested resolution on the FULL (all-clone) filtered object, and colours each
#' clustree node by the fraction of its cells that were removed during hypoxia
#' splitting. A companion UMAP panel shows the removed cells directly.
#'
#' @param seu_path Path to a filtered Seurat `.rds` (all clones).
#' @param scna_of_interest SCNA token used only in the output basename label.
#' @param resolutions Resolution vector (default `seq(0.2, 1.4, by = 0.2)`).
#' @param high_hypoxia_paths Character vector of `*_hypoxia_high_seu.rds` paths
#'   (e.g. the `seus_high_hypoxia` target). The one whose SRX id matches this
#'   sample supplies the removed-cell barcodes; if none matches, nothing is marked
#'   removed and the collage is still produced.
#' @param nb_paths,clone_simplifications Passed through to
#'   [plot_seu_marker_heatmap()] for the clone/segment tree panels.
#' @param assay Assay whose `<assay>_snn_res.<r>` columns are swept (default "SCT").
#' @param bar_signif,bar_signif_min_cells Passed through to
#'   [plot_seu_marker_heatmap()].
#' @param score_annotations Continuous per-cell metadata columns to draw as
#'   heatmap column annotations beside the cell-cycle scores, e.g.
#'   `c("hypoxia_score", "mito_score")`. Missing columns are computed via
#'   [.ensure_hypoxia_scores()]. `NULL` (default) keeps the cell-cycle-only
#'   annotation.
#' @return Character vector of written PDF paths, one per resolution (`NA` for a
#'   skipped resolution). Never errors -- best-effort, like the two-clone sibling.
#' @export
plot_scna_all_clone_res_collages <- function(seu_path,
                                             scna_of_interest,
                                             resolutions = seq(0.2, 1.4, by = 0.2),
                                             high_hypoxia_paths = NULL,
                                             nb_paths = NULL,
                                             clone_simplifications = NULL,
                                             assay = "SCT",
                                             bar_signif = TRUE,
                                             bar_signif_min_cells = 20,
                                             score_annotations = NULL) {
  if (is.null(seu_path) || length(seu_path) == 0 ||
      is.na(seu_path) || !file.exists(seu_path)) {
    return(NA_character_)
  }
  resolutions <- sort(unique(as.numeric(resolutions)))
  sample_id <- stringr::str_extract(seu_path, "SR[RX][0-9]+")

  seu <- readRDS(seu_path)
  if (assay %in% names(seu@assays)) Seurat::DefaultAssay(seu) <- assay
  if (!"clone_opt" %in% colnames(seu@meta.data)) {
    message(sample_id, ": no clone_opt column -> skip all-clone collages")
    return(NA_character_)
  }
  seu <- .ensure_cc_scores(seu)
  # Hypoxia / mitochondrial scores for the column annotations, computed here
  # because the *_filtered_seu.rds objects carry neither (only the hypoxia-split
  # objects have been through add_hypoxia_score()). Cached in cell_scores keyed
  # on seu_path -- the two-clone builder reads the same object, so whichever of
  # the two runs first pays for the scoring and the rest reuse it, on one scale.
  if (length(score_annotations) > 0)
    seu <- .ensure_hypoxia_scores(seu, cache_path = seu_path)

  # --- per-cell "removed during hypoxia splitting" flag (0/1) --------------
  # Match this sample's high-hypoxia object and mark its barcodes as removed.
  # Zero overlap between two non-empty barcode sets means a barcode-format
  # mismatch, not "nothing removed" -- warn loudly rather than mark all clean.
  seu@meta.data$hypoxia_removed <- 0
  hi <- NULL
  if (!is.null(high_hypoxia_paths)) {
    hi_paths <- unlist(high_hypoxia_paths, use.names = FALSE)
    hi_paths <- hi_paths[!is.na(hi_paths)]
    hi <- hi_paths[stringr::str_detect(hi_paths, stringr::fixed(sample_id))]
    hi <- hi[file.exists(hi)]
  }
  if (length(hi) >= 1) {
    removed_bc <- tryCatch(colnames(readRDS(hi[[1]])),
                           error = function(e) {
                             warning(sample_id, ": could not read high object ",
                                     hi[[1]], ": ", conditionMessage(e))
                             character(0)
                           })
    ov <- intersect(colnames(seu), removed_bc)
    if (length(removed_bc) > 0 && length(ov) == 0) {
      warning(sample_id, ": ", length(removed_bc), " removed barcodes but 0 ",
              "overlap this object's ", ncol(seu), " cells -- barcode mismatch; ",
              "clustree will show no removal.")
    }
    seu@meta.data[ov, "hypoxia_removed"] <- 1
    message(sample_id, " ", scna_of_interest, ": marked ", length(ov), "/",
            ncol(seu), " cells removed during hypoxia splitting")
  } else {
    message(sample_id, " ", scna_of_interest,
            ": no matching high-hypoxia object -> nothing marked removed")
  }

  # --- fresh PCA + SNN + full resolution sweep on ALL cells ----------------
  # Mirrors the guarded recompute in plot_scna_two_clone_res_collages(): cluster
  # every requested resolution here so the clustree spans the whole sweep, and
  # stash a clustree_res.* copy because each per-resolution collage below
  # overwrites SCT_snn_res.0.6 with its own resolution. On failure fall back to
  # whatever persisted columns exist rather than losing the collages.
  snn_name <- glue::glue("{assay}_snn")
  npcs     <- max(2L, min(30L, ncol(seu) - 1L))
  k_param  <- max(2L, min(20L, ncol(seu) - 1L))
  recomputed <- tryCatch({
    seu <- Seurat::RunPCA(seu, assay = assay, npcs = npcs, verbose = FALSE)
    seu <- Seurat::FindNeighbors(seu, dims = 1:npcs, reduction = "pca",
                                 k.param = k_param,
                                 graph.name = paste0(assay, c("_nn", "_snn")),
                                 verbose = FALSE)
    if (!"umap" %in% names(seu@reductions)) {
      seu <- Seurat::RunUMAP(seu, dims = 1:npcs, reduction = "pca",
                             verbose = FALSE)
    }
    for (r in resolutions) {
      seu <- Seurat::FindClusters(seu, graph.name = snn_name, resolution = r,
                                  verbose = FALSE)
      seu@meta.data[[glue::glue("{assay}_snn_res.{r}")]] <- seu$seurat_clusters
      seu@meta.data[[glue::glue("clustree_res.{r}")]] <- seu$seurat_clusters
    }
    message(sample_id, " ", scna_of_interest, ": recomputed PCA on ", ncol(seu),
            " cells; clustered resolutions ", paste(resolutions, collapse = " "))
    TRUE
  }, error = function(e) {
    message("!! ", sample_id, " ", scna_of_interest, ": PCA recompute failed (",
            conditionMessage(e), "); falling back to persisted clusters")
    FALSE
  })
  if (!isTRUE(recomputed)) seu <- .stash_clustree_sweep(seu, assay)

  # --- one collage per resolution ------------------------------------------
  purrr::map_chr(resolutions, function(res) {
    tryCatch({
      col <- glue::glue("{assay}_snn_res.{res}")
      s   <- seu
      if (!col %in% colnames(s@meta.data)) {
        message(sample_id, " ", scna_of_interest, ": no ", col,
                " -> skip res ", res)
        return(NA_character_)
      }
      # plot_seu_marker_heatmap keys cluster order / group.by on SCT_snn_res.0.6.
      s@meta.data[["SCT_snn_res.0.6"]] <- factor(as.character(s@meta.data[[col]]))

      tmp <- file.path(tempdir(), basename(seu_path))
      saveRDS(s, tmp)

      out <- plot_seu_marker_heatmap(
        tmp, cluster_order = NULL, nb_paths = nb_paths,
        clone_simplifications = clone_simplifications,
        bar_var = "clone",
        bar_signif = bar_signif, bar_signif_min_cells = bar_signif_min_cells,
        mark_removed_col = "hypoxia_removed",
        score_annotations = score_annotations,
        label = glue::glue("_allclone_{scna_of_interest}_res{res}_"))

      if (length(out) != 1 || is.na(out) || !file.exists(out)) {
        message("!! all-clone collage produced NO pdf (got '",
                paste(out, collapse = ", "), "')")
        return(NA_character_)
      }
      message("wrote all-clone scna collage ", out)
      out
    }, error = function(e) {
      message("!! all-clone collage failed for ", sample_id, " ",
              scna_of_interest, " res ", res, ": ", conditionMessage(e))
      NA_character_
    })
  })
}
