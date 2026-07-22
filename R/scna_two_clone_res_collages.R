# Two-clone (acquiring vs immediately-preceding) marker-heatmap collages of an
# SCNA subclone across the full clustering-resolution sweep.
#
# For a given SCNA of interest (e.g. "1q", "2p", "16q"), large_clone_comparisons
# names the comparison `<N>_v_<M>_<scna>` where numbat clone N ACQUIRED the SCNA
# and clone M is the immediately preceding clone WITHOUT it. This restricts the
# object to exactly those two clones and emits one marker-heatmap + phase-scatter
# collage per requested resolution -- the same plot as
# plot_hypoxia_low_res_collages / plot_seu_marker_heatmap, but showing only the
# two compared clones so the SCNA effect can be read against its own background.
#
# Design mirrors plot_hypoxia_low_res_collages(): reuse the proven
# plot_seu_marker_heatmap() (which already displays clones via clone_opt, drops
# the clone annotation gracefully when only one clone survives filtering, and
# guards degenerate subsets) rather than plot_seu_marker_heatmap_by_scna(), which
# assumes >= 2 surviving `scna` levels and dies "subscript out of bounds" on
# degenerate two-clone subsets.

#' Two-clone SCNA marker-heatmap collages across a resolution sweep
#'
#' Restricts a Seurat object to the two clones of an SCNA comparison from
#' `large_clone_comparisons` (the clone that acquired `scna_of_interest` and its
#' immediately preceding clone) and emits a marker-heatmap + phase-scatter collage
#' at each requested resolution.
#'
#' The object is restricted to the two compared clones FIRST, and the PCA, SNN
#' graph, and clusters are then recomputed on that subset -- so the clustering
#' reflects structure among the cells actually under analysis rather than
#' structure of the whole low-hypoxia population. Cluster order, phase labels, and
#' marker-gene rows are likewise derived on the subset inside
#' [plot_seu_marker_heatmap()].
#'
#' Consequence: cluster IDs are NOT comparable with the full-population collage,
#' nor across SCNAs for the same sample -- each two-clone collage is its own
#' re-analysis. (This replaced an earlier "lock to full" behaviour, where the
#' persisted `<assay>_snn_res.<r>` columns were read as-is and the panels were
#' merely restricted to the two clones at the end, making the collage a true zoom
#' of the full one. If you need that view back, [plot_seu_marker_heatmap()] still
#' accepts `display_cells`, which is what implemented it.)
#'
#' Because the clusters are computed here, any requested resolution works --
#' nothing depends on which `<assay>_snn_res.<r>` columns the hypoxia split
#' happened to persist. The default sweep (0.2..0.8) still mirrors the full
#' low-hypoxia sweep so every two-clone panel has a same-resolution twin.
#'
#' Output naming follows [plot_seu_marker_heatmap()]'s `label`, e.g.
#' `results/<SRX>_hypoxia_low_seu.rds__1q_res0.2_heatmap_phase_scatter_patchwork.pdf`.
#'
#' @param seu_path Path to a Seurat `.rds` carrying `clone_opt` metadata.
#' @param scna_of_interest SCNA tag matched against comparison names, e.g.
#'   "1q", "2p", "16q" (matches "3_v_2_1q+", "2_v_1_16q-", ...).
#' @param large_clone_comparisons The parsed `large_clone_comparisons.yaml` list.
#' @param resolutions Resolution vector (default `seq(0.2, 1.2, 0.2)`).
#' @param nb_paths Numbat RDS paths for the clone-tree panel; `NULL` (default)
#'   omits it -- a two-clone tree is uninformative and the tree is the fragile /
#'   expensive panel.
#' @param clone_simplifications Passed through to [plot_seu_marker_heatmap()].
#' @param assay Assay whose `<assay>_snn_res.<r>` columns are used (default "SCT").
#' @param bar_signif Mark each cluster bar with the significance of its
#'   SCNA-status composition versus all other cells (Fisher's exact, BH-adjusted
#'   within the panel). Default `TRUE` here -- "is this cluster enriched for the
#'   clone that acquired the SCNA" is the question these collages exist to answer.
#'   The per-cluster q-values are also written to
#'   `results/<basename>_<scna>_res<r>_bar_enrichment.csv`.
#' @param bar_signif_min_cells Clusters smaller than this are not tested; their
#'   axis label reads `(n<min)` rather than silently going unmarked (default 20).
#' @return Character vector of PDF paths written (`NA_character_` per failed or
#'   skipped resolution). Never errors -- best-effort, like the sibling collage
#'   builders.
#' @export
plot_scna_two_clone_res_collages <- function(seu_path,
                                             scna_of_interest,
                                             large_clone_comparisons,
                                             resolutions = seq(0.2, 0.8, by = 0.2),
                                             nb_paths = NULL,
                                             clone_simplifications = NULL,
                                             assay = "SCT",
                                             bar_signif = TRUE,
                                             bar_signif_min_cells = 20) {
  if (is.null(seu_path) || length(seu_path) == 0 ||
      is.na(seu_path) || !file.exists(seu_path)) {
    return(NA_character_)
  }
  resolutions <- as.numeric(resolutions)
  sample_id <- stringr::str_extract(seu_path, "SR[RX][0-9]+")

  # Resolve the config key the same way plot_seu_marker_heatmap_by_scna() does:
  # try the filename slug, else fall back to the SRX tumor id.
  slug <- stringr::str_remove(fs::path_file(seu_path), "_filtered_seu.*")
  key  <- if (slug %in% names(large_clone_comparisons)) slug else sample_id
  comps <- names(large_clone_comparisons[[key]])
  comp  <- comps[stringr::str_detect(comps, stringr::fixed(scna_of_interest))]
  if (length(comp) == 0) {
    message(sample_id, ": no ", scna_of_interest,
            " clone comparison in large_clone_comparisons -> skip")
    return(NA_character_)
  }

  # Acquiring clone N and preceding clone M from the `<N>_v_<M>` prefix, unioned
  # across every comparison matching this SCNA (usually one).
  retained_clones <- comp |>
    stringr::str_extract("[0-9]+_v_[0-9]+") |>
    stringr::str_split("_v_", simplify = TRUE) |>
    as.vector()
  retained_clones <- unique(retained_clones[!is.na(retained_clones) &
                                              retained_clones != ""])
  if (length(retained_clones) < 1) {
    message(sample_id, " ", scna_of_interest,
            ": could not parse retained clones from '",
            paste(comp, collapse = ", "), "' -> skip")
    return(NA_character_)
  }

  # SCNA-of-interest status label per clone for the stacked-bar panel: the
  # acquiring clone N gets the signed SCNA token from the comparison name
  # (e.g. "1q+", "16q-"), the immediately-preceding clone M gets "preceding",
  # each tagged with its numbat clone id -> e.g. "1q+ (clone 2)" /
  # "preceding (clone 1)". Parsed per matching comparison (usually one).
  scna_label_map <- character(0)
  for (cmp in comp) {
    pair <- stringr::str_split(stringr::str_extract(cmp, "[0-9]+_v_[0-9]+"),
                               "_v_", simplify = TRUE)
    if (length(pair) < 2) next
    signed <- stringr::str_remove(cmp, "^[0-9]+_v_[0-9]+_")  # "1q+", "16q-", ...
    scna_label_map[pair[1]] <- glue::glue("{signed} (clone {pair[1]})")
    scna_label_map[pair[2]] <- glue::glue("preceding (clone {pair[2]})")
  }

  seu <- readRDS(seu_path)
  if (assay %in% names(seu@assays)) Seurat::DefaultAssay(seu) <- assay

  if (!"clone_opt" %in% colnames(seu@meta.data)) {
    message(sample_id, ": no clone_opt column -> skip")
    return(NA_character_)
  }

  # Attach the SCNA-of-interest status label the stacked-bar panel groups by
  # (bar_var = "scna_status"). Non-retained clones -> NA (they are not displayed).
  # Levels order the acquiring (signed) label before "preceding".
  scna_lvls <- unique(unname(scna_label_map[order(grepl("^preceding", scna_label_map))]))
  seu@meta.data$scna_status <- factor(
    unname(scna_label_map[as.character(seu@meta.data$clone_opt)]),
    levels = scna_lvls)

  # Restrict to the two compared clones BEFORE any dimensionality reduction, so
  # every downstream step -- PCA, SNN graph, clusters, marker selection, cluster
  # order, phase labels -- is computed on the cells actually under analysis.
  keep <- as.character(seu@meta.data$clone_opt) %in% retained_clones
  n_keep <- sum(keep, na.rm = TRUE)
  if (n_keep < 20) {
    message(sample_id, " ", scna_of_interest, ": only ", n_keep,
            " cells in clones ", paste(retained_clones, collapse = "/"),
            " -> skip")
    return(NA_character_)
  }
  seu <- seu[, colnames(seu)[which(keep)]]
  message(sample_id, " ", scna_of_interest, " two-clone collages (clones ",
          paste(retained_clones, collapse = " vs "), ", ", n_keep,
          " cells) at res ", paste(resolutions, collapse = " "))

  # Fresh PCA + SNN graph on the two-clone subset. npcs and k.param must stay
  # below the cell count or Seurat errors outright on the small subsets these
  # comparisons routinely produce, so both are clamped. On any failure fall back to
  # the persisted (full-population) columns rather than losing the collage --
  # mirrors the guarded recompute in plot_hypoxia_low_res_collages().
  snn_name <- glue::glue("{assay}_snn")
  npcs     <- max(2L, min(30L, ncol(seu) - 1L))
  k_param  <- max(2L, min(20L, ncol(seu) - 1L))
  recomputed <- tryCatch({
    seu <- Seurat::RunPCA(seu, assay = assay, npcs = npcs, verbose = FALSE)
    seu <- Seurat::FindNeighbors(seu, dims = 1:npcs, reduction = "pca",
                                 k.param = k_param,
                                 graph.name = paste0(assay, c("_nn", "_snn")),
                                 verbose = FALSE)
    # Cluster the WHOLE sweep here, not one resolution at a time inside the map:
    # the clustree panel needs every resolution to have been clustered on this
    # same two-clone graph. The persisted <assay>_snn_res.* columns carried over
    # from the parent object describe the full population, so a tree built from
    # them would not describe these cells. Stashed under clustree_res.* because
    # each collage below overwrites SCT_snn_res.0.6 with its own resolution.
    for (r in resolutions) {
      seu <- Seurat::FindClusters(seu, graph.name = snn_name, resolution = r,
                                  verbose = FALSE)
      seu@meta.data[[glue::glue("{assay}_snn_res.{r}")]] <- seu$seurat_clusters
      seu@meta.data[[glue::glue("clustree_res.{r}")]] <- seu$seurat_clusters
    }
    message(sample_id, " ", scna_of_interest, ": recomputed PCA on ", ncol(seu),
            " two-clone cells (npcs ", npcs, ", k ", k_param,
            "); clustered resolutions ", paste(resolutions, collapse = " "),
            " on it")
    TRUE
  }, error = function(e) {
    message("!! ", sample_id, " ", scna_of_interest, ": PCA recompute failed (",
            conditionMessage(e), "); falling back to persisted clusters")
    FALSE
  })

  if (!isTRUE(recomputed)) seu <- .stash_clustree_sweep(seu, assay)

  purrr::map_chr(resolutions, function(res) {
    tryCatch({
      col <- glue::glue("{assay}_snn_res.{res}")
      s   <- seu
      if (!col %in% colnames(s@meta.data)) {
        message(sample_id, " ", scna_of_interest, ": no ", col,
                " and no usable PCA -> skip res ", res)
        return(NA_character_)
      }
      # dummy_cluster_order()/group.by key on SCT_snn_res.0.6, so copy the
      # requested resolution's clusters into that column.
      s@meta.data[["SCT_snn_res.0.6"]] <- factor(as.character(s@meta.data[[col]]))

      # Temp copy keeps the source basename so plot_seu_marker_heatmap() derives a
      # predictable output name; the SCNA + resolution ride in `label`. The object
      # written here is already the two-clone subset, so no display_cells is
      # needed -- markers, order, and phase are derived from these cells alone.
      tmp <- file.path(tempdir(), basename(seu_path))
      saveRDS(s, tmp)

      plot_seu_marker_heatmap(
        tmp, cluster_order = NULL, nb_paths = nb_paths,
        clone_simplifications = clone_simplifications,
        bar_var = "scna_status",
        bar_signif = bar_signif, bar_signif_min_cells = bar_signif_min_cells,
        label = glue::glue("_{scna_of_interest}_res{res}_"))

      expected <- glue::glue(
        "results/{basename(seu_path)}__{scna_of_interest}_res{res}_heatmap_phase_scatter_patchwork.pdf")
      if (!file.exists(expected)) {
        message("!! two-clone scna collage produced NO pdf at ", expected)
        return(NA_character_)
      }
      message("wrote two-clone scna collage ", expected)
      as.character(expected)
    }, error = function(e) {
      message("!! two-clone scna collage FAILED for ", sample_id, " ",
              scna_of_interest, " res ", res, ": ", conditionMessage(e))
      NA_character_
    })
  })
}
