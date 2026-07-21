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
#' No PCA or clustering is performed here: the persisted `<assay>_snn_res.<r>`
#' columns (written by the hypoxia split, which recomputes PCA on the low-hypoxia
#' population) are read as-is. Moreover the collage is a true ZOOM of the
#' full-population collage: cluster order, phase labels, and marker-gene rows are
#' all derived on the FULL low object inside [plot_seu_marker_heatmap()], which is
#' then handed `display_cells` and restricts every panel to the two compared clones
#' only at the end. So the two-clone panels carry exactly the same cluster IDs,
#' marker rows, and ordering as the full collage -- just fewer cells. The default
#' resolution sweep (0.2..0.8) mirrors the full low-hypoxia sweep so every
#' two-clone panel has a same-resolution full-population twin. Resolutions the
#' split did not persist are skipped, never clustered on demand.
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
                                             assay = "SCT") {
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

  # Use ONLY the cluster columns the hypoxia split already persisted on this low
  # object (its PCA was recomputed on the low-hypoxia population upstream). No PCA
  # and no FindClusters happen here, so the two-clone collages show exactly the
  # clusters of the full-population low-hypoxia collages -- just restricted to the
  # two compared clones. Any requested resolution the split did not persist is
  # skipped, never clustered on demand.
  have <- purrr::map_lgl(resolutions, function(res)
    glue::glue("{assay}_snn_res.{res}") %in% colnames(seu@meta.data))
  if (any(!have)) {
    message(sample_id, " ", scna_of_interest, ": no persisted ", assay,
            "_snn_res column for res ", paste(resolutions[!have], collapse = " "),
            " -> skipping those (not clustering on demand)")
  }
  resolutions <- resolutions[have]
  if (length(resolutions) == 0) {
    message(sample_id, " ", scna_of_interest,
            ": none of the requested resolutions are persisted -> skip")
    return(NA_character_)
  }

  # Identify the two clones' cells to DISPLAY. We do NOT subset the object here:
  # plot_seu_marker_heatmap() derives the cluster order, phase labels, and
  # marker-gene rows on the FULL low object (identical to the full-population
  # collage) and only then restricts every panel to `display_cells`. So the
  # two-clone collage is a true zoom of the full collage -- same cluster IDs, same
  # marker rows, same order -- rather than a re-analysis of the two-clone slice.
  keep <- as.character(seu@meta.data$clone_opt) %in% retained_clones
  n_keep <- sum(keep, na.rm = TRUE)
  if (n_keep < 20) {
    message(sample_id, " ", scna_of_interest, ": only ", n_keep,
            " cells in clones ", paste(retained_clones, collapse = "/"),
            " -> skip")
    return(NA_character_)
  }
  display_cells <- colnames(seu)[which(keep)]
  message(sample_id, " ", scna_of_interest, " two-clone collages (clones ",
          paste(retained_clones, collapse = " vs "), ", ", n_keep,
          " cells) at res ", paste(resolutions, collapse = " "))

  purrr::map_chr(resolutions, function(res) {
    tryCatch({
      col <- glue::glue("{assay}_snn_res.{res}")
      s   <- seu
      # dummy_cluster_order()/group.by key on SCT_snn_res.0.6, so copy the
      # requested resolution's clusters into that column ON THE FULL object.
      s@meta.data[["SCT_snn_res.0.6"]] <- factor(as.character(s@meta.data[[col]]))

      # Temp copy keeps the source basename so plot_seu_marker_heatmap() derives a
      # predictable output name; the SCNA + resolution ride in `label`. The FULL
      # object is written -- markers/order/phase are computed on it, then the panels
      # are restricted to `display_cells` inside plot_seu_marker_heatmap().
      tmp <- file.path(tempdir(), basename(seu_path))
      saveRDS(s, tmp)

      plot_seu_marker_heatmap(
        tmp, cluster_order = NULL, nb_paths = nb_paths,
        clone_simplifications = clone_simplifications,
        display_cells = display_cells, bar_var = "scna_status",
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
