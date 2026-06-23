# Cluster-based iterative hypoxia split.
#
# Alternative to the per-cell hypoxia_score threshold used by
# subset_seu_by_expression(): cluster the cells, identify clusters driven by
# hypoxia marker genes, and assign WHOLE clusters to the high-hypoxia subset.
# Iterated (default twice) so hypoxia substructure that is masked by the
# non-hypoxia majority in round 1 surfaces once those cells are peeled off and
# the residual is re-clustered.
#
# Design mirrors the marker-overlap rule of classify_cluster_abbreviation()
# (cluster_dictionary_auto.R): a cluster is "hypoxia" if enough of its top
# markers are hypoxia genes OR (corroboration) its mean hypoxia_score is high.

# snoRNA-host genes (the ZFAS1/GAS5 family that drives these clusters) plus
# classic canonical-hypoxia genes.
.hypoxia_marker_set <- c(
  "ZFAS1", "GAS5", "SNHG5", "SNHG6", "SNHG7", "SNHG8", "SNHG15", "SNHG29",
  "VEGFA", "NDRG1", "BNIP3", "BNIP3L", "SLC2A1", "CA9", "PGK1", "LDHA",
  "ENO1", "ANKRD37", "ADM", "P4HA1", "PLOD2"
)

#' Identify hypoxia-driven clusters at one grouping
#'
#' A cluster is flagged hypoxia if >= `min_markers` of its top-`top_n` markers
#' (by logFC, padj < 0.5) are in `hypoxia_genes`, OR its mean `hypoxia_score`
#' is at/above the `score_quantile` quantile of the per-cell scores.
#'
#' @param seu Seurat object with the grouping column `grp` and (for the score
#'   corroboration) a `hypoxia_score` metadata column.
#' @param grp Name of the cluster grouping column in `seu@meta.data`.
#' @param hypoxia_genes Character vector of marker genes (default
#'   [.hypoxia_marker_set]).
#' @param min_markers,top_n Marker-overlap rule parameters.
#' @param score_quantile Quantile of per-cell `hypoxia_score` above which a
#'   cluster's mean is considered elevated.
#' @return Character vector of hypoxia cluster ids (possibly empty). Never
#'   errors — returns `character(0)` on failure.
#' @export
identify_hypoxia_clusters <- function(seu, grp,
                                      hypoxia_genes = .hypoxia_marker_set,
                                      min_markers = 2, top_n = 20,
                                      score_quantile = 0.75) {
  tryCatch({
    da <- Seurat::DefaultAssay(seu)
    # marker-overlap hits
    marker_hits <- character(0)
    n_grp <- length(unique(as.character(seu@meta.data[[grp]])))
    if (n_grp >= 2) {
      m <- presto::wilcoxauc(seu, grp, seurat_assay = da)
      marker_hits <- m |>
        dplyr::filter(.data$padj < 0.5) |>
        dplyr::group_by(.data$group) |>
        dplyr::arrange(dplyr::desc(.data$logFC), .by_group = TRUE) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::summarise(n_hyp = sum(.data$feature %in% hypoxia_genes),
                         .groups = "drop") |>
        dplyr::filter(.data$n_hyp >= min_markers) |>
        dplyr::pull(.data$group) |>
        as.character()
    }
    # score-corroboration hits
    score_hits <- character(0)
    if ("hypoxia_score" %in% colnames(seu@meta.data)) {
      sc   <- seu$hypoxia_score
      grpv <- as.character(seu@meta.data[[grp]])
      thr  <- stats::quantile(sc, score_quantile, na.rm = TRUE)
      cl_mean <- tapply(sc, grpv, mean, na.rm = TRUE)
      score_hits <- names(cl_mean)[!is.na(cl_mean) & cl_mean >= thr]
    }
    as.character(union(marker_hits, score_hits))
  }, error = function(e) {
    warning("identify_hypoxia_clusters failed for ", grp, ": ",
            conditionMessage(e))
    character(0)
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
#' Iteratively clusters cells, flags hypoxia-driven clusters (see
#' [identify_hypoxia_clusters()]), and peels their cells into the high subset.
#' A cell is high if it lands in a hypoxia cluster at ANY resolution in
#' `resolutions`. After labelling, the established [subset_seu_by_expression()]
#' path produces the `*_hypoxia_low_seu.rds` / `*_hypoxia_high_seu.rds` outputs
#' (re-clustering, markers, UMAP, hash metadata, barcode DB) so every
#' downstream target is unchanged.
#'
#' @param hypoxia_seu_path Path to a `*_seu_hypoxia.rds` (has `hypoxia_score`,
#'   a `pca` reduction, and ideally a `umap`).
#' @param resolutions Louvain resolutions to cluster each round at.
#' @param n_iter Number of peeling rounds (default 2).
#' @param hypoxia_genes,min_markers,top_n,score_quantile Passed to
#'   [identify_hypoxia_clusters()].
#' @param split_assay Assay used for the partition clustering (default "gene").
#' @param low_assay,high_assay Assays passed to [subset_seu_by_expression()]
#'   for the low/high subset re-clustering (defaults match the prior pipeline:
#'   gene for low, SCT for high).
#' @param write_qc Whether to emit a per-sample QC PDF.
#' @return Named character vector `c(low = ..., high = ...)` of output RDS
#'   paths, or `NA_character_` for a missing input.
#' @export
split_hypoxia_by_clusters <- function(hypoxia_seu_path,
                                      resolutions = c(0.2, 0.4, 0.6),
                                      n_iter = 2,
                                      hypoxia_genes = .hypoxia_marker_set,
                                      min_markers = 2, top_n = 20,
                                      score_quantile = 0.75,
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

  seu$hypoxia_round <- NA_integer_
  pool <- colnames(seu)  # current low pool

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
        hyp_cl <- identify_hypoxia_clusters(
          sub, grp, hypoxia_genes = hypoxia_genes, min_markers = min_markers,
          top_n = top_n, score_quantile = score_quantile)
        if (length(hyp_cl) > 0) {
          cells_res <- colnames(sub)[as.character(sub@meta.data[[grp]]) %in% hyp_cl]
          hyp_cells <- union(hyp_cells, cells_res)
        }
      }
      if (length(hyp_cells) == 0) break
      seu$hypoxia_round[colnames(seu) %in% hyp_cells] <- round
      pool <- setdiff(pool, hyp_cells)
    }
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

  low_path <- subset_seu_by_expression(
    labeled_path, run_hypoxia_clustering = TRUE,
    hypoxia_expr = "hypoxia_partition == 'low'",
    slug = "hypoxia_low", assay = low_assay)
  high_path <- subset_seu_by_expression(
    labeled_path, run_hypoxia_clustering = TRUE,
    hypoxia_expr = "hypoxia_partition == 'high'",
    slug = "hypoxia_high", assay = high_assay)

  c(low = low_path, high = high_path)
}
