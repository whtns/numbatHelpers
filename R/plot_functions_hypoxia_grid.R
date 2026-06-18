# Hypoxia-threshold x clustering-resolution grid search.
#
# Adapted from plot_effect_of_filtering(): instead of comparing filtering
# stages, this sweeps a grid of hypoxia thresholds (which define the low- and
# high-hypoxia subsets) crossed with clustering resolutions, and reports for
# each (threshold, subset, resolution):
#   * a marker-gene dotplot annotated with the cell count,
#   * a UMAP coloured by cluster at that resolution,
#   * a UMAP coloured by diploid-vs-aneuploid status,
# plus a tidy cell-count / cluster-count summary table across the whole grid.

#' Top-marker dotplot for one Seurat (sub)object at a given grouping.
#'
#' Computes markers with find_all_markers() (presto), takes the top genes per
#' cluster, and draws a Seurat DotPlot. Title is annotated with the cell count.
#' Returns a placeholder panel (never errors) when markers can't be produced.
#'
#' @keywords internal
hypoxia_grid_marker_dotplot <- function(sub, grp, title, n_top = 5, max_cells = 50000L) {
  blank <- function(msg) {
    ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = msg, size = 3) +
      ggplot2::labs(title = title)
  }
  n <- ncol(sub)
  if (n < 20) return(blank(glue::glue("n = {n}\n(too few cells)")))
  if (n > max_cells) return(blank(glue::glue("n = {format(n, big.mark = ',')}\n(too many cells)")))
  n_grp <- length(unique(as.character(sub@meta.data[[grp]])))
  if (n_grp < 2) return(blank(glue::glue("n = {format(n, big.mark = ',')}\n(1 cluster)")))

  sub <- tryCatch(
    if (exists("find_all_markers", mode = "function")) find_all_markers(sub, grp) else sub,
    error = function(e) sub
  )
  presto <- tryCatch(sub@misc$markers[[grp]]$presto, error = function(e) NULL)
  if (is.null(presto)) return(blank(glue::glue("n = {format(n, big.mark = ',')}\n(no markers)")))

  top <- presto |>
    dplyr::group_by(Cluster) |>
    dplyr::slice_head(n = n_top) |>
    dplyr::ungroup() |>
    dplyr::arrange(Cluster, Gene.Name)
  feats <- unique(top$Gene.Name)
  feats <- feats[feats %in% rownames(sub)]
  if (length(feats) == 0) return(blank(glue::glue("n = {format(n, big.mark = ',')}\n(no valid genes)")))

  tryCatch(
    Seurat::DotPlot(sub, features = feats, group.by = grp, dot.scale = 5) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0)) +
      ggplot2::labs(title = title, x = NULL, y = NULL),
    error = function(e) blank(glue::glue("n = {format(n, big.mark = ',')}\n(dotplot failed)"))
  )
}

#' Assess hypoxia threshold x clustering resolution for one sample.
#'
#' @param hypoxia_seu_path Path to a *_seu_hypoxia.rds file (must contain a
#'   hypoxia_score column, and ideally a precomputed SNN graph + UMAP).
#' @param thresholds Numeric hypoxia thresholds; low = score <= t, high = > t.
#' @param resolutions Numeric Louvain resolutions to cluster each subset at.
#' @param assay Assay whose SNN graph is used for clustering (default "gene").
#' @param plot_path Output PDF path (auto-derived from sample id if NULL).
#' @param csv_path Output CSV summary path (auto-derived if NULL).
#' @return list(plot = <pdf path>, summary = <csv path>), or NULL on bad input.
#' @export
assess_hypoxia_clustering_grid <- function(hypoxia_seu_path,
                                           thresholds = c(0.5, 0.4, 0.3),
                                           resolutions = c(0.2, 0.4),
                                           assay = "gene",
                                           plot_path = NULL,
                                           csv_path = NULL) {
  if (is.null(hypoxia_seu_path) || length(hypoxia_seu_path) == 0 || is.na(hypoxia_seu_path)) {
    return(NULL)
  }
  sample_id <- stringr::str_extract(hypoxia_seu_path, "SR[RX][0-9]+")
  fs::dir_create("results/hypoxia_clustering_grid")
  if (is.null(plot_path)) {
    plot_path <- glue::glue("results/hypoxia_clustering_grid/{sample_id}_hypoxia_grid.pdf")
  }
  if (is.null(csv_path)) {
    csv_path <- glue::glue("results/hypoxia_clustering_grid/{sample_id}_hypoxia_grid.csv")
  }

  seu <- readRDS(hypoxia_seu_path)
  if (!"hypoxia_score" %in% colnames(seu@meta.data)) {
    warning("No hypoxia_score in ", sample_id, "; skipping grid.")
    return(NULL)
  }
  da <- if (assay %in% names(seu@assays)) assay else Seurat::DefaultAssay(seu)
  Seurat::DefaultAssay(seu) <- da
  graph_name <- paste0(da, "_snn")
  has_umap <- "umap" %in% names(seu@reductions)

  # diploid-vs-aneuploid label from scna (fall back to GT_opt)
  diploid_label <- function(s) {
    scna <- if ("scna" %in% colnames(s@meta.data)) as.character(s@meta.data$scna) else rep(NA_character_, ncol(s))
    gt   <- if ("GT_opt" %in% colnames(s@meta.data)) as.character(s@meta.data$GT_opt) else rep(NA_character_, ncol(s))
    val <- ifelse(!is.na(scna) & scna != "", scna, gt)
    factor(ifelse(is.na(val) | val == "" | tolower(val) == "diploid", "diploid", "aneuploid"),
           levels = c("diploid", "aneuploid"))
  }
  diploid_cols <- c(diploid = "grey70", aneuploid = "#d73027")

  panels <- list()
  summary_rows <- list()

  for (t in thresholds) {
    subsets <- list(
      low  = seu[, !is.na(seu$hypoxia_score) & seu$hypoxia_score <= t],
      high = seu[, !is.na(seu$hypoxia_score) & seu$hypoxia_score >  t]
    )
    for (state in names(subsets)) {
      sub <- subsets[[state]]
      n_cells <- ncol(sub)
      sub$diploid_status <- diploid_label(sub)
      n_aneu <- sum(sub$diploid_status == "aneuploid")

      for (res in resolutions) {
        grp <- glue::glue("gridclust_res.{res}")
        n_clusters <- NA_integer_
        # (re)cluster the subset at this resolution on the inherited SNN graph
        if (n_cells >= 20 && graph_name %in% names(sub@graphs)) {
          sub <- tryCatch(
            Seurat::FindClusters(sub, graph.name = graph_name, resolution = res, verbose = FALSE),
            error = function(e) sub
          )
          if ("seurat_clusters" %in% colnames(sub@meta.data)) {
            sub@meta.data[[grp]] <- sub$seurat_clusters
          }
        }
        # fall back to an existing resolution column if clustering didn't run
        if (!grp %in% colnames(sub@meta.data)) {
          fb <- paste0(da, "_snn_res.", res)
          if (fb %in% colnames(sub@meta.data)) sub@meta.data[[grp]] <- sub@meta.data[[fb]]
        }
        if (grp %in% colnames(sub@meta.data)) {
          n_clusters <- length(unique(as.character(sub@meta.data[[grp]])))
        }

        title_base <- glue::glue("{state} hyp≤{t} · res {res} · n = {format(n_cells, big.mark = ',')}")

        # marker dotplot (panel 1)
        p_markers <- hypoxia_grid_marker_dotplot(sub, grp, glue::glue("markers · {title_base}"))

        # UMAP by cluster (panel 2)
        p_clusters <- if (has_umap && grp %in% colnames(sub@meta.data)) {
          tryCatch(
            Seurat::DimPlot(sub, group.by = grp, reduction = "umap") +
              ggplot2::theme_minimal(base_size = 8) +
              ggplot2::labs(title = glue::glue("clusters · {title_base}")),
            error = function(e) NULL
          )
        } else NULL

        # UMAP by diploid status (panel 3)
        p_diploid <- if (has_umap) {
          tryCatch(
            Seurat::DimPlot(sub, group.by = "diploid_status", reduction = "umap",
                            cols = diploid_cols) +
              ggplot2::theme_minimal(base_size = 8) +
              ggplot2::labs(title = glue::glue("diploid · {title_base}")),
            error = function(e) NULL
          )
        } else NULL

        row_panels <- purrr::compact(list(p_clusters, p_diploid, p_markers))
        if (length(row_panels) > 0) {
          panels[[length(panels) + 1]] <- patchwork::wrap_plots(row_panels, nrow = 1)
        }

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          sample_id   = sample_id,
          threshold   = t,
          subset      = state,
          resolution  = res,
          n_cells     = n_cells,
          n_clusters  = n_clusters,
          n_aneuploid = n_aneu,
          n_diploid   = n_cells - n_aneu,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  summary_df <- dplyr::bind_rows(summary_rows)
  readr::write_csv(summary_df, csv_path)

  if (length(panels) > 0) {
    grid <- patchwork::wrap_plots(panels, ncol = 1) +
      patchwork::plot_annotation(
        title = glue::glue("{sample_id} — hypoxia threshold × resolution grid")
      )
    # height scales with number of (threshold x subset x resolution) rows
    tryCatch(
      ggplot2::ggsave(plot_path, grid, height = 3.2 * length(panels), width = 15,
                      units = "in", limitsize = FALSE),
      error = function(e) warning("Failed to save hypoxia grid PDF: ", conditionMessage(e))
    )
  }

  list(plot = plot_path, summary = csv_path)
}
