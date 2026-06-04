# Plot Functions (115)

#' Perform clustering analysis
#'
#' @param seu Seurat object
#' @param seu_name Parameter for seu name
#' @param clusters Cluster information
#' @return List object
#' @export
table_distribution_of_clones_across_clusters <- function(seu, seu_name, clusters = "SCT_snn_res.0.4") {
  meta <- seu@meta.data %>%
    dplyr::mutate(cluster = factor(.data[[clusters]])) %>%
    identity()

  cluster_per_scna <-
    janitor::tabyl(meta, cluster, scna) %>%
    janitor::adorn_percentages("col") %>%
    # tidyr::pivot_longer(-cluster, names_to = "scna", values_to = "percent_scna") %>%
    dplyr::mutate(sample_id = seu_name) %>%
    identity()

  scna_per_cluster <-
    janitor::tabyl(meta, scna, cluster) %>%
    janitor::adorn_percentages("col") %>%
    # tidyr::pivot_longer(-scna, names_to = "cluster", values_to = "percent_cluster") %>%
    dplyr::mutate(sample_id = seu_name) %>%
    identity()

  return(list("cluster_per_scna" = cluster_per_scna, "scna_per_cluster" = scna_per_cluster))
}

#' Perform clustering analysis
#'
#' @param seu Seurat object
#' @param assay Character string (default: "SCT")
#' @return Function result
#' @export
table_cluster_markers <- function(seu, assay = "SCT") {
  cluster_names <- glue("{assay}_snn_res.{seq(0.2, 1.0, by = 0.2)}") %>%
    set_names(.)

  map(cluster_names, ~ (seu@misc$markers[[.x]]$presto))
}

#' Create a numbat-related plot visualization
#'
#' @param seu_path File path
#' @param numbat_rds_files File path
#' @param cluster_dictionary Cluster information
#' @param filter_expressions Parameter for filter expressions
#' @param clone_simplifications Parameter for clone simplifications
#' @param extension Character string (default: "")
#' @return ggplot2 plot object
#' @export
make_numbat_plot_files <- function(seu_path, numbat_rds_files, cluster_dictionary, filter_expressions = NULL, clone_simplifications = NULL, extension = "") {

  sample_id <- str_extract(seu_path, "SR[RX][0-9]+")

  tryCatch({

  names(numbat_rds_files) <- str_extract(numbat_rds_files, "SR[RX][0-9]+")
  numbat_rds_file <- numbat_rds_files[[sample_id]]

  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(seu_path)

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[, !is.na(seu$clone_opt)]

  # Compute scna per cell from GT_opt + clone_simplifications (not stored in *_seu.rds)
  rb_scnas_lookup <- tibble::enframe(clone_simplifications[[sample_id]], "scna", "seg") %>%
    tidyr::unnest(seg) %>%
    dplyr::mutate(seg = as.character(seg))

  scna_labels <- vapply(
    seu@meta.data$GT_opt,
    function(gt_opt) {
      wrapped <- wrap_scna_labels(simplify_gt_col(gt_opt, rb_scnas_lookup))
      if (length(wrapped) == 0) {
        return(NA_character_)
      }
      as.character(wrapped[[1]])
    },
    FUN.VALUE = character(1)
  )

  scna_per_cell <- data.frame(
    scna = scna_labels,
    row.names = rownames(seu@meta.data),
    stringsAsFactors = FALSE
  )

  scna_per_cell <- scna_per_cell[colnames(seu), , drop = FALSE]

  seu <- Seurat::AddMetaData(seu, scna_per_cell)

  # Drop unused factor levels from seurat_clusters after cell filtering to avoid
  # wilcoxauc dimension mismatch ("number of columns of X does not match length of y")
  seu$seurat_clusters <- droplevels(seu$seurat_clusters)

  if (length(levels(seu$seurat_clusters)) > 1 && ncol(seu) > 0) {
    plot_markers(seu, metavar = "seurat_clusters", marker_method = "presto", return_plotly = FALSE, hide_technical = "all", seurat_assay = "gene") +
      ggplot2::scale_y_discrete(position = "left") +
      labs(title = sample_id)

    ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_sample_marker{extension}.pdf"))
  } else {
    warning("Skipping marker plot for ", sample_id, ": insufficient clusters or cells")
  }

  # output_plots[["merged_marker"]] + output_plots[["sample_marker"]]
  # ggsave(glue("results/{numbat_dir}/{sample_id}_combined_marker{extension}.pdf"))

  # seu <- annotate_seu_with_rb_subtype_gene_expression(seu)

  DimPlot(seu, group.by = c("abbreviation", "clone_opt", "Phase")) +
    plot_annotation(title = sample_id)
  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_dimplot{extension}.pdf"), width = 10)


  ## clone distribution ------------------------------
  plot_distribution_of_clones_across_clusters(seu, sample_id, var_y = "seurat_clusters")
  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_clone_distribution{extension}.pdf"), width = 8, height = 4)

  ## clone tree ------------------------------

  nclones <- length(unique(mynb$clone_post$clone_opt))

  mypal <- scales::hue_pal()(nclones) %>%
    set_names(1:nclones)

  rb_scnas <- clone_simplifications[[sample_id]]

  mynb <- simplify_gt(mynb, rb_scnas)

  mynb$plot_mut_history(pal = mypal, legend = FALSE, horizontal = FALSE) +
    labs(title = sample_id)

  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_tree{extension}.pdf"), width = 2, height = 5)

  # plot types ------------------------------
  plot_types <- c("dimplot", "sample_marker", "clone_distribution", "tree")

  plot_files <- glue("results/{numbat_dir}/{sample_id}/{sample_id}_{plot_types}{extension}.pdf") %>%
    set_names(plot_types)

  return(plot_files)

  }, error = function(e) {
    stop(paste0("[", sample_id, "] ", conditionMessage(e)), call. = FALSE)
  })
}