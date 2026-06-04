# Plot Functions (11)

#' Filter data based on specified criteria
#'
#' @param sample_id Parameter for sample id
#' @param myseus Parameter for myseus
#' @param mynbs Parameter for mynbs
#' @param merged_metadata Parameter for merged metadata
#' @param myexpressions Parameter for myexpressions
#' @return Filtered data
#' @export
filter_numbat_cells <- function(sample_id, myseus, mynbs, merged_metadata, myexpressions) {
  seu <- readRDS(myseus[[sample_id]])

  seu_meta <- seu@meta.data %>%
    tibble::rownames_to_column("cell")

  merged_metadata_transfer <-
    merged_metadata %>%
    dplyr::filter(sample_id == {{ sample_id }}) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)

  mynb <- readRDS(mynbs[[sample_id]])

  nb_meta <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  myannot <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, GT_opt, clone_opt, nCount_gene, nFeature_gene) %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-"))

  if (!"" %in% unlist(myexpressions)) {
    filtered_nb <- filter_phylo_plot(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9, expressions = myexpressions)
  }

  returned_meta <- dplyr::left_join(
    filtered_nb[[3]][["data"]],
    myannot,
    by = "cell"
  ) %>%
    dplyr::mutate(cell = str_replace(cell, "-", ".")) %>%
    dplyr::filter(!is.na(cell)) %>%
    identity()

  return(returned_meta)
}

#' Create a plot visualization
#'
#' @param seu Seurat object
#' @param seu_name Parameter for seu name
#' @param subtype_hallmarks Parameter for subtype hallmarks
#' @return ggplot2 plot object
#' @export
plot_rb_subtype_expression <- function(seu, seu_name, subtype_hallmarks) {
  myfeatures <- c("exprs_gp1", "exprs_gp2", c(subtype_hallmarks))
  featureplots <- FeaturePlot(seu, myfeatures, combine = FALSE)

  featureplots <- map(featureplots, ~ (.x + labs(subtitle = seu_name))) %>%
    set_names(myfeatures)
}

#' Create a plot visualization
#'
#' @param seu Seurat object
#' @param checked_cluster_markers Cluster information
#' @return ggplot2 plot object
#' @export
plot_cluster_markers_by_cell_type <- function(seu, checked_cluster_markers) {
  cluster_plots <- map(checked_cluster_markers, ~ VlnPlot(seu, features = .x, group.by = "Phase"))

  cluster_plots <- map2(cluster_plots, names(checked_cluster_markers), ~ (.x + labs(subtitle = .y)))

  return(cluster_plots)
}
#' Generate a feature plot
#'
#' @param seu Seurat object
#' @param feature Parameter for feature
#' @return ggplot2 plot object
#' @export
compplot_feature_and_clusters <- function(seu, feature) {
  fp <- FeaturePlot(seu, feature)

  cp <- DimPlot(seu, group.by = "Phase")

  dp1 <- DimPlot(seu, group.by = c("gene_snn_res.0.15", "merged_leiden"))

  mypatch <- wrap_plots(fp, cp, nrow = 1) / dp1

  return(mypatch)
}

