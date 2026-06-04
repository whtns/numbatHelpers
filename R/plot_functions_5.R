# Feature plotting and differential expression functions (6)

#' Generate a feature plot
#'
#' @param seu Seurat object
#' @param seu_title Plot title
#' @param myvar Parameter for myvar
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - multiple_joins: Combine multiple joins into single join operation where possible

plot_feature_across_seus <- function(seu, seu_title, myvar) {
  varplot <- seu %>%
    FeaturePlot(features = myvar)
  dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
  phase_plot <- DimPlot(seu, group.by = "Phase")
  (varplot | (phase_plot / dimplot)) +
    plot_annotation(title = seu_title)
}
#' Create a plot visualization
#'
#' @param nb Numbat object
#' @param myseu Seurat object
#' @param myannot Parameter for myannot
#' @param mytitle Plot title
#' @param expressions Parameter for expressions
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
filter_phylo_plot <- function(nb, myseu, myannot, mytitle, expressions, ...) {
  celltypes <-
    myseu@meta.data["type"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    identity()
  myannot <- dplyr::left_join(myannot, celltypes, by = "cell")
  mypal <- c("1" = "gray", "2" = "#377EB8", "3" = "#4DAF4A", "4" = "#984EA3")
  initial_phylo_heatmap <- nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = myannot,
    show_phylo = FALSE,
    sort_by = "GT_opt",
    ...
  ) +
    labs(title = mytitle)
  phylo_heatmap_data <- initial_phylo_heatmap$data %>%
    dplyr::left_join(myannot, by = "cell")
}