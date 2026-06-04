# Plot Functions (126)
#' Create a plot visualization
#'
#' @param diffex_path File path
#' @return ggplot2 plot object
#' @export
gse_plot_from_clone_diffex <- function(diffex_path) {
  sample_id <- str_extract(diffex_path, "SR[RX][0-9]+")

  numbat_dir <- path_split(diffex_path)[[1]][[2]]

  location <- str_extract(diffex_path, "(?<=diffex_).*_segment")

  diffex <-
    diffex_path %>%
    read_csv() %>%
    split(.$clone_comparison) %>%
    identity()

  annotable_cols <- colnames(annotables::grch38)
  annotable_cols <- annotable_cols[!annotable_cols == "symbol"]
  g
eu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[, !is.na(seu$clone_opt)]

  Seurat::FeaturePlot(seu, ...)
}