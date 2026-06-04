# Plot Functions (118)

#' Merge or join datasets
#'
#' @param meta_path File path
#' @return Data frame
#' @export
get_merged_metadata <- function(meta_path) {
  metadata <- read_csv(meta_path) %>%
    set_names(c("cell", "sample_id", "merged_leiden")) %>%
    dplyr::mutate(sample_id = str_remove(sample_id, ".h5ad")) %>%
    dplyr::mutate(cell = str_replace(cell, "-", ".")) %>%
    identity()

  return(metadata)
}
#' Perform differential expression analysis
#'
#' @param sample_id Parameter for sample id
#' @param myseus Parameter for myseus
#' @param celldf Cell identifiers or information
#' @param ... Additional arguments passed to other functions
#' @return Differential expression results
#' @export
diffex_groups_old <- function(sample_id, myseus, celldf, ...) {
  seu <- readRDS(myseus[[sample_id]])

  celldf <-
    celldf %>%
    dplyr::distinct(cell, .keep_all = TRUE) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, celldf)

  Seurat::FindMarkers(seu, ...)
}