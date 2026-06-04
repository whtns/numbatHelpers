# Data retrieval and differential expression functions (7)
#' Perform differential expression analysis
#'
#' @param numbat_rds_file File path
#' @param filter_expressions Parameter for filter expressions
#' @param idents Cell identities or groups
#' @return ggplot2 plot object
#' @export
diffex_groups <- function(numbat_rds_file, filter_expressions, idents = NULL) {
  sample_id <- str_extract(numbat_rds_file, "SR[RX][0-9]+")
  ident.1 <- idents[[sample_id]]
  filter_expressions <- filter_expressions[[sample_id]]
  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]
  dir_create(glue("results/{numbat_dir}"))
  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()
  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))
  mynb <- readRDS(numbat_rds_file)
  joined_nb_meta <-
    mynb$clone_post %>%
    dplyr::select(-p_cnv) %>%
    dplyr::left_join(mynb$joint_post, by = "cell")
  excluded_cells <- map(filter_expressions, pull_cells_matching_expression, joined_nb_meta) %>%
    unlist()
  nb_meta <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    dplyr::filter(!cell %in% excluded_cells) %>%
    tibble::column_to_rownames("cell")
  seu <- Seurat::AddMetaData(seu, nb_meta)
  seu <- seu[, !is.na(seu$clone_opt)]
  diffex <- Seurat::FindMarkers(seu, ident.1)
  gse_plot_path <- glue("results/{numbat_dir}/{sample_id}_gsea.pdf")
  gse_plot <- enrichment_analysis(diffex) +
    labs(title = glue("{sample_id}"))
  ggsave(gse_plot_path)
  diffex_path <- glue("results/{numbat_dir}/{sample_id}_diffex.csv")
  diffex <-
    diffex %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::distinct(ensgene, .keep_all = TRUE)
  write_csv(diffex, diffex_path)
  return(list(sample_id, diffex_path, gse_plot_path))
}

#' Perform differential expression analysis
#'
#' @param myseus Parameter for myseus
#' @param cells.1 Cell identifiers or information
#' @param cells.2 Cell identifiers or information
#' @param ... Additional arguments passed to other functions
#' @return Differential expression results
#' @export
diffex_cells <- function(myseus, cells.1, cells.2, ...) {
  seu <- readRDS(myseus[[sample_id]])
  Seurat::FindMarkers(seu$gene, cells.1 = cells.1, cells.2 = cells.2)
}

#' Perform split by cnv operation
#'
#' @param returned_meta Parameter for returned meta
#' @param myseg Character string (default: "2a")
#' @return Function result
#' @export
split_by_cnv <- function(returned_meta, myseg = "2a") {
  test0 <-
    returned_meta %>%
    dplyr::filter(!is.na(cell)) %>%
    dplyr::filter(seg == myseg) %>%
    dplyr::mutate(cnv_status = ifelse(p_cnv > 0.5, "present", "absent")) %>%
    dplyr::group_by(cnv_status) %>%
    dplyr::group_split() %>%
    map(pull, cell) %>%
    identity()
}

#' Create a plot visualization
#'
#' @param tbl Parameter for tbl
#' @param sample_id Parameter for sample id
#' @return ggplot2 plot object
#' @export
plot_pcnv_by_nsnp <- function(tbl, sample_id) {
  ggplot(tbl, aes(x = p_cnv, y = n_snp)) +
    geom_point(size = 0.1, alpha = 0.1) +
    facet_wrap(~seg, ncol = 1)
  ggsave(glue("results/{sample_id}_pcnv_by_nsnp.pdf"))
  return(glue("results/{sample_id}_pcnv_by_nsnp.pdf"))
}