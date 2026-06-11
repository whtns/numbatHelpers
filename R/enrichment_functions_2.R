# Enrichment and cluster analysis functions (9)

#' Perform differential expression analysis
#'
#' @param sample_id Parameter for sample id
#' @param myseus Parameter for myseus
#' @param celldf Cell identifiers or information
#' @param ... Additional arguments passed to other functions
#' @return Data frame
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O

enrich_diffex_by_cluster <- function(sample_id, myseus, celldf, ...) {
  
  seu <- myseus[[sample_id]]
  celldf <-
    celldf %>%
    dplyr::distinct(cell, .keep_all = TRUE) %>%
    tibble::column_to_rownames("cell") %>%
    identity()
  seu <- Seurat::AddMetaData(seu, celldf)
  clusters <- unique(seu$gene_snn_res.0.2) %>%
    set_names(.)
  
  clusters <- janitor::tabyl(seu@meta.data, gene_snn_res.0.2, clone_opt) %>%
    rowwise() %>%
    mutate(empty_clone = min(c_across(!any_of("gene_snn_res.0.2")))) %>%
    dplyr::filter(!empty_clone < 3) %>%
    dplyr::pull(`gene_snn_res.0.2`) %>%
    set_names(.) %>%
    identity()
  
  split_seu <- map(clusters, filter_seu_to_cluster, seu)
  possible_FindMarkers <- purrr::possibly(FindMarkers, otherwise = NA_real_)
  cluster_diffex <- map(split_seu, possible_FindMarkers, ...) %>%
    identity()
  cluster_diffex <- cluster_diffex[!is.na(cluster_diffex)]
  safe_enrichment_analysis <- purrr::safely(enrichment_analysis, otherwise = NA_real_)
  enrich_plots <- map(cluster_diffex, safe_enrichment_analysis) %>%
    map("result") %>%
    identity()
  enrich_plots <- enrich_plots[!is.na(enrich_plots)] %>%
    compact()
  enrich_plots <- compact(enrich_plots) %>%
    imap(~ (.x + labs(title = .y))) %>%
    identity()
  pdf_path <- tempfile(tmpdir = "results", fileext = ".pdf")
  pdf(pdf_path, width = 10)
  print(enrich_plots)
  dev.off()
  return(pdf_path)
}

#' Perform enrichment analysis
#'
#' @param seu_path File path
#' @return ggplot2 plot object
#' @export
enrich_by_cluster <- function(seu_path) {
  
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
  seu <- readRDS(seu_path)
  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))
  cluster_diffex <- seu@misc$markers$clusters$presto %>%
    split(.$Cluster) %>%
    map(tibble::column_to_rownames, "Gene.Name")
  
  
#' Perform drop cc genes operation
#'
#' @param df Input data frame or dataset
#' @param cc.genes Gene names or identifiers
#' @return Function result
#' @export
drop_cc_genes <- function(df, cc.genes = Seurat::cc.genes) {
  
  
    cc_genes <- unlist(cc.genes)
    dplyr::filter(df, !rownames(df) %in% cc_genes)
  }
  
  cluster_diffex <-
    cluster_diffex %>%
    identity()
  
  safe_enrichment_analysis <- purrr::safely(enrichment_analysis, otherwise = NA_real_)
  hallmark_enrich_results <- map(cluster_diffex, safe_enrichment_analysis, fold_change_col = "Average.Log.Fold.Change", gene_set = "hallmark") %>%
    map("result") %>%
    identity()
  gobp_enrich_results <- map(cluster_diffex, safe_enrichment_analysis, fold_change_col = "Average.Log.Fold.Change", gene_set = "gobp") %>%
    map("result") %>%
    identity()
  enrich_results <- c(rbind(hallmark_enrich_results, gobp_enrich_results))
  names(enrich_results) <- c(rbind(names(hallmark_enrich_results), names(gobp_enrich_results)))
  enrich_plots <- enrich_results |>
    map(plot_enrichment, p_val_cutoff = 1, signed = FALSE) |>
    compact() %>%
    imap(~ (.x + labs(title = sample_id, subtitle = .y))) %>%
    identity()
  gse_plot_path <- tempfile(tmpdir = "results", fileext = ".pdf")
  pdf(gse_plot_path, width = 8, height = 10)
  print(enrich_plots)
  dev.off()
  return(gse_plot_path)
}

#' Create a plot visualization
#'
#' @param tbl Parameter for tbl
#' @param sample_id Parameter for sample id
#' @return ggplot2 plot object
#' @export
plot_pcnv_by_reads <- function(tbl, sample_id) {
  
  
  ggplot(tbl, aes(x = p_cnv, y = nCount_gene)) +
    geom_point(size = 0.1, alpha = 0.1) +
    facet_wrap(~seg, ncol = 1)
  ggsave(glue("results/{sample_id}_pcnv_by_reads.pdf"))
  return(glue("results/{sample_id}_pcnv_by_reads.pdf"))
}

#' Perform numbat-related analysis
#'
#' @param numbat_rds_file File path
#' @return Data frame
#' @export
convert_numbat_pngs <- function(numbat_rds_file, n_sample = 6000, lim = 0.8) {

  numbat_output_dir <- str_remove(numbat_rds_file, "_numbat.*")
  sample_id <- str_extract(numbat_rds_file, "SR[RX][0-9]+")
  numbat_pngs <- dir_ls(numbat_output_dir, glob = "*.png") %>%
    set_names(path_file(.))
  numbat_pdfs <- stringr::str_replace(path_file(numbat_pngs), ".png", ".pdf")
  numbat_pdfs <- glue("results/numbat_plots/{sample_id}/{numbat_pdfs}")
  dir_create(glue("results/numbat_plots/{sample_id}"))

  # exp_roll_clust: regenerate as a true vector PDF using nb$plot_exp_roll() so
  # the many thin per-cell strokes produce the characteristic hazy look.
  exp_roll_idx <- which(names(numbat_pngs) == "exp_roll_clust.png")
  if (length(exp_roll_idx) > 0) {
    nb <- readRDS(numbat_rds_file)
    k  <- length(unique(nb$clone_post$clone_opt))
    pdf(numbat_pdfs[exp_roll_idx], width = 14, height = 5)
    # plot_exp_roll returns the patchwork object; must print() to render to device.
    print(nb$plot_exp_roll(k = k, n_sample = n_sample, lim = lim))
    dev.off()
  }

  # All other PNGs: convert via magick as before.
  other_idx    <- setdiff(seq_along(numbat_pngs), exp_roll_idx)
  other_images <- purrr::map(numbat_pngs[other_idx], magick::image_read) %>%
    imap(~ magick::image_annotate(.x, sample_id, size = 50))
  map2(other_images, numbat_pdfs[other_idx], ~ magick::image_write(.x, format = "pdf", .y))

  return(numbat_pdfs)
}

#' Perform compare infercnv operation
#'
#' @param myseus Parameter for myseus
#' @param sample_id Parameter for sample id
#' @return Function result
#' @export
compare_infercnv <- function(myseus, sample_id) {
  
  
  seu <- myseus[[sample_id]]
}