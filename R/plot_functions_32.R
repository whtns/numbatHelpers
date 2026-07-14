# Plot Functions (131)

#' Perform differential expression analysis
#'
#' @param cluster_diffex_clones Cluster information
#' @param cluster_xlsx Cluster information
#' @param cluster_by_chr_xlsx Cluster information
#' @param total_diffex_clones Parameter for total diffex clones
#' @param total_xlsx Character string (default: "results/straight_diffex_bw_clones_large.xlsx")
#' @param total_by_chr_xlsx Character string (default: "results/straight_diffex_bw_clones_large_by_chr.xlsx")
#' @return Differential expression results
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

tabulate_diffex_clones <- function(cluster_diffex_clones,
                                   cluster_xlsx = "results/diffex_bw_clones_per_cluster_large.xlsx",
                                   cluster_by_chr_xlsx = "results/diffex_bw_clones_per_cluster_large_by_chr.xlsx",
                                   total_diffex_clones,
                                   total_xlsx = "results/straight_diffex_bw_clones_large.xlsx",
                                   total_by_chr_xlsx = "results/straight_diffex_bw_clones_large_by_chr.xlsx") {
  # Drop skipped samples (NA_character_) and unwritten files before anything reads
  # them as filenames -- see .drop_missing_paths().
  cluster_diffex_clones <- .drop_missing_paths(cluster_diffex_clones)
  total_diffex_clones   <- .drop_missing_paths(total_diffex_clones)
  if (length(cluster_diffex_clones) == 0 && length(total_diffex_clones) == 0) {
    message("tabulate_diffex_clones: no usable diffex CSVs, returning NULL")
    return(invisible(NULL))
  }

  kooi_candidates <- read_csv("data/kooi_candidates.csv")

  cc_genes <- Seurat::cc.genes %>%
    tibble::enframe("phase_of_gene", "symbol") %>%
    tidyr::unnest(symbol)

  sample_ids <- str_extract(cluster_diffex_clones, "SR[RX][0-9]+")

  
#' Perform differential expression analysis
#'
#' @param diffex Parameter for diffex
#' @return Differential expression results
#' @export
annotate_percent_segment_diffex <- function(diffex) {
    if ("genes_in_segment" %in% colnames(diffex)) {
      diffex <-
        diffex %>%
        dplyr::group_by(seg, genes_in_segment) %>%
        dplyr::select(symbol, genes_in_segment) %>%
        dplyr::summarize(diffex_genes = list(symbol)) %>%
        dplyr::mutate(genes_in_segment = str_split(genes_in_segment, ", ")) %>%
        mutate(diffex_genes_in_segment = map2(diffex_genes, genes_in_segment, safely(base::intersect))) %>%
        mutate(diffex_genes_in_segment = map(diffex_genes_in_segment, "result")) %>%
        dplyr::mutate(ratio_genes_diffex_in_segment = map2_dbl(diffex_genes_in_segment, genes_in_segment, ~ (length(.x) / length(.y)))) %>%
        dplyr::select(seg, ratio_genes_diffex_in_segment) %>%
        dplyr::left_join(diffex, by = "seg") %>%
        dplyr::select(-c("genes_in_segment")) %>%
        identity()
    }

    return(diffex)
  }

  total_diffex_clones <-
    total_diffex_clones %>%
    set_names(sample_ids) %>%
    map(bind_rows, .id = "clone_comparison") %>%
    purrr::keep(~ ncol(.x) > 0) %>%
  	purrr::compact() |>
    map(dplyr::left_join, cc_genes, by = "symbol") %>%
    map(group_by, clone_comparison) %>%
    map(dplyr::filter, p_val_adj < 1) %>%
    purrr::keep(~ nrow(.x) > 0) %>%
    map(dplyr::arrange, clone_comparison, p_val_adj) %>%
    map(dplyr::select, clone_comparison, chr, symbol, description, everything()) %>%
    map(dplyr::left_join, kooi_candidates, by = "symbol") %>%
    map(annotate_percent_segment_diffex) %>%
    identity()

  # cluster ------------------------------
  cluster_diffex_clones <-
    cluster_diffex_clones %>%
    set_names(str_extract(., "SR[RX][0-9]+")) %>%
    map(read_csv) %>%
    purrr::keep(~ nrow(.x) > 0) %>%
    map(dplyr::filter, p_val_adj < 1) %>%
    purrr::keep(~ nrow(.x) > 0) %>%
    map(dplyr::left_join, cc_genes, by = "symbol") %>%
    map(dplyr::arrange, clone_comparison, cluster, p_val_adj) %>%
    map(dplyr::select, clone_comparison, cluster, chr, symbol, description, everything()) %>%
    map(dplyr::left_join, kooi_candidates, by = "symbol") %>%
    map(annotate_percent_segment_diffex) %>%
    identity()

  cluster_diffex_clones <-
    cluster_diffex_clones %>%
    # map2(total_diffex_clones, annotate_cluster_membership, "is_total_clone_diffex") %>%
    map(dplyr::distinct, clone_comparison, cluster, chr, symbol, .keep_all = TRUE)

  write_xlsx(cluster_diffex_clones, cluster_xlsx)

  cluster_diffex_clones_by_chr <-
    if (length(cluster_diffex_clones) > 0) {
      cluster_diffex_clones %>%
        dplyr::bind_rows(.id = "sample_id") %>%
        dplyr::distinct(sample_id, clone_comparison, cluster, chr, symbol, .keep_all = TRUE) %>%
        dplyr::group_by(symbol) %>%
        dplyr::mutate(num_samples = length(unique(sample_id))) %>%
        dplyr::arrange(desc(num_samples), symbol) %>%
        dplyr::filter(num_samples > 1) %>%
        dplyr::filter(!str_detect(chr, "CHR_")) %>%
        split(.$chr)
    } else {
      list()
    }

  write_xlsx(cluster_diffex_clones_by_chr, cluster_by_chr_xlsx)

  # total ------------------------------

  total_diffex_clones <-
    total_diffex_clones %>%
    # map2(total_diffex_clones, cluster_diffex_clones, annotate_cluster_membership, "is_cluster_clone_diffex") %>%
    map(dplyr::distinct, clone_comparison, chr, symbol, .keep_all = TRUE)

  write_xlsx(total_diffex_clones, total_xlsx)

  total_diffex_clones_by_chr <-
    if (length(total_diffex_clones) > 0) {
      total_diffex_clones %>%
        dplyr::bind_rows(.id = "sample_id") %>%
        dplyr::distinct(sample_id, clone_comparison, chr, symbol, .keep_all = TRUE) %>%
        dplyr::group_by(symbol) %>%
        dplyr::mutate(num_samples = length(unique(sample_id))) %>%
        dplyr::arrange(desc(num_samples), symbol) %>%
        dplyr::filter(num_samples > 1) %>%
        dplyr::filter(!str_detect(chr, "CHR_")) %>%
        split(.$chr)
    } else {
      list()
    }

  write_xlsx(total_diffex_clones_by_chr, total_by_chr_xlsx)

  return(list("total" = c(total_xlsx, total_by_chr_xlsx), "cluster" = c(cluster_xlsx, cluster_by_chr_xlsx)))
}
#' Generate a volcano plot for differential expression data
#'
#' @param clone_diffex Parameter for clone diffex
#' @param sample_id Parameter for sample id
#' @return ggplot2 plot object
#' @export
volcano_plot_clone_clusters <- function(clone_diffex, sample_id) {
  #
  myres <- clone_diffex %>%
    group_by(clone_comparison, cluster)

  myres <-
    myres %>%
    group_split() %>%
    map(tibble::column_to_rownames, "symbol") %>%
    identity()

  names(myres) <- myres %>%
    map_chr(~ glue("{unique(.x$clone_comparison)} {unique(.x$cluster)}"))

  volcano_plots <- myres %>%
    map(dplyr::mutate, diffex_comparison = str_replace(str_extract(clone_comparison, "[0-9]_v_[0-9]"), "_v_", "_")) %>%
    imap(make_volcano_plots, sample_id)

  return(volcano_plots)
}

#' Generate a volcano plot for differential expression data
#'
#' @param clone_diffex Parameter for clone diffex
#' @param sample_id Parameter for sample id
#' @return ggplot2 plot object
#' @export
volcano_plot_clones <- function(clone_diffex, sample_id) {
  myres <- clone_diffex %>%
    group_by(clone_comparison)

  myres <-
    myres %>%
    group_split() %>%
    map(tibble::column_to_rownames, "symbol") %>%
    identity()

  names(myres) <-
    myres %>%
    map_chr(~ glue("{unique(.x$clone_comparison)}"))

  volcano_plots <-
    myres %>%
    map(dplyr::mutate, diffex_comparison = str_replace(str_extract(clone_comparison, "[0-9]_v_[0-9]"), "_v_", "_")) %>%
    imap(make_volcano_plots, sample_id)

  return(volcano_plots)
}

#' Perform differential expression analysis
#'
#' @param cluster_diffex_clones Cluster information
#' @param cluster_pdf Cluster information
#' @param total_diffex_clones Parameter for total diffex clones
#' @param total_pdf Character string (default: "results/straight_diffex_bw_clones_large.pdf")
#' @return Data frame
#' @export
make_volcano_diffex_clones <- function(cluster_diffex_clones,
                                       cluster_pdf = "results/diffex_bw_clones_per_cluster_large.pdf",
                                       total_diffex_clones,
                                       total_pdf = "results/straight_diffex_bw_clones_large.pdf") {
  # Drop skipped samples (NA_character_) and unwritten files before anything reads
  # them as filenames -- see .drop_missing_paths().
  cluster_diffex_clones <- .drop_missing_paths(cluster_diffex_clones)
  total_diffex_clones   <- .drop_missing_paths(total_diffex_clones)
  if (length(cluster_diffex_clones) == 0 && length(total_diffex_clones) == 0) {
    message("make_volcano_diffex_clones: no usable diffex CSVs, returning NULL")
    return(invisible(NULL))
  }

  kooi_candidates <- read_csv("data/kooi_candidates.csv")

  cc_genes <- Seurat::cc.genes %>%
    tibble::enframe("phase_of_gene", "symbol") %>%
    tidyr::unnest(symbol)

  sample_ids <- str_extract(cluster_diffex_clones, "SR[RX][0-9]+")

  # total ------------------------------
  total_diffex_clones <-
    total_diffex_clones %>%
    set_names(sample_ids) %>%
    map(bind_rows, .id = "clone_comparison") %>%
    purrr::discard(~ nrow(.x) < 1) %>%
    map(dplyr::left_join, cc_genes, by = "symbol") %>%
    map(group_by, clone_comparison) %>%
    # map(dplyr::filter, p_val_adj < 0.05) %>%
    map(dplyr::arrange, clone_comparison, p_val_adj) %>%
    map(dplyr::select, clone_comparison, chr, symbol, description, everything()) %>%
    map(dplyr::left_join, kooi_candidates, by = "symbol") %>%
    identity()

  # cluster------------------------------
  cluster_diffex_clones <-
    cluster_diffex_clones %>%
    set_names(str_extract(., "SR[RX][0-9]+")) %>%
    map(read_csv) %>%
    purrr::keep(~ nrow(.x) > 0) %>%
    map(dplyr::left_join, cc_genes, by = "symbol") %>%
    # map(dplyr::filter, p_val_adj < 0.05) %>%
    map(dplyr::arrange, clone_comparison, cluster, p_val_adj) %>%
    map(dplyr::select, clone_comparison, cluster, chr, symbol, description, everything()) %>%
    map(dplyr::left_join, kooi_candidates, by = "symbol") %>%
    identity()

  cluster_diffex_clones <-
    cluster_diffex_clones %>%
    # map2(total_diffex_clones, annotate_cluster_membership, "is_total_clone_diffex") %>%
    map(dplyr::distinct, clone_comparison, cluster, chr, symbol, .keep_all = TRUE)

  total_diffex_clones <-
    total_diffex_clones %>%
    # map2(cluster_diffex_clones, annotate_cluster_membership, "is_cluster_clone_diffex") %>%
    map(dplyr::distinct, clone_comparison, chr, symbol, .keep_all = TRUE)

  # cluster ------------------------------
  clone_cluster_comparison_volcanos <- imap(cluster_diffex_clones, volcano_plot_clone_clusters)

  pdf(cluster_pdf)
  print(clone_cluster_comparison_volcanos)
  dev.off()

  # total ------------------------------

  clone_comparison_volcanos <- imap(total_diffex_clones, volcano_plot_clones)

  pdf(total_pdf)
  print(clone_comparison_volcanos)
  dev.off()

  return(list("cluster" = cluster_pdf, "total" = total_pdf))
}

