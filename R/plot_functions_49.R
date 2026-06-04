# Plot Functions (148)

#' Perform differential expression analysis
#'
#' @param clustree_diffexes Parameter for clustree diffexes
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

find_candidate_cis_in_clustree_diffexes <- function(clustree_diffexes) {
  sample_id_list <- clustree_diffexes %>%
    map("diffex") %>%
    map("sample_id") %>%
    map_chr(unique)

  to_clust_list <- clustree_diffexes %>%
    map("diffex") %>%
    map("to_clust") %>%
    map_chr(unique) %>%
    identity()

  table_names <- glue("{sample_id_list}_{to_clust_list}")

  test1 <-
    clustree_diffexes %>%
    map("diffex") %>%
    map(function(mytable) {
      test0 <-
        mytable %>%
        dplyr::arrange(to_clust, avg_log2FC) %>%
        dplyr::filter(p_val_adj < 0.1, location == "cis") %>%
        dplyr::slice_max(abs(avg_log2FC), n = 5) %>%
        identity()
    }) %>%
    set_names(table_names)

  table_path <- writexl::write_xlsx(test1, "results/cis_divergent_cluster_diffex.xlsx")

  clustree_plots <- clustree_diffexes %>%
    map("plot")

  fs::dir_create("results/cis_clustree_diffex")

  cis_clustree_plots <- clustree_plots %>%
    str_replace(".pdf", "_cis.pdf") %>%
    str_replace("results", "results/cis_clustree_diffex")

  cis_pages_to_extract <-
    clustree_plots %>%
    map(qpdf::pdf_length) %>%
    map(~ seq(2, .x, 3))

  pmap(list(clustree_plots, cis_pages_to_extract, cis_clustree_plots), qpdf::pdf_subset)

  plot_path <- glue("results/divergent_cluster_diffex_cis.pdf")

  qpdf::pdf_combine(cis_clustree_plots, plot_path)

  return(list("table" = table_path, "plot" = plot_path))
}

#' Perform differential expression analysis
#'
#' @param clustree_diffexes Parameter for clustree diffexes
#' @param gene_location Gene names or identifiers
#' @return ggplot2 plot object
#' @export
find_candidate_trans_in_clustree_diffexes <- function(clustree_diffexes, gene_location = "trans") {
  sample_id_list <- clustree_diffexes %>%
    map("diffex") %>%
    map("sample_id") %>%
    map_chr(unique)

  to_clust_list <- clustree_diffexes %>%
    map("diffex") %>%
    map("to_clust") %>%
    map_chr(unique) %>%
    identity()

  table_names <- glue("{sample_id_list}_{to_clust_list}")

  test1 <-
    clustree_diffexes %>%
    map("diffex") %>%
    map(function(mytable) {
      test0 <-
        mytable %>%
        dplyr::arrange(to_clust, avg_log2FC) %>%
        dplyr::filter(p_val_adj < 0.1, location == gene_location) %>%
        dplyr::slice_max(abs(avg_log2FC), n = 5) %>%
        identity()
    }) %>%
    set_names(table_names)

  table_path <- writexl::write_xlsx(test1, "results/trans_divergent_cluster_diffex.xlsx")

  clustree_plots <- clustree_diffexes %>%
    map("plot")

  fs::dir_create("results/trans_clustree_diffex")

  trans_clustree_plots <- clustree_plots %>%
    str_replace(".pdf", "_trans.pdf") %>%
    str_replace("results", "results/trans_clustree_diffex")

  trans_pages_to_extract <-
    clustree_plots %>%
    map(qpdf::pdf_length) %>%
    map(~ seq(3, .x, 3))

  pmap(list(clustree_plots, trans_pages_to_extract, trans_clustree_plots), qpdf::pdf_subset)

  plot_path <- glue("results/divergent_cluster_diffex_trans.pdf")

  qpdf::pdf_combine(trans_clustree_plots, plot_path)

  return(list("table" = table_path, "plot" = plot_path))
}

#' Perform differential expression analysis
#'
#' @param clustree_diffexes Parameter for clustree diffexes
#' @param gene_location Gene names or identifiers
#' @return ggplot2 plot object
#' @export
find_candidate_all_in_clustree_diffexes <- function(clustree_diffexes, gene_location = "all") {
  sample_id_list <- clustree_diffexes %>%
    map("diffex") %>%
    map("sample_id") %>%
    map_chr(unique)

  to_clust_list <- clustree_diffexes %>%
    map("diffex") %>%
    map("to_clust") %>%
    map_chr(unique) %>%
    identity()

  table_names <- glue("{sample_id_list}_{to_clust_list}")

  test1 <-
    clustree_diffexes %>%
    map("diffex") %>%
    map(function(mytable) {
      test0 <-
        mytable %>%
        dplyr::arrange(to_clust, avg_log2FC) %>%
        dplyr::filter(p_val_adj < 0.1, location == gene_location) %>%
        dplyr::slice_max(abs(avg_log2FC), n = 5) %>%
        identity()
    }) %>%
    set_names(table_names)

  table_path <- writexl::write_xlsx(test1, "results/all_divergent_cluster_diffex.xlsx")

  clustree_plots <- clustree_diffexes %>%
    map("plot")

  fs::dir_create("results/all_clustree_diffex")

  all_clustree_plots <- clustree_plots %>%
    str_replace(".pdf", "_all.pdf") %>%
    str_replace("results", "results/all_clustree_diffex")

  all_pages_to_extract <-
    clustree_plots %>%
    map(qpdf::pdf_length) %>%
    map(~ seq(1, .x, 3))

  pmap(list(clustree_plots, all_pages_to_extract, all_clustree_plots), qpdf::pdf_subset)

  plot_path <- glue("results/divergent_cluster_diffex_all.pdf")

  qpdf::pdf_combine(all_clustree_plots, plot_path)

  return(list("table" = table_path, "plot" = plot_path))
}

#' Perform differential expression analysis
#'
#' @param oncoprint_input_by_scna_unfiltered Parameter for oncoprint input by scna unfiltered
#' @return Differential expression results
#' @export
select_genes_from_arbitrary_diffex <- function(oncoprint_input_by_scna_unfiltered) {
  samples_1q_without_preceding_16q <- c("SRX11133594", "SRX11133593", "SRX11133592")

  
#' Perform arrange by recurrence operation
#'
#' @param df Input data frame or dataset
#' @return Function result
#' @export
arrange_by_recurrence <- function(df) {
    df %>%
      dplyr::group_by(symbol) %>%
      dplyr::mutate(abs_mean_FC = abs(mean(avg_log2FC))) %>%
      dplyr::mutate(recurrence = dplyr::n()) %>%
      dplyr::arrange(desc(recurrence), desc(abs_mean_FC))
  }


  test0 <-
    oncoprint_input_by_scna_unfiltered$cis$`1q+` %>%
    dplyr::filter(sample_id %in% samples_1q_without_preceding_16q) %>%
    arrange_by_recurrence() %>%
    identity()
}

