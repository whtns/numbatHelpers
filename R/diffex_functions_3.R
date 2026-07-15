# Diffex Functions (3)

#' Filter data based on specified criteria
#'
#' @param df Input data frame or dataset
#' @param segment_region Parameter for segment region
#' @return Filtered data
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

filter_input_by_region <- function(df, segment_region, oncoprint_settings) {
  #

  region_settings <-
    dplyr::filter(oncoprint_settings, region == segment_region)

  df <- df[[segment_region]]

  test0 <- pmap(list(df, region_settings$scna, region_settings$p_val, region_settings$fc, region_settings$recurrence), filter_input_by_scna)
}

#' Filter data based on specified criteria
#'
#' @param unfiltered_oncoprint_input_by_scna Parameter for unfiltered oncoprint input by scna
#' @param oncoprint_settings Parameter for oncoprint settings
#' @return Filtered data
#' @export
filter_oncoprint_diffex <- function(unfiltered_oncoprint_input_by_scna, oncoprint_settings) {
  oncoprint_input_by_scna <- map(c("cis" = "cis", "trans" = "trans", "all" = "all"), ~ filter_input_by_region(unfiltered_oncoprint_input_by_scna, .x, oncoprint_settings))

  return(oncoprint_input_by_scna)
}

#' Perform differential expression analysis
#'
#' @param large_filter_expressions Parameter for large filter expressions
#' @param cluster_dictionary Cluster information
#' @param interesting_samples Parameter for interesting samples
#' @param cis_diffex_clones Parameter for cis diffex clones
#' @param trans_diffex_clones Parameter for trans diffex clones
#' @param all_diffex_clones Parameter for all diffex clones
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @param rb_scna_samples Parameter for rb scna samples
#' @param by_cluster Cluster information
#' @param ... Additional arguments passed to other functions
#' @return Differential expression results
#' @export
make_oncoprint_diffex <- function(large_filter_expressions, cluster_dictionary, interesting_samples, cis_diffex_clones, trans_diffex_clones, all_diffex_clones, large_clone_comparisons, rb_scna_samples, by_cluster = FALSE, ...) {
  
#' Filter data based on specified criteria
#'
#' @param df Input data frame or dataset
#' @param n_slice Parameter for n slice
#' @return Filtered data
#' @export
filter_diffex <- function(df, n_slice = 10) {
    #
  required_cols <- c("p_val_adj", "symbol", "avg_log2FC")
  if (nrow(df) == 0 || !all(required_cols %in% colnames(df))) return(df)

    test0 <-
      df %>%
      dplyr::arrange(.data$p_val_adj) %>%
      dplyr::group_by(.data$symbol) %>%
      dplyr::mutate(neg_log_p_val_adj = -log(.data$p_val_adj, base = 10)) %>%
      dplyr::mutate(abs_log2FC = abs(.data$avg_log2FC)) %>%
      identity()

    test1 <-
      test0 %>%
      dplyr::summarize(mean_FC = mean(abs(avg_log2FC))) %>%
      dplyr::inner_join(test0, by = "symbol")

    return(test1)
  }

  cis_diffex_clones   <- .name_by_sample(cis_diffex_clones,   interesting_samples, "cis_diffex_clones")
  trans_diffex_clones <- .name_by_sample(trans_diffex_clones, interesting_samples, "trans_diffex_clones")
  all_diffex_clones   <- .name_by_sample(all_diffex_clones,   interesting_samples, "all_diffex_clones")

  
#' Perform differential expression analysis
#'
#' @param diffex_clones Parameter for diffex clones
#' @param selected_scna Character string (default: "1q+")
#' @param rb_scna_samples Parameter for rb scna samples
#' @param ... Additional arguments passed to other functions
#' @return Differential expression results
#' @export
select_scna_diffex <- function(diffex_clones, selected_scna = "1q+", rb_scna_samples, ...) {
    #

    segment <- str_extract(selected_scna, "[0-9]*[a-z]")

    sign <- str_extract(selected_scna, "[+,-]")

    comparisons <-
      diffex_clones[str_extract(names(diffex_clones), "SR[RX][0-9]+") %in% rb_scna_samples[[segment]]] %>%
      map(~ .x[str_detect(names(.x), glue("{segment}\\{sign}$"))]) %>%
      compact() %>%
      map(dplyr::bind_rows, .id = "clone_comparison") %>%
      dplyr::bind_rows(.id = "sample_id") %>%
      # dplyr::filter(seqnames == str_extract(segment, "[0-9]*")) %>%
      filter_diffex(...) %>%
      identity()
  }

  if (by_cluster) {
    cis_diffex_clones <- map(cis_diffex_clones, read_csv) %>%
      purrr::compact() %>%
      map(~ split(.x, .x$clone_comparison))
    trans_diffex_clones <- map(trans_diffex_clones, read_csv) %>%
      purrr::compact() %>%
      map(~ split(.x, .x$clone_comparison))
    all_diffex_clones <- map(all_diffex_clones, read_csv) %>%
      purrr::compact() %>%
      map(~ split(.x, .x$clone_comparison))
  }

  # cis ------------------------------

  comparisons_of_1q <- select_scna_diffex(cis_diffex_clones, "1q+", rb_scna_samples)

  comparisons_of_2p <- select_scna_diffex(cis_diffex_clones, "2p+", rb_scna_samples)

  comparisons_of_6p <- select_scna_diffex(cis_diffex_clones, "6p+", rb_scna_samples)

  comparisons_of_16q <- select_scna_diffex(cis_diffex_clones, "16q-", rb_scna_samples)

  cis_comps <- list(
    comparisons_of_1q,
    comparisons_of_2p,
    comparisons_of_6p,
    comparisons_of_16q
  ) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  inspect_cis_comps <-
    cis_comps %>%
    map(~ dplyr::distinct(.x, symbol, description)) %>%
    dplyr::bind_rows(.id = "clone_comparison")

  # trans ------------------------------
  #
  comparisons_of_1q <- select_scna_diffex(trans_diffex_clones, "1q+", rb_scna_samples)

  comparisons_of_2p <- select_scna_diffex(trans_diffex_clones, "2p+", rb_scna_samples)

  comparisons_of_6p <- select_scna_diffex(trans_diffex_clones, "6p+", rb_scna_samples)

  comparisons_of_16q <- select_scna_diffex(trans_diffex_clones, "16q-", rb_scna_samples)

  trans_comps <- list(
    comparisons_of_1q,
    comparisons_of_2p,
    comparisons_of_6p,
    comparisons_of_16q
  ) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  # all ------------------------------

  comparisons_of_1q <- select_scna_diffex(all_diffex_clones, "1q+", rb_scna_samples)

  comparisons_of_2p <- select_scna_diffex(all_diffex_clones, "2p+", rb_scna_samples)

  comparisons_of_6p <- select_scna_diffex(all_diffex_clones, "6p+", rb_scna_samples)

  comparisons_of_16q <- select_scna_diffex(all_diffex_clones, "16q-", rb_scna_samples)

  all_comps <- list(
    comparisons_of_1q,
    comparisons_of_2p,
    comparisons_of_6p,
    comparisons_of_16q
  ) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  if (by_cluster) {
    cis_comps <-
      cis_comps %>%
      map(tidyr::unite, "sample_id", sample_id, cluster, sep = "-")

    trans_comps <-
      trans_comps %>%
      map(tidyr::unite, "sample_id", sample_id, cluster, sep = "-")

    all_comps <-
      all_comps %>%
      map(tidyr::unite, "sample_id", sample_id, cluster, sep = "-")
  }

  return(list("cis" = cis_comps, "trans" = trans_comps, "all" = all_comps))
}