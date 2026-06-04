# Plot Functions (134)

#' Perform enrichment analysis
#'
#' @param large_filter_expressions Parameter for large filter expressions
#' @param cluster_dictionary Cluster information
#' @param interesting_samples Parameter for interesting samples
#' @param cis_diffex_clones Parameter for cis diffex clones
#' @param trans_diffex_clones Parameter for trans diffex clones
#' @param all_diffex_clones Parameter for all diffex clones
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @param ... Additional arguments passed to other functions
#' @return List object
#' @export
enrich_oncoprints <- function(large_filter_expressions, cluster_dictionary, interesting_samples, cis_diffex_clones, trans_diffex_clones, all_diffex_clones, large_clone_comparisons, ...) {
  filter_diffex_for_recurrence <- function(df, num_recur = 2, n_slice = 10) {
    #

    test0 <-
      df %>%
      # dplyr::arrange(symbol, sample_id) %>%
      dplyr::arrange(.data$p_val_adj) %>%
      dplyr::group_by(.data$symbol) %>%
      dplyr::mutate(neg_log_p_val_adj = -log(.data$p_val_adj, base = 10)) %>%
      dplyr::mutate(abs_log2FC = abs(.data$avg_log2FC)) %>%
      dplyr::filter(.data$p_val_adj < 0.05) %>%
      dplyr::filter(n_distinct(.data$sample_id) >= num_recur) %>%
      identity()

    test1 <-
      test0 %>%
      dplyr::summarize(mean_FC = mean(abs(.data$avg_log2FC))) %>%
      # dplyr::slice_max(abs(mean_FC), n = n_slice) %>%
      dplyr::inner_join(test0, by = "symbol")

    return(test1)
  }

  # cis ------------------------------
  names(cis_diffex_clones) <- interesting_samples

  add_clone_comparison_column <- function(df, mycomparison) {
    dplyr::mutate(df, clone_comparison = mycomparison)
  }

  cis_diffex_clones <- purrr::map(cis_diffex_clones, ~ purrr::imap(.x, add_clone_comparison_column))

  clone_comparisons <- map(cis_diffex_clones, names)


  comparisons_of_1q <- map(clone_comparisons, str_detect, "1q\\+$") %>%
    map2(cis_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(~ map(.x, dplyr::filter, seqnames == "1")) %>%
    identity()

  comparisons_of_2p <- map(clone_comparisons, str_detect, "2p\\+$") %>%
    map2(cis_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(~ map(.x, dplyr::filter, seqnames == "2")) %>%
    identity()

  comparisons_of_6p <- map(clone_comparisons, str_detect, "6p\\+$") %>%
    map2(cis_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(~ map(.x, dplyr::filter, seqnames == "6")) %>%
    identity()

  comparisons_of_16q <- map(clone_comparisons, str_detect, "[0-9]_v_[0-9]_16q\\-$") %>%
    map2(cis_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(~ map(.x, dplyr::filter, seqnames == "16")) %>%
    identity()

  cis_comps <- list(
    comparisons_of_1q,
    comparisons_of_2p,
    comparisons_of_6p,
    comparisons_of_16q
  ) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  # trans ------------------------------
  names(trans_diffex_clones) <- interesting_samples

  trans_diffex_clones <- purrr::map(trans_diffex_clones, ~ purrr::imap(.x, add_clone_comparison_column))

  clone_comparisons <- map(trans_diffex_clones, names)

  comparisons_of_1q <- map(clone_comparisons, str_detect, "1q\\+$") %>%
    map2(trans_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    identity()

  comparisons_of_2p <- map(clone_comparisons, str_detect, "2p\\+$") %>%
    map2(trans_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    identity()

  comparisons_of_6p <- map(clone_comparisons, str_detect, "6p\\+$") %>%
    map2(trans_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    identity()

  comparisons_of_16q <- map(clone_comparisons, str_detect, "[0-9]_v_[0-9]_16q\\-$") %>%
    map2(trans_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    identity()

  trans_comps <- list(
    comparisons_of_1q,
    comparisons_of_2p,
    comparisons_of_6p,
    comparisons_of_16q
  ) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  # all ------------------------------
  names(all_diffex_clones) <- interesting_samples

  all_diffex_clones <- purrr::map(all_diffex_clones, ~ purrr::imap(.x, add_clone_comparison_column))

  clone_comparisons <- map(all_diffex_clones, names)

  comparisons_of_1q <- map(clone_comparisons, str_detect, "1q\\+$") %>%
    map2(all_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    identity()

  comparisons_of_2p <- map(clone_comparisons, str_detect, "2p\\+$") %>%
    map2(all_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    identity()

  comparisons_of_6p <- map(clone_comparisons, str_detect, "6p\\+$") %>%
    map2(all_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    identity()

  comparisons_of_16q <- map(clone_comparisons, str_detect, "[0-9]_v_[0-9]_16q\\-$") %>%
    map2(all_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    identity()

  all_comps <- list(
    comparisons_of_1q,
    comparisons_of_2p,
    comparisons_of_6p,
    comparisons_of_16q
  ) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  # proceed ------------------------------

  possible_prep_for_enrichment <- purrr::possibly(prep_for_enrichment, otherwise = NA_real_)

  cis_enrichment_tables <- modify_depth(cis_comps, 3, possible_prep_for_enrichment, .ragged = TRUE, ...)

  trans_enrichment_tables <- modify_depth(trans_comps, 3, possible_prep_for_enrichment, .ragged = TRUE, ...)

  all_enrichment_tables <- modify_depth(all_comps, 3, possible_prep_for_enrichment, .ragged = TRUE, ...)


  return(list("cis" = cis_enrichment_tables, "trans" = trans_enrichment_tables, "all" = all_enrichment_tables))
}

#' Perform enrichment analysis
#'
#' @param large_filter_expressions Parameter for large filter expressions
#' @param cluster_dictionary Cluster information
#' @param interesting_samples Parameter for interesting samples
#' @param cis_diffex_clones_for_each_cluster Cluster information
#' @param trans_diffex_clones_for_each_cluster Cluster information
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @param ... Additional arguments passed to other functions
#' @return List object
#' @export
enrich_oncoprints_clusters <- function(large_filter_expressions, cluster_dictionary, interesting_samples, cis_diffex_clones_for_each_cluster, trans_diffex_clones_for_each_cluster, large_clone_comparisons, ...) {
  cis_diffex_clones_for_each_cluster <- map(cis_diffex_clones_for_each_cluster, read_csv)

  trans_diffex_clones_for_each_cluster <- map(trans_diffex_clones_for_each_cluster, read_csv)

  filter_diffex_for_recurrence <- function(df, num_recur = 2, n_slice = 10) {
    #

    test0 <-
      df %>%
      # dplyr::arrange(symbol, sample_id) %>%
      dplyr::arrange(.data$p_val_adj) %>%
      dplyr::group_by(.data$symbol) %>%
      dplyr::mutate(neg_log_p_val_adj = -log(.data$p_val_adj, base = 10)) %>%
      dplyr::mutate(abs_log2FC = abs(.data$avg_log2FC)) %>%
      dplyr::filter(.data$p_val_adj < 0.05) %>%
      dplyr::filter(n_distinct(.data$sample_id) >= num_recur) %>%
      identity()

    test1 <-
      test0 %>%
      dplyr::summarize(mean_FC = mean(abs(.data$avg_log2FC))) %>%
      # dplyr::slice_max(abs(mean_FC), n = n_slice) %>%
      dplyr::inner_join(test0, by = "symbol")

    return(test1)
  }

  # cis ------------------------------
  names(cis_diffex_clones_for_each_cluster) <- interesting_samples

  cis_diffex_clones_for_each_cluster <-
    cis_diffex_clones_for_each_cluster %>%
    keep(~ nrow(.x) > 0)

  clone_comparisons <- map(cis_diffex_clones_for_each_cluster, ~ unique(.x[["clone_comparison"]]))

  shape_cluster_comparison <- function(mylist) {
    mylist %>%
      purrr::discard(~ nrow(.x) == 0) %>%
      map(~ split(.x, .x$cluster)) %>%
      identity()
  }

  comparisons_of_1q <- map(cis_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "1q\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_2p <- map(cis_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "2p\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_6p <- map(cis_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "6p\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_16q <- map(cis_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "[0-9]_v_[0-9]_16q\\-$")) %>%
    shape_cluster_comparison()

  cis_comps <- list(
    comparisons_of_1q,
    comparisons_of_2p,
    comparisons_of_6p,
    comparisons_of_16q
  ) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  # trans ------------------------------

  names(trans_diffex_clones_for_each_cluster) <- interesting_samples

  trans_diffex_clones_for_each_cluster <-
    trans_diffex_clones_for_each_cluster %>%
    keep(~ nrow(.x) > 0)

  clone_comparisons <- map(trans_diffex_clones_for_each_cluster, ~ unique(.x[["clone_comparison"]]))

  comparisons_of_1q <- map(trans_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "1q\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_2p <- map(trans_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "2p\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_6p <- map(trans_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "6p\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_16q <- map(trans_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "[0-9]_v_[0-9]_16q\\-$")) %>%
    shape_cluster_comparison()

  trans_comps <- list(
    comparisons_of_1q,
    comparisons_of_2p,
    comparisons_of_6p,
    comparisons_of_16q
  ) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  # proceed ------------------------------

  possible_prep_for_enrichment <- purrr::possibly(prep_for_enrichment, otherwise = NA_real_)

  cis_enrichment_tables <- modify_depth(cis_comps, 3, possible_prep_for_enrichment, .ragged = TRUE, ...)

  trans_enrichment_tables <- modify_depth(trans_comps, 3, possible_prep_for_enrichment, .ragged = TRUE, ...)

  # cis_enrichment_plots <- modify_depth(cis_enrichment_tables, 3, plot_enrichment, .ragged = TRUE) %>%
  #   purrr::list_flatten() %>%
  #   purrr::list_flatten() %>%
  #   imap(~{.x + labs(title = .y)}) %>%
  #   identity()
  #
  # pdf("cis_cluster.pdf")
  # cis_enrichment_plots
  # dev.off()
  #
  # trans_enrichment_plots <- modify_depth(trans_enrichment_tables, 3, plot_enrichment, .ragged = TRUE) %>%
  #   purrr::list_flatten() %>%
  #   purrr::list_flatten() %>%
  #   imap(~{.x + labs(title = .y)}) %>%
  #   identity()
  #
  # pdf("trans_cluster.pdf")
  # trans_enrichment_plots
  # dev.off()

  return(list("cis" = cis_enrichment_tables, "trans" = trans_enrichment_tables))
}