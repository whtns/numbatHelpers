# Plot Functions (137)

#' Perform append clone nums operation
#'
#' @param diffex Parameter for diffex
#' @param clone_comparison Parameter for clone comparison
#' @param seu Seurat object
#' @return Function result
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging
# - multiple_joins: Combine multiple joins into single join operation where possible

append_clone_nums <- function(diffex, clone_comparison, seu) {
  if (nrow(diffex) == 0) return(diffex)

  idents <-
    clone_comparison %>%
    str_extract("[0-9]_v_[0-9]") %>%
    str_split(pattern = "_v_") %>%
    set_names(clone_comparison) %>%
    identity()

  clone_nums <- map_chr(idents, ~ {
    paste(table(seu@meta.data[["clone_opt"]])[.x], collapse = "_v_")
  }) %>%
    identity()

  diffex0 <-
    diffex %>%
    dplyr::mutate(clone_nums = clone_nums)

  return(diffex0)
}

#' Perform enrichment analysis
#'
#' @param df Input data frame or dataset
#' @param ... Additional arguments passed to other functions
#' @return Enrichment analysis results
#' @export
prep_for_enrichment <- function(df, ...) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    return(NULL)
  }

  if (!("symbol" %in% colnames(df))) {
    return(NULL)
  }

  enrich_table <-
    df %>%
    dplyr::select(any_of(c("symbol", "avg_log2FC", "p_val", "clone_nums", "clone_comparison"))) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    tibble::column_to_rownames("symbol") %>%
    enrichment_analysis(...)

  return(enrich_table)
}

#' Create a plot visualization
#'
#' @param enrichment_table_list Parameter for enrichment table list
#' @param scna Parameter for scna
#' @param input_plot_file File path
#' @return ggplot2 plot object
#' @export
plot_enrichment_per_scna <- function(enrichment_table_list, scna, input_plot_file) {
  scna_plot_file <- str_replace(input_plot_file, ".pdf", glue("_{scna}.pdf"))

  pdf(scna_plot_file)
  enrichment_table_list %>%
    purrr::list_flatten() %>%
    purrr::compact() %>%
    purrr::discard(~ nrow(.x) == 0) %>%
    map(clusterProfiler::dotplot) %>%
    imap(~ (.x + labs(title = .y))) %>%
    map(print) %>%
    identity()
  dev.off()

  return(scna_plot_file)
}

#' Create a plot visualization
#'
#' @param enrichment_tables Parameter for enrichment tables
#' @param num_recur Parameter for num recur
#' @param mytitle Plot title
#' @param n_slice Parameter for n slice
#' @param by_cluster Cluster information
#' @param pvalueCutoff Threshold value for filtering
#' @return ggplot2 plot object
#' @export
plot_enrichment_recurrence <- function(enrichment_tables, num_recur = 2, mytitle = "", n_slice = 10, by_cluster = FALSE, pvalueCutoff = 0.3) {
  #

  names(enrichment_tables) <- str_replace_all(names(enrichment_tables), "_", "-")

  for (sample in names(enrichment_tables)) {
    names(enrichment_tables[[sample]]) <- str_replace_all(names(enrichment_tables[[sample]]), "_", "-")
  }

  df <-
    enrichment_tables %>%
    purrr::list_flatten() %>%
    purrr::discard(is.na) %>%
    map(~ {
      .x@result
    }) %>%
    dplyr::bind_rows(.id = "clone_comparison") %>%
    identity()

  test0 <-
    df %>%
    # dplyr::arrange(symbol, sample_id) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::filter(p.adjust <= pvalueCutoff) %>%
    dplyr::group_by(Description) %>%
    dplyr::mutate(neg_log10_p_val_adj = -log(p.adjust, base = 10)) %>%
    dplyr::mutate(n_samples = n_distinct(clone_comparison)) %>%
    dplyr::arrange(desc(n_samples), Description)

  test0 <-
    test0 %>%
    dplyr::select(Description, core_enrichment) %>%
    dplyr::mutate(core_enrichment = str_split(core_enrichment, pattern = "/")) %>%
    tidyr::unnest(core_enrichment) %>%
    dplyr::mutate(core_enrichment = as.integer(core_enrichment)) %>%
    dplyr::left_join(annotables::grch38[c("symbol", "entrez")], by = c("core_enrichment" = "entrez")) %>%
    dplyr::select(Description, symbol) %>%
    dplyr::group_by(Description) %>%
    dplyr::summarize(genes = list(symbol), set_size = n_distinct(symbol)) %>%
    dplyr::mutate(genes = purrr::map_chr(genes, ~ paste(unique(.x), collapse = ","))) %>%
    dplyr::left_join(test0, by = "Description") %>%
    dplyr::select(-c("set_size", "setSize")) %>%
    identity()


  test1 <-
    test0 %>%
    dplyr::filter(p.adjust < 0.5) %>%
    dplyr::group_by(Description) %>%
    # dplyr::filter(n_distinct(clone_comparison) >= num_recur) %>%
    dplyr::summarize(mean_NES = mean(NES)) %>%
    # dplyr::slice_max(abs(mean_NES), n = n_slice) %>%
    dplyr::inner_join(test0, by = c("Description")) %>%
    dplyr::mutate(comparison = str_remove(clone_comparison, "SR[RX][0-9]+_")) %>%
    dplyr::mutate(comparison = factor(str_remove(comparison, "_.*"))) %>%
    # dplyr::mutate(clone_comparison = str_replace_all(clone_comparison, "_", "\n")) %>%
    dplyr::mutate(clone_comparison = factor(clone_comparison)) %>%
    dplyr::mutate(clone_comparison = fct_reorder(clone_comparison, as.integer(comparison))) %>%
    dplyr::mutate(genes = str_split(genes, ",")) %>%
    # dplyr::mutate(phase = str_extract(clone_comparison, "(?<=_).*(?=-[0-9]*)")) %>%
    # dplyr::mutate(sample_id = str_extract(clone_comparison, ".*(?=_)")) %>%
    # dplyr::mutate(sample_id = factor(sample_id)) %>%
    identity()

  color_scale_lim <- ceiling(max(abs(test1$NES)))

  color_scale_lim <- ifelse(is.infinite(color_scale_lim), 2, color_scale_lim)

  enrich_plot <-
    test1 %>%
    ggplot(aes(x = clone_comparison, y = Description, size = neg_log10_p_val_adj, color = NES)) + # by_cluster
    scale_color_gradient2(breaks = seq((-color_scale_lim), color_scale_lim, 1), limits = c((-color_scale_lim - 0.5), (color_scale_lim + 0.5))) +
    # scale_color_gradient2() +
    geom_point() +
    labs(title = mytitle, size = "-log10 p.adj") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(str_replace_all(x, "_", " "), width = 20)) +
    # scale_y_discrete(labels = function(x) str_wrap(x, width = 4, whitespace_only = FALSE)) +
    labs(color = "enrichment") +
    # facet_wrap(~comparison) +
    NULL

  # upset_plot <-
  #   test1 %>%
  #   ggplot(aes(x=genes, y = Description)) +
  #   geom_point() +
  #   scale_x_upset() +
  #   theme_minimal() +
  #   labs(title = glue("{mytitle} Gene set overlap"))

  return(list("table" = test0, "enrich_plot" = enrich_plot))
}

#' Create a plot visualization
#'
#' @param enrichment_tables Parameter for enrichment tables
#' @param num_recur Parameter for num recur
#' @param mytitle Plot title
#' @param n_slice Parameter for n slice
#' @param by_cluster Cluster information
#' @param pvalueCutoff Threshold value for filtering
#' @return ggplot2 plot object
#' @export
plot_enrichment_recurrence_by_cluster <- function(enrichment_tables, num_recur = 2, mytitle = "", n_slice = 10, by_cluster = TRUE, pvalueCutoff = 0.3) {
  #

  names(enrichment_tables) <- str_replace_all(names(enrichment_tables), "_", "-")

  for (sample in names(enrichment_tables)) {
    names(enrichment_tables[[sample]]) <- str_replace_all(names(enrichment_tables[[sample]]), "_", "-")
  }

  df <-
    enrichment_tables %>%
    purrr::list_flatten() %>%
    purrr::discard(is.na) %>%
    map(~ {
      .x@result
    }) %>%
    dplyr::bind_rows(.id = "clone_comparison") %>%
    identity()

  test0 <-
    df %>%
    # dplyr::arrange(symbol, sample_id) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::filter(p.adjust <= pvalueCutoff) %>%
    dplyr::group_by(Description) %>%
    dplyr::mutate(neg_log10_p_val_adj = -log(p.adjust, base = 10)) %>%
    dplyr::mutate(n_samples = n_distinct(clone_comparison)) %>%
    dplyr::arrange(desc(n_samples), Description)

  test0 <-
    test0 %>%
    dplyr::select(Description, core_enrichment) %>%
    dplyr::mutate(core_enrichment = str_split(core_enrichment, pattern = "/")) %>%
    tidyr::unnest(core_enrichment) %>%
    dplyr::mutate(core_enrichment = as.integer(core_enrichment)) %>%
    dplyr::left_join(annotables::grch38[c("symbol", "entrez")], by = c("core_enrichment" = "entrez")) %>%
    dplyr::select(Description, symbol) %>%
    dplyr::group_by(Description) %>%
    dplyr::summarize(genes = list(symbol), set_size = n_distinct(symbol)) %>%
    dplyr::mutate(genes = purrr::map_chr(genes, ~ paste(unique(.x), collapse = ","))) %>%
    dplyr::left_join(test0, by = "Description") %>%
    dplyr::select(-c("set_size", "setSize")) %>%
    identity()


  test1 <-
    test0 %>%
    dplyr::filter(p.adjust < 0.5) %>%
    dplyr::group_by(Description) %>%
    # dplyr::filter(n_distinct(clone_comparison) >= num_recur) %>%
    dplyr::summarize(mean_NES = mean(NES)) %>%
    # dplyr::slice_max(abs(mean_NES), n = n_slice) %>%
    dplyr::inner_join(test0, by = c("Description")) %>%
    dplyr::mutate(comparison = str_remove(clone_comparison, "SR[RX][0-9]+_")) %>%
    dplyr::mutate(comparison = factor(str_remove(comparison, "_.*"))) %>%
    # dplyr::mutate(clone_comparison = str_replace_all(clone_comparison, "_", "\n")) %>%
    dplyr::mutate(clone_comparison = factor(clone_comparison)) %>%
    dplyr::mutate(clone_comparison = fct_reorder(clone_comparison, as.integer(comparison))) %>%
    dplyr::mutate(genes = str_split(genes, ",")) %>%
    dplyr::mutate(phase = str_extract(clone_comparison, "(?<=_).*(?=-[0-9]*)")) %>%
    dplyr::mutate(sample_id = str_extract(clone_comparison, ".*(?=_)")) %>%
    dplyr::mutate(sample_id = factor(sample_id)) %>%
    identity()

  phase_levels <- c("pm", "g1", "g1-s", "s", "s-g2", "g2", "g2-m", "hsp", "hypoxia", "other", "s-2")

  phase_levels <- phase_levels[phase_levels %in% unique(test1$phase)]

  color_scale_lim <- ceiling(max(abs(test1$NES)))

  color_scale_lim <- ifelse(is.infinite(color_scale_lim), 2, color_scale_lim)

  enrich_plot_by_phase <-
    test1 %>%
    dplyr::mutate(phase = factor(phase, levels = phase_levels)) %>% # by_cluster
    dplyr::mutate(phase = as.numeric(phase)) %>% # by_cluster
    ggplot(aes(x = fct_reorder(clone_comparison, phase), y = Description, size = neg_log10_p_val_adj, color = NES)) + # by_cluster
    scale_color_gradient2(breaks = seq((-color_scale_lim), color_scale_lim, 1), limits = c((-color_scale_lim - 0.5), (color_scale_lim + 0.5))) +
    # scale_color_gradient2() +
    geom_point() +
    labs(title = mytitle, size = "-log10 p.adj") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    labs(color = "enrichment") +
    # facet_wrap(~comparison) +
    NULL

  enrich_plot_by_sample <-
    test1 %>%
    dplyr::mutate(sample_id = factor(sample_id)) %>% # by_cluster
    dplyr::mutate(sample_id = as.numeric(sample_id)) %>% # by_cluster
    ggplot(aes(x = fct_reorder(clone_comparison, sample_id), y = Description, size = neg_log10_p_val_adj, color = NES)) + # by_cluster
    scale_color_gradient2(breaks = seq((-color_scale_lim), color_scale_lim, 1), limits = c((-color_scale_lim - 0.5), (color_scale_lim + 0.5))) +
    # scale_color_gradient2() +
    geom_point() +
    labs(title = mytitle, size = "-log10 p.adj") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    labs(color = "enrichment") +
    # facet_wrap(~comparison) +
    NULL

  # upset_plot <-
  #   test1 %>%
  #   ggplot(aes(x=genes, y = Description)) +
  #   geom_point() +
  #   scale_x_upset() +
  #   theme_minimal() +
  #   labs(title = glue("{mytitle} Gene set overlap"))

  return(list("table" = test0, "enrich_plot_by_phase" = enrich_plot_by_phase, "enrich_plot_by_sample" = enrich_plot_by_sample))
}

