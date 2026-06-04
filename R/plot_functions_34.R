# Plot Functions (133)

#' Filter data based on specified criteria
#'
#' @param large_filter_expressions Parameter for large filter expressions
#' @param cluster_dictionary Cluster information
#' @param interesting_samples Parameter for interesting samples
#' @param cis_diffex_clones Parameter for cis diffex clones
#' @param trans_diffex_clones Parameter for trans diffex clones
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @param ... Additional arguments passed to other functions
#' @return List object
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

make_oncoprint_diffex_unfiltered <- function(large_filter_expressions, cluster_dictionary, interesting_samples, cis_diffex_clones, trans_diffex_clones, large_clone_comparisons, ...) {
  # cis ------------------------------
  names(cis_diffex_clones) <- interesting_samples

  names(trans_diffex_clones) <- interesting_samples

  clone_comparisons <- map(cis_diffex_clones, names)


  comparisons_of_1q <- map(clone_comparisons, str_detect, "1q\\+$") %>%
    map2(cis_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(seqnames == "1") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_2p <- map(clone_comparisons, str_detect, "2p\\+$") %>%
    map2(cis_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(seqnames == "2") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_6p <- map(clone_comparisons, str_detect, "6p\\+$") %>%
    map2(cis_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(seqnames == "6") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_16q <- map(clone_comparisons, str_detect, "[0-9]_v_[0-9]_16q\\-$") %>%
    map2(cis_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(seqnames == "16") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

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
  names(trans_diffex_clones) <- interesting_samples

  names(trans_diffex_clones) <- interesting_samples

  clone_comparisons <- map(trans_diffex_clones, names)

  comparisons_of_1q <- map(clone_comparisons, str_detect, "1q\\+$") %>%
    map2(trans_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(abs(avg_log2FC) > 0.25) %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_2p <- map(clone_comparisons, str_detect, "2p\\+$") %>%
    map2(trans_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_6p <- map(clone_comparisons, str_detect, "6p\\+$") %>%
    map2(trans_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_16q <- map(clone_comparisons, str_detect, "[0-9]_v_[0-9]_16q\\-$") %>%
    map2(trans_diffex_clones, ~ {
      .y[.x]
    }) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  trans_comps <- list(
    comparisons_of_1q,
    comparisons_of_2p,
    comparisons_of_6p,
    comparisons_of_16q
  ) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  return(list("cis" = cis_comps, "trans" = trans_comps))
}
#' Perform inspect oncoprints operation
#'
#' @param cis_comps Parameter for cis comps
#' @param trans_comps Parameter for trans comps
#' @param all_comps Parameter for all comps
#' @return List object
#' @export
inspect_oncoprints <- function(cis_comps, trans_comps, all_comps) {
  inspect_trans_comps <-
    trans_comps %>%
    map(~ dplyr::distinct(.x, symbol, description, .keep_all = TRUE)) %>%
    dplyr::bind_rows(.id = "clone_comparison") %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(mean_FC = mean(avg_log2FC)) %>%
    dplyr::arrange(clone_comparison, desc(mean_FC)) %>%
    dplyr::select(symbol, description, everything())

  inspect_cis_comps <-
    cis_comps %>%
    map(~ dplyr::distinct(.x, symbol, description, .keep_all = TRUE)) %>%
    dplyr::bind_rows(.id = "clone_comparison") %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(mean_FC = mean(avg_log2FC)) %>%
    dplyr::arrange(clone_comparison, desc(mean_FC)) %>%
    dplyr::select(symbol, description, everything())

  inspect_all_comps <-
    all_comps %>%
    map(~ dplyr::distinct(.x, symbol, description, .keep_all = TRUE)) %>%
    dplyr::bind_rows(.id = "clone_comparison") %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(mean_FC = mean(avg_log2FC)) %>%
    dplyr::arrange(clone_comparison, desc(mean_FC)) %>%
    dplyr::select(symbol, description, everything())

  return(list("cis" = inspect_cis_comps, "trans" = inspect_trans_comps, "all" = inspect_all_comps))
}

#' Create a plot visualization
#'
#' @param diffex_input Parameter for diffex input
#' @param scna_of_interest Parameter for scna of interest
#' @param segment_region Character string (default: "cis")
#' @param oncoprint_settings Parameter for oncoprint settings
#' @param clone_trees Parameter for clone trees
#' @param n_genes Gene names or identifiers
#' @return ggplot2 plot object
#' @export
plot_recurrence <- function(diffex_input, scna_of_interest, segment_region = "cis", oncoprint_settings, clone_trees, n_genes = 20) {
  #

  required_cols <- c("symbol", "avg_log2FC", "p_val_adj", "sample_id", "description")
  if (is.null(diffex_input) || !is.data.frame(diffex_input) || nrow(diffex_input) == 0 || !all(required_cols %in% colnames(diffex_input))) {
    empty_plot <- ggplot() + theme_void() + labs(title = scna_of_interest, subtitle = "No valid recurrence data")
    empty_table <- tibble::tibble(note = "No valid recurrence data")
    return(list("table" = empty_table, "plot" = empty_plot))
  }

  region_settings <-
    dplyr::filter(oncoprint_settings, region == segment_region) |>
    dplyr::filter(scna == scna_of_interest)

  clone_trees <-
    clone_trees %>%
    set_names(map(., c("labels", "title")))

  mytitle <- scna_of_interest
  mysubtitle <- glue("
  									recurrence: {region_settings$recurrence}
  									p_val <= {region_settings$p_val}
  									fold-change >= {region_settings$fc}
  									")

  n_recurrence <- region_settings$recurrence

  diffex_table <-
    diffex_input %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(mean_FC = mean(avg_log2FC)) %>%
    dplyr::ungroup() %>%
    # dplyr::arrange(desc(mean_FC)) %>%
  	dplyr::arrange(p_val_adj) %>%
    dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) %>%
    dplyr::mutate(neg_log_p_val_adj = -log10(p_val_adj)) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    dplyr::group_by(symbol) %>%
    dplyr::filter(n_distinct(sample_id) >= n_recurrence) %>%
    dplyr::ungroup() %>%
    # dplyr::slice_head(n= 10) %>%
    identity()

  plot_input <-
    diffex_table %>%
    dplyr::filter(p_val_adj <= region_settings$p_val) %>%
    dplyr::filter(abs_log2FC >= region_settings$fc) %>%
    # dplyr::filter(abs_log2FC >= 0.5) %>%
    dplyr::select(sample_id, symbol, description, neg_log_p_val_adj, abs_log2FC, avg_log2FC, mean_FC) %>%
    group_by(symbol) %>%
    dplyr::mutate(num_positive = sum(avg_log2FC > 0)) %>%
    dplyr::mutate(num_negative = sum(avg_log2FC < 0)) %>%
    dplyr::mutate(same_sign = all_same_sign(avg_log2FC)) %>%
    dplyr::mutate(major_sign = abs(num_positive - num_negative)) %>%
    dplyr::filter(same_sign > 0 | major_sign > 1) %>%
    dplyr::mutate(minor_sign = min(num_positive, num_negative)) %>%
    dplyr::filter(major_sign >= n_recurrence) %>%
    # dplyr::filter(minor_sign < 2) %>%
    identity()

  top_genes <- unique(plot_input$symbol)[1:n_genes]
  top_genes <- top_genes[!is.na(top_genes)]


  plot_input <-
    plot_input %>%
    dplyr::filter(symbol %in% top_genes) %>%
    dplyr::filter(minor_sign < 2) %>%
    dplyr::mutate(comparison_sign = ifelse(num_positive > num_negative, 1, -1)) %>%
    dplyr::filter(sign(avg_log2FC) == sign(comparison_sign)) %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(n_samples = n_distinct(sample_id)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(n_samples, desc(mean_FC)) %>%
    dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) %>%
    identity()

  # clone_trees

  diffex_plot <- ggplot(plot_input, aes(x = sample_id, y = symbol)) +
    geom_point(aes(color = neg_log_p_val_adj, size = abs_log2FC)) +
    # geom_tile(fill = NA, color = "black", linewidth = 0.5) +
    labs(title = mytitle, subtitle = mysubtitle, color = "-log10 \n p_adj") +
    theme_bw() +
    # scale_size_continuous(
    #   limits = c(0.1, 0.9)
    # ) +
    scale_color_continuous(
      limits = c(1, 150),
      oob = squish
    ) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
      # legend.position="none"
    )

  sub_clone_trees <-
    clone_trees[names(clone_trees) %in% sort(unique(str_remove(plot_input$sample_id, "-.*")))] %>%
    map(~ {
      .x + theme(plot.title = element_blank())
    }) %>%
    wrap_plots(nrow = 1)

  diffex_plot <- wrap_plots(diffex_plot, sub_clone_trees, ncol = 1)


  return(list("table" = diffex_table, "plot" = diffex_plot))
}

#' Create a plot visualization
#'
#' @param comps Parameter for comps
#' @param clone_trees Parameter for clone trees
#' @param oncoprint_settings Parameter for oncoprint settings
#' @param label Character string (default: "_by_clone")
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
make_oncoprint_plots <- function(comps, clone_trees, oncoprint_settings, label = "_by_clone", ...) {
  #

  comps_res <- list()
  # cis ------------------------------
  for (region in names(comps)) {
    comps_res[[region]] <- imap(comps[[region]], plot_recurrence, region, oncoprint_settings, clone_trees, ...)
  }


  cis_comps_plots <- map(comps_res[["cis"]], "plot")

  # wrap_plots(cis_comps_plots) +  plot_layout(guides = 'collect') +
  #   plot_annotation(title = "top 10 DE genes in segment")

  in_segment_plot_path <- glue("results/diffex_oncoprints{label}_in_segment.pdf")

  pdf(in_segment_plot_path, height = 8, width = 6)
  cis_comps_plots %>%
    walk(print)
  dev.off()

  in_segment_table_path <- glue("results/diffex_oncoprints{label}_in_segment.xlsx")

  cis_comps_tables <-
    comps_res[["cis"]] %>%
    map("table") %>%
    map(~ dplyr::select(.x, -dplyr::any_of("genes_in_segment"))) %>%
    writexl::write_xlsx(in_segment_table_path)

  # trans ------------------------------

  trans_comps_plots <- map(comps_res[["trans"]], "plot")

  # wrap_plots(trans_comps_plots) +  plot_layout(guides = 'collect') +
  #   plot_annotation(title = "top 10 DE genes out of segment")

  out_segment_plot_path <- glue("results/diffex_oncoprints{label}_out_of_segment.pdf")

  pdf(out_segment_plot_path, height = 8, width = 6)
  trans_comps_plots %>%
    walk(print)
  dev.off()

  out_segment_table_path <- glue("results/diffex_oncoprints{label}_out_of_segment.xlsx")

  comps_res[["trans"]] %>%
    map("table") %>%
    writexl::write_xlsx(out_segment_table_path)

  # all ------------------------------

  all_comps_plots <- map(comps_res[["all"]], "plot")

  all_plot_path <- glue("results/diffex_oncoprints{label}_all.pdf")


  pdf(all_plot_path, height = 8, width = 6)
  all_comps_plots %>%
    walk(print)
  dev.off()

  all_table_path <- glue("results/diffex_oncoprints{label}_all.xlsx")

  comps_res[["all"]] %>%
    map("table") %>%
    writexl::write_xlsx(all_table_path)

  return(
    list(
      "cis" =
        list(
          "plot" = in_segment_plot_path,
          "table" = in_segment_table_path
        ),
      "trans" =
        list(
          "plot" = out_segment_plot_path,
          "table" = out_segment_table_path
        ),
      "all" =
        list(
          "plot" = all_plot_path,
          "table" = all_table_path
        )
    )
  )
}

