# Plot Functions (143)

#' Create a plot visualization
#'
#' @param all_cells_meta Cell identifiers or information
#' @param scna_meta Parameter for scna meta
#' @param qc_meta Parameter for qc meta
#' @param cell_type_meta Cell identifiers or information
#' @param sample_id Parameter for sample id
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

plot_filtering_timeline <- function(all_cells_meta, scna_meta, qc_meta, cell_type_meta, sample_id) {
  #

  n_cells <- nrow(all_cells_meta)

  all_cells_meta$sample_id <- sample_id

  all_cells_bar <-
    all_cells_meta %>%
    tibble::rownames_to_column("cell") %>%
    ggplot(fill = "gray", aes(x = sample_id)) +
    geom_bar(position = "stack", width = 0.1) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt"),
      legend.position = "bottom"
    ) +
    ylim(0, n_cells) +
    labs(title = "all cells") +
    scale_x_discrete(expand = c(0, 0)) +
    # scale_y_log10() +
    NULL

  # scna_meta <- dplyr::mutate(scna_meta, scna = ifelse(scna == "", "none", scna))
  #
  scna_levels <- levels(factor(scna_meta$scna))

  scna_pal <- scales::hue_pal()(length(scna_levels)) %>%
    set_names(scna_levels) %>%
    identity()

  scna_bar <-
    scna_meta %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(scna = wrap_scna_labels(scna)) %>%
    ggplot(aes(fill = scna, x = sample_id)) +
    geom_bar(position = "stack", width = 0.1) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt"),
      legend.position = "bottom"
    ) +
    ylim(0, n_cells) +
    # labs(title = "scna") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = scna_pal) +
    # scale_y_log10() +
    NULL

  qc_bar <-
    qc_meta %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(scna = wrap_scna_labels(scna)) %>%
    ggplot(aes(fill = scna, x = sample_id)) +
    geom_bar(position = "stack", width = 0.1) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt"),
      legend.position = "bottom"
    ) +
    ylim(0, n_cells) +
    # labs(title = "qc") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = scna_pal) +
    # scale_y_log10() +
    NULL

  cell_type_bar <-
    cell_type_meta %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(scna = wrap_scna_labels(scna)) %>%
    ggplot(aes(fill = .data[["scna"]], x = sample_id)) +
    geom_bar(position = "stack", width = 0.2) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt"),
      # panel.grid = element_blank(),
      # panel.border = element_blank(),
      legend.position = "bottom"
    ) +
    ylim(0, n_cells) +
    # labs(title = "bad cells") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = scna_pal) +
    theme(legend.position = "none") +
    # scale_y_log10() +
    NULL

  wrap_plots(scna_bar, qc_bar, cell_type_bar, nrow = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(title = sample_id)
}

#' Extract or pull specific data elements
#'
#' @param filtered_seus Parameter for filtered seus
#' @param gene_lists Gene names or identifiers
#' @return Extracted data elements
#' @export
pull_common_markers <- function(filtered_seus, gene_lists) {
  annotated_genes <-
    gene_lists %>%
    tibble::enframe("mp", "symbol") %>%
    tidyr::unnest(symbol) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    dplyr::arrange(mp) %>%
    identity()

  names(filtered_seus) <- str_extract(filtered_seus, "SR[RX][0-9]+")

  # load filtered_seus ------------------------------
  # find cluster markers of every seu, compare with zinovyev markers
  # find common markers for clusters
  
#' Perform clustering analysis
#'
#' @param seu_path File path
#' @return Extracted data elements
#' @export
pull_cluster_markers <- function(seu_path) {
    seu <- readRDS(seu_path)

    table_cluster_markers(seu)
  }

  my_cluster_markers <- map(filtered_seus, pull_cluster_markers)

  common_genes <-
    map(my_cluster_markers, "SCT_snn_res.0.4") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(abs(Average.Log.Fold.Change) > 0.5) %>%
    dplyr::arrange(Gene.Name) %>%
    dplyr::group_by(Gene.Name) %>%
    dplyr::filter(dplyr::n() > 3) %>%
    dplyr::filter(!str_detect(Gene.Name, pattern = "^RP.*")) %>%
    dplyr::filter(!str_detect(Gene.Name, pattern = "^MT.*")) %>%
    dplyr::distinct(Gene.Name) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(annotated_genes, by = c("Gene.Name" = "symbol")) %>%
    # dplyr::slice_sample(n =50) %>%
    # dplyr::pull(Gene.Name) %>%
    # unique() %>%
    # sample(50) %>%
    identity()

  return(common_genes)
}
#' Perform wrap scna labels operation
#'
#' @param scna_labels Parameter for scna labels
#' @return Function result
#' @export
wrap_scna_labels <- function(scna_labels) {
  str_wrap(str_replace_all(scna_labels, ",", " "), width = 10, whitespace_only = FALSE)
}

#' Create a plot visualization
#'
#' @param clone_set Parameter for clone set
#' @param clone_comparison Parameter for clone comparison
#' @param seu Seurat object
#' @param sample_id Parameter for sample id
#' @param var_y Character string (default: "clusters")
#' @return ggplot2 plot object
#' @export
make_pairwise_plots <- function(clone_set, clone_comparison, seu, sample_id, var_y = "clusters") {
  #
  pair_seu <- seu[, seu$scna %in% clone_set]

  clone_ratio <- janitor::tabyl(as.character(pair_seu$scna))$percent[[2]]

  comparison_scna <-
    janitor::tabyl(as.character(pair_seu$scna))[2, 1]

  myplot <- plot_distribution_of_clones_across_clusters(
    pair_seu,
    seu_name = clone_comparison, var_x = "scna", var_y = var_y, avg_line = clone_ratio, signif = TRUE, plot_type = "clone"
  )

  cluster_values <-
    pair_seu@meta.data %>%
    dplyr::group_by(.data[[var_y]], scna) %>%
    dplyr::summarize(value = dplyr::n())

  all_plot_table <-
    cluster_values %>%
    dplyr::group_by(scna) %>%
    dplyr::summarize(value = sum(value)) %>%
    mutate(percent = value / sum(value)) %>%
    dplyr::select(-value) %>%
    dplyr::mutate({{ var_y }} := "all") %>%
    identity()

  mytable <-
    cluster_values %>%
    mutate(percent = value / sum(value)) %>%
    dplyr::select(-value) %>%
    dplyr::bind_rows(all_plot_table) %>%
    tidyr::pivot_wider(names_from = "scna", values_from = "percent") %>%
    dplyr::mutate(sample_id = sample_id) %>%
    identity()

  mytable <-
    mytable %>%
    dplyr::mutate(up = ifelse(.data[[as.character(comparison_scna)]] > (clone_ratio + 0.03), 1, 0)) %>%
    dplyr::mutate(down = ifelse(.data[[as.character(comparison_scna)]] < (clone_ratio - 0.03), 1, 0)) %>%
    identity()

  return(
    list(
      "plot" = myplot,
      "table" = mytable
    )
  )
}

