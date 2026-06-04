#' plot distribution of clones across clusters
#'
#' @param seu
#' @param seu_name
#' @param clusters
#'
#' @return
#' @export
#'
#' @examples
plot_distribution_of_clones_across_clusters <- function(seu, seu_name, var_x = "scna", var_y = "SCT_snn_res.0.6", plot_type = c("both", "clone", "cluster"), avg_line = NULL, signif = FALSE, integrated = FALSE) {
  plot_type <- match.arg(plot_type)

  seu_meta <- seu@meta.data
  summarized_clones <- dplyr::mutate(dplyr::select(seu_meta, .data[[var_x]], .data[[var_y]]), scna = "all")
  cluster_plot <- ggplot(seu_meta) +
    geom_bar(position = "fill", aes(x = .data[[var_x]], fill = .data[[var_y]])) +
    geom_bar(data = summarized_clones, position = "fill", aes(x = .data[[var_x]], fill = .data[[var_y]])) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20), limits = rev) +
    coord_flip()

  summarized_clusters <- dplyr::mutate(dplyr::select(seu_meta, .data[[var_x]], .data[[var_y]]), clusters = "all")
  cluster_levels <- c("all", levels(seu_meta$clusters))

  if (signif) {
    summarized_clusters_tabyl <- summarized_clusters %>%
      dplyr::mutate(scna = as.character(scna)) %>%
      janitor::tabyl(clusters, scna) %>%
      tibble::column_to_rownames("clusters")

    clusters_tabyl <- seu_meta %>%
      dplyr::select(any_of(c("clusters", "scna"))) %>%
      janitor::tabyl(clusters, scna)

    fisher_results <- clusters_tabyl %>%
      tibble::column_to_rownames("clusters") %>%
      dplyr::rowwise() %>%
      dplyr::group_split("clusters") %>%
      lapply(function(x) dplyr::bind_rows(x, summarized_clusters_tabyl)) %>%
      lapply(stats::fisher.test, simulate.p.value = TRUE) %>%
      map_dfr(broom::tidy)

    fisher_results <- clusters_tabyl %>%
      cbind(fisher_results) %>%
      dplyr::mutate(p.adjust = p.adjust(p.value)) %>%
      dplyr::mutate(signif = symnum(p.adjust, corr = FALSE,
                                    cutpoints = c(0, .001, .01, .05, .1, 1),
                                    symbols = c("***", "**", "*", ".", " "))) %>%
      dplyr::select(clusters, signif) %>%
      dplyr::mutate(clusters = factor(clusters, levels = cluster_levels))

    cluster_levels <- tidyr::unite(fisher_results, "clusters", clusters, signif, sep = " ") %>%
      dplyr::pull(clusters)

    clone_input <- seu_meta %>%
      dplyr::bind_rows(summarized_clusters) %>%
      dplyr::left_join(fisher_results, by = "clusters") %>%
      tidyr::unite("clusters", clusters, signif, sep = " ") %>%
      dplyr::mutate(clusters = str_replace(clusters, "all NA", "all")) %>%
      dplyr::mutate(clusters = factor(clusters, levels = c("all", cluster_levels)))

    clone_plot <- ggplot(clone_input) +
      geom_bar(position = "fill", aes(x = .data[[var_y]], fill = .data[[var_x]])) +
      scale_x_discrete(limits = rev) +
      coord_flip() +
      theme_minimal()
  } else {
    clone_plot <- ggplot(seu_meta) +
      geom_bar(position = "fill", aes(x = .data[[var_y]], fill = .data[[var_x]])) +
      geom_bar(data = summarized_clusters, position = "fill", aes(x = .data[[var_y]], fill = .data[[var_x]])) +
      scale_x_discrete(limits = rev) +
      coord_flip() +
      theme_minimal()
  }

  if (!is.null(avg_line)) {
    clone_plot <- clone_plot + geom_hline(yintercept = avg_line)
  }

  plot_return <- switch(plot_type,
    clone = clone_plot,
    cluster = cluster_plot,
    both = (clone_plot / cluster_plot) +
      plot_layout(ncol = 1) +
      plot_annotation(title = seu_name)
  )

  return(plot_return)
}
