# Plot Functions (146)

#' Perform collate clone distribution tables operation
#'
#' @param clone_distribution_plots Parameter for clone distribution plots
#' @param interesting_samples Parameter for interesting samples
#' @return List object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - map_bind_rows: Use map_dfr() instead of map() %>% bind_rows() for better performance

collate_clone_distribution_tables <- function(clone_distribution_plots, interesting_samples) {
  
  
  names(clone_distribution_plots) <- interesting_samples

  cluster_by_clone_tables <- clone_distribution_plots %>%
    map("table") %>%
    map("cluster_by_clone") %>%
    map(dplyr::select, -value) %>%
    map(dplyr::mutate, scna = na_if(scna, "")) %>%
    map(dplyr::mutate, scna = replace_na(scna, "diploid")) %>%
    map(tidyr::pivot_wider, names_from = "scna", values_from = "percent") %>%
    writexl::write_xlsx("results/cluster_by_clone_tables.xlsx") %>%
    identity()

  clone_by_cluster_tables <- clone_distribution_plots %>%
    map("table") %>%
    map("clone_by_cluster") %>%
    map(dplyr::select, -value) %>%
    map(dplyr::mutate, scna = na_if(scna, "")) %>%
    map(dplyr::mutate, scna = replace_na(scna, "diploid")) %>%
    map(tidyr::pivot_wider, names_from = "scna", values_from = "percent") %>%
    writexl::write_xlsx("results/clone_by_cluster_tables.xlsx")

  return(list(
    "cluster_by_clone" = "results/cluster_by_clone_tables.xlsx",
    "clone_by_cluster" = "results/clone_by_cluster_tables.xlsx",
  ))
}

#' Calculate scores for the given data
#'
#' @param clone_distribution_plots Parameter for clone distribution plots
#' @param interesting_samples Parameter for interesting samples
#' @param csv_file File path
#' @return Numeric scores or scored data
#' @export
score_clusters_up_down <- function(clone_distribution_plots, csv_file = "results/cluster_up_down_scorecard.csv") {

  sample_ids <- map_chr(clone_distribution_plots, "sample_id")

  test0 <-
    clone_distribution_plots %>%
    map("table") %>%
    set_names(sample_ids) %>%
    identity()

  
#' Perform clustering analysis
#'
#' @param df Input data frame or dataset
#' @return Function result
#' @export
summarize_clusters <- function(df) {
  
  
    df %>%
      dplyr::select(clusters, up, down) %>%
      tidyr::pivot_longer(-clusters, names_to = "direction", values_to = "whether") %>%
      dplyr::filter(whether == 1) %>%
      dplyr::select(-whether) %>%
      identity()
  }

  test1 <- map(
    test0,
    ~ map(.x, summarize_clusters)
  ) %>%
    map(dplyr::bind_rows, .id = "comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    write_csv(csv_file) %>%
    identity()
}
#' Create a plot visualization
#'
#' @param seu_path File path
#' @param nb_paths File path
#' @param clone_simplifications Parameter for clone simplifications
#' @param label Character string (default: "_clone_tree")
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_clone_tree_from_path <- function(seu_path, nb_paths, clone_simplifications, label = "_clone_tree", ...) {


  seu <- readRDS(seu_path)
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  nb_paths <- nb_paths %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  nb_path <- nb_paths[[tumor_id]]

  clone_tree <- plot_clone_tree(seu, tumor_id = tumor_id, nb_path, clone_simplifications, sample_id = sample_id, ...)

  return(clone_tree)
}

#' Save or write data to file
#'
#' @param seu_path File path
#' @param nb_paths File path
#' @param clone_simplifications Parameter for clone simplifications
#' @param label Character string (default: "_clone_tree")
#' @param ... Additional arguments passed to other functions
#' @return path to a ggplot pdf file
#' @export
save_clone_tree_from_path <- function(seu_path, nb_paths, clone_simplifications, label = "_clone_tree", ...) {

  # Handle NA inputs (e.g., when clone_post is NULL for a sample)
  if (is.na(seu_path)) {
    return(NA_character_)
  }

  seu <- readRDS(seu_path)
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  nb_paths <- nb_paths %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  nb_path <- nb_paths[[tumor_id]]

  p <- plot_clone_tree(seu, tumor_id = tumor_id, nb_path, clone_simplifications, sample_id = sample_id, ...)

  plot_path <- ggsave(glue("results/{sample_id}{label}.pdf"), plot = p, width = 4, height = 4)
  return(plot_path)
}

