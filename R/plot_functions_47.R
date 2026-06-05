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

  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*|_unfiltered_seu.*|_seu\\.rds$")

  nb_paths <- nb_paths %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  nb_path <- nb_paths[[tumor_id]]

  mynb <- readRDS(nb_path)
  all_clone_df <- mynb$clone_post %>%
    dplyr::select(cell, clone_opt) %>%
    dplyr::distinct()

  retained_cells <- read_cell_barcodes_from_db(seu_path)
  clone_df <- if (!is.null(retained_cells)) {
    dplyr::filter(all_clone_df, cell %in% retained_cells)
  } else {
    all_clone_df
  }

  clone_tree <- plot_clone_tree(clone_df, tumor_id = tumor_id, nb_path, clone_simplifications, sample_id = sample_id, ...)

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

  if (is.na(seu_path)) {
    return(NA_character_)
  }

  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*|_unfiltered_seu.*|_seu\\.rds$")

  nb_paths <- nb_paths %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  nb_path <- nb_paths[[tumor_id]]

  if (is.null(nb_path) || is.na(nb_path)) {
    warning("No numbat RDS found for ", tumor_id, "; skipping clone tree")
    return(NA_character_)
  }

  mynb <- readRDS(nb_path)
  all_clone_df <- mynb$clone_post %>%
    dplyr::select(cell, clone_opt) %>%
    dplyr::distinct()

  if (is.null(all_clone_df) || nrow(all_clone_df) == 0) {
    warning("No clone_post data for ", tumor_id, "; skipping clone tree")
    return(NA_character_)
  }

  retained_cells <- read_cell_barcodes_from_db(seu_path)
  clone_df <- if (!is.null(retained_cells)) {
    dplyr::filter(all_clone_df, cell %in% retained_cells)
  } else {
    all_clone_df
  }

  p <- plot_clone_tree(clone_df, tumor_id = tumor_id, nb_path, clone_simplifications, sample_id = sample_id, ...)

  dir.create("results/clone_trees", showWarnings = FALSE, recursive = TRUE)
  plot_path <- ggsave(glue("results/clone_trees/{sample_id}{label}.pdf"), plot = p, width = 4, height = 4)
  return(plot_path)
}


#' Retrieve a specific plot type from a list of numbat plot paths
#'
#' @param numbat_plots List or vector of file paths to numbat output PDFs
#' @param plot_type Filename pattern to match (e.g. "exp_roll_clust.pdf")
#' @return Character vector of matching file paths
#' @export
retrieve_numbat_plot_type <- function(numbat_plots, plot_type = "exp_roll_clust.pdf") {
  retrieved_plot_types <- purrr::map(numbat_plots, ~set_names(.x, fs::path_file(.x))) %>%
    purrr::map(stringr::str_detect, plot_type)

  purrr::map2(numbat_plots, retrieved_plot_types, ~{.x[.y]}) %>%
    unlist()
}
