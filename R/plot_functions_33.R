# Plot Functions (132)
#' Generate a volcano plot for differential expression data
#'
#' @param myplot Plot object (ggplot2)
#' @param out_html Parameter for out html
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

convert_volcano_to_plotly <- function(myplot, out_html) {
  myplotly <- ggplotly(myplot + aes(x = avg_log2FC, y = -log10(p_val_adj), symbol = symbol), tooltip = "symbol")

  saveWidget(myplotly, out_html, selfcontained = T, libdir = "lib")
}

#' Perform tabulate clone comparisons operation
#'
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @return Function result
#' @export
tabulate_clone_comparisons <- function(large_clone_comparisons) {
  clone_comparison_table <-
    large_clone_comparisons %>%
    map(~ {
      .x %>%
        tibble::enframe("comparison", "segment") %>%
        tidyr::unnest(segment)
    }) %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    identity()

  return(clone_comparison_table)
}

#' Perform remove non tumor cells operation
#'
#' @param seu Seurat object
#' @return Function result
#' @export
remove_non_tumor_cells <- function(seu) {
  # identify retinal cell type clusters
  # identify cells with no SCNAs

  # filter out cells that are
  # 1) not cones or RB cells
  # 2) have no SCNAs
}

#' Filter data based on specified criteria
#'
#' @param df Input data frame or dataset
#' @param scna_of_interest Parameter for scna of interest
#' @param p_val_threshold Threshold value for filtering
#' @param fc_threshold Threshold value for filtering
#' @param n_recur Parameter for n recur
#' @return Filtered data
#' @export
filter_input_by_scna <- function(df, scna_of_interest, p_val_threshold, fc_threshold, n_recur) {
  #
  required_cols <- c("symbol", "p_val_adj", "avg_log2FC", "sample_id")
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || !all(required_cols %in% colnames(df))) {
    return(df)
  }

  filtered_df <-
    df %>%
    dplyr::group_by(symbol) %>%
    dplyr::filter(p_val_adj <= p_val_threshold) %>%
    dplyr::filter(abs(avg_log2FC) >= fc_threshold) %>%
    dplyr::filter(n_distinct(sample_id) >= n_recur) %>%
    identity()
}

