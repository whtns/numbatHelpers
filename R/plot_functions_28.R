read_cells_to_remove <- function(cell_remove_file = "data/cells_to_remove_final.xlsx") {
  mysheets <- excel_sheets(cell_remove_file) %>%
    set_names(.)

  cells_to_remove <-
    mysheets %>%
    map(~ readxl::read_xlsx(cell_remove_file, .x))

  return(cells_to_remove)
}

# Plot Functions (127)
#' Perform make pdf montages operation
#'
#' @param plot_files File path
#' @param heatmaps Parameter for heatmaps
#' @param tile Character string (default: "6")
#' @return Function result
#' @export
make_pdf_montages <- function(plot_files, heatmaps, tile = "6") {
  numbat_dir <-
    path_split(plot_files[[1]][[1]])[[1]][[2]] %>%
    identity()

  plot_files <- map2(plot_files, heatmaps, c)

  # plot_files <- map2(plot_files, expression_files, c)

  sample_ids <-
    plot_files %>%
    map(1) %>%
    map(fs::path_split) %>%
    map(c(1, 3))

  names(plot_files) <- sample_ids

  montage_paths <- imap(plot_files, montage_images, tile = tile)

  return(montage_paths)
}

#' Perform make expression heatmap comparison operation
#'
#' @param large_numbat_pdfs Parameter for large numbat pdfs
#' @param heatmaps Parameter for heatmaps
#' @param montage_path File path
#' @return Function result
#' @export
make_expression_heatmap_comparison <- function(large_numbat_pdfs, heatmaps, montage_path = "results/expr_heatmap_comparison.pdf") {
  numbat_dir <-
    path_split(large_numbat_pdfs[[1]][[1]])[[1]][[2]] %>%
    identity()

  expression_files <- map(large_numbat_pdfs, 6)
  names(expression_files) <- str_extract(expression_files, "SR[RX][0-9]+")

  heatmaps <- map(heatmaps, 1)
  names(heatmaps) <- str_extract(heatmaps, "SR[RX][0-9]+")
  
  heatmaps <- heatmaps[names(heatmaps) %in% names(expression_files)]
  
  expression_files <- expression_files[names(expression_files) %in% names(heatmaps)]

  plot_files <- map2(expression_files, heatmaps, c)

  sample_ids <-
    plot_files %>%
    map(1) %>%
    map(fs::path_split) %>%
    map(c(1, 3))

  names(plot_files) <- names(sample_ids)
  montage_paths <- imap(plot_files, montage_images)

  qpdf::pdf_combine(montage_paths, montage_path)	
  
  return(montage_path)
}