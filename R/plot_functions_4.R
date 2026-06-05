# Numbat heatmap and plotting functions (5)

#' Perform numbat-related analysis
#'
#' @param seu_path File path
#' @param numbat_rds_files File path
#' @param p_min Parameter for p min
#' @param line_width Parameter for line width
#' @param extension Character string (default: "")
#' @param midline_threshold Threshold value for filtering
#' @param show_segment_names_on_x Logical; if TRUE, keep segment labels visible on the numbat heatmap x-axis.
#' @return ggplot2 plot object
#' @export
make_numbat_heatmaps <- function(seu_path, numbat_rds_files, p_min = 0.9, line_width = 0.1, extension = "", midline_threshold = 0.4, show_segment_names_on_x = FALSE, numbat_rds_filtered_files = NULL, filter_midline = TRUE) {
  sample_id <- str_extract(seu_path, "SR[RX][0-9]+")
  names(numbat_rds_files) <- str_extract(numbat_rds_files, "SR[RX][0-9]+")
  if (!is.null(numbat_rds_filtered_files) && length(numbat_rds_filtered_files) > 0) {
    names(numbat_rds_filtered_files) <- str_extract(numbat_rds_filtered_files, "SR[RX][0-9]+")
    filt_idx <- which(names(numbat_rds_filtered_files) == sample_id)
    if (length(filt_idx) > 0) {
      numbat_rds_files[[sample_id]] <- numbat_rds_filtered_files[[filt_idx[[1]]]]
    }
  }
  match_idx <- which(names(numbat_rds_files) == sample_id)
  if (length(match_idx) == 0) {
    warning("No numbat RDS file found for sample: ", sample_id)
    return(c(NA_character_, NA_character_))
  }
  numbat_rds_file <- numbat_rds_files[[match_idx[[1]]]]
  numbat_dir <- basename(dirname(dirname(numbat_rds_file)))
  dir_create(glue("results/{numbat_dir}/"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))
  mynb <- readRDS(numbat_rds_file)
  if (filter_midline) {
    retained_segs <- mynb$joint_post |>
      dplyr::mutate(at_midline = dplyr::case_when(
        dplyr::between(p_cnv, 0.3, 0.7) ~ 1,
        .default = 0
      )) |>
      group_by(seg) |>
      dplyr::summarise(percent_at_midline = sum(at_midline)/dplyr::n()) |>
      dplyr::filter(percent_at_midline <= midline_threshold) |>
      dplyr::arrange(desc(percent_at_midline)) |>
      dplyr::pull(seg) |>
      identity()
    mynb$joint_post <- mynb$joint_post[mynb$joint_post$seg %in% retained_segs,]
  }
  retained_cells <- unique(mynb$joint_post$cell)
  clone_annot <- mynb$clone_post %>%
    dplyr::select(cell, clone_opt) %>%
    dplyr::distinct() %>%
    dplyr::filter(cell %in% retained_cells)
  plot_result <- safe_plot_numbat(
    mynb,
    NULL,
    clone_annot,
    sample_id,
    clone_bar = FALSE,
    p_min = p_min,
    line_width = line_width,
    show_segment_names_on_x = show_segment_names_on_x
  )
  if (!is.null(plot_result$error)) {
    cat("WARNING: plot_numbat failed:", conditionMessage(plot_result$error), "\n")
  }
  numbat_heatmap <- plot_result[["result"]]
  if (!is.null(numbat_heatmap) && !identical(numbat_heatmap, NA_real_)) {
    heatmap_no_phylo_path <- tempfile(fileext = ".pdf")
    ggsave(heatmap_no_phylo_path, plot = numbat_heatmap, w = 10, h = 5)
    if (length(numbat_heatmap) >= 3) {
      scna_variability_plot <-
        numbat_heatmap[[3]][["data"]] |>
        dplyr::left_join(clone_annot, by = "cell") |>
        plot_variability_at_SCNA(p_min = p_min)
      scna_var_path <- tempfile(fileext = ".pdf")
      ggsave(scna_var_path, plot = scna_variability_plot, w = 10, h = 5)
    } else {
      scna_var_path <- NULL
    }
  } else {
    heatmap_no_phylo_path <- NULL
    scna_var_path <- NULL
  }
  numbat_heatmap_w_phylo <- safe_plot_numbat_w_phylo(
    mynb,
    NULL,
    clone_annot,
    sample_id,
    clone_bar = FALSE,
    p_min = p_min,
    show_segment_names_on_x = show_segment_names_on_x
  )[["result"]]
  heatmap_w_phylo_path <- tempfile(fileext = ".pdf")
  if (!is.null(numbat_heatmap_w_phylo) && !identical(numbat_heatmap_w_phylo, NA_real_)) {
    tryCatch(
      ggsave(heatmap_w_phylo_path, plot = numbat_heatmap_w_phylo, w = 10, h = 5),
      error = function(e) {
        warning("ggsave failed for w_phylo heatmap (", sample_id, "): ", conditionMessage(e))
        heatmap_w_phylo_path <<- NULL
      }
    )
  } else {
    heatmap_w_phylo_path <- NULL
  }
  plot_path <- glue("results/{numbat_dir}/{sample_id}/{sample_id}{extension}.pdf")
  pdf_inputs <- unlist(purrr::compact(list(heatmap_no_phylo_path, heatmap_w_phylo_path)))
  if (length(pdf_inputs) > 0) {
    qpdf::pdf_combine(pdf_inputs, plot_path)
  } else {
    plot_path <- NA_character_
  }
  scna_var_path_final <- if (!is.null(scna_var_path)) {
    scna_var_final <- glue("results/{numbat_dir}/{sample_id}/{sample_id}{extension}_scna_var.pdf")
    file.copy(scna_var_path, scna_var_final, overwrite = TRUE)
    scna_var_final
  } else {
    NA_character_
  }
  return(c(plot_path, scna_var_path_final))
}
