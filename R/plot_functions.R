# Numbat plotting and analysis helper functions

#' Plot numbat heatmap
#'
#' @param nb Numbat object
#' @param myseu Seurat object
#' @param myannot Annotation data frame
#' @param mytitle Plot title
#' @param sort_by Sorting variable
#' @param ... Additional arguments
#'
#' @return ggplot object
#' @export
plot_numbat <- function(nb, myseu, myannot, mytitle, sort_by = "scna", ...) {
  num_cols <- length(unique(nb$clone_post$clone_opt))
  mypal <- scales::hue_pal()(num_cols) %>%
    set_names(seq(num_cols))

  myheatmap <- nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = myannot,
    show_phylo = FALSE,
    sort_by = sort_by,
    annot_bar_width = 1,
    raster = FALSE,
    ...
  ) + ggplot2::labs(title = mytitle)

  return(myheatmap)
}

#' Plot numbat heatmap with phylogeny
#'
#' @param nb Numbat object
#' @param myseu Seurat object
#' @param myannot Annotation data frame
#' @param mytitle Plot title
#' @param ... Additional arguments
#'
#' @return ggplot object
#' @export
plot_numbat_w_phylo <- function(nb, myseu, myannot, mytitle, ...) {
  celltypes <- myseu@meta.data["type"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(cell = stringr::str_replace(cell, "\\.", "-")) %>%
    identity()

  myannot <- dplyr::left_join(myannot, celltypes, by = "cell")
  mypal <- c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

  nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = myannot,
    show_phylo = TRUE,
    annot_bar_width = 1,
    ...
  ) + ggplot2::labs(title = mytitle)
}

safe_plot_numbat <- purrr::safely(plot_numbat, otherwise = NA_real_)
safe_plot_numbat_w_phylo <- purrr::safely(plot_numbat_w_phylo, otherwise = NA_real_)

#' Plot variability at SCNA
#'
#' @param phylo_plot_output Data from phylo plot
#' @param chrom Chromosome
#' @param p_min Minimum p value
#'
#' @return ggplot object
#' @export
plot_variability_at_SCNA <- function(phylo_plot_output, chrom = "1", p_min = 0.9){
  test0 <- phylo_plot_output %>%
    dplyr::mutate(seg = factor(seg, levels = stringr::str_sort(unique(seg), numeric = TRUE))) %>%
    identity()

  p_cnv_plot <- ggplot2::ggplot(test0, ggplot2::aes(x = cell_index, y = p_cnv, color = cnv_state)) +
    ggplot2::geom_point(size = 0.1, alpha = 0.1) +
    ggplot2::scale_x_reverse() +
    ggplot2::scale_color_manual(values = c("amp" = "#7f180f",
                                  "bamp" = "pink",
                                  "del" = "#010185",
                                  "loh" = "#387229")) +
    ggplot2::facet_wrap(~seg)

  return(p_cnv_plot)
}

#' Make numbat heatmaps from Seurat and numbat objects
#'
#' @param seu Seurat object
#' @param numbat_rds_files Numbat RDS file paths
#' @param p_min Minimum p value
#' @param line_width Line width for plots
#' @param extension File extension
#' @param midline_threshold Threshold for midline filtering
#' @param show_segment_names_on_x Show segment names on x-axis
#' @param numbat_rds_filtered_files Optional filtered numbat files
#'
#' @return Character vector with plot paths
#' @export
make_numbat_heatmaps <- function(seu, numbat_rds_files, p_min = 0.9, line_width = 0.1, extension = "", midline_threshold = 0.4, show_segment_names_on_x = FALSE, numbat_rds_filtered_files = NULL, filter_midline = TRUE) {
  # seu may be a Seurat object, a file path, or a list containing a path
  if (!inherits(seu, "Seurat")) {
    seu_path <- unlist(seu, use.names = FALSE)[1]
    sample_id <- stringr::str_extract(seu_path, "SRX[0-9]+")
    seu <- readRDS(seu_path)
  } else {
    sample_id <- stringr::str_extract(colnames(seu)[1], "SRX[0-9]+")
  }
  names(numbat_rds_files) <- stringr::str_extract(numbat_rds_files, "SRX[0-9]+")
  if (!is.null(numbat_rds_filtered_files) && length(numbat_rds_filtered_files) > 0) {
    names(numbat_rds_filtered_files) <- stringr::str_extract(numbat_rds_filtered_files, "SRX[0-9]+")
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
  dir.create(glue::glue("results/{sample_id}"), showWarnings = FALSE, recursive = TRUE)

  seu <- Seurat::RenameCells(seu, new.names = stringr::str_replace(colnames(seu), "\\.", "-"))
  mynb <- readRDS(numbat_rds_file)
  if (filter_midline) {
    retained_segs <- mynb$joint_post |>
      dplyr::mutate(at_midline = dplyr::case_when(
        dplyr::between(p_cnv, 0.3, 0.7) ~ 1,
        .default = 0
      )) |>
      dplyr::group_by(seg) |>
      dplyr::summarise(percent_at_midline = sum(at_midline)/dplyr::n()) |>
      dplyr::filter(percent_at_midline <= midline_threshold) |>
      dplyr::arrange(dplyr::desc(percent_at_midline)) |>
      dplyr::pull(seg) |>
      identity()
    mynb$joint_post <- mynb$joint_post[mynb$joint_post$seg %in% retained_segs,]
  }
  cell_names <- mynb$joint_post$cell
  seu <- seu[, colnames(seu) %in% cell_names]
  myannot <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, scna) %>%
    identity()
  myannot$scna[myannot$scna == ""] <- ".diploid"
  clone_annot <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, clone_opt) %>%
    identity()
  numbat_heatmap <- safe_plot_numbat(
    mynb,
    seu,
    clone_annot,
    sample_id,
    clone_bar = FALSE,
    p_min = p_min,
    line_width = line_width,
    show_segment_names_on_x = show_segment_names_on_x
  )[["result"]]
  if (!is.null(numbat_heatmap) && !identical(numbat_heatmap, NA_real_)) {
    heatmap_no_phylo_path <- tempfile(fileext = ".pdf")
    ggplot2::ggsave(heatmap_no_phylo_path, plot = numbat_heatmap, w = 10, h = 5)
    if (length(numbat_heatmap) >= 3) {
      scna_variability_plot <- numbat_heatmap[[3]][["data"]] |>
        dplyr::left_join(clone_annot, by = "cell") |>
        plot_variability_at_SCNA(p_min = p_min)
      scna_var_path <- tempfile(fileext = ".pdf")
      ggplot2::ggsave(scna_var_path, plot = scna_variability_plot, w = 10, h = 5)
    } else {
      scna_var_path <- NULL
    }
  } else {
    heatmap_no_phylo_path <- NULL
    scna_var_path <- NULL
  }
  numbat_heatmap_w_phylo <- safe_plot_numbat_w_phylo(
    mynb,
    seu,
    clone_annot,
    sample_id,
    clone_bar = FALSE,
    p_min = p_min,
    show_segment_names_on_x = show_segment_names_on_x
  )[["result"]]
  heatmap_w_phylo_path <- tempfile(fileext = ".pdf")
  if (!is.null(numbat_heatmap_w_phylo) && !identical(numbat_heatmap_w_phylo, NA_real_)) {
    tryCatch(
      ggplot2::ggsave(heatmap_w_phylo_path, plot = numbat_heatmap_w_phylo, w = 10, h = 5),
      error = function(e) {
        warning("ggsave failed for w_phylo heatmap (", sample_id, "): ", conditionMessage(e))
        heatmap_w_phylo_path <<- NULL
      }
    )
  } else {
    heatmap_w_phylo_path <- NULL
  }
  plot_path <- glue::glue("results/{sample_id}/{sample_id}{extension}.pdf")
  pdf_inputs <- unlist(purrr::compact(list(heatmap_no_phylo_path, heatmap_w_phylo_path)))
  if (length(pdf_inputs) > 0) {
    qpdf::pdf_combine(pdf_inputs, plot_path)
  } else {
    plot_path <- NA_character_
  }
  scna_var_path_final <- if (!is.null(scna_var_path)) {
    scna_var_final <- glue::glue("results/{sample_id}/{sample_id}{extension}_scna_var.pdf")
    file.copy(scna_var_path, scna_var_final, overwrite = TRUE)
    scna_var_final
  } else {
    NA_character_
  }
  return(c(plot_path, scna_var_path_final))
}

#' Score chromosomal instability
#'
#' @param seu Seurat object
#'
#' @return Numeric vector of instability scores
#' @export
score_chrom_instability <- function(seu) {
  # Placeholder for chromosome instability scoring
  # Would need more information about implementation
  rep(0, ncol(seu))
}
