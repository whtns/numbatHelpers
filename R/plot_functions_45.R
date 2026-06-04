# Plot Functions (144)

#' Perform tidy eval arrange operation
#'
#' @param .data Parameter for .data
#' @param ... Additional arguments passed to other functions
#' @return Function result
#' @export
tidy_eval_arrange <- function(.data, ...) {
				.data %>%
					arrange(...)
			}

#' Perform run speckle operation
#'
#' @param seu_path File path
#' @param group.by Character string (default: "SCT_snn_res.0.6")
#' @return ggplot2 plot object
#' @export
run_speckle <- function(seu_path = "output/seurat/SRX11133594_filtered_seu.rds", group.by = "SCT_snn_res.0.6") {
  #
  seu <- readRDS(seu_path)
  sample_id <- str_extract(seu_path, "SR[RX][0-9]+")

  if (!"sample_id" %in% colnames(seu@meta.data)) {
    seu <- AddMetaData(seu, sample(2, ncol(seu), replace = TRUE), "sample_id")
  }

  seu_meta <-
    seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(clusters = .data[[group.by]], group = scna, sample_id, `orig.ident`) %>%
    # dplyr::mutate(group = as.numeric(factor(group))) %>%
    dplyr::mutate(sample_group = paste0(sample_id, "_", group)) %>%
    identity()

  # Run propeller testing for cell type proportion differences between the two
  # groups
  mytable <- propeller(
    clusters = seu_meta$clusters, sample = seu_meta$sample_group,
    group = seu_meta$group
  ) %>%
    tibble::rownames_to_column("cluster")

  seu_meta$clusters <- factor(seu_meta$clusters, levels = mytable$cluster)

  cell_prop_plot <- plotCellTypeProps(sample = seu_meta$clusters, clusters = seu_meta$group) +
    labs(title = sample_id)

  return(list("table" = mytable, "plot" = cell_prop_plot))
}

#' Perform run speckle for set operation
#'
#' @param seu_paths File path
#' @param group.by Character string (default: "SCT_snn_res.0.6")
#' @return Function result
#' @export
run_speckle_for_set <- function(seu_paths, group.by = "SCT_snn_res.0.6") {
  #
  dir_create("results/speckle")
  test0 <- map(seu_paths, safe_run_speckle, group.by = group.by)

  plotlist <- compact(map(test0, c("result", "plot")))

  pdf(glue("results/speckle/{group.by}.pdf"))
  print(plotlist)
  dev.off()

  tables <- compact(map(test0, c("result", "table")))

  write_xlsx(tables, glue("results/speckle/{group.by}.xlsx"))
}

#' Create a plot visualization
#'
#' @param mylabel Parameter for mylabel
#' @param clustree_plot Parameter for clustree plot
#' @param speckle_proportions Parameter for speckle proportions
#' @param sample_id Parameter for sample id
#' @return ggplot2 plot object
#' @export
plot_clustree_per_comparison <- function(mylabel, clustree_plot, speckle_proportions, sample_id) {
  #

  comparison_clones <- str_split(mylabel, pattern = "_v_") %>%
    unlist()

  comparison_clones[comparison_clones == "x"] <- "diploid"

  clone_names <- colnames(speckle_proportions)[!colnames(speckle_proportions) %in% c("samples", "clusters")]

  brewer_palettes <- c("Reds", "Greens", "Blues", "Purples", "Oranges")

  clone_colors <- brewer_palettes[1:length(clone_names)] %>%
    set_names(clone_names)

  # clone_colors <- clone_colors[names(clone_colors) %in% comparison_clones]

  clustree_res <- imap(clone_colors, ~ color_clustree_by_clone(clustree_plot, .x, .y, mylabel = mylabel, sample_id))

  return(
    list(
      "plot" = map(clustree_res, "plot"),
      "table" = map(clustree_res, "table")
    )
  )
}

#' Perform color clustree by clone operation
#'
#' @param clustree_plot Parameter for clustree plot
#' @param mycolor Color specification
#' @param myclone Parameter for myclone
#' @param mylabel Character string (default: "asdf")
#' @param sample_id Parameter for sample id
#' @return ggplot2 plot object
#' @export
color_clustree_by_clone <- function(clustree_plot, mycolor, myclone, mylabel = "asdf", sample_id) {
  #

  label_data <-
    clustree_plot$data %>%
    dplyr::filter(!!sym(myclone) < 0.05) %>%
    # dplyr::filter(x12p_16q_1q_2p_v_x12p_16q_1q_2p_11p_8p < 0.05) %>%
    identity()

  myplot <- clustree_plot +
    ggraph::geom_node_point(aes(colour = .data[[myclone]], size = size)) +
    ggraph::geom_node_text(aes(label = cluster)) +
    ggraph::geom_node_label(data = label_data, aes(label = cluster)) +
    labs(title = mylabel, subtitle = myclone, caption = sample_id, colour = "clone %") +
    scale_color_distiller(palette = mycolor, direction = 1) +
    NULL

  return(list("plot" = myplot, "table" = clustree_plot$data))
}

