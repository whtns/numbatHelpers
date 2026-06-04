# Plot Functions (140)

#' Perform ora effect of regression operation
#'
#' @param filtered_seu_path File path
#' @param regressed_seu_path File path
#' @param resolution Character string (default: "0.4")
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O

ora_effect_of_regression <- function(filtered_seu_path, regressed_seu_path, resolution = "0.4") {
  
  
  #

  plot_list <- list()

  sample_id <- str_extract(filtered_seu_path, "SR[RX][0-9]+")

  filtered_seu <- readRDS(filtered_seu_path)

  regressed_seu <- readRDS(regressed_seu_path)

  fs::dir_create("results/effect_of_regression")

  # cluster ora_analysis ------------------------------
  fs::dir_create(glue("results/effect_of_regression/enrichment"))
  fs::dir_create(glue("results/effect_of_regression/enrichment/filtered"))

  filtered_ora_output <- ora_analysis(filtered_seu, "SCT_snn_res.0.4")

  filtered_ora_tables <- filtered_ora_output[["tables"]]

  table_path <- glue("results/effect_of_regression/enrichment/filtered/{sample_id}_filtered_cluster_enrichment.xlsx")
  plot_list["filtered_enrichment_tables"] <- table_path
  writexl::write_xlsx(filtered_ora_tables, table_path)

  filtered_ora_plots <- filtered_ora_output[["plots"]] %>%
    wrap_plots() +
    plot_annotation(title = sample_id)

  plot_path <- glue("results/effect_of_regression/enrichment/filtered/{sample_id}_filtered_cluster_enrichment.pdf")
  plot_list["filtered_enrichment"] <- plot_path
  ggsave(plot_path, plot = filtered_ora_plots, height = 32, width = 40)

  fs::dir_create(glue("results/effect_of_regression/enrichment"))
  fs::dir_create(glue("results/effect_of_regression/enrichment/regressed"))

  regressed_ora_output <- ora_analysis(regressed_seu, "SCT_snn_res.0.4")
  regressed_ora_tables <- regressed_ora_output[["tables"]]

  table_path <- glue("results/effect_of_regression/enrichment/regressed/{sample_id}_regressed_cluster_enrichment.xlsx")
  plot_list["regressed_enrichment_tables"] <- table_path
  writexl::write_xlsx(regressed_ora_tables, table_path)

  regressed_ora_plots <- regressed_ora_output[["plots"]] %>%
    wrap_plots() +
    plot_annotation(title = sample_id)

  plot_path <- glue("results/effect_of_regression/enrichment/regressed/{sample_id}_regressed_cluster_enrichment.pdf")
  plot_list["regressed_enrichment"] <- plot_path
  ggsave(plot_path, plot = regressed_ora_plots, height = 32, width = 40)

  # diffex ora_analysis ------------------------------
  fs::dir_create(glue("results/effect_of_regression/enrichment"))
  fs::dir_create(glue("results/effect_of_regression/enrichment/filtered"))

  filtered_ora_output <- ora_analysis(filtered_seu, "SCT_snn_res.0.4")

  filtered_ora_tables <- filtered_ora_output[["tables"]]

  table_path <- glue("results/effect_of_regression/enrichment/filtered/{sample_id}_filtered_diffex_enrichment.xlsx")
  plot_list["filtered_enrichment_tables"] <- table_path
  writexl::write_xlsx(filtered_ora_tables, table_path)

  filtered_ora_plots <- filtered_ora_output[["plots"]] %>%
    wrap_plots() +
    plot_annotation(title = sample_id)

  plot_path <- glue("results/effect_of_regression/enrichment/filtered/{sample_id}_filtered_diffex_enrichment.pdf")
  plot_list["filtered_enrichment"] <- plot_path
  ggsave(plot_path, plot = filtered_ora_plots, height = 32, width = 40)

  fs::dir_create(glue("results/effect_of_regression/enrichment"))
  fs::dir_create(glue("results/effect_of_regression/enrichment/regressed"))

  regressed_ora_output <- ora_analysis(regressed_seu, "SCT_snn_res.0.4")
  regressed_ora_tables <- regressed_ora_output[["tables"]]

  table_path <- glue("results/effect_of_regression/enrichment/regressed/{sample_id}_regressed_diffex_enrichment.xlsx")
  plot_list["regressed_enrichment_tables"] <- table_path
  writexl::write_xlsx(regressed_ora_tables, table_path)

  regressed_ora_plots <- regressed_ora_output[["plots"]] %>%
    wrap_plots() +
    plot_annotation(title = sample_id)

  plot_path <- glue("results/effect_of_regression/enrichment/regressed/{sample_id}_regressed_diffex_enrichment.pdf")
  plot_list["regressed_enrichment"] <- plot_path
  ggsave(plot_path, plot = regressed_ora_plots, height = 32, width = 40)


  return(plot_list)
}

#' Perform split label line operation
#'
#' @param label Parameter for label
#' @param n_comma_values Parameter for n comma values
#' @return Function result
#' @export
split_label_line <- function(label, n_comma_values = 2) {
  if (is.na(label)) return(NA_character_)

  label_vec <- label %>% stringr::str_split_1(",")

  label_groups <- ceiling(seq_along(label_vec) / n_comma_values)

  split_label <- split(label_vec, label_groups) %>%
    purrr::map(paste, collapse = ", ") %>%
    paste(collapse = "\n")
}

vec_split_label_line <- Vectorize(split_label_line)

#' Calculate scores for the given data
#'
#' @param unfiltered_seu_path File path
#' @param filtered_seu_path File path
#' @param celltype_markers Cell identifiers or information
#' @return ggplot2 plot object
#' @export
score_samples_for_celltype_enrichment <- function(unfiltered_seu_path, filtered_seu_path, celltype_markers) {
  
  
  #
  sample_id <- str_extract(unfiltered_seu_path, "SR[RX][0-9]+")

  unfiltered_seu <- readRDS(unfiltered_seu_path) %>%
    Seurat::AddModuleScore(features = celltype_markers, name = "celltype")

  filtered_seu <- filtered_seu_path %>%
    Seurat::AddModuleScore(features = celltype_markers, name = "celltype")

  module_names <- paste0("celltype", seq(length(celltype_markers)))

  unfiltered_seu@meta.data[names(celltype_markers)] <- unfiltered_seu@meta.data[module_names]
  unfiltered_seu@meta.data[module_names] <- NULL

  filtered_seu@meta.data[names(celltype_markers)] <- filtered_seu@meta.data[module_names]
  filtered_seu@meta.data[module_names] <- NULL

  # unfiltered featureplot ------------------------------
  (FeaturePlot(unfiltered_seu, names(celltype_markers)) +
    plot_annotation(subtitle = "unfiltered"))

  file_tag <- str_extract(unfiltered_seu_path, "SR[RX][0-9]+_[a-z]*")
  dir_create("results/celltype_plots")
  unfiltered_celltype_plot_path <- glue("results/celltype_plots/{file_tag}_featureplot.pdf")
  ggsave(unfiltered_celltype_plot_path, height = 8, width = 8)

  # filtered featureplot ------------------------------
  (FeaturePlot(filtered_seu, names(celltype_markers)) +
    plot_annotation(subtitle = "filtered"))

  file_tag <- str_extract(filtered_seu_path, "SR[RX][0-9]+_[a-z]*")
  dir_create("results/celltype_plots")
  filtered_celltype_plot_path <- glue("results/celltype_plots/{file_tag}_featureplot.pdf")
  ggsave(filtered_celltype_plot_path, height = 8, width = 8)

  unfiltered_seu <- score_binary_celltype_markers(unfiltered_seu, celltype_markers)
  filtered_seu <- score_binary_celltype_markers(filtered_seu, celltype_markers)

  # unfiltered dimplot ------------------------------
  (DimPlot(unfiltered_seu, group.by = paste0(names(celltype_markers), "_id")) +
    plot_annotation(subtitle = "unfiltered"))

  file_tag <- str_extract(unfiltered_seu_path, "SR[RX][0-9]+_[a-z]*")
  dir_create("results/celltype_plots")
  unfiltered_celltype_plot_path <- glue("results/celltype_plots/{file_tag}_dimplot.pdf")
  ggsave(unfiltered_celltype_plot_path, height = 8, width = 8)

  # filtered dimplot ------------------------------
  (DimPlot(filtered_seu, group.by = paste0(names(celltype_markers), "_id")) +
    plot_annotation(subtitle = "filtered"))

  file_tag <- str_extract(filtered_seu_path, "SR[RX][0-9]+_[a-z]*")
  dir_create("results/celltype_plots")
  filtered_celltype_plot_path <- glue("results/celltype_plots/{file_tag}_dimplot.pdf")
  ggsave(filtered_celltype_plot_path, height = 8, width = 8)

  return(list(filtered_celltype_plot_path, unfiltered_celltype_plot_path))
}

#' Calculate scores for the given data
#'
#' @param seu Seurat object
#' @param celltype_markers Cell identifiers or information
#' @return Modified Seurat object
#' @export
score_binary_celltype_markers <- function(seu, celltype_markers) {
  
  
  #
  upper_quartiles <-
    names(celltype_markers) %>%
    set_names(.) %>%
    map(~ quantile(seu@meta.data[[.x]], 0.75))

  for (celltype in names(celltype_markers)) {
    seu@meta.data[[paste0(celltype, "_id")]] <- seu@meta.data[[celltype]] > upper_quartiles[[celltype]]
  }

  return(seu)
}

#' Calculate scores for the given data
#'
#' @param seu Seurat object
#' @param celltype_markers Cell identifiers or information
#' @return Modified Seurat object
#' @export
score_binary_celltype_clusters <- function(seu, celltype_markers) {
  
  
  upper_quartiles <-
    names(celltype_markers) %>%
    set_names(.) %>%
    map(~ quantile(seu@meta.data[[.x]], 0.95))

  seu <- seurat_cluster(seu, resolution = 2.0)

  
#' Perform clustering analysis
#'
#' @param upper_quartile Parameter for upper quartile
#' @param cluster_name Cluster information
#' @param seu Seurat object
#' @return Function result
#' @export
pick_max_cluster <- function(upper_quartile, cluster_name, seu) {
  
  
    test0 <-
      seu@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::group_by(`gene_snn_res.2`) %>%
      dplyr::filter(dplyr::n() < 30) %>%
      dplyr::summarize(score = mean(.data[[cluster_name]])) %>%
      # dplyr::filter(score > upper_quartile) %>%
      # dplyr::slice_max(score) %>%
      # dplyr::pull(`gene_snn_res.2`) %>%
      identity()
  }

  test0 <- imap(upper_quartiles, pick_max_cluster, seu)

  for (celltype in names(celltype_markers)) {
    seu@meta.data[[paste0(celltype, "_id")]] <- seu@meta.data[["gene_snn_res.2"]] %in% test0[[celltype]]
  }

  return(seu)
}


simplify_gt <- function(mynb, rb_scnas = c("1" = "1q", "2" = "2p", "6" = "6p", "8" = "8p", "11" = "11p", "13" = "13q", "16" = "16q")) {
  mynb$mut_graph

  if (all(diff(igraph::vertex_attr(mynb$mut_graph)$clone) <= 0)) {
    mynb$mut_graph <- reverse_edges(mynb$mut_graph)
  }

  to_labels <- igraph::edge_attr(mynb$mut_graph, "to_label") %>%
    str_split(",") %>%
    map(~ (paste(names(rb_scnas[rb_scnas %in% .x]), collapse = ","))) %>%
    map_chr(str_pad, side = "left", width = 2) %>%
    identity()

  from_labels <- igraph::edge_attr(mynb$mut_graph, "from_label") %>%
    str_split(",") %>%
    map(~ (paste(names(rb_scnas[rb_scnas %in% .x]), collapse = ","))) %>%
    map_chr(str_pad, side = "left", width = 2) %>%
    identity()

  mynb$mut_graph <- mynb$mut_graph %>%
    igraph::set_edge_attr("to_label", value = to_labels) %>%
    igraph::set_edge_attr("from_label", value = from_labels) %>%
    identity()

  return(mynb)
}

seu_integrate_rbs <- function(numbat_dir = "output/numbat_sridhar", kept_samples = c("SRX11133594", "SRX11133593", "SRX11133592"), cluster_dictionary = cluster_dictionary) {
  numbat_rds_files <- retrieve_numbat_rds_files(numbat_dir, kept_samples)

  seus <- map(numbat_rds_files, retrieve_numbat_seurat, cluster_dictionary)

  integrated_seu <- seuratTools::seurat_integration_pipeline(seus, resolution = c(0.2, 0.4))

  sample_slug <- paste(kept_samples, collapse = "_")

  seu_path <- glue("output/seurat/{sample_slug}_seu.rds")

  add_batch_hash_metadata(seu = integrated_seu, filepath = seu_path)

  return(seu_path)
}


score_pca <- function(seu) {
  mygroups <- seu$clone_opt %>%
    tibble::enframe("cell", "group")

  mypca <- seu@reductions$pca@cell.embeddings %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(-c("cell"), names_to = "PC", values_to = "value") %>%
    dplyr::left_join(mygroups, by = "cell")

  test0 <-
    mypca %>%
    dplyr::group_by(PC, group) %>%
    dplyr::summarize(value = mean(value)) %>%
    tidyr::pivot_wider(names_from = "PC", values_from = "value") %>%
    tibble::column_to_rownames("group") %>%
    as.matrix() %>%
    t() %>%
    cor() %>%
    # dist() %>%
    identity()

  test1 <- dist(1 - test0)

  return(test1)

  # # aggregate values within categories using 'mean'
  # mean_df = rep_df.groupby(level=0).mean()
  #
  # import scipy.cluster.hierarchy as sch
  # from scipy.spatial import distance
  #
  # corr_matrix = mean_df.T.corr(method=cor_method)
  # corr_condensed = distance.squareform(1 - corr_matrix)
  # z_var = sch.linkage(
  #   corr_condensed, method=linkage_method, optimal_ordering=optimal_ordering
  # )
  # dendro_info = sch.dendrogram(z_var, labels=list(categories), no_plot=True)
}
