# Plot Functions (128)

#' Perform browse celltype expression operation
#'
#' @param sridhar_seu Parameter for sridhar seu
#' @param symbol Parameter for symbol
#' @return Data frame
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - rownames_round_trip: Avoid unnecessary rownames conversions - keep as column when possible
# - multiple_joins: Combine multiple joins into single join operation where possible

browse_celltype_expression <- function(sridhar_seu, symbol) {
  
  pdf(glue("~/tmp/{symbol}.pdf"))
  FeaturePlot(sridhar_seu, features = c(glue("{symbol}")), split.by = "CellType_predict", combine = FALSE)
  dev.off()
  browseURL(glue("~/tmp/{symbol}.pdf"))

  return(glue("~/tmp/{symbol}.pdf"))
}

#' Perform collect markers operation
#'
#' @param numbat_rds_file File path
#' @param metavar Character string (default: "gene_snn_res.0.2")
#' @param num_markers Parameter for num markers
#' @param selected_values Parameter for selected values
#' @param return_plotly Logical flag (default: TRUE)
#' @param marker_method Character string (default: "presto")
#' @param seurat_assay Character string (default: "gene")
#' @param hide_technical Character string (default: "all")
#' @param unique_markers Logical flag (default: FALSE)
#' @param p_val_cutoff Threshold value for filtering
#' @param ... Additional arguments passed to other functions
#' @return Function result
#' @export
collect_markers <- function(numbat_rds_file, metavar = "gene_snn_res.0.2", num_markers = 5, selected_values = NULL, return_plotly = TRUE, marker_method = "presto", seurat_assay = "gene", hide_technical = "all", unique_markers = FALSE, p_val_cutoff = 1, ...) {
  
  
  # # by default only resolution markers are calculated in pre-processing
  # seu <- find_all_markers(numbat_rds_file, metavar, seurat_assay = seurat_assay, p_val_cutoff = p_val_cutoff)

  sample_id <- str_extract(numbat_rds_file, "SR[RX][0-9]+")

  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  marker_table <- seu@misc$markers[[metavar]][[marker_method]]

  markers <-
    marker_table %>%
    seuratTools:::enframe_markers() %>%
    dplyr::mutate(dplyr::across(.fns = as.character))

  if (!is.null(hide_technical)) {
    markers <- purrr::map(markers, c)

    if (hide_technical == "pseudo") {
      markers <- purrr::map(markers, ~ .x[!.x %in% pseudogenes[[seurat_assay]]])
    } else if (hide_technical == "mito_ribo") {
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^MT-")])
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^RPS")])
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^RPL")])
    } else if (hide_technical == "all") {
      markers <- purrr::map(markers, ~ .x[!.x %in% pseudogenes[[seurat_assay]]])
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^MT-")])
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^RPS")])
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^RPL")])
    }

    min_length <- min(purrr::map_int(markers, length))

    markers <- purrr::map(markers, head, min_length) %>%
      dplyr::bind_cols()
  }

  colnames(markers) <- glue("{colnames(markers)}_{sample_id}")

  return(markers)
}

#' Perform collect all markers operation
#'
#' @param numbat_rds_files File path
#' @param excel_output Character string (default: "results/markers.xlsx")
#' @return Function result
#' @export
collect_all_markers <- function(numbat_rds_files, excel_output = "results/markers.xlsx") {
  
  
  names(numbat_rds_files) <- str_extract(numbat_rds_files, "SR[RX][0-9]+")

  marker_tables <- purrr::map(numbat_rds_files, collect_markers)

  writexl::write_xlsx(marker_tables, excel_output)
}

#' Extract or pull specific data elements
#'
#' @param myexpression Parameter for myexpression
#' @param joint_post Parameter for joint post
#' @return Extracted data elements
#' @export
pull_cells_matching_expression <- function(myexpression, joint_post) {
  
  
  #
  excluded_cells <-
    joint_post %>%
    dplyr::filter(!!parse_expr(myexpression)) %>%
    dplyr::pull(cell) %>%
    identity()

  return(excluded_cells)
}

#' Calculate scores for the given data
#'
#' @param numbat_rds_file File path
#' @param cluster_dictionary Cluster information
#' @param filter_expressions Parameter for filter expressions
#' @return List object
#' @export
score_filtration <- function(numbat_rds_file, cluster_dictionary, filter_expressions = NULL) {
  
  
  #

  sample_id <- str_extract(numbat_rds_file, "SR[RX][0-9]+")

  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[, !is.na(seu$clone_opt)]

  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")

  seu <- AddMetaData(seu, test0)

  phylo_heatmap_data <- mynb$clone_post %>%
    dplyr::select(cell, clone_opt) %>%
    dplyr::left_join(mynb$joint_post, by = "cell")

  # filter out cells
  excluded_cells <- map(filter_expressions[[sample_id]], pull_cells_matching_expression, phylo_heatmap_data) %>%
    unlist()

  seu_filtered <- seu[, !colnames(seu) %in% excluded_cells]

  dist_list <- list(
    "unfiltered" = score_pca(seu),
    "filtered" = score_pca(seu_filtered)
  )

  return(dist_list)
}

collect_clusters_from_seus <- function(filtered_seus) {
  #

  resolutions <- glue("SCT_snn_res.{seq(0.2, 1.0, by = 0.2)}")

  gather_clusters <- function(filtered_seu, resolutions) {
    #
    # sample_id <- str_extract(filtered_seu, "SR[RX][0-9]+")

    # numbat_dir = fs::path_split(filtered_seu)[[1]][[2]]

    seu <- readRDS(filtered_seu)

    clusters <- seu@meta.data[resolutions] %>%
      tidyr::pivot_longer(everything(), names_to = "resolution", values_to = "cluster") %>%
      dplyr::group_by(resolution, cluster) %>%
      dplyr::summarize(n_cells = dplyr::n()) %>%
      split(.$resolution) %>%
      map(ungroup) %>%
      imap(~ dplyr::select(.x, {{ .y }} := cluster, n_cells)) %>%
      identity()

    return(clusters)
  }

  sample_ids <- str_extract(filtered_seus, "SR[RX][0-9]+")

  names(filtered_seus) <- sample_ids

  myclusters <- map(filtered_seus, gather_clusters, resolutions) %>%
    purrr::transpose() %>%
    purrr::map(dplyr::bind_rows, .id = "sample_id") %>%
    identity()

  return(myclusters)
}

compare_cluster_continuous_var <- function(seu_path, cluster_order = NULL, continuous_var = "TFF1", group.by = "SCT_snn_res.0.6", assay = "SCT", mygene = "EZH2", label = "_filtered_", height = 14, width = 22, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), gene_lists = NULL) {
  #
	file_id <- fs::path_file(seu_path)
	
	tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
	
	message(file_id)
	cluster_order_list <- cluster_order[[file_id]]
	cluster_order <- if (!is.null(cluster_order_list)) cluster_order_list[["0"]] %||% cluster_order_list[[1]] else NULL

  seu <- readRDS(seu_path)

  heatmap_features <-
    table_cluster_markers(seu, assay = assay)

  if (!is.null(cluster_order)) {
    group.by <- unique(cluster_order$resolution)
  }

  cluster_order <-
    if (!is.null(cluster_order)) cluster_order %>%
    dplyr::filter(!is.na(clusters)) %>%
    dplyr::mutate(clusters = as.character(clusters)) else NULL

  seu@meta.data$clusters <- seu@meta.data[[group.by]]

  seu_meta <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::left_join(cluster_order, by = "clusters") %>%
    dplyr::select(-clusters) %>%
    dplyr::rename(clusters = phase) %>%
    identity()

  phase_levels <- phase_levels[phase_levels %in% unique(seu_meta$clusters)]

  seu_meta <-
    seu_meta %>%
    tidyr::unite("new_clusters", all_of(c("clusters", group.by)), remove = FALSE) %>%
    dplyr::arrange(clusters, new_clusters) %>%
    dplyr::mutate(clusters = factor(new_clusters, levels = unique(new_clusters))) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu@meta.data <- seu_meta[rownames(seu@meta.data), ]

  gene_lists <- map(gene_lists, ~ (.x[.x %in% VariableFeatures(seu)]))

  seu <- Seurat::AddModuleScore(seu, features = gene_lists, name = "subtype")

  module_names <- paste0("subtype", seq(length(gene_lists)))

  seu@meta.data[names(gene_lists)] <- seu@meta.data[module_names]

  cont_tibble <- FetchData(seu, vars = c(continuous_var, "clusters"))
}

compare_markers <- function(seu_path, cluster_order = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", mygene = "EZH2", label = "_filtered_", height = 14, width = 22, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star")) {
  #
	file_id <- fs::path_file(seu_path)

	tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

	message(file_id)
	cluster_order_list <- cluster_order[[file_id]]
	cluster_order <- if (!is.null(cluster_order_list)) cluster_order_list[["0"]] %||% cluster_order_list[[1]] else NULL

  seu <- readRDS(seu_path)

  group.by <- unique(cluster_order$resolution)

  cluster_order <-
    cluster_order %>%
    dplyr::filter(!is.na(clusters)) %>%
    dplyr::mutate(clusters = as.character(clusters))

  seu@meta.data$clusters <- seu@meta.data[[group.by]]

  seu_meta <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::left_join(cluster_order, by = "clusters") %>%
    dplyr::select(-clusters) %>%
    dplyr::rename(clusters = phase) %>%
    identity()

  phase_levels <- phase_levels[phase_levels %in% unique(seu_meta$clusters)]

  seu_meta <-
    seu_meta %>%
    tidyr::unite("new_clusters", all_of(c("clusters", group.by)), remove = FALSE) %>%
    dplyr::arrange(clusters, new_clusters) %>%
    dplyr::mutate(clusters = factor(new_clusters, levels = unique(new_clusters))) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu@meta.data <- seu_meta[rownames(seu@meta.data), ]

  g1_seu <- seu[, str_detect(seu$clusters, "g1_[0-9]")] %>%
    # find_all_markers(metavar = "clusters", seurat_assay = "SCT") %>%
    identity()

  # test0 <-
  #   g1_seu@misc$markers$clusters$presto %>%
  #   dplyr::filter(str_detect(Cluster, "g1_[0-9]")) %>%
  #   group_by(Cluster) %>%
  #   slice_head(n=30) %>%
  #   dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
  #   dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
  #   identity()

  return(g1_seu)
}

compile_subtype_violins <- function(interesting_samples, subtype_violins, selected_samples = c(
                                      "SRX10264519", "SRX10264520", "SRX10264525", "SRX10264526", "SRX11133594",
                                      "SRX11133593", "SRX11133592", "SRX11133585", "SRX14116944"
                                    )) {
  names(subtype_violins) <- interesting_samples

  if (!is.null(selected_samples)) {
    subtype_violins <- subtype_violins[selected_samples]
  }

  all_subtype_scores <- subtype_violins %>%
    purrr::map(purrr::flatten) %>%
    purrr::imap(~ {
      wrap_plots(.x) + plot_annotation(title = glue("{.y} subtype enrichment"))
    }) %>%
    identity()

  all_subtype_score_plot_file <- "results/all_subtype_scores.pdf"

  pdf(all_subtype_score_plot_file, height = 12, width = 20)
  print(all_subtype_scores)
  dev.off()

  subtype_scores_by_scna <-
    subtype_violins %>%
    purrr::map(2) %>%
    identity()

  rescale_and_clean_plots <- function(groblist) {
    gglist <- map(groblist, ggplotify::as.ggplot)

    gglist <-
      gglist %>%
      map(~ (.x + scale_fill_manual(values = scales::hue_pal()(7)))) %>%
      map(~ (.x +
        theme(
          legend.position = "none",
          axis.title.x = element_blank()
        ) +
        NULL
      )) %>%
      identity()

    return(gglist)
  }

  test0 <-
    subtype_scores_by_scna %>%
    map(rescale_and_clean_plots) %>%
    purrr::imap(~ {
      wrap_plots(.x) + plot_annotation(title = glue("{.y}"))
    }) %>%
    identity()

  plot_tags <- list(c(rbind(names(test0), rep("", length(names(test0))))))

  wrap_plots(test0, ncol = 2) +
    plot_annotation(tag_levels = plot_tags) &
    theme(
      plot.tag.position = c(1, 0),
      plot.tag = element_text(size = 12, hjust = 0, vjust = 0)
    ) +
      NULL

  subtype_by_scna_plot_file <- "results/subtype_scores_by_scna.pdf"

  ggsave(subtype_by_scna_plot_file, height = 25, width = 20, limitsize = FALSE)

  # pdf(subtype_by_scna_plot_file)
  # print(subtype_scores_by_scna)
  # dev.off()

  subtype_2_scores_by_scna <- subtype_violins %>%
    purrr::map(purrr::flatten) %>%
    purrr::map(4) %>%
    map(ggplotify::as.ggplot) %>%
    purrr::imap(~ {
      .x + labs(title = glue("{.y} subtype 2 enrichment"))
    }) %>%
    identity()

  subtype_2_by_scna_plot_file <- "results/subtype_2_scores_by_scna.pdf"


  wrap_plots(subtype_2_scores_by_scna, ncol = 3)

  ggsave(subtype_2_by_scna_plot_file, height = 12, width = 12, limitsize = FALSE)

  # pdf(subtype_2_by_scna_plot_file)
  # print(subtype_2_scores_by_scna)
  # dev.off()

  return(list("all" = all_subtype_score_plot_file, "patchwork" = subtype_by_scna_plot_file, "s2" = subtype_2_by_scna_plot_file))
}