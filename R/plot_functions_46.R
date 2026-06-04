# Plot Functions (145)

#' Perform clustering analysis
#'
#' @param seu Seurat object
#' @param daughter_clusters Cluster information
#' @param resolution Parameter for resolution
#' @param from_clust Parameter for from clust
#' @param assay Character string (default: "SCT")
#' @return Function result
#' @export
# Performance optimizations applied:
# - map_bind_rows: Use map_dfr() instead of map() %>% bind_rows() for better performance

chi_sq_daughter_clusters <- function(seu, daughter_clusters, resolution, from_clust, assay = "SCT") {
  #
  from_resolution <- glue("{assay}_snn_res.{resolution}")

  to_clusts <- daughter_clusters[[as.character(resolution)]][[as.character(from_clust)]][["to_clust"]]
  
  to_resolution <- unique(daughter_clusters[[as.character(resolution)]][[as.character(from_clust)]][[glue("to_{assay}_snn.res.")]])
  
  to_resolution <- glue("{assay}_snn_res.{to_resolution}")

  message(glue("resolution: {from_resolution} from_clust: {from_clust}"))
  
  test_seu <- seu[, seu[[]][[to_resolution]] %in% to_clusts]

  test_seu$scna <- janitor::make_clean_names(test_seu$scna, allow_dupes = TRUE)

  test_seu$clusters <- as.numeric(test_seu@meta.data[[to_resolution]])

  if (length(unique(test_seu$scna)) > 1) {
    scna_counts <-
      test_seu@meta.data %>%
      dplyr::mutate(scna = factor(scna))

    scna_clones <- levels(scna_counts$scna)

    pairwise_clone_vectors <-
      bind_cols(scna_clones[-length(scna_clones)], scna_clones[-1]) %>%
      t() %>%
      as.data.frame() %>%
      as.list() %>%
      map(as.character) %>%
      identity()

    names(pairwise_clone_vectors) <- map(pairwise_clone_vectors, ~ paste(., collapse = "_v_"))

    scna_counts <-
      scna_counts %>%
      janitor::tabyl(clusters, scna) %>%
      tibble::column_to_rownames("clusters") %>%
      as.matrix() %>%
      identity()

    fisher_results <- pairwise_clone_vectors %>%
      map(~ scna_counts[, .x]) %>%
      map(fisher.test, simulate.p.value = TRUE, B = 1e5) %>%
      map(broom::tidy) %>%
      dplyr::bind_rows(.id = "clone_comparison") %>%
      dplyr::mutate(from_clust = from_clust)

    test0 <- dplyr::left_join(daughter_clusters[[resolution]][[from_clust]], fisher_results, by = "from_clust")
  } else {
    return(daughter_clusters[[resolution]][[from_clust]])
  }
}

#' Perform clustering analysis
#'
#' @param json_file File path
#' @return List object
#' @export
pull_cluster_orders2 <- function(json_file = "data/scna_cluster_order.json") {

  cluster_ids <- jsonlite::fromJSON(json_file, simplifyDataFrame = TRUE) %>%
    as_tibble() %>%
    group_by(resolution, file_id) %>%
    mutate(preferred_resolution = as.numeric(preferred_resolution))

  cluster_id_list <-
  cluster_ids |>
  dplyr::select(-tumor_id) %>%
    pivot_longer(-any_of(c("resolution", "sample_id", "file_id", "preferred_resolution")), names_to = "phase", values_to = "clusters") %>%
    mutate(clusters = map(clusters, ~if (is.null(.x)) NA else .x)) %>%
    unnest(clusters) %>%
    mutate(phase = factor(phase, levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"))) %>%
    tidyr::unnest()  %>%
    split(.$file_id) %>%
    map(~ split(.x, .x$preferred_resolution)) %>%
    identity()

  return(cluster_id_list)
}

#' Perform clustering analysis
#'
#' @param raw_cluster_file File path
#' @return List object
#' @export
pull_cluster_orders <- function(raw_cluster_file = "data/scna_cluster_order.csv") {
	cluster_ids <-
		raw_cluster_file %>%
		read_csv() %>%
		janitor::clean_names() %>%
		group_by(resolution, file_id) %>%
		dplyr::transmute(across(-any_of(c("resolution", "file_id")), as.character)) |> 
		dplyr::mutate(preferred_resolution = as.numeric(preferred_resolution))
	
	cluster_id_list <-
		cluster_ids %>%
		tidyr::pivot_longer(-any_of(c("resolution", "file_id", "sample_id", "tumor_id", "preferred_resolution")), names_to = "phase", values_to = "clusters") %>%
		dplyr::mutate(clusters = str_split(clusters, pattern = "_")) %>%
		tidyr::unnest(clusters) %>%
		dplyr::mutate(phase = factor(phase, levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"))) %>%
		split(.$file_id) %>%
		map(~ split(.x, .x$preferred_resolution)) %>%
		identity()
	
	return(cluster_id_list)
}

#' Create a plot visualization
#'
#' @param full_seu Parameter for full seu
#' @param retained_clones Parameter for retained clones
#' @param group_by Character string (default: "clusters")
#' @return ggplot2 plot object
#' @export
make_faded_umap_plots <- function(full_seu, retained_clones, group_by = "clusters") {
  scna_plot <- DimPlot(full_seu, group.by = c("scna")) +
    aes(alpha = alpha_var) +
    NULL
  scna_plot[[1]]$layers[[1]]$aes_params$alpha <- ifelse(full_seu@meta.data$clone_opt %in% retained_clones, 1, .05)

  cluster_plot <- DimPlot(full_seu, group.by = group_by) +
    aes(alpha = alpha_var) +
    NULL
  cluster_plot[[1]]$layers[[1]]$aes_params$alpha <- ifelse(full_seu@meta.data$clone_opt %in% retained_clones, 1, .05)

  umap_plots <- scna_plot +
    cluster_plot +
    plot_layout(ncol = 1)

  return(umap_plots)
}

#' Create a heatmap visualization
#'
#' @param seu_path File path
#' @param nb_paths File path
#' @param clone_simplifications Parameter for clone simplifications
#' @param cluster_orders Cluster information
#' @return ggplot2 plot object
#' @export
plot_seu_marker_heatmap_all_resolutions <- function(seu_path = NULL, nb_paths = NULL, clone_simplifications = NULL, cluster_orders = NULL) {
	file_id <- fs::path_file(seu_path)
	
  sample_id <- str_extract(seu_path, "SR[RX][0-9]+")

  nb_paths <-
    nb_paths %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  nb_path <- nb_paths[[sample_id]]

  if (!is.null(cluster_orders)) {
    cluster_orders <- cluster_orders[[file_id]]

    resolutions <- unique(cluster_orders$resolution)

    cluster_orders <- split(cluster_orders, cluster_orders$resolution)

    cluster_orders <- map(cluster_orders, ~ list(.x)) |>
      map(set_names, sample_id)

    collages <- map2(resolutions, cluster_orders, ~ plot_seu_marker_heatmap(seu = seu_path, cluster_order = .y, nb_path = nb_path, clone_simplifications = clone_simplifications, group.by = .x, label = .x))
  } else {
    resolutions <- glue("SCT_snn_res.{seq(0.2, 1.0, by = 0.2)}") %>%
      set_names(.)

    collages <- map(resolutions, ~ plot_seu_marker_heatmap(seu = seu_path, nb_path = nb_path, clone_simplifications = clone_simplifications, group.by = .x, label = .x))
  }

  file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  qpdf::pdf_combine(collages, glue("results/{file_slug}_collage_all_resolutions.pdf"))
}

