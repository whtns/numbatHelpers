# Plotting and annotation functions (1)

read_cluster_dictionary <- function(cluster_dictionary_path = "data/cluster_dictionary.csv") {
  cluster_dictionary <- read_tsv(cluster_dictionary_path) %>%
    split(.$sample_id)

  return(cluster_dictionary)
}

assign_phase_cluster_at_resolution <- function(seu_path = NULL, cluster_order = NULL, assay = "SCT", resolution = 1, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star")) {
  #

	file_id <- fs::path_file(seu_path)
	
	tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
	
	message(file_id)
	cluster_order <- cluster_order[[file_id]]

  seu <- readRDS(seu_path)

  # start loop ------------------------------


  if (!is.null(cluster_order)) {
    single_cluster_order <- cluster_order[[resolution]]

    group.by <- unique(single_cluster_order$resolution)
  }

  if (!is.null(single_cluster_order)) {
    single_cluster_order <-
      single_cluster_order |>
      dplyr::mutate(order = dplyr::row_number()) %>%
      dplyr::filter(!is.na(clusters)) %>%
      dplyr::mutate(clusters = as.character(clusters))

    group.by <- unique(single_cluster_order$resolution)

    seu@meta.data$clusters <- seu@meta.data[[group.by]]

    seu_meta <- seu@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::select(-any_of(c("phase_level", "order"))) %>%
      dplyr::left_join(single_cluster_order, by = "clusters") %>%
      dplyr::select(-clusters) %>%
      dplyr::rename(phase_level = phase) %>%
      identity()

    phase_levels <- phase_levels[phase_levels %in% unique(seu_meta$phase_level)]

    seu_meta <-
      seu_meta %>%
      tidyr::unite("clusters", all_of(c("phase_level", group.by)), remove = FALSE) %>%
      dplyr::arrange(phase_level, order) %>%
      dplyr::mutate(clusters = factor(clusters, levels = unique(clusters))) %>%
      tibble::column_to_rownames("cell") %>%
      identity()

    seu@meta.data <- seu_meta[rownames(seu@meta.data), ]

    seu <- tryCatch(
      find_all_markers(seu, metavar = "clusters", seurat_assay = "SCT"),
      error = function(e) {
        if (grepl("JoinLayers", conditionMessage(e), fixed = TRUE)) {
          warning("SCT marker JoinLayers failed; using stash_marker_features fallback.")
          seu@misc$markers[["clusters"]] <- seuratTools:::stash_marker_features("clusters", seu, seurat_assay = "SCT")
          return(seu)
        }
        stop(e)
      }
    )

    seu@meta.data$clusters <- forcats::fct_drop(seu@meta.data$clusters)

    # mysec ------------------------------

    heatmap_features <-
      seu@misc$markers[["clusters"]][["presto"]] %>%
      dplyr::filter(Gene.Name %in% VariableFeatures(seu))

    tidy_eval_arrange <- function(.data, ...) {
      .data %>%
        arrange(...)
    }
  }

  add_hash_metadata(seu = seu, filepath = seu_path)

  return(seu_path)
}

calculate_clone_distribution <- function(seu_path = NULL, cluster_order = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", height = 5, width = 9, equalize_scna_clones = FALSE, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL, large_clone_comparisons = NULL, scna_of_interest = "1q", min_cells_per_cluster = 50, return_plots = FALSE, split_columns = "clusters", pairwise = TRUE) {
  if (is.list(seu_path)) {
    seu_path <- unlist(seu_path, use.names = FALSE)[[1]]
  }

  kept_phases <- kept_phases %||% phase_levels

  file_id <- fs::path_file(seu_path)
  
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
  
  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
  
  message(file_id)
  if (!is.null(cluster_order) && file_id %in% names(cluster_order)) {
    cluster_order <- cluster_order[[file_id]]
  } else {
    cluster_order <- NULL
  }

  full_seu <- readRDS(seu_path)

  full_seu$scna[full_seu$scna == ""] <- ".diploid"
  full_seu$scna <- factor(full_seu$scna)

  # subset by retained clones ------------------------------
  retained_clones <- sort(unique(full_seu$clone_opt))
  if (!is.null(large_clone_comparisons) && sample_id %in% names(large_clone_comparisons)) {
    clone_comparisons <- names(large_clone_comparisons[[sample_id]])
    clone_comparison <- clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
    if (length(clone_comparison) > 0) {
      retained_clones <- clone_comparison %>%
        str_extract("[0-9]_v_[0-9]") %>%
        str_split("_v_", simplify = TRUE)
    }
  }

  file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
  plot_path <- tempfile(tmpdir = "results/numbat_sridhar", fileext = ".pdf")
  table_path <- tempfile(tmpdir = "results/numbat_sridhar", fileext = ".xlsx")

  pdf(plot_path, height = height, width = width)

  pairwise_seu_tables <- list()
  processed_any_resolution <- FALSE
  for (resolution in names(cluster_order)) {
    processed_any_resolution <- TRUE
    # start loop ------------------------------

    seu <- full_seu[, full_seu$clone_opt %in% retained_clones]

    if (!is.null(cluster_order)) {
      single_cluster_order <- cluster_order[[resolution]]

      group.by <- unique(single_cluster_order$resolution)
    }


    if (!is.null(single_cluster_order)) {
      single_cluster_order <-
        single_cluster_order |>
        dplyr::mutate(order = dplyr::row_number()) %>%
        dplyr::filter(!is.na(clusters)) %>%
        dplyr::mutate(clusters = as.character(clusters))

      group.by <- unique(single_cluster_order$resolution)

      seu@meta.data$clusters <- seu@meta.data[[group.by]]

      seu_meta <- seu@meta.data %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::select(-any_of(c("phase_level", "order"))) %>%
        dplyr::left_join(single_cluster_order, by = "clusters") %>%
        dplyr::select(-clusters) %>%
        dplyr::rename(phase_level = phase) %>%
        identity()

      phase_levels <- phase_levels[phase_levels %in% unique(seu_meta$phase_level)]

      seu_meta <-
        seu_meta %>%
        tidyr::unite("clusters", all_of(c("phase_level", group.by)), remove = FALSE) %>%
        dplyr::arrange(phase_level, order) %>%
        dplyr::mutate(clusters = factor(clusters, levels = unique(clusters))) %>%
        tibble::column_to_rownames("cell") %>%
        identity()

      seu@meta.data <- seu_meta[rownames(seu@meta.data), ]

      seu <- seu[, seu$phase_level %in% kept_phases]
      seu <- tryCatch(
        find_all_markers(seu, metavar = "clusters", seurat_assay = "SCT"),
        error = function(e) {
          if (grepl("JoinLayers", conditionMessage(e), fixed = TRUE)) {
            warning("SCT marker JoinLayers failed; using stash_marker_features fallback.")
            seu@misc$markers[["clusters"]] <- seuratTools:::stash_marker_features("clusters", seu, seurat_assay = "SCT")
            return(seu)
          }
          stop(e)
        }
      )

      seu@meta.data$clusters <- forcats::fct_drop(seu@meta.data$clusters)

      # mysec ------------------------------

      heatmap_features <-
        seu@misc$markers[["clusters"]][["presto"]] %>%
        dplyr::filter(Gene.Name %in% VariableFeatures(seu))

      tidy_eval_arrange <- function(.data, ...) {
        .data %>%
          arrange(...)
      }
      #
      single_cluster_order_vec <-
        seu@meta.data %>%
        dplyr::select(clusters, !!group.by) %>%
        dplyr::arrange(clusters, !!sym(group.by)) %>%
        dplyr::select(clusters, !!group.by) |>
        dplyr::distinct(.data[[group.by]], .keep_all = TRUE) |>
        dplyr::mutate(!!group.by := as.character(.data[[group.by]])) |>
        tibble::deframe() |>
        identity()

      heatmap_features[["Cluster"]] <-
        factor(heatmap_features[["Cluster"]], levels = levels(seu_meta$clusters))

      heatmap_features <-
        heatmap_features %>%
        dplyr::group_by(Gene.Name) |>
        dplyr::slice_max(order_by = Average.Log.Fold.Change, n = 1) |>
        dplyr::ungroup() |>
        dplyr::arrange(Cluster, desc(Average.Log.Fold.Change)) |>
        group_by(Cluster) %>%
        slice_head(n = 5) %>%
        dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(Gene.Name, .keep_all = TRUE) |>
        identity()
    } else {
      heatmap_features <-
        seu@misc$markers[[group.by]][["presto"]]

      single_cluster_order <- levels(seu@meta.data[[group.by]]) %>%
        set_names(.)

      seu@meta.data[[group.by]] <-
        factor(seu@meta.data[[group.by]], levels = single_cluster_order)

      group_by_clusters <- seu@meta.data[[group.by]]

      seu@meta.data$clusters <- names(single_cluster_order[group_by_clusters])

      seu@meta.data$clusters <- factor(seu@meta.data$clusters, levels = unique(setNames(names(single_cluster_order), single_cluster_order)[levels(seu@meta.data[[group.by]])]))

      heatmap_features <-
        heatmap_features %>%
        dplyr::arrange(Cluster) %>%
        group_by(Cluster) %>%
        slice_head(n = 6) %>%
        dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
        dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
        identity()
    }

    all_seu_plot <- plot_distribution_of_clones_across_clusters(seu, seu_name = sample_id, var_x = "scna", var_y = "clusters")

    if (pairwise) {
      pairwise_seu_plots <- list()

      scna_clones <- unique(sort(as.factor(seu@meta.data$scna)))

      pairwise_clone_vectors <-
        bind_cols(scna_clones[-length(scna_clones)], scna_clones[-1]) %>%
        t() %>%
        as.data.frame() %>%
        as.list() %>%
        map(as.character) %>%
        identity()

      names(pairwise_clone_vectors) <- map(pairwise_clone_vectors, ~ paste(., collapse = "_v_"))


      pairwise_res <- imap(pairwise_clone_vectors, make_pairwise_plots, seu, sample_id, var_y = "phase_level")

      pairwise_seu_plots <- map(pairwise_res, "plot") |>
        map(~ {
          .x + labs(title = sample_id, subtitle = resolution)
        })

      pairwise_seu_tables[resolution] <- map(pairwise_res, "table")

      print(pairwise_seu_plots)
    }
  }

  if (!processed_any_resolution) {
    dev.off()
    return(list("table" = NA_character_, "plot" = NA_character_, "sample_id" = sample_id))
  }

  dev.off()

  table_path <- writexl::write_xlsx(pairwise_seu_tables, table_path)

  return(list("table" = table_path, "plot" = plot_path, "sample_id" = sample_id))
}