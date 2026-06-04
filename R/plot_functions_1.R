# Plotting and annotation functions (1)

#' Create a numbat-related plot visualization
#'
#' @param nb Numbat object
#' @param myseu Seurat object
#' @param myannot Parameter for myannot
#' @param mytitle Plot title
#' @param sort_by Character string (default: "scna")
#' @param show_segment_names_on_x Logical; if TRUE, keep x-axis labels visible on the numbat heatmap panel.
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O

plot_numbat <- function(nb, myseu, myannot, mytitle, sort_by = "clone_opt", show_segment_names_on_x = FALSE, ...) {

  # Ensure clone_opt is a factor to avoid continuous/discrete scale errors
  # Must do this BEFORE accessing clone levels for palette
  if ("clone_opt" %in% colnames(nb$clone_post)) {
    nb$clone_post$clone_opt <- as.character(nb$clone_post$clone_opt)
  }
  if ("clone_opt" %in% colnames(myannot)) {
    myannot$clone_opt <- as.character(myannot$clone_opt)
  }
  
  all_clones <- na.omit(unique(c(as.character(myannot$clone_opt), as.character(nb$clone_post$clone_opt))))
  nclones <- max(as.integer(all_clones), na.rm = TRUE)
  clone_pal <- scales::hue_pal()(nclones) %>% set_names(as.character(1:nclones))
  myheatmap <- nb$plot_phylo_heatmap(
    pal_clone = clone_pal,
    pal_annot = clone_pal,
    annot = myannot,
    show_phylo = FALSE,
    sort_by = sort_by,
    annot_bar_width = 1,
    raster = FALSE,
    ...
  ) +
    labs(title = mytitle) +
    theme(legend.position = "none")


  .add_segment_labels_to_heatmap <- function(heatmap_obj, show_labels) {
    if (!isTRUE(show_labels)) {
      if (length(heatmap_obj) >= 2) {
        heatmap_obj[[2]] <- heatmap_obj[[2]] + theme(legend.position = "none", axis.text.x = element_blank())
      }
      return(heatmap_obj)
    }

    heatmap_obj <- heatmap_obj &
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
      )

    if (length(heatmap_obj) >= 3 && !is.null(heatmap_obj[[3]]$data) && all(c("seg", "CHROM", "seg_start", "seg_end") %in% colnames(heatmap_obj[[3]]$data))) {
      seg_labels <- heatmap_obj[[3]]$data %>%
        dplyr::distinct(CHROM, seg, seg_start, seg_end) %>%
        dplyr::arrange(CHROM, seg) %>%
        dplyr::mutate(seg_mid = (seg_start + seg_end) / 2)

      heatmap_obj[[3]] <- heatmap_obj[[3]] +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
        ) +
        geom_text(
          data = seg_labels,
          aes(x = seg_mid, y = -0.03, label = seg),
          inherit.aes = FALSE,
          size = 2.8,
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          check_overlap = TRUE
        ) +
        theme(
          strip.text = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(5.5, 5.5, 24, 5.5)
        ) +
        coord_cartesian(clip = "off")
    }

    heatmap_obj
  }
  myheatmap <- .add_segment_labels_to_heatmap(myheatmap, show_segment_names_on_x)
  return(myheatmap)
}

#' Create a numbat-related plot visualization
#'
#' @param nb Numbat object
#' @param myseu Seurat object
#' @param myannot Parameter for myannot
#' @param mytitle Plot title
#' @param show_segment_names_on_x Logical; if TRUE, keep segment labels on applicable genomic regions.
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_numbat_w_phylo <- function(nb, myseu, myannot, mytitle, show_segment_names_on_x = FALSE, ...) {
  
  if ("clone_opt" %in% colnames(nb$clone_post)) {
    nb$clone_post$clone_opt <- as.character(nb$clone_post$clone_opt)
  }
  if ("clone_opt" %in% colnames(myannot)) {
    myannot$clone_opt <- as.character(myannot$clone_opt)
  }
  all_clones <- na.omit(unique(c(as.character(myannot$clone_opt), as.character(nb$clone_post$clone_opt))))
  nclones <- max(as.integer(all_clones), na.rm = TRUE)
  mypal <- scales::hue_pal()(nclones) %>% set_names(as.character(1:nclones))
  myheatmap <- nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = myannot,
    show_phylo = TRUE,
    annot_bar_width = 1,
    ...
  ) + labs(title = mytitle)

  if (isTRUE(show_segment_names_on_x)) {
    myheatmap <- myheatmap &
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
      )

    if (length(myheatmap) >= 3 && !is.null(myheatmap[[3]]$data) && all(c("seg", "CHROM", "seg_start", "seg_end") %in% colnames(myheatmap[[3]]$data))) {
      seg_labels <- myheatmap[[3]]$data %>%
        dplyr::distinct(CHROM, seg, seg_start, seg_end) %>%
        dplyr::arrange(CHROM, seg) %>%
        dplyr::mutate(seg_mid = (seg_start + seg_end) / 2)

      myheatmap[[3]] <- myheatmap[[3]] +
        geom_text(
          data = seg_labels,
          aes(x = seg_mid, y = -0.03, label = seg),
          inherit.aes = FALSE,
          size = 2.8,
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          check_overlap = TRUE
        ) +
        theme(
          strip.text = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(5.5, 5.5, 24, 5.5)
        ) +
        coord_cartesian(clip = "off")
    }
  }

  return(myheatmap)
}

safe_plot_numbat <- purrr::safely(plot_numbat, otherwise = NA_real_)
safe_plot_numbat_w_phylo <- purrr::safely(plot_numbat_w_phylo, otherwise = NA_real_)

#' Create a plot visualization
#'
#' @param phylo_plot_output Parameter for phylo plot output
#' @param chrom Character string (default: "1")
#' @param p_min Parameter for p min
#' @return ggplot2 plot object
#' @export
plot_variability_at_SCNA <- function(phylo_plot_output, chrom = "1", p_min = 0.9) {
  
  plot_input <- phylo_plot_output
  plot_input$seg <- factor(plot_input$seg, levels = str_sort(unique(plot_input$seg), numeric = TRUE))
  plot_input <- plot_input[order(plot_input$seg), ]
  plot_input$cnv_state <- dplyr::case_when(
    plot_input$cnv_state == "amp" ~ "gain",
    plot_input$cnv_state == "bamp" ~ "balanced gain",
    plot_input$cnv_state == "del" ~ "loss",
    plot_input$cnv_state == "loh" ~ "CNLoH",
    TRUE ~ plot_input$cnv_state
  )
  p_cnv_plot <- plot_input %>%
    dplyr::group_by(seg) %>%
    dplyr::filter(!all(is.na(LLR))) %>%
    ggplot(aes(x = cell_index, y = p_cnv)) +
    geom_point(aes(color = cnv_state), size = 0.1, alpha = 0.1) +
    geom_hline(aes(yintercept = p_min), linetype = 'dashed', color = "grey") +
    scale_x_reverse() +
    scale_color_manual(values = c(
      "gain" = "#7f180f",
      "balanced gain" = "pink",
      "loss" = "#010185",
      "CNLoH" = "#387229"
    )) +
    labs(color = "SCNA state", fill = "Clone", y = "Probability of SCNA") +
    facet_wrap(~seg) +
    geom_tile(aes(y = -0.2, height = 0.1, fill = factor(clone_opt))) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))

  return(p_cnv_plot)
}


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