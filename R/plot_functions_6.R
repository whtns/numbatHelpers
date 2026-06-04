
#' Create a numbat-related plot visualization
#'
#' @param output_file File path
#' @param sample_id Parameter for sample id
#' @param myseus Parameter for myseus
#' @param mynbs Parameter for mynbs
#' @param merged_metadata Parameter for merged metadata
#' @param myexpressions Parameter for myexpressions
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging
# - multiple_joins: Combine multiple joins into single join operation where possible

make_filtered_numbat_plots <- function(output_file, sample_id, myseus, mynbs, merged_metadata, myexpressions) {
  seu <- readRDS(myseus[[sample_id]])

  merged_metadata_transfer <-
    merged_metadata %>%
    dplyr::filter(sample_id == {{ sample_id }}) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)

  mynb <- readRDS(mynbs[[sample_id]])

  nb_meta <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  myannot <- mynb$clone_post[, c("cell", "GT_opt", "clone_opt")]

  if (!"" %in% unlist(myexpressions)) {
    final_phylo_heatmap <- filter_phylo_plot(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9, expressions = myexpressions)

    test0 <- final_phylo_heatmap[[3]][["data"]] %>%
      plot_variability_at_SCNA()

    patchwork::wrap_plots(final_phylo_heatmap / test0)

    ggsave(output_file)
  }

  return(output_file)
}


#' Title
#' @export
retrieve_numbat_seurat <- function(numbat_rds_file, cluster_dictionary) {
  #
  sample_id <- str_extract(numbat_rds_file, "SR[RX][0-9]+")

  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)


  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")

  seu <- AddMetaData(seu, test0)

  return(seu)
}

plot_study_cell_stats <- function(study_cell_stats, plot_path = tempfile(fileext = ".pdf"), umi_threshold = 1e3, genes_threshold = 1e3, mito_threshold = 5, mito_expansion = 0.8, bandwidth = 0.35, plot_height = NULL, ...) {
  #

  study_cell_stats <-
    study_cell_stats |>
    dplyr::mutate(study = dplyr::case_when(
      study == "collin" ~ "Collin et al. 2021",
      study == "wu" ~ "Wu et al. 2022",
      study == "yang" ~ "Yang et al. 2021",
      study == "field" ~ "Field et al. 2022",
      study == "liu" ~ "Liu et al. 2024"
    )) |>
    dplyr::mutate(study = factor(study, levels = c("Collin et al. 2021", "Wu et al. 2022", "Yang et al. 2021", "Field et al. 2022", "Liu et al. 2024")))

  mypal <- scales::hue_pal()(5) %>%
    set_names(levels(study_cell_stats$study))

  umis_per_cell <- ggplot(study_cell_stats, aes(x = nCount_gene, y = sample_id, fill = study, group = sample_id)) +
    geom_density_ridges(scale = 1, panel_scaling = FALSE, bandwidth = bandwidth) +
    scale_x_log10() +
    scale_y_discrete(limits = rev) +
    scale_fill_manual(values = mypal) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), # remove y axis labels,
      # axis.text.y=element_blank(),  #remove y axis labels
      # axis.ticks.y=element_blank()  #remove y axis ticks
    ) +
    geom_vline(xintercept = umi_threshold, linetype = "dotted") +
    labs(title = "UMIs/cell")


  genes_per_cell <- ggplot(study_cell_stats, aes(x = nFeature_gene, y = sample_id, fill = study, group = sample_id)) +
    geom_density_ridges(scale = 1, panel_scaling = FALSE, bandwidth = bandwidth) +
    scale_x_log10() +
    scale_y_discrete(limits = rev) +
    scale_fill_manual(values = mypal) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), # remove y axis labels,
      # axis.text.y=element_blank(),  #remove y axis labels
      # axis.ticks.y=element_blank()  #remove y axis ticks
    ) +
    geom_vline(xintercept = genes_threshold, linetype = "dotted") +
    labs(title = "genes/cell")

  percent_mito_per_cell <- ggplot(study_cell_stats, aes(x = `percent.mt`, y = sample_id, fill = study, group = sample_id)) +
    geom_density_ridges(scale = 1, panel_scaling = FALSE, bandwidth = bandwidth) +
    scale_y_discrete(limits = rev, expand = expansion(add = c(0.55, mito_expansion))) +
    # scale_y_discrete(limits=rev) +
    scale_fill_manual(values = mypal) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), # remove y axis labels,
      # axis.text.y=element_blank(),  #remove y axis labels
      # axis.ticks.y=element_blank()  #remove y axis ticks
    ) +
    geom_vline(xintercept = mito_threshold, linetype = "dotted") +
    labs(title = "% mito/cell") +
    xlim(0, 25) +
    NULL

  retention_labels <- study_cell_stats %>%
    dplyr::select(study, sample_id, exclusion_criteria) %>%
    dplyr::distinct() %>%
    dplyr::mutate(exclusion_criteria = replace_na(exclusion_criteria, "retained")) %>%
    dplyr::mutate(sample_id = factor(sample_id)) %>%
    identity()

  retention_labels_plot <-
    retention_labels %>%
    ggplot(aes(y = sample_id, x = 1, fill = exclusion_criteria)) +
    geom_tile(color = "black") +
    scale_y_discrete(limits = rev) +
    theme_void()

  list(umis_per_cell, genes_per_cell, percent_mito_per_cell, retention_labels_plot) %>%
    wrap_plots() +
    plot_layout(axes = "collect", guides = "collect", widths = c(3, 3, 3, 0.2))

  plot_height <- ifelse(is.null(plot_height), 0.33 * n_distinct(study_cell_stats$sample_id), plot_height)

  ggsave(plot_path, height = plot_height, ...)
}

#' Title
#' @export
plot_several_diffex_clones_in_phase <- function(enrichments, scna_of_interest = "2p", phase = "g1") {
  enrichment_plot <-
    enrichments |>
    purrr::list_flatten() |>
    clusterProfiler::merge_result() |>
    plot_enrichment(p_val_cutoff = 1, result_slot = "compareClusterResult")

  pdf_path <- ggsave(glue("results/{scna_of_interest}_{phase}_enrichment.pdf"), h = 10, w = 8)
}

plot_seu_marker_heatmap_integrated <- function(seu_path = NULL, cluster_order = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", label = "_filtered_", height = 10, width = 18, equalize_scna_clones = FALSE, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL) {
  if (is.list(seu_path)) {
    seu_path <- unlist(seu_path, use.names = FALSE)[[1]]
  }

  kept_phases <- kept_phases %||% phase_levels

  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  message(sample_id)

  seu <- readRDS(seu_path)


  if (!is.null(cluster_order)) {
    group.by <- unique(cluster_order$resolution)
  }

  if (equalize_scna_clones) {
    seu_meta <- seu@meta.data %>%
      tibble::rownames_to_column("cell")

    clones <- table(seu_meta$scna)

    min_clone_num <- clones[which.min(clones)]

    selected_cells <-
      seu_meta %>%
      dplyr::group_by(scna) %>%
      slice_sample(n = min_clone_num) %>%
      pull(cell)

    seu <- seu[, selected_cells]
  }


  if (!is.null(cluster_order)) {
    group.by <- unique(cluster_order$resolution)

    cluster_order <-
      cluster_order %>%
      dplyr::mutate(order = dplyr::row_number()) %>%
      dplyr::filter(!is.na(clusters)) %>%
      dplyr::mutate(clusters = as.character(clusters))

    seu@meta.data$clusters <- seu@meta.data[[group.by]]

    seu_meta <- seu@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::left_join(cluster_order, by = "clusters") %>%
      dplyr::select(-clusters) %>%
      dplyr::select(-any_of(c("phase_level"))) %>%
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

    heatmap_features <-
      seu@misc$markers[["clusters"]][["presto"]] %>%
      dplyr::filter(Gene.Name %in% VariableFeatures(seu))

    tidy_eval_arrange <- function(.data, ...) {
      .data %>%
        arrange(...)
    }

    cluster_order_vec <-
      seu@meta.data %>%
      dplyr::select(clusters, !!group.by) %>%
      dplyr::arrange(clusters, !!sym(group.by)) %>%
      dplyr::pull(!!group.by) %>%
      unique() %>%
      as.character() %>%
      identity()

    heatmap_features[["Cluster"]] <-
      factor(heatmap_features[["Cluster"]], levels = levels(seu_meta$clusters))

    heatmap_features <-
      heatmap_features %>%
      dplyr::arrange(Cluster) %>%
      dplyr::group_by(Cluster) %>%
      slice_max(Average.Log.Fold.Change, n = 5) %>%
      identity()
  } else {

    seu <- find_all_markers(seu, metavar = group.by)

    heatmap_features <-
      seu@misc$markers[[group.by]][["presto"]]

    if (!is.ordered(seu@meta.data[[group.by]])) {
      seu@meta.data[[group.by]] <- factor(as.numeric(seu@meta.data[[group.by]]))
    }

    cluster_order <- levels(seu@meta.data[[group.by]]) %>%
      set_names(.)

    seu@meta.data[[group.by]] <-
      factor(seu@meta.data[[group.by]], levels = cluster_order)

    group_by_clusters <- seu@meta.data[[group.by]]

    seu@meta.data$clusters <- names(cluster_order[group_by_clusters])

    seu@meta.data$clusters <- factor(seu@meta.data$clusters, levels = unique(setNames(names(cluster_order), cluster_order)[levels(seu@meta.data[[group.by]])]))

    heatmap_features <-
      heatmap_features %>%
      dplyr::arrange(Cluster) %>%
      group_by(Cluster) %>%
      slice_head(n = 6) %>%
      dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
      dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
      identity()
  }

  seu$scna <- factor(seu$scna)
  levels(seu$scna)[1] <- "none"

  giotti_genes <- read_giotti_genes()

  heatmap_features <-
    heatmap_features %>%
    dplyr::ungroup() %>%
    left_join(giotti_genes, by = c("Gene.Name" = "symbol")) %>%
    # select(Gene.Name, term) %>%
    dplyr::mutate(term = replace_na(term, "")) %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE)

  row_ha <- ComplexHeatmap::rowAnnotation(term = rev(heatmap_features$term))

  seu_heatmap <- ggplotify::as.ggplot(
    seu_complex_heatmap(seu,
      features = heatmap_features$Gene.Name,
      group.by = c("G2M.Score", "S.Score", "scna", "clusters"),
      col_arrangement = c("clusters", "scna"),
      cluster_rows = FALSE,
      column_split = as.character(sort(as.character(seu@meta.data$clusters))),
      row_split = as.character(rev(as.character(heatmap_features$Cluster))),
      row_title_rot = 0,
      # row_split = sort(seu@meta.data$clusters)
    )
  ) +
    labs(title = sample_id) +
    theme()


  #
  labels <- data.frame(clusters = unique(seu[[]][["clusters"]]), label = unique(seu[[]][["clusters"]])) %>%
    # dplyr::rename({{group.by}} := cluster) %>%
    identity()

  cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))

  centroid_data <-
    cc_data %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
    dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
    dplyr::mutate(centroid = "centroids") %>%
    identity()

  centroid_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
    geom_point(size = 0.1) +
    theme_light() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
    NULL


  facet_cell_cycle_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
    geom_point(size = 0.1) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
    facet_wrap(~ .data[["clusters"]], ncol = 2) +
    theme_light() +
    geom_label(
      data = labels,
      aes(label = label),
      # x = Inf,
      # y = -Inf,
      x = max(cc_data$S.Score) + 0.05,
      y = max(cc_data$G2M.Score) - 0.1,
      hjust = 1,
      vjust = 1,
      inherit.aes = FALSE
    ) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    # guides(color = "none") +
    NULL

  appender <- function(string) str_wrap(string, width = 40)

  labels <- data.frame(scna = unique(seu$scna), label = str_replace(unique(seu$scna), "^$", "diploid"))

  clone_distribution_plot <-
    plot_distribution_of_clones_across_clusters(seu, tumor_id, var_x = "scna", var_y = "clusters")

  umap_plots <- DimPlot(seu, group.by = c("scna", "clusters"), combine = FALSE) %>%
    # map(~(.x + theme(legend.position = "bottom"))) %>%
    wrap_plots(ncol = 1)

  layout <- "
            AAAAAAAAAABBBBCCCC
            AAAAAAAAAABBBBCCCC
            AAAAAAAAAABBBBDDDD
            AAAAAAAAAABBBBDDDD
            AAAAAAAAAABBBBDDDD
    "

  collage_plots <- list(seu_heatmap, facet_cell_cycle_plot, centroid_plot, clone_distribution_plot)

  mypatch <- wrap_plots(collage_plots) +
    # plot_layout(widths = c(16, 4)) +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = "A") +
    NULL


  file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  plot_path <- glue("results/{file_slug}_{label}heatmap_phase_scatter_patchwork.pdf")
  ggsave(plot_path, plot = mypatch, height = height, width = width)
  return(plot_path)
}

plot_seu_marker_heatmap_by_scna_ara <- function(seu_path = NULL, cluster_order = NULL, nb_paths = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", label = "_filtered_", height = 10, width = 18, equalize_scna_clones = FALSE, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL, rb_scna_samples, large_clone_comparisons, scna_of_interest = "1q", min_cells_per_cluster = 50, return_plots = FALSE, column_split = "clusters") {
  kept_phases <- kept_phases %||% phase_levels

  #

  file_id <- fs::path_file(seu_path)
  
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
  
  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
  
  message(file_id)
  cluster_order <- cluster_order[[file_id]]

  full_seu <- readRDS(seu_path)

  # subset by retained clones ------------------------------
  clone_comparisons <- names(large_clone_comparisons[[sample_id]])
  clone_comparison <- clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
  retained_clones <- clone_comparison %>%
    str_extract("[0-9]_v_[0-9]") %>%
    str_split("_v_", simplify = TRUE)
  seu <- full_seu[, full_seu$clone_opt %in% retained_clones]


  if (!is.null(cluster_order)) {
    heatmap_features <-
      seu@misc$markers[["clusters"]][["presto"]] %>%
      dplyr::filter(Gene.Name %in% VariableFeatures(seu))

    heatmap_features[["Cluster"]] <-
      factor(heatmap_features[["Cluster"]], levels = levels(seu$clusters))

    heatmap_features <-
      heatmap_features %>%
      dplyr::arrange(Cluster) %>%
      dplyr::group_by(Cluster) %>%
      slice_max(Average.Log.Fold.Change, n = 5) %>%
      identity()
  } else {
    heatmap_features <-
      heatmap_features %>%
      dplyr::arrange(Cluster) %>%
      group_by(Cluster) %>%
      slice_head(n = 6) %>%
      dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
      dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
      identity()
  }

  large_enough_clusters <-
    seu@meta.data %>%
    dplyr::group_by(clusters) %>%
    dplyr::count()

  large_enough_clusters <-
    large_enough_clusters %>%
    dplyr::filter(n >= min_cells_per_cluster) %>%
    dplyr::pull(clusters)

  seu <- seu[, seu$clusters %in% large_enough_clusters]

  seu$scna[seu$scna == ""] <- ".diploid"
  seu$scna <- factor(seu$scna)
  # levels(seu$scna)[1] <- "none"

  heatmap_features <-
    heatmap_features %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE)

  if (!is.null(column_split)) {
    column_split <- sort(seu@meta.data[[column_split]])
    column_title <- unique(column_split)
  } else {
    column_title <- NULL
  }

  #
  labels <- data.frame(clusters = unique(seu[[]][["clusters"]]), label = unique(seu[[]][["clusters"]])) %>%
    # dplyr::rename({{group.by}} := cluster) %>%
    identity()

  cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))

  centroid_data <-
    cc_data %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
    dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
    dplyr::mutate(centroid = "centroids") %>%
    identity()

  centroid_plot_by_cluster <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["clusters"]])) +
    geom_point(size = 0.1) +
    theme_light() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.text = element_text(size = 14), # increase legend text size
      legend.title = element_text(size = 16)
    ) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
    NULL

  centroid_plot_by_scna <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
    geom_point(size = 0.1) +
    theme_light() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.text = element_text(size = 14), # increase legend text size
      legend.title = element_text(size = 16)
    ) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
    guides(
      colour = guide_legend(override.aes = list(size = 10)),
      fill = FALSE
    ) +
    NULL


  facet_cell_cycle_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
    geom_point(size = 0.1) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
    facet_wrap(~ .data[["clusters"]], ncol = 2) +
    theme_light() +
    geom_label(
      data = labels,
      aes(label = label),
      # x = Inf,
      # y = -Inf,
      x = max(cc_data$S.Score) + 0.05,
      y = max(cc_data$G2M.Score) - 0.1,
      hjust = 1,
      vjust = 1,
      inherit.aes = FALSE
    ) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.text = element_text(size = 14), # increase legend text size
      legend.title = element_text(size = 16),
      legend.position = "none"
    ) +
    guides(colour = guide_legend(override.aes = list(size = 10))) +
    # guides(color = "none") +
    NULL

  comparison_scna <-
    janitor::tabyl(as.character(seu$scna))[2, 1]

  clone_distribution_plot <- plot_distribution_of_clones_across_clusters(
    seu,
    seu_name = glue("{tumor_id} {comparison_scna}"), var_x = "scna", var_y = "clusters", signif = TRUE, plot_type = "clone"
  )

  collage_plots <- list(
    "centroid_plot_by_cluster" = centroid_plot_by_cluster,
    "centroid_plot_by_scna" = centroid_plot_by_scna,
    "facet_cell_cycle_plot" = facet_cell_cycle_plot,
    plot_spacer(),
    "clone_distribution_plot" = clone_distribution_plot,
    plot_spacer()
  )

  layout <- "
		EECCCDDDDDD
		EECCCDDDDDD
		AACCCDDDDDD
		AACCCDDDDDD
		BBCCCDDDDDD
		BBCCCDDDDDD
		"

  plot_collage <- wrap_plots(collage_plots) +
    # plot_layout(widths = c(16, 4)) +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = "A") +
    NULL

  file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  if (is.character(return_plots)) {
    plot_path <- ggsave(glue("results/{file_slug}_{scna_of_interest}_{return_plots}_ara.pdf"), collage_plots[[return_plots]], height = height, width = width)
  } else {
    plot_path <- ggsave(glue("results/{file_slug}_{scna_of_interest}_heatmap_phase_scatter_patchwork_ara.pdf"), plot_collage, height = height, width = width)
  }
}