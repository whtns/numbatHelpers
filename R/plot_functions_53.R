# Plot Functions (152)

#' Create a plot visualization
#'
#' @param seu_path File path
#' @param pdf_path File path
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O

plot_phase_wo_arm <- function(seu_path, pdf_path = NULL) {
  
  
  pdf_path <- pdf_path %||% str_replace("seu_path", ".rds", "_cc_wo_1q.pdf")

  seu <- seu_path

  ccplot1 <- DimPlot(seu, group.by = "Phase") +
    labs(subtitle = paste(glue("{names(table(seu$Phase))}: {table(seu$Phase)}"), collapse = "; "))

  seu <- annotate_cell_cycle_without_1q(seu)

  ccplot2 <- DimPlot(seu, group.by = "Phase") +
    labs(subtitle = paste(glue("{names(table(seu$Phase))}: {table(seu$Phase)}"), collapse = "; "))

  ccplot1 + ccplot2

  ggsave(pdf_path)
}

#' Create a plot visualization
#'
#' @param seu_path File path
#' @param facet Logical flag (default: FALSE)
#' @param group_by Character string (default: "clusters")
#' @param color_by Color specification
#' @return ggplot2 plot object
#' @export
plot_cc_space_plot <- function(seu_path = "output/seurat/SRX11133594_filtered_seu.rds", facet = FALSE, group_by = "clusters", color_by = "clusters") {
  
  
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  seu <- seu_path

  cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))

  centroid_data <-
    cc_data %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
    dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
    dplyr::mutate(centroid = "centroids") %>%
    identity()

  if (!facet) {
    centroid_plot <-
      cc_data %>%
      ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[group_by]], color = .data[[color_by]])) +
      geom_point(size = 0.1) +
      theme_light() +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
      ) +
      geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
      labs(title = sample_id)
  } else {
    facet_cell_cycle_plot <-
      cc_data %>%
      ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[group_by]], color = .data[[color_by]])) +
      geom_point(size = 0.1) +
      geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[group_by]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
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
      labs(title = sample_id) +
      # guides(color = "none") +
      NULL
  }

  plot_path <- ggsave(glue("results/{sample_id}_cc_space_plot.pdf"))

  return(plot_path)
}

#' Perform seu complex heatmap operation
#'
#' @param seu Seurat object
#' @param features Parameter for features
#' @param group.by Character string (default: "ident")
#' @param cells Cell identifiers or information
#' @param layer Character string (default: "scale.data")
#' @param assay Parameter for assay
#' @param group.bar.height Parameter for group.bar.height
#' @param column_split Color specification
#' @param col_arrangement Color specification
#' @param mm_col_dend Color specification
#' @param embedding Character string (default: "pca")
#' @param ... Additional arguments passed to other functions
#' @return Function result
#' @export
seu_complex_heatmap <- function(seu, features = NULL, group.by = "ident", cells = NULL,
                                layer = "scale.data", assay = NULL, group.bar.height = 0.01,
                                column_split = NULL, col_arrangement = "ward.D2", mm_col_dend = 30,
                                embedding = "pca", ...) {
  
  
  if (length(GetAssayData(seu, layer = "scale.data")) == 0) {
    message("seurat object has not been scaled. Please run `Seurat::ScaleData` to view a scaled heatmap; showing unscaled expression data")
    layer <- "data"
  }
  cells <- cells %||% colnames(x = seu)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = seu)[cells]
  }
  assay <- assay %||% Seurat::DefaultAssay(object = seu)
  Seurat::DefaultAssay(object = seu) <- assay
  features <- features %||% VariableFeatures(object = seu)
  features <- rev(x = unique(x = features))
  possible.features <- rownames(x = GetAssayData(
    object = seu,
    layer = layer
  ))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop(
        "No requested features found in the ", layer,
        " layer for the ", assay, " assay."
      )
    }
    warning(
      "The following features were omitted as they were not found in the ",
      layer, " layer for the ", assay, " assay: ", paste(bad.features,
        collapse = ", "
      )
    )
  }
  data <- as.data.frame(x = t(x = as.matrix(x = GetAssayData(
    object = seu,
    layer = layer
  )[features, cells, drop = FALSE])))
  seu <- suppressMessages(expr = StashIdent(
    object = seu,
    save.name = "ident"
  ))
  if (any(col_arrangement %in% c(
    "ward.D", "single", "complete",
    "average", "mcquitty", "median", "centroid", "ward.D2"
  ))) {
    if ("pca" %in% Seurat::Reductions(seu)) {
      cluster_columns <- Seurat::Embeddings(seu, embedding) %>%
        dist() %>%
        hclust(col_arrangement)
    } else {
      message(glue("{embedding} not computed for this dataset; cells will be clustered by displayed features"))
    }
  } else {
  	cells <- seu %>%
  		Seurat::FetchData(vars = col_arrangement) %>%
  		sample_frac() %>%
  		dplyr::arrange(across(all_of(col_arrangement))) %>%
  		rownames()
    data <- data[cells, ]
    group.by <- base::union(group.by, col_arrangement)
    cluster_columns <- FALSE
  }
  group.by <- group.by %||% "ident"
  groups.use <- seu[[group.by]][cells, , drop = FALSE]
  groups.use <- groups.use %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::mutate(across(where(is.character), ~ str_wrap(str_replace_all(
      .x,
      ",", " "
    ), 10))) %>%
    dplyr::mutate(across(
      where(is.character),
      as.factor
    )) %>%
    data.frame(row.names = 1) %>%
    identity()
  groups.use.factor <- groups.use[sapply(groups.use, is.factor)]
  ha_cols.factor <- NULL
  if (length(groups.use.factor) > 0) {
    ha_col_names.factor <- lapply(groups.use.factor, levels)
    ha_cols.factor <- purrr::map(ha_col_names.factor, ~ (scales::hue_pal())(length(.x))) %>%
      purrr::map2(ha_col_names.factor, purrr::set_names)
  }
  groups.use.numeric <- groups.use[sapply(groups.use, is.numeric)]
  ha_cols.numeric <- NULL
  if (length(groups.use.numeric) > 0) {
    ha_col_names.numeric <- names(groups.use.numeric)
    ha_col_hues.numeric <- (scales::hue_pal())(length(ha_col_names.numeric))
    ha_cols.numeric <- purrr::map2(
      groups.use[ha_col_names.numeric],
      ha_col_hues.numeric, numeric_col_fun
    )
  }
  ha_cols <- c(ha_cols.factor, ha_cols.numeric)
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    df = groups.use,
    height = grid::unit(group.bar.height, "points"), col = ha_cols
  )
  hm <- ComplexHeatmap::Heatmap(t(data),
    name = "log expression",
    top_annotation = column_ha, cluster_columns = cluster_columns,
    show_column_names = FALSE, column_dend_height = grid::unit(
      mm_col_dend,
      "mm"
    ), column_split = column_split,
    ...
  )
  return(hm)
}

#' Perform seu gene heatmap operation
#'
#' @param seu Seurat object
#' @param marker_col Color specification
#' @param group.by Character string (default: "ident")
#' @param cells Cell identifiers or information
#' @param layer Character string (default: "scale.data")
#' @param assay Parameter for assay
#' @param group.bar.height Parameter for group.bar.height
#' @param column_split Color specification
#' @param row_split Parameter for row split
#' @param col_arrangement Color specification
#' @param mm_col_dend Color specification
#' @param embedding Character string (default: "pca")
#' @param factor_cols Color specification
#' @param summarize Logical flag (default: FALSE)
#' @param hide_legends Parameter for hide legends
#' @param heatmap_features Parameter for heatmap features
#' @param seurat_assay Character string (default: "SCT")
#' @param features Parameter for features
#' @param cluster_rows Cluster information
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
seu_gene_heatmap <- function(seu, marker_col = "clusters", group.by = "ident", cells = NULL,
                             layer = "scale.data", assay = NULL, group.bar.height = 0.01,
                             column_split = NULL, row_split = NULL, col_arrangement = "ward.D2", mm_col_dend = 30,
                             embedding = "pca", factor_cols = NULL, summarize = FALSE, hide_legends = NULL, heatmap_features = NULL, seurat_assay = "SCT", features = NULL, cluster_rows = FALSE, ...) {
  
  
	
	DefaultAssay(seu) <- "gene"
	
  seu$scna <- factor(seu$scna)
  levels(seu$scna)[1] <- "none"

  if (is.null(heatmap_features)) {
    cluster_order <- levels(seu@meta.data[[marker_col]]) %>%
      set_names(.)

    heatmap_features <-
      seu@misc$markers[[marker_col]][["presto"]] %>%
      # dplyr::filter(Gene.Name %in% rownames(seu@assays[[seurat_assay]]@scale.data)) %>%
      # dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
      dplyr::group_by(Gene.Name) |>
      dplyr::slice_max(order_by = Average.Log.Fold.Change, n = 1) |>
      dplyr::ungroup() |>
      dplyr::arrange(Cluster, Average.Log.Fold.Change) |>
      dplyr::mutate(Cluster = factor(Cluster, levels = cluster_order)) %>%
      dplyr::arrange(Cluster, desc(Average.Log.Fold.Change)) %>%
      group_by(Cluster) %>%
      slice_head(n = 5) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(Gene.Name, .keep_all = TRUE) |>
      identity()
  } else {
    heatmap_features <-
      heatmap_features %>%
      dplyr::arrange(Cluster) %>%
      group_by(Cluster) %>%
      slice_head(n = 6) %>%
      dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
      dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(Gene.Name, .keep_all = TRUE)
  }

  # row_ha = ComplexHeatmap::rowAnnotation(mp = rev(test0$mp))

  if (!is.null(column_split)) {
    column_split <- sort(seu@meta.data[[column_split]])
    column_title <- unique(column_split)
  } else {
    column_title <- NULL
  }

  
  if(is.null(features)){
  	features <- heatmap_features$Gene.Name
  	
  	if (!is.null(row_split)) {
  		row_split <- rev(heatmap_features$Cluster)
  	} 
  }

  myplot <- ggplotify::as.ggplot(seu_complex_heatmap2(seu,
    features = features, group.by = group.by,
    col_arrangement = col_arrangement, cluster_rows = cluster_rows,
    column_split = column_split, column_title = column_title,
    column_title_rot = 90, row_split = row_split,
    row_title_rot = 0, hide_legends = hide_legends, use_raster = TRUE
  ))

  return(myplot)
}

#' Perform seu complex heatmap2 operation
#'
#' @param seu Seurat object
#' @param features Parameter for features
#' @param group.by Character string (default: "ident")
#' @param cells Cell identifiers or information
#' @param layer Character string (default: "scale.data")
#' @param assay Parameter for assay
#' @param group.bar.height Parameter for group.bar.height
#' @param column_split Color specification
#' @param col_arrangement Color specification
#' @param mm_col_dend Color specification
#' @param embedding Character string (default: "pca")
#' @param hide_legends Parameter for hide legends
#' @param ... Additional arguments passed to other functions
#' @return Function result
#' @export
seu_complex_heatmap2 <- function(seu, features = NULL, group.by = "ident", cells = NULL,
                                 layer = "scale.data", assay = NULL, group.bar.height = 0.01,
                                 column_split = NULL, col_arrangement = "ward.D2", mm_col_dend = 30,
                                 embedding = "pca", hide_legends = NULL, ...) {
  if (length(GetAssayData(seu, layer = "scale.data")) == 0) {
    message("seurat object has not been scaled. Please run `Seurat::ScaleData` to view a scaled heatmap; showing unscaled expression data")
    layer <- "data"
  }
  cells <- cells %||% colnames(x = seu)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = seu)[cells]
  }
  assay <- assay %||% Seurat::DefaultAssay(object = seu)
  Seurat::DefaultAssay(object = seu) <- assay
  features <- features %||% VariableFeatures(object = seu)
  features <- rev(x = unique(x = features))
  possible.features <- rownames(x = GetAssayData(
    object = seu,
    layer = layer
  ))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop(
        "No requested features found in the ", layer,
        " layer for the ", assay, " assay."
      )
    }
    warning(
      "The following features were omitted as they were not found in the ",
      layer, " layer for the ", assay, " assay: ", paste(bad.features,
        collapse = ", "
      )
    )
  }
  data <- as.data.frame(x = t(x = as.matrix(x = GetAssayData(
    object = seu,
    layer = layer
  )[features, cells, drop = FALSE])))
  seu <- suppressMessages(expr = StashIdent(
    object = seu,
    save.name = "ident"
  ))
  if (any(col_arrangement %in% c(
    "ward.D", "single", "complete",
    "average", "mcquitty", "median", "centroid", "ward.D2"
  ))) {
    if ("pca" %in% Seurat::Reductions(seu)) {
      cluster_columns <- Seurat::Embeddings(seu, embedding) %>%
        dist() %>%
        hclust(col_arrangement)
    } else {
      message(glue("{embedding} not computed for this dataset; cells will be clustered by displayed features"))
      cluster_columns <- function(m) {
        as.dendrogram(cluster::agnes(m),
          method = col_arrangement
        )
      }
    }
  } else {
    cells <- seu %>%
      Seurat::FetchData(vars = col_arrangement) %>%
      sample_frac() %>%
      dplyr::arrange(across(all_of(col_arrangement))) %>%
      rownames()
    data <- data[cells, ]
    group.by <- base::union(group.by, col_arrangement)
    cluster_columns <- FALSE
  }
  group.by <- group.by %||% "ident"
  groups.use <- seu[[group.by]][cells, , drop = FALSE]
  groups.use <- groups.use %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::mutate(across(where(is.character), ~ str_wrap(str_replace_all(
      .x,
      ",", " "
    ), 10))) %>%
    dplyr::mutate(across(
      where(is.character),
      as.factor
    )) %>%
    data.frame(row.names = 1) %>%
    identity()
  groups.use.factor <- groups.use[sapply(groups.use, is.factor)]
  ha_cols.factor <- NULL
  if (length(groups.use.factor) > 0) {
    ha_col_names.factor <- lapply(groups.use.factor, levels)
    ha_cols.factor <- purrr::map(ha_col_names.factor, ~ (scales::hue_pal())(length(.x))) %>%
      purrr::map2(ha_col_names.factor, purrr::set_names)
  }
  groups.use.numeric <- groups.use[sapply(groups.use, is.numeric)]
  ha_cols.numeric <- NULL
  if (length(groups.use.numeric) > 0) {
    numeric_col_fun <- function(myvec, color) {
    	# 
      circlize::colorRamp2(c(min(myvec), quantile(myvec, probs = seq(0,1, 0.1))[10]), c("white", color))
    }
    ha_col_names.numeric <- names(groups.use.numeric)
    ha_col_hues.numeric <- (scales::hue_pal())(length(ha_col_names.numeric))
    ha_cols.numeric <- purrr::map2(
      groups.use[ha_col_names.numeric],
      ha_col_hues.numeric, numeric_col_fun
    )
  }
  ha_cols <- c(ha_cols.factor, ha_cols.numeric)

  if (!is.null(hide_legends)) {
    show_legend <- !colnames(groups.use) %in% hide_legends
  } else {
    show_legend <- TRUE
  }

  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    df = groups.use,
    height = grid::unit(group.bar.height, "points"), col = ha_cols, show_legend = show_legend
  )


  q <- quantile(abs(t(data)), 0.99)
  col_fun <- circlize::colorRamp2(c(-q, 0, q), c("blue", "white", "red"))

  hm <- ComplexHeatmap::Heatmap(t(data),
    name = "log expression",
    top_annotation = column_ha, cluster_columns = cluster_columns,
    show_column_names = FALSE, column_dend_height = grid::unit(
      mm_col_dend,
      "mm"
    ), column_split = column_split,
    col = col_fun,
    ...
  )
  return(hm)
}