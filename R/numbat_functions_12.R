# Numbat Functions (12)

#' Convert data from one format to another
#'
#' @param seu_v5 Parameter for seu v5
#' @return Modified Seurat object
#' @export
convert_v5_to_v3 <- function(seu_v5) {
  meta <- seu_v5@meta.data

  seu_v3 <- CreateSeuratObject(counts = seu_v5$gene@counts, data = seu_v5$gene@data, assay = "gene", meta.data = meta)

  # transcript_assay.v5 <- CreateAssay5Object(counts = seu_v3$transcript@counts, data = seu_v3$transcript@data)
  # seu_v5$transcript <- transcript_assay.v5

  seu_v3$gene <- seurat_preprocess(seu_v5$gene, normalize = FALSE)

  # seu_v5 <- clustering_workflow(seu_v5)
  seu_v3@reductions <- seu_v5@reductions
  seu_v3@graphs <- seu_v5@graphs
  seu_v3@neighbors <- seu_v5@neighbors

  seu_v3@misc <- seu_v5@misc

  Idents(seu_v3) <- Idents(seu_v5)

  return(seu_v3)
}

#' Perform seu factor heatmap operation
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
#' @param factor_cols Color specification
#' @param ... Additional arguments passed to other functions
#' @return Function result
#' @export
seu_factor_heatmap <- function(seu, features = NULL, group.by = "ident", cells = NULL,
                               layer = "scale.data", assay = NULL, group.bar.height = 0.01,
                               column_split = NULL, col_arrangement = "ward.D2", mm_col_dend = 30,
                               embedding = "pca", factor_cols = NULL, ...) {
  #
  if (!is.null(factor_cols)) {
    data <- seu@meta.data[, factor_cols]
  } else {
    data <- seu@meta.data[, str_detect(colnames(seu@meta.data), pattern = "^k[0-9]*_[0-9]*$")]
  }

  if (any(col_arrangement %in% c(
    "ward.D", "single", "complete",
    "average", "mcquitty", "median", "centroid", "ward.D2"
  ))) {
    
#' Perform clustering analysis
#'
#' @param m Parameter for m
#' @return Function result
#' @export
cluster_columns <- function(m) {
      as.dendrogram(cluster::agnes(m),
        method = col_arrangement
      )
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
  # groups.use <- seu[[group.by]]
  groups.use <- groups.use %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::mutate(across(where(is.character), ~ str_wrap(str_replace_all(.x, ",", " "), 10))) %>%
    dplyr::mutate(across(where(is.character), as.factor)) %>%
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
    name = "normalized usage",
    top_annotation = column_ha, cluster_columns = cluster_columns,
    cluster_rows = FALSE,
    show_column_names = FALSE, column_dend_height = grid::unit(
      mm_col_dend,
      "mm"
    ), column_split = sort(seu@meta.data$clusters), column_title = NULL,
    ...
  )
  return(hm)
}
#' Add annotations to data
#'
#' @param seu Seurat object
#' @param organism Character string (default: "human")
#' @param ... Additional arguments passed to other functions
#' @return Function result
#' @export
annotate_cell_cycle_without_1q <- function(seu, organism = "human", ...) {
  cc_genes_by_arm <- find_cc_genes_by_arm() |>
    dplyr::distinct(symbol, .keep_all = TRUE) |>
    group_by(phase_of_gene) |>
    identity()

  cc_wo_1q <-
    cc_genes_by_arm |>
    dplyr::filter(!(seqnames == "01" & arm == "q")) |>
    identity()

  cc_wo_1q <- split(cc_wo_1q, cc_wo_1q$phase_of_gene) |>
    map(pull, symbol)

  s_genes <- cc_wo_1q$s.genes
  g2m_genes <- cc_wo_1q$g2m.genes
  if (organism == "mouse") {
    s_genes <- dplyr::filter(human_to_mouse_homologs, HGNC.symbol %in%
      s_genes) %>% dplyr::pull(MGI.symbol)
    g2m_genes <- dplyr::filter(human_to_mouse_homologs, HGNC.symbol %in%
      g2m_genes) %>% dplyr::pull(MGI.symbol)
  }
  seu <- CellCycleScoring(seu,
    s.features = s_genes, g2m.features = g2m_genes,
    set.ident = FALSE
  )
}

