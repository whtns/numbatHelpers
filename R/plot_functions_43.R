# Plot Functions (142)

#' Perform generate plae ref operation
#'
#' @param plae_seu_path File path
#' @return Function result
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

generate_plae_ref <- function(plae_seu_path = "data/plae_human_fetal_seu.rds") {
  
  
  plae_human_fetal_seu <- readRDS(plae_seu_path)

  plae_ref <- AggregateExpression(plae_human_fetal_seu, group.by = "CellType_predict") %>%
    pluck("RNA")

  return(plae_ref)
}

#' Create a plot visualization
#'
#' @param seu_path File path
#' @param sample_id Parameter for sample id
#' @param plae_ref Parameter for plae ref
#' @param group.by Character string (default: "SCT_snn_res.0.4")
#' @param query_genes Gene names or identifiers
#' @return ggplot2 plot object
#' @export
plot_celltype_predictions <- function(seu_path, sample_id, plae_ref = NULL, group.by = "SCT_snn_res.0.4", query_genes = NULL) {
  
  
  #

  seu <- seu_path

  celltypes <- c(
    "Amacrine Cells", "Bipolar Cells", "Cones",
    "Early RPCs", "Horizontal Cells", "Late RPCs",
    "Neurogenic Cells", "Photoreceptor Precursors",
    "Retinal Ganglion Cells", "Rods", "RPCs", "RPE"
  )

  if (is.null(query_genes)) {
    query_genes <- VariableFeatures(seu)
  }

  for (assay_name in SeuratObject::Assays(seu)) {
    if (inherits(seu[[assay_name]], "Assay5"))
      seu[[assay_name]] <- SeuratObject::JoinLayers(seu[[assay_name]])
  }
  if (inherits(seu[[DefaultAssay(seu)]], "Assay5")) {
    seu_mat <- GetAssayData(seu, layer = "data")
  } else {
    seu_mat <- GetAssayData(seu, layer = "data")
  }

  sub_annotable <-
    annotables::grch38 %>%
    dplyr::filter(symbol %in% query_genes)

  if (is.null(plae_ref)) {
    plae_ref <-
      "data/plae_pseudobulk_counts.csv" %>%
      read_csv() %>%
      dplyr::inner_join(sub_annotable, by = c("Gene" = "ensgene"), relationship = "many-to-many") %>%
      dplyr::distinct(study, type, symbol, .keep_all = TRUE) %>%
      dplyr::filter(type %in% str_remove(celltypes, "s$")) %>%
      dplyr::group_by(symbol, type) %>%
      dplyr::summarize(total_counts = mean(counts)) %>%
      tidyr::pivot_wider(names_from = "type", values_from = "total_counts") %>%
      tibble::column_to_rownames("symbol") %>%
      as.matrix() %>%
      identity()
  } else {
    plae_ref <- plae_ref[rownames(plae_ref) %in% query_genes, colnames(plae_ref) %in% celltypes]
  }

  res <- clustify(
    input = seu_mat,
    metadata = seu[[group.by]][[1]],
    ref_mat = plae_ref,
    query_genes = query_genes
  )

  cor_to_call(res)

  res2 <- cor_to_call(
    cor_mat = res, # matrix correlation coefficients
    cluster_col = group.by # name of column in meta.data containing cell clusters
  )

  res3 <- cor_to_call_rank(
    cor_mat = res, # matrix correlation coefficients
    cluster_col = group.by # name of column in meta.data containing cell clusters
  ) %>%
    dplyr::mutate(cluster = .data[[group.by]]) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(rank)

  dir_create("results/clustify")
  table_path <- glue("results/clustify/{sample_id}_clustifyr.csv")
  write_csv(res3, table_path)


  # Insert into original metadata as "type" column
  seu@meta.data <- call_to_metadata(
    res = res2, # data.frame of called cell type for each cluster
    metadata = seu@meta.data, # original meta.data table containing cell clusters
    cluster_col = group.by # name of column in meta.data containing cell clusters
  )

  neurogenic_table <- seu@meta.data %>%
    dplyr::mutate(cluster = .data[[group.by]]) %>%
    janitor::tabyl(cluster, type) %>%
    dplyr::mutate(sample = sample_id) %>%
    # dplyr::mutate(percent_neurogenic = `Neurogenic Cells`/(Cones+`Neurogenic Cells`)) %>%
    identity()

  allcells_dimplot <- DimPlot(
    seu,
    group.by = "type",
    split.by = group.by
  ) +
    plot_annotation(title = sample_id) +
    NULL

  dir_create("results/clustify")
  plot_path <- glue("results/clustify/{sample_id}_clustifyr.pdf")
  ggsave(plot_path)

  return(list("plot" = plot_path, "table" = neurogenic_table, "seu" = seu))
}

#' Load or read data from file
#'
#' @param zinovyev_file File path
#' @return Loaded data object
#' @export
read_zinovyev_genes <- function(zinovyev_file = "data/zinovyev_cc_genes.tsv") {
  
  
  zinovyev_cc_genes <-
    read_tsv(zinovyev_file) %>%
    dplyr::group_by(term) %>%
    tidyr::nest(data = symbol) %>%
    # split(.$term) %>%
    tibble::deframe() %>%
    map(tibble::deframe) %>%
    identity()
}

#' Load or read data from file
#'
#' @param cc_file File path
#' @return Loaded data object
#' @export
read_giotti_genes <- function(cc_file = "data/giotti_cc_genes.tsv") {
  
  
  giotti_cc_genes <-
    read_tsv(cc_file) %>%
    dplyr::filter(!(term %in% c("Function known but link to cell division not well established", "Uncharacterized"))) %>%
    dplyr::mutate(term = str_wrap(term, width = 10)) %>%
    identity()
}

#' Create a plot visualization
#'
#' @param seu_path File path
#' @param seu_name Parameter for seu name
#' @return ggplot2 plot object
#' @export
plot_phase_distribution_by_scna <- function(seu_path, seu_name) {


    seu <- seu_path
    seu$Phase <- factor(seu$Phase, levels = c("G1", "S", "G2M"))
    plot_distribution_of_clones_across_clusters(seu, seu_name, var_x = "scna", var_y = "Phase", both_ways = FALSE)
  }

#' Plot UMAP and feature plots for a diploid Seurat object and save to PDF
#'
#' @param diploid_seu Path to diploid Seurat RDS file
#' @param celltype_markers Named list of celltype -> gene vectors
#' @param out_path Output PDF path
#' @param use_integrated_clusters Logical (default: FALSE). If TRUE, plot integrated_snn_res.0.2 instead of gene_snn_res.0.2
#' @param include_batch_series Logical (default: TRUE). If TRUE, create serial batch plots where each batch is colored and others grayed
#' @return Path to the saved PDF
#' @export
plot_diploid_seu_umaps <- function(
    diploid_seu,
    celltype_markers,
    out_path = NULL,
    use_integrated_clusters = FALSE,
    include_batch_series = TRUE) {
  if (is.null(out_path)) {
    stem <- tools::file_path_sans_ext(basename(diploid_seu))
    out_path <- file.path("results", paste0(stem, "_umap_plots.pdf"))
  }
  seu <- readRDS(diploid_seu)
  cluster_var <- if (use_integrated_clusters) "integrated_snn_res.0.2" else "gene_snn_res.0.2"
  meta_vars <- c(cluster_var, "Phase", "type", "clone_opt", "batch")
  meta_plots <- purrr::map(meta_vars, function(var) {
    if (!var %in% colnames(seu@meta.data)) return(NULL)
    Seurat::DimPlot(seu, group.by = var, label = TRUE, repel = TRUE) +
      ggplot2::ggtitle(var) +
      ggplot2::theme(legend.position = "bottom")
  }) |> purrr::compact()

  # Serial batch plots: each batch colored, all others grayed
  batch_series_plots <- list()
  if (include_batch_series && "batch" %in% colnames(seu@meta.data)) {
    unique_batches <- unique(seu$batch)
    # Get colors matching Seurat's default DimPlot color scheme
    batch_colors <- scales::hue_pal()(length(unique_batches))
    names(batch_colors) <- sort(as.character(unique_batches))
    
    batch_series_plots <- purrr::map(unique_batches, function(batch_id) {
      # Create a copy to avoid modifying original seu
      seu_temp <- seu
      # Create display column: target batch name or "Other"
      seu_temp@meta.data$batch_display <- factor(
        ifelse(seu_temp$batch == batch_id, as.character(batch_id), "Other"),
        levels = c(as.character(batch_id), "Other")
      )
      
      # Get the color for this batch from the default palette
      target_color <- batch_colors[as.character(batch_id)]
      
      # Create color vector with proper naming
      color_vals <- c(target_color, "lightgray")
      names(color_vals) <- c(as.character(batch_id), "Other")
      
      p <- Seurat::DimPlot(seu_temp, group.by = "batch_display", label = FALSE, repel = FALSE) +
        ggplot2::scale_color_manual(
          values = color_vals,
          labels = c(as.character(batch_id), "Other batches")
        ) +
        ggplot2::ggtitle(paste("Batch:", batch_id)) +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::guides(color = ggplot2::guide_legend(title = "Batch"))
      
      p
    })
  }

  marker_celltypes <- c("cones", "rods", "microglia", "bipolar cells")
  feature_plots <- purrr::imap(
    celltype_markers[marker_celltypes],
    function(genes, celltype) {
      genes <- genes[genes %in% rownames(seu)]
      purrr::map(genes, function(gene) {
        Seurat::FeaturePlot(seu, features = gene) +
          ggplot2::ggtitle(sprintf("%s (%s)", gene, celltype))
      })
    }
  ) |> purrr::flatten()

  all_plots <- c(meta_plots, batch_series_plots, feature_plots)
  pdf(out_path, width = 10, height = 8)
  purrr::walk(all_plots, print)
  dev.off()
  out_path
}

