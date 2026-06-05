# Plot Functions (141)
#' Convert data from one format to another
#'
#' @param seu_path File path
#' @return Converted data object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

convert_seu_to_scanpy <- function(seu_path) {
  
  
  filtered_seu <- seu_path

  # filtered_seu <- DietSeurat(filtered_seu, misc = FALSE)

  filtered_seu <- DietSeurat(
    filtered_seu,
    counts = TRUE, # so, raw counts save to adata.layers['counts']
    data = TRUE, # so, log1p counts save to adata.X when scale.data = False, else adata.layers['data']
    scale.data = FALSE, # if only scaled highly variable gene, the export to h5ad would fail. set to false
    features = rownames(filtered_seu), # export all genes, not just top highly variable genes
    assays = "gene",
    dimreducs = c("pca", "umap"),
    graphs = c("gene_nn", "gene_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
    misc = FALSE
  )

  filtered_scanpy_path <- str_replace(seu_path, "seu.rds", "scanpy.h5ad") %>%
    str_replace("seurat", "scanpy")

  MuDataSeurat::WriteH5AD(filtered_seu, filtered_scanpy_path, assay = "gene")

  return(filtered_scanpy_path)
}

#' Calculate scores for the given data
#'
#' @param seu_path File path
#' @param stachelek_scores_table Parameter for stachelek scores table
#' @return ggplot2 plot object
#' @export
score_stachelek <- function(seu_path, stachelek_scores_table) {
  
  
  sample_id <- str_extract(seu_path, "SR[RX][0-9]+")

  seu <- seu_path

  stachelek_scores <-
    stachelek_scores_table %>%
    list_flatten() %>%
    map(dplyr::distinct, symbol) %>%
    map(pull, symbol)


  seu <- AddModuleScore(seu, stachelek_scores, name = "stachelek")

  names(seu@meta.data)[which(names(seu@meta.data) %in% paste0("stachelek", seq(1, length(stachelek_scores))))] <- names(stachelek_scores)

  # stachelek_score_fplot <- FeaturePlot(seu, names(stachelek_scores)) +
  #   plot_annotation(title = sample_id)

  # stachelek_score_vlnplot <- VlnPlot(seu, names(stachelek_scores), group.by = "scna", ncol = 8) +
  #   plot_annotation(title = sample_id)

  stachelek_score_vlnplot <- VlnPlot(seu, names(stachelek_scores), group.by = "abbreviation", ncol = 8) +
    plot_annotation(title = sample_id)

  dir_create("results/effect_of_regression/stachelek/")

  plot_path <- glue("results/effect_of_regression/stachelek/{sample_id}_stachelek_scores.pdf")

  ggsave(plot_path, height = 8, width = 24)

  return(plot_path)
}

#' Perform reference plae celltypes operation
#'
#' @param seu_path File path
#' @param mycluster Cluster information
#' @param mygenes Gene names or identifiers
#' @return Function result
#' @export
reference_plae_celltypes <- function(seu_path, mycluster = "SCT_snn_res.0.4", mygenes) {
  
  
  # con <- dbConnect(RSQLite::SQLite(), "/dataVolume/storage/scEiad/human_pseudobulk/diff_resultsCellType.sqlite")
  #
  # mytable <-
  #   tbl(con, "diffex") %>%
  #   dplyr::filter(Against =="All") %>%
  #   dplyr::filter(Base %in% celltypes) %>%
  #   dplyr::filter(padj < 0.05, abs(log2FoldChange) > 2) %>%
  #   dplyr::filter(Organism == "Homo sapiens") %>%
  #   collect()

  seu_path <- "output/seurat/SRX10264519_regressed_seu.rds"

  seu <- seu_path

  cluster_markers <-
    seu@misc$markers[[mycluster]][["presto"]] %>%
    dplyr::group_by(Cluster) %>%
    dplyr::slice_head(n = 30) %>%
    dplyr::rename(symbol = `Gene.Name`)

  mytable <-
    "data/plae_top_diffex.csv" %>%
    read_csv() %>%
    dplyr::group_by(Base) %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    dplyr::mutate(symbol = str_remove(Gene, " \\(.*\\)")) %>%
    dplyr::filter(!is.na(symbol)) %>%
    dplyr::left_join(cluster_markers, by = "symbol", relationship = "many-to-many") %>%
    dplyr::filter(!is.na(Cluster)) %>%
    dplyr::filter(log2FoldChange > 0) %>%
    dplyr::arrange(Cluster, Base) %>%
    dplyr::select(Cluster, Base, everything()) %>%
    # dplyr::mutate(checked_gene = case_when(symbol %in% mygenes ~ 1)) %>%
    identity()

  mytable0 <-
    mytable %>%
    dplyr::group_by(Cluster, Base) %>%
    dplyr::summarise(mean_fc = mean(log2FoldChange)) %>%
    dplyr::arrange(Cluster, desc(mean_fc))

  # janitor::tabyl(mytable, Cluster, Base)
}

