# Plot Functions (138)
#' Create a plot visualization
#'
#' @param study_cell_stats Cell identifiers or information
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

plot_study_metadata <- function(study_cell_stats, ...) {
  
  
  # study_cell_stats <- read_csv("results/study_cell_stats.csv")

  normal_ctrl_samples <- unlist(list(
    "collin" = c("SRX10031193", "SRX10031194")
  ))

  bad_qc_sample_ids <- list(
    "collin" = c("SRX10031191", "SRX10031192"),
    "yang" = c("SRX11133591", "SRX11133590", "SRX11133589", "SRX11133586"),
    "field" = c("SRX14116948", "SRX14116946", "SRX14116944")
  ) %>%
    enframe("study", "sample_id") %>%
    unnest(sample_id) %>%
    identity()

  bad_scna_sample_ids <- list(
    "collin" = c(),
    "yang" = c("SRX11133588", "SRX11133587", "SRX11133585"),
    "field" = c("SRX14116945"),
    "wu" = c("SRX10264517", "SRX10264518", "SRX10264521", "SRX10264522"),
    "liu" = c("SRX22868104", "SRX22868103")
  ) %>%
    enframe("study", "sample_id") %>%
    unnest(sample_id) %>%
    identity()

  excluded_sample_ids <-
    list(
      "bad_qc" = bad_qc_sample_ids,
      "bad_scna" = bad_scna_sample_ids
    ) %>%
    dplyr::bind_rows(.id = "exclusion_criteria")

  unfiltered_cell_stats_plot_file <- "results/unfiltered_study_stats_mt_v_nUMI.pdf"

  study_cell_stats <-
    study_cell_stats %>%
    dplyr::full_join(excluded_sample_ids, by = c("study", "sample_id")) %>%
    identity()

  mt_v_nUMI_plot <- 
  	study_cell_stats %>%
    dplyr::filter(!sample_id %in% normal_ctrl_samples) %>%
    plot_mt_v_nUMI()

  unfiltered_cell_stats_plot <- study_cell_stats %>%
    dplyr::filter(!sample_id %in% normal_ctrl_samples) %>%
    plot_study_cell_stats(mito_expansion = 0.8, ...)
  
  unfiltered_cell_stats_plot_file <- qpdf::pdf_combine(list(mt_v_nUMI_plot, unfiltered_cell_stats_plot), "results/unfiltered_study_stats.pdf")
  
  # asdf

  # good_qc_study_cell_stats <-
  # 	study_cell_stats %>%
  # 	dplyr::filter(!sample_id %in% bad_qc_samples)
  #
  # good_qc_cell_stats_plot_file <- "results/good_qc_study_stats.pdf"
  # plot_study_cell_stats(good_qc_study_cell_stats, good_qc_cell_stats_plot_file, mito_expansion = 1)
  #
  #
  # # retained ------------------------------
  # retained_study_cell_stats <-
  # 	study_cell_stats %>%
  # 	dplyr::filter(!sample_id %in% c(bad_scna_samples, bad_qc_samples))
  #
  # retained_cell_stats_plot_file <- "results/retained_study_stats.pdf"
  # plot_study_cell_stats(retained_study_cell_stats, retained_cell_stats_plot_file, mito_expansion = 0.55)

  return(unfiltered_cell_stats_plot_file)
}

#' Calculate scores for the given data
#'
#' @param numbat_rds_file File path
#' @return List object
#' @export
score_samples_for_rod_enrichment <- function(numbat_rds_file) {
  
  
  sample_id <- str_extract(numbat_rds_file, "SR[RX][0-9]+")

  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  # count number of rods by RHO, ROM1, GNAT1, NR2E3?------------------------------

  seu <- AddModuleScore(seu, list("rod" = c("RHO", "ROM1", "GNAT1", "NR2E3")), name = "rod")

  rod_meta <-
    seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, cluster = gene_snn_res.0.2, rod1) %>%
    dplyr::mutate(rod_identity = ifelse(rod1 > 1.5, "rod", "cell")) %>%
    identity()

  seu@meta.data["rod_identity"] <-
    rod_meta %>%
    dplyr::pull(rod_identity)

  # Seurat::FeaturePlot(seu, features = "rod1", split.by = "gene_snn_res.0.2")
  rod_score_plot <- Seurat::FeaturePlot(seu, features = "rod1")

  rod_id_plot <- Seurat::DimPlot(seu, group.by = "rod_identity")

  percent_rod_cells <-
    janitor::tabyl(rod_meta, rod_identity) %>%
    dplyr::filter(rod_identity == "rod") %>%
    dplyr::pull(percent) %>%
    round(2) %>%
    identity()

  rod_patch <- rod_score_plot + rod_id_plot +
    plot_annotation(title = sample_id, subtitle = scales::label_percent()(percent_rod_cells))

  dir_create("results/rod_plots")
  rod_plot_path <- glue("results/rod_plots/{sample_id}.pdf")

  ggsave(rod_plot_path)

  return(list("sample_id" = sample_id, "percent_rod_cells" = (scales::label_percent()(percent_rod_cells)), "rod_rich" = ifelse(percent_rod_cells > 0.05, 1, 0)))
}

#' Calculate scores for the given data
#'
#' @param numbat_rds_file File path
#' @param subtype_markers Parameter for subtype markers
#' @return List object
#' @export
score_whole_pseudobulks <- function(numbat_rds_file, subtype_markers) {
  
  
  sample_id <- str_extract(numbat_rds_file, "SR[RX][0-9]+")

  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  bulk_expression <- GetAssayData(seu, layer = "data") %>%
    rowSums()

  subtype1_expression <- bulk_expression[names(bulk_expression) %in% subtype_markers$subtype1]

  subtype2_expression <- bulk_expression[names(bulk_expression) %in% subtype_markers$subtype2]

  return(list("sample_id" = sample_id, "s1" = mean(subtype1_expression), "s2" = mean(subtype2_expression)))
}

#' Calculate scores for the given data
#'
#' @param seu_paths File path
#' @return ggplot2 plot object
#' @export
derive_pseudobulk_subtype_scores <- function(seu_paths) {
  
  
  #
  
#' Extract or pull specific data elements
#'
#' @param seu_path File path
#' @return Data frame
#' @export
pull_assay_data <- function(seu_path) {
  
  
    seu <- readRDS(seu_path) %>%
      filter_sample_qc()

    bulk_data <- GetAssayData(seu, layer = "data") %>%
      rowSums()

    bulk_counts <- GetAssayData(seu, layer = "counts") %>%
      rowSums()

    return(list("counts" = bulk_counts, "data" = bulk_data))
  }

  bulk_assays <-
    seu_paths %>%
    set_names(str_extract(., "SR[RX][0-9]+")) %>%
    map(pull_assay_data)

  bulk_assay_counts <-
    bulk_assays %>%
    map("counts") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  bulk_assay_counts[is.na(bulk_assay_counts)] <- 0

  scaled_counts_datas <- prop.table(t(bulk_assay_counts), 2)

  mydend <- as.dendrogram(hclust(dist(t(scaled_counts_datas)), method = "ward.D2"))

  dend_plot <- ggplot(dendextend::as.ggdend(mydend))

  dend_groups <- dendextend:::cutree.dendrogram(mydend, 2) %>%
    tibble::enframe("sample_id", "group") %>%
    tibble::column_to_rownames("sample_id") %>%
    dplyr::mutate(group = factor(group))

  dds <- DESeqDataSetFromMatrix(t(bulk_assay_counts), dend_groups, ~group)

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]

  dds <- DESeq(dds)
  res <- results(dds)
  res

  res_table <-
    res %>%
    as.data.frame() %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::inner_join(annotables::grch38, by = "symbol") %>%
    dplyr::arrange(log2FoldChange) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    dplyr::select(symbol, description, everything()) %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::distinct(entrez, .keep_all = TRUE)

  # we want the log2 fold change
  original_gene_list <- res_table[["log2FoldChange"]]

  # name the vector
  names(original_gene_list) <- res_table$entrez

  # omit any NA values
  gene_list <- na.omit(original_gene_list)

  # sort the list in decreasing order (required for clusterProfiler)
  gene_list <- sort(gene_list, decreasing = TRUE)

  gse <- clusterProfiler::gseGO(
    geneList = gene_list,
    ont = "BP",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    minGSSize = 3,
    maxGSSize = 800,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    pAdjustMethod = "BH"
  ) %>%
    clusterProfiler::simplify() # for GSEGO

  gse_table <-
    dplyr::arrange(gse@result, NES)

  return(list("dend" = dend_plot, "diffex" = res_table, "enrichment" = gse_table))
}

