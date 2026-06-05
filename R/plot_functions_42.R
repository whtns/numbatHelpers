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

#' Score chromosomal instability using CIN gene signatures
#'
#' @param seu_path File path to Seurat RDS, or list containing a path
#' @return File path of output PDF
#' @export
score_chrom_instability <- function(seu_path) {
  if (is.list(seu_path)) seu_path <- unlist(seu_path, use.names = FALSE)[1]
  sample_id <- stringr::str_extract(seu_path, "SR[RX][0-9]+")
  seu <- readRDS(seu_path)

  cin_scores <- list(
    "cin70" = c("TPX2", "PRC1", "FOXM1", "CDC2", "TGIF2", "MCM2", "H2AFZ", "TOP2A", "PCNA", "UBE2C", "MELK", "TRIP13", "CNAP1", "MCM7", "RNASEH2A", "RAD51AP1", "KIF20A", "CDC45L", "MAD2L1", "ESPL1", "CCNB2", "FEN1", "TTK", "CCT5", "RFC4", "ATAD2", "ch-TOG", "NUP205", "CDC20", "CKS2", "RRM2", "ELAVL1", "CCNB1", "RRM1", "AURKB", "MSH6", "EZH2", "CTPS", "DKC1", "OIP5", "CDCA8", "PTTG1", "CEP55", "H2AFX", "CMAS", "BRRN1", "MCM10", "LSM4", "MTB", "ASF1B", "ZWINT", "TOPK", "FLJ10036", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2", "RAD21", "ACTL6A", "GPIandMGC13096", "SFRS2", "HDGF", "NXT1", "NEK2", "DHCR7", "STK6", "NDUFAB1", "KIAA0286", "KIF4A"),
    "pos_tri70" = c("ANXA7", "ATG7", "BDNF", "CDKN2B", "CHFR", "CNDP2", "ENO3", "F3", "GIPC2", "GJB5", "HSPB7", "IMPACT", "P2RY14", "PCDH7", "PKD1", "PLCG2", "SNCA", "SNCG", "TMEM140", "TMEM40"),
    "neg_tri70" = c("AURKA", "BCL11B", "BIRC5", "BLMH", "BRD8", "BUB1B", "CCNA2", "CDC5L", "CDK1", "CENPE", "CENPN", "CTH", "DLGAP5", "HMGB2", "IDH2", "ISOC1", "KIAA0101", "KIF22", "LIG1", "LSM2", "MCM2", "MCM5", "MCM7", "MYBL2", "NASP", "NCAPD2", "NCAPH", "NMI", "NUDT21", "PCNA", "PIGO", "PLK1", "PLK4", "POLD2", "RACGAP1", "RAD51", "RFC2", "RFC3", "RFC5", "RPA2", "SMAD4", "SMC4", "SSRP1", "TAB2", "TCOF1", "TIMELESS", "TIPIN", "TOP2A", "UBE2C", "USP1"),
    "het70" = c("AHCYL1", "AKT3", "ANO10", "ANTXR1", "ATP6V0E1", "ATXN1", "B4GALT2", "BASP1", "BHLHE40", "BLVRA", "CALU", "CAP1", "CAST", "CAV1", "CLIC4", "CTSL1", "CYB5R3", "ELOVL1", "EMP3", "FKBP14", "FN1", "FST", "GNA12", "GOLT1B", "HECTD3", "HEG1", "HOMER3", "IGFBP3", "IL6ST", "ITCH", "LEPRE1", "LEPREL1", "LEPROT", "LGALS1", "LIMA1", "LPP", "MED8", "MMP2", "MUL1", "MYO10", "NAGK", "NR1D2", "NRIP3", "P4HA2", "PKIG", "PLOD2", "PMP22", "POFUT2", "POMGNT1", "PRKAR2A", "RAGE", "RHOC", "RRAGC", "SEC22B", "SERPINB8", "SPAG9", "SQSTM1", "TIMP2", "TMEM111", "TRIM16", "TRIO", "TUBB2A", "VEGFC", "VIM", "WASL", "YIPF5", "YKT6", "ZBTB38", "ZCCHC24", "ZMPSTE24"),
    "bucc_up" = c("HPDL", "L1CAM", "SMPD4", "ELF4", "IRAK1", "ISG20L2", "ZNF512B", "TMEM164", "RCE1", "AMPD2", "SCXB", "HSF1", "PLXNA3", "RBP4", "MRPL11", "NSUN5", "PNMA5", "FADS3", "TRAIP", "IGSF9", "TFAP4", "BOP1", "FAM131A", "EIF1AD", "FAM122C", "SLC10A3", "CCNK", "ELK1", "DGAT1", "DHCR7", "GNAZ", "SCAMP3", "ZNF7", "CCDC86", "SLC35A2", "CPSF4", "USP39", "FTSJD2", "DBF4B", "MEPCE", "ZBTB2", "SEMA4D", "LRRC14", "YKT6", "METTL7B", "SCAMP5", "IP6K1", "MTCH1", "SLC6A9", "SMCR7L", "SAMD4B", "CD3EAP", "PPP1R16A", "TMCO6", "U2AF2", "POM121", "TSR2", "UNKL", "IFRD2", "CYHR1", "MAGEA9B", "TBC1D25", "ZNF275", "TJAP1", "CLDN4", "SPNS1", "TMEM120A", "CPSF1", "FABP3", "FAM58A"),
    "bucc_down" = c("CYTIP", "HLA-DMB", "FYCO1", "NUBP1", "AGPS", "CADPS2", "C3orf52", "FGL2", "PDIA3P", "IL18", "APOBEC3G", "CA8", "HNRNPA1L2", "CD44", "TMEM2", "ELL2", "FCGBP", "TAPT1", "AHR", "HLCS", "AEN", "AGR3", "STEAP1", "TLR4", "HNRNPA1", "HLA-DRB5", "SLITRK6", "LPXN", "HLA-DRB1", "BACE2", "SYTL1", "KYNU", "RPS27L", "CD96", "HLA-DPA1", "TRA2A", "FBXO25", "HFE", "DRAM1", "DYRK4", "FBXL5", "GLUL", "HLA-DRA", "KCNMA1", "IL18R1", "SIPA1L2", "STEAP2", "TP53I3", "LIMA1", "BMX", "HLA-DMA", "C17orf97", "RPH3AL", "SLC37A1", "PHLDA3", "TK2", "NHS", "TMEM229B", "CYFIP1", "B3GALT5", "FAM107B", "ADH6", "CCND1", "LCK", "RABL2B", "SPA17", "GADD45A", "FHL2", "RPL36AL", "CCDC60", "MUC5B", "PHLDA1", "ZG16B", "SIDT1", "BTD", "TTC9", "GSTZ1", "ITGA6", "C10orf32", "FAM46C", "SPDEF", "ASRGL1", "TCN1", "NUDT2", "CDKN1A", "POLR3GL", "C4BPA", "LAMC2", "SLC6A14")
  )

  seu <- Seurat::AddModuleScore(seu, cin_scores, name = "cin")
  names(seu@meta.data)[which(names(seu@meta.data) %in% paste0("cin", seq(1, length(cin_scores))))] <- names(cin_scores)

  cin_score_fplot <- Seurat::FeaturePlot(seu, names(cin_scores)) +
    patchwork::plot_annotation(title = sample_id)

  fs::dir_create("results/effect_of_regression/cin_scores/")
  plot_path <- glue::glue("results/effect_of_regression/cin_scores/{sample_id}_cin_scores.pdf")
  ggplot2::ggsave(plot_path, height = 8, width = 8)

  return(plot_path)
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

