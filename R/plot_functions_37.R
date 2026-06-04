# Plot Functions (136)
#' Calculate scores for the given data
#'
#' @param seu_path File path
#' @param gene_lists Gene names or identifiers
#' @param group.by Character string (default: "SCT_snn_res.0.4")
#' @param leiden_cluster_file File path
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

score_and_heatmap_seu <- function(seu_path, gene_lists, group.by = "SCT_snn_res.0.4", leiden_cluster_file = NULL) {
  #

  mysample <- str_extract(seu_path, "SR[RX][0-9]+")

  # for(geneset in names(gene_lists)){
  #   # seu <- Seurat::MetaFeature(seu, features = gene_lists[[geneset]], meta.name = geneset)
  # }

  seu <- readRDS(seu_path)

  seu <- Seurat::AddModuleScore(seu, features = gene_lists, name = "subtype")

  module_names <- paste0("subtype", seq(length(gene_lists)))

  seu@meta.data[names(gene_lists)] <- seu@meta.data[module_names]

  cluster_mps <-
    seu@meta.data[, c(names(gene_lists), group.by)] %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(-c(group.by, "cell"), names_to = "mp", values_to = "score") %>%
    group_by(.data[[group.by]], mp) %>%
    dplyr::summarize(score = mean(score)) %>%
    dplyr::mutate({{ group.by }} := as.factor(.data[[group.by]])) %>%
    dplyr::group_by(mp) %>%
    dplyr::filter(any(score > 0.1))

  cluster_heatmap <-
    cluster_mps %>%
    ggplot(aes(x = .data[[group.by]], y = mp, fill = score)) +
    geom_tile() +
    theme_minimal() +
    scale_fill_gradient(low = "red", high = "yellow", na.value = NA)

  clone_heatmap <- seu@meta.data[, c(names(gene_lists), "scna")] %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(-c("scna", "cell"), names_to = "mp", values_to = "score") %>%
    group_by(scna, mp) %>%
    dplyr::summarize(score = mean(score)) %>%
    dplyr::mutate(scna = as.factor(scna)) %>%
    dplyr::group_by(mp) %>%
    dplyr::filter(any(score > 0.1)) %>%
    ggplot(aes(x = scna, y = mp, fill = score)) +
    geom_tile() +
    theme_minimal()

  cluster_heatmap + clone_heatmap
  mp_score_heatmap_file <- glue("results/{mysample}_mp_score_heatmaps.pdf")
  ggsave(mp_score_heatmap_file)

  cluster_mp_genes <- gene_lists[unique(cluster_mps$mp)] %>%
    tibble::enframe("mp", "symbol") %>%
    tidyr::unnest(symbol) %>%
    dplyr::filter(symbol %in% VariableFeatures(seu)) %>%
    dplyr::filter(symbol %in% rownames(seu$gene@scale.data)) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    dplyr::slice_sample(n = 50) %>%
    dplyr::arrange(mp) %>%
    # dplyr::bind_rows(.id = "mp") %>%
    # sample(100) %>%
    identity()

  row_ha <- ComplexHeatmap::rowAnnotation(mp = rev(cluster_mp_genes$mp))

  ggplotify::as.ggplot(seu_complex_heatmap(seu, features = cluster_mp_genes$symbol, group.by = c("SCT_snn_res.0.4", "scna"), col_arrangement = c("SCT_snn_res.0.4", "scna"), cluster_rows = FALSE, right_annotation = row_ha, use_raster = TRUE)) +
    labs(title = mysample)

  cluster_mp_gene_heatmap_file <- glue("results/{mysample}_mp_gene_heatmaps.pdf")
  ggsave(cluster_mp_gene_heatmap_file, height = 8, width = 10)

  return(list("mp_score_heatmap" = mp_score_heatmap_file, "mp_gene_heatmap" = cluster_mp_gene_heatmap_file))
}

#' Extract or pull specific data elements
#'
#' @param supp_excel Character string (default: "data/Liu_Radvanyi_2022_supp_data/41467_2021_25792_MOESM6_ESM.xlsx")
#' @param subtype_chr_preference_file File path
#' @return Extracted data elements
#' @export
pull_subtype_genes <- function(supp_excel = "data/Liu_Radvanyi_2022_supp_data/41467_2021_25792_MOESM6_ESM.xlsx", subtype_chr_preference_file = "results/chromosome_distribution_of_subtype_genes.xlsx") {
  genes_diff_expressed <-
    supp_excel %>%
    readxl::read_excel(sheet = 1, skip = 2) %>%
    janitor::clean_names() %>%
    dplyr::rename(symbol = gene) %>%
    dplyr::filter(gene_cluster %in% c("1.2", "2")) %>%
    dplyr::mutate(gene_cluster = paste0("subtype", str_remove(gene_cluster, "\\..*"))) %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::arrange(chr) %>%
    dplyr::filter(abs(log_fc_subtype_2_vs_subtype_1) > 0.5) %>%
    dplyr::filter(adjusted_p_value < 0.05) %>%
    dplyr::filter(chr %in% c(1:22, "X")) %>%
    identity()

  list(
    "chr_distribution" = janitor::tabyl(genes_diff_expressed, chr, gene_cluster),
    "all_subtype_genes" = genes_diff_expressed
  ) %>%
    writexl::write_xlsx(subtype_chr_preference_file)

  genes_diff_expressed <-
    genes_diff_expressed %>%
    dplyr::select(gene_cluster, symbol) %>%
    split(.$gene_cluster) %>%
    map(tibble::deframe)
}

#' Create a plot visualization
#'
#' @param groblist Parameter for groblist
#' @return ggplot2 plot object
#' @export
rescale_and_clean_plots <- function(groblist) {
    gglist <- map(groblist, ggplotify::as.ggplot)

    gglist <-
      gglist %>%
      map(~ (.x + scale_fill_manual(values = scales::hue_pal()(7)))) %>%
      map(~ (.x +
        theme(
          legend.position = "none",
          axis.title.x = element_blank()
        ) +
        NULL
      )) %>%
      identity()

    return(gglist)
  }

#' Extract or pull specific data elements
#'
#' @return List object
#' @export
pull_stem_cell_markers <- function() {
  smith_markers <- list(
    adult = c(
      "VSNL1", "AKR1B1", "NOTCH4", "TMEM237", "SAMD5",
      "PKD2", "NAP1L1", "PTTG1", "CDK6", "CDCA7", "ACSL4", "HELLS",
      "IKBIP", "PLTP", "TMEM201", "CACHD1", "ILF3", "DNMT1", "USP31",
      "FAM216A", "SLC41A1", "PFKM", "KANK1", "SUPT16H", "ADCY3", "FGD1",
      "PTPN14", "C200rf27", "LGR6", "SLC16A7", "JAM3", "FBL", "NASP",
      "RANBP1", "PRNP", "DSE", "GPX7", "KDELC1", "FCHSD2", "SLCO3A1",
      "CONB11P1", "LOC284023", "NOL9", "NKRE", "NUP107", "RCC2", "ARHGAP25",
      "DDX46", "TCOF1", "GMPS"
    ),
    naive = c(
      "DNMT3L", "ALPPL2", "NLRP7",
      "SLC25A16", "DPPA5", "ATAD3B", "OLAH", "SAMHD1", "PYGB", "TFCP2L1",
      "CP1", "NEFH", "RAB15", "SUSD2", "VSIG10", "PRSS12", "PTPRU",
      "ASRGL1", "A4GALT", "DNAJC15", "CBFA2T2", "KHDC1L", "SLC8B1",
      "SLC35F6", "AARS2", "CDHR1", "SLC25A44", "SLC7A7", "REEP1", "PINK1",
      "GAS7", "KLHL18", "HYAL4", "DCAF4", "TGFBR3", "SLC23A2", "TUBB4A",
      "VAV2", "SLC16A10", "IL6R", "ARPC1B", "MYBL2", "TNS3", "CACNA2D2",
      "ITGAM", "GALNT6", "NDUFAB1", "KIE5", "UPP1", "DACT2"
    ),
    primed = c(
      "SALL2", "DUSP6", "LRRN1", "P3H2", "CDH2", "VRTN", "UCHL1", "STC1", "FGFBP3",
      "NELL2", "ANOS1", "FZD7", "THY1", "PHLDA1", "USP44", "NAP1L3",
      "SPRY4", "PTPRZ1", "EDNRB", "ADAMTS19", "FREM2", "PCYT1B", "TAGLN",
      "CAV1", "COLZA1", "CYP26A1", "MAP7", "PODXL", "GRPR", "NTS",
      "PLCH1", "COL18A1", "PCDH18", "CRABP1", "EPHA2", "VIM", "NECTIN3",
      "GI12", "FAM13A", "DPYSL2", "ATP8A1", "PMEL", "ZIC2", "GPC6",
      "FKBP10", "SEMA3A", "SALL1", "ROR1", "НЕРН", "МЕХЗ"
    )
  )

  return(list("smith" = smith_markers))
}

