# Plot Functions (135)

#' Create a plot visualization
#'
#' @param . Parameter for .
#' @param readRDS Parameter for readRDS
#' @param FeaturePlot Parameter for FeaturePlot
#' @param .x Parameter for .x
#' @param .y Parameter for .y
#' @param VlnPlot Parameter for VlnPlot
#' @return ggplot2 plot object
#' @export
plot_scna_violins <- function(., readRDS, FeaturePlot, .x, .y, VlnPlot) {
  # plot drivers ------------------------------

  samples_1q <- c(
    "SRX11133594", "SRX11133592", "SRX11133593", "SRX10264526",
    "SRX14116944"
  )

  samples_2p <- c("SRX10264526", "SRX14116944", "SRX14116947", "SRX10264524")

  samples_6p <- c("SRX14116944", "SRX10264524", "SRX14116947")

  samples_16q <- c("SRX11133588", "SRX10264519", "SRX11133593", "SRX11133594", "SRX11133585", "SRX11133592", "SRX10264520")

  # samples_16q <- c("SRX11133588", "SRX11133585", "SRX11133592")

  
#' Generate a feature plot
#'
#' @param sample_ids Parameter for sample ids
#' @param features Parameter for features
#' @return ggplot2 plot object
#' @export
plot_markers_featureplot <- function(sample_ids, features) {
    #

    seus <- glue("output/seurat/{sample_ids}_filtered_seu.rds") %>%
      set_names(str_extract(., "SR[RX][0-9]+")) %>%
      map(readRDS) %>%
      identity()

    markerplots <- map(seus, FeaturePlot, features = features, combine = TRUE)

    markerplots0 <- imap(markerplots, ~ {
      .x + labs(subtitle = .y)
    })

    vlnplots <- map(seus, VlnPlot, features = features, group.by = "GT_opt", combine = TRUE, pt.size = 0)

    vlnplots0 <- imap(vlnplots, ~ {
      .x +
        labs(title = .y) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
        )
    })

    # vln_wrapped = wrap_plots(vlnplots0, ncol = 2) + plot_annotation(title = features)

    return(vlnplots0)
  }

  vlns_1q <- plot_markers_featureplot(samples_1q, "CENPF")

  vlns_2p <- plot_markers_featureplot(samples_2p, "SOX11")

  vlns_6p <- plot_markers_featureplot(samples_6p, "DEK")

  vlns_16q <- plot_markers_featureplot(samples_16q, "CDT1")


  vlns_1q +
    wrap_plots(nrow = 1) + plot_annotation(title = "CENPF")
  ggsave("results/vln_plots_from_oncoprint_1q.pdf", width = 6, height = 8)
  vlns_2p +
    wrap_plots(nrow = 1) + plot_annotation(title = "SOX11")
  ggsave("results/vln_plots_from_oncoprint_2p.pdf", width = 6, height = 6)
  vlns_6p +
    wrap_plots(nrow = 1) + plot_annotation(title = "DEK")
  ggsave("results/vln_plots_from_oncoprint_6p.pdf", width = 6, height = 6)
  vlns_16q +
    wrap_plots(nrow = 1) + plot_annotation(title = "CDT1")
  ggsave("results/vln_plots_from_oncoprint_16q.pdf", width = 6, height = 6)
}
#' Perform enrichment analysis
#'
#' @param enrichment_tables Parameter for enrichment tables
#' @param cis_plot_file File path
#' @param trans_plot_file File path
#' @param cis_table_file File path
#' @param trans_table_file File path
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
compile_cis_trans_enrichment_recurrence <- function(enrichment_tables,
                                                    cis_plot_file = "results/cis_enrichment_plots",
                                                    trans_plot_file = "results/trans_enrichment_plots",
                                                    cis_table_file = "results/cis_enrichment_tables.xlsx",
                                                    trans_table_file = "results/trans_enrichment_tables.xlsx",
                                                    ...) {
  # cis ------------------------------
  recurrences <- c(3, 1, 1, 1)

  titles <- c(
    "1q+ cis enrichment",
    "2p+ cis enrichment",
    "6p+ cis enrichment",
    "16q- cis enrichment"
  )

  cis_enrich_results <- pmap(list(enrichment_tables$cis, recurrences, titles), plot_enrichment_recurrence, ...)

  # trans ------------------------------

  recurrences <- c(2, 2, 2, 3)

  titles <- c(
    "1q+ trans enrichment",
    "2p+ trans enrichment",
    "6p+ trans enrichment",
    "16q- trans enrichment"
  )

  trans_enrich_results <- pmap(list(enrichment_tables$trans, recurrences, titles), plot_enrichment_recurrence, ...)

  # cis ------------------------------
  # enrich plot
  cis_enrich_plots <- map(cis_enrich_results, "enrich_plot") %>%
    purrr::discard(~ nrow(.x$data) < 1)

  names(cis_enrich_plots) <- str_remove_all(names(cis_enrich_plots), "[+]|-")

  enrich_plot_widths <- map_int(cis_enrich_plots, ~ n_distinct(.x$data$clone_comparison)) * 1 + 1

  fs::dir_create(fs::path_ext_remove(cis_plot_file))
  cis_pdfs <- map2(names(cis_enrich_plots), enrich_plot_widths, ~ ggsave(glue("{fs::path_ext_remove(cis_plot_file)}/{.x}_enrich.pdf"), cis_enrich_plots[[.x]], width = .y, height = 7))

  qpdf::pdf_combine(cis_pdfs, cis_plot_file)

  per_scna_enrichment_cis <- imap(enrichment_tables$cis, plot_enrichment_per_scna, cis_plot_file)

  # trans ------------------------------
  # enrich plot
  trans_enrich_plots <- map(trans_enrich_results, "enrich_plot") %>%
    purrr::discard(~ nrow(.x$data) < 1)

  names(trans_enrich_plots) <- str_remove_all(names(trans_enrich_plots), "[+]|-")

  enrich_plot_widths <- map_int(trans_enrich_plots, ~ n_distinct(.x$data$clone_comparison)) * 1 + 1

  fs::dir_create(fs::path_ext_remove(trans_plot_file))
  trans_pdfs <- map2(names(trans_enrich_plots), enrich_plot_widths, ~ ggsave(glue("{fs::path_ext_remove(trans_plot_file)}/{.x}_enrich.pdf"), trans_enrich_plots[[.x]], width = .y, height = 7, limitsize = FALSE))

  qpdf::pdf_combine(trans_pdfs, trans_plot_file)

  per_scna_enrichment_trans <- imap(enrichment_tables$trans, plot_enrichment_per_scna, trans_plot_file)

  # cis tables
  cis_tables <- map(cis_enrich_results, "table") %>%
    purrr::discard(~ nrow(.x) < 1) %>%
    purrr::map(dplyr::mutate, comparison = purrr::pluck(str_split_1(clone_comparison, "_"), 2)) %>%
    purrr::map(dplyr::arrange, comparison, p.adjust)

  writexl::write_xlsx(cis_tables, cis_table_file)

  # trans tables
  trans_tables <- map(trans_enrich_results, "table") %>%
    purrr::discard(~ nrow(.x) < 1) %>%
    purrr::map(dplyr::mutate, comparison = purrr::pluck(str_split_1(clone_comparison, "_"), 2)) %>%
    purrr::map(dplyr::arrange, comparison, p.adjust)

  writexl::write_xlsx(trans_tables, trans_table_file)

  return(
    list(
      "cis_enrich_compiled" = cis_plot_file,
      "trans_enrich_compiled" = trans_plot_file,
      "cis_table" = cis_table_file,
      "trans_table" = trans_table_file,
      "cis_per_sample" = per_scna_enrichment_cis,
      "trans_per_sample" = per_scna_enrichment_trans
    )
  )
}

#' Perform enrichment analysis
#'
#' @param enrichment_tables Parameter for enrichment tables
#' @param cis_plot_file File path
#' @param trans_plot_file File path
#' @param cis_table_file File path
#' @param trans_table_file File path
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
compile_cis_trans_enrichment_recurrence_by_cluster <- function(enrichment_tables,
                                                               cis_plot_file = "results/cis_enrichment_plots",
                                                               trans_plot_file = "results/trans_enrichment_plots",
                                                               cis_table_file = "results/cis_enrichment_tables.xlsx",
                                                               trans_table_file = "results/trans_enrichment_tables.xlsx",
                                                               ...) {
  # cis ------------------------------
  recurrences <- c(3, 1, 1, 1)

  titles <- c(
    "1q+ cis enrichment",
    "2p+ cis enrichment",
    "6p+ cis enrichment",
    "16q- cis enrichment"
  )

  cis_enrich_results <- pmap(list(enrichment_tables$cis, recurrences, titles), plot_enrichment_recurrence_by_cluster, ...)

  # trans ------------------------------

  recurrences <- c(2, 2, 2, 3)

  titles <- c(
    "1q+ trans enrichment",
    "2p+ trans enrichment",
    "6p+ trans enrichment",
    "16q- trans enrichment"
  )

  trans_enrich_results <- pmap(list(enrichment_tables$trans, recurrences, titles), plot_enrichment_recurrence_by_cluster, ...)

  # cis ------------------------------
  # enrich plot
  cis_enrich_plots <- map(cis_enrich_results, "enrich_plot_by_phase") %>%
    purrr::discard(~ nrow(.x$data) < 1)

  names(cis_enrich_plots) <- str_remove_all(names(cis_enrich_plots), "[+]|-")

  enrich_plot_widths <- map_int(cis_enrich_plots, ~ n_distinct(.x$data$clone_comparison)) * 0.3 + 1

  fs::dir_create(fs::path_ext_remove(cis_plot_file))
  cis_pdfs <- map2(names(cis_enrich_plots), enrich_plot_widths, ~ ggsave(glue("{fs::path_ext_remove(cis_plot_file)}/{.x}_enrich.pdf"), cis_enrich_plots[[.x]], width = .y, height = 10))

  qpdf::pdf_combine(cis_pdfs, cis_plot_file)

  per_scna_enrichment_cis <- imap(enrichment_tables$cis, plot_enrichment_per_scna, cis_plot_file)

  # trans ------------------------------
  # enrich plot
  trans_enrich_plots <- map(trans_enrich_results, "enrich_plot_by_phase") %>%
    purrr::discard(~ nrow(.x$data) < 1)

  names(trans_enrich_plots) <- str_remove_all(names(trans_enrich_plots), "[+]|-")

  enrich_plot_widths <- map_int(trans_enrich_plots, ~ n_distinct(.x$data$clone_comparison)) * 0.3 + 1

  fs::dir_create(fs::path_ext_remove(trans_plot_file))
  trans_pdfs <- map2(names(trans_enrich_plots), enrich_plot_widths, ~ ggsave(glue("{fs::path_ext_remove(trans_plot_file)}/{.x}_enrich.pdf"), trans_enrich_plots[[.x]], width = .y, height = 10, limitsize = FALSE))

  qpdf::pdf_combine(trans_pdfs, trans_plot_file)

  per_scna_enrichment_trans <- imap(enrichment_tables$trans, plot_enrichment_per_scna, trans_plot_file)

  # cis tables
  cis_tables <- map(cis_enrich_results, "table") %>%
    purrr::discard(~ nrow(.x) < 1) %>%
    purrr::map(dplyr::mutate, comparison = purrr::pluck(str_split_1(clone_comparison, "_"), 2)) %>%
    purrr::map(dplyr::arrange, comparison, p.adjust)

  writexl::write_xlsx(cis_tables, cis_table_file)

  # trans tables
  trans_tables <- map(trans_enrich_results, "table") %>%
    purrr::discard(~ nrow(.x) < 1) %>%
    purrr::map(dplyr::mutate, comparison = purrr::pluck(str_split_1(clone_comparison, "_"), 2)) %>%
    purrr::map(dplyr::arrange, comparison, p.adjust)

  writexl::write_xlsx(trans_tables, trans_table_file)

  return(
    list(
      "cis_enrich_compiled" = cis_plot_file,
      "trans_enrich_compiled" = trans_plot_file,
      "cis_table" = cis_table_file,
      "trans_table" = trans_table_file,
      "cis_per_sample" = per_scna_enrichment_cis,
      "trans_per_sample" = per_scna_enrichment_trans
    )
  )
}

#' Create a plot visualization
#'
#' @param seu_path File path
#' @param nb_path File path
#' @param clone_simplifications Parameter for clone simplifications
#' @param gene_lists Gene names or identifiers
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
score_and_vlnplot_seu <- function(seu_path, nb_path, clone_simplifications, gene_lists, ...) {
  #

  mysample <- str_extract(seu_path, "SR[RX][0-9]+")

  seu <- readRDS(seu_path)

  # gene_lists <- map(gene_lists, ~(.x[.x %in% rownames(seu$SCT)]))

  # DefaultAssay(seu) <- "gene"
  nbin <- floor(length(VariableFeatures(seu)) / 100)

  seu <- Seurat::AddModuleScore(seu, features = gene_lists, name = "subtype", nbin = nbin, ctrl = 100)

  module_names <- paste0("subtype", seq(length(gene_lists)))

  seu@meta.data[names(gene_lists)] <- seu@meta.data[module_names]

  # cluster_seu <- seu[,!is.na(seu$leiden)]

  
#' Perform pub violin operation
#'
#' @param score.by Character string (default: "subtype2")
#' @param group.by Character string (default: "scna")
#' @param seu Seurat object
#' @param y_lim Parameter for y lim
#' @param step Parameter for step
#' @return Function result
#' @export
pub_violin <- function(score.by = "subtype2", group.by = "scna", seu, y_lim = 0.4, step = 0.1) {
    mydf <- seu@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::mutate(scna = fct_reorder(scna, clone_opt)) %>%
      identity()

    if (is.factor(mydf[[group.by]])) {
      mycomparisons <- seq(1, length(levels(mydf[[group.by]])))
    } else {
      mycomparisons <- rev(unique(mydf[[group.by]]))
    }

    mycomparisons <- map2(head(mycomparisons, -1), mycomparisons[-1], c)

    ggpubr::ggviolin(mydf, x = group.by, y = score.by, fill = group.by, add = "boxplot", add.params = list(fill = "white")) +
      # stat_compare_means(label.y = y_lim) +
      stat_compare_means(comparisons = mycomparisons, label = "p.signif") +
      scale_y_continuous(limits = c(0, y_lim), breaks = seq(0, y_lim, step))
  }

  # cluster_vln = map(names(gene_lists), pub_violin, group.by = "abbreviation", seu = seu, ...)
  #
  # cluster_vln <- map(cluster_vln, ~{.x + guides(fill="none") + scale_x_discrete(labels = function(x) wrap_scna_labels(x))})
  #
  # names(cluster_vln) <- names(gene_lists)

  scna_vln <- map(names(gene_lists), pub_violin, group.by = "clone_opt", seu = seu, ...)

  names(scna_vln) <- paste0(names(gene_lists), " genes")

  scna_vln <- imap(scna_vln, ~ {
    .x + guides(fill = "none") + scale_x_discrete(labels = function(x) wrap_scna_labels(x)) + labs(x = "clone", y = .y)
  })


  myclonetree <- plot_clone_tree(seu, mysample, nb_path, clone_simplifications)

  design <-
    "AAAA
  AAAA
  BBCC
  BBCC
  BBCC"

  test0 <- list(myclonetree, scna_vln[[1]], scna_vln[[2]]) %>%
    wrap_plots() +
    plot_layout(design = design)

  scna_vln_plot_path <- ggsave(glue("results/{mysample}_scna_vln.pdf"), height = 6, width = 6)

  return(scna_vln_plot_path)
}

make_all_numbat_plots <- function(numbat_dir, num_iter = 2, min_LLR = 2, genome = "hg38", init_k = 3, gtf = gtf_hg38, overwrite = FALSE) {
  sample_id <- path_file(numbat_dir)

  print(numbat_dir)

  for (i in seq(num_iter)) {
    bulk_clone_path <- glue("{numbat_dir}/bulk_clones_{i}.png")
    if (!file.exists(bulk_clone_path) | (overwrite <- TRUE)) {
      # Plot bulk clones
      bulk_clones <- read_tsv(glue("{numbat_dir}/bulk_clones_{i}.tsv.gz"), col_types = cols())
      p <- plot_bulks(bulk_clones,
        min_LLR = min_LLR, use_pos = TRUE,
        genome = genome
      ) +
        labs(title = sample_id)

      ggsave(bulk_clone_path, p,
        width = 13, height = 2 * length(unique(bulk_clones$sample)),
        dpi = 250
      )
      print(glue("plotted {numbat_dir}/bulk_clones_{i}.png"))
    }


    bulk_subtrees_path <- glue("{numbat_dir}/bulk_subtrees_{i}.png")
    if (!file.exists(bulk_subtrees_path) | (overwrite <- TRUE)) {
      # Plot bulk subtrees
      bulk_subtrees <- read_tsv(glue("{numbat_dir}/bulk_subtrees_{i}.tsv.gz"), col_types = cols())
      p <- plot_bulks(bulk_subtrees,
        min_LLR = min_LLR,
        use_pos = TRUE, genome = genome
      ) +
        labs(title = sample_id)

      ggsave(glue("{numbat_dir}/bulk_subtrees_{i}.png"),
        p,
        width = 13, height = 2 * length(unique(bulk_subtrees$sample)),
        dpi = 250
      )
    }
  }

  final_bulk_clones_path <- glue("{numbat_dir}/bulk_clones_final.png")
  if (!file.exists(final_bulk_clones_path) | (overwrite <- TRUE)) {
    bulk_clones <- read_tsv(glue("{numbat_dir}/bulk_clones_final.tsv.gz"), col_types = cols())
    p <- plot_bulks(bulk_clones,
      min_LLR = min_LLR, use_pos = TRUE,
      genome = genome
    ) +
      labs(title = sample_id)
    ggsave(final_bulk_clones_path, p,
      width = 13,
      height = 2 * length(unique(bulk_clones$sample)),
      dpi = 250
    )
  }

  # exp_clust_path = glue("{numbat_dir}/exp_roll_clust.png")
  # if(!file.exists(exp_clust_path)){
  #
  #   # # Plot single-cell smoothed expression magnitude heatmap
  #   gexp_roll_wide <- read_tsv(glue("{numbat_dir}/gexp_roll_wide.tsv.gz"), col_types = cols()) %>%
  #     column_to_rownames("cell")
  #   hc <- readRDS(glue("{numbat_dir}/hc.rds"))
  #   p = plot_exp_roll(gexp_roll_wide = gexp_roll_wide,
  #                     hc = hc, k = init_k, gtf = gtf, n_sample = 10000)
  #   labs(title = sample_id)
  #   ggsave(exp_clust_path, p,
  #          width = 8, height = 4, dpi = 200)
  #
  # }

  # phylo_heatmap_path = glue("{numbat_dir}/phylo_heatmap.png")
  #
  #   if(!file.exists(phylo_heatmap_path)){
  #
  #     # # Plot single-cell CNV calls along with the clonal phylogeny
  #     nb <- readRDS(glue("{numbat_dir}_numbat.rds"))
  #     mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")
  #
  #     phylo_heatmap <- nb$plot_phylo_heatmap(
  #       pal_clone = mypal,
  #       show_phylo = TRUE
  #     ) +
  #       labs(title = sample_id)
  #     ggsave(phylo_heatmap_path, width = 13,
  #            height = 10,
  #            dpi = 250)
  #
  #   }


  return("success!")
}

heatmap_marker_genes_debug <- function(seu, gene_lists, label = "", marker_col = "clusters", group.by = c("clusters", "scna", "hypoxia1", "hypoxia2"), col_arrangement = c("clusters", "hypoxia1"), split_columns = "clusters", ...) {
    
    test0 <- 
        gene_lists |> 
        tibble::enframe("group", "Gene.Name") |> 
        tidyr::unnest(Gene.Name) |> 
        dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
        # dplyr::mutate(mp = replace_na(mp, "")) %>%
        identity()
    
    mymarkers <-
        test0 %>%
        dplyr::pull(Gene.Name)
    
    seu$scna <- factor(seu$scna)
    levels(seu$scna)[1] <- "none"
    
    # row_ha <- ComplexHeatmap::rowAnnotation(mp = rev(test0$clusters))
    
    if (!is.null(split_columns)) {
        column_split <- sort(seu@meta.data[[split_columns]])
        column_title <- unique(column_split)
    } else {
        column_split <- split_columns
        column_title <- NULL
    }
    
    ggplotify::as.ggplot(seu_complex_heatmap(seu, features = mymarkers, group.by = group.by, col_arrangement = col_arrangement, cluster_rows = FALSE, column_split = column_split)) +
        labs(...) +
        NULL
    
    # seu_heatmap <- Seurat::DoHeatmap(seu, features = mymarkers, group.by = "SCT_snn_res.0.4")
    
    heatmap_file <- tempfile(tmpdir = "results", fileext = ".pdf")
    
    ggsave(heatmap_file, ...)
    
    return(heatmap_file)
}

heatmap_marker_genes <- function(seu_path, common_genes, gene_lists, label = "", sample_id = NULL, marker_col = "SCT_snn_res.0.4", group.by = c("SCT_snn_res.0.4", "scna", "subtype1", "subtype2"), col_arrangement = c("SCT_snn_res.0.4", "scna")) {
  #
  if (is.null(sample_id)) {
    sample_id <- str_extract(seu_path, "SR[RX][0-9]+")
  }

  seu <- readRDS(seu_path)

  seu <- Seurat::AddModuleScore(seu, features = gene_lists, name = "subtype")

  module_names <- paste0("subtype", seq(length(gene_lists)))

  seu@meta.data[names(gene_lists)] <- seu@meta.data[module_names]

  test0 <- seu@misc$markers[[marker_col]]$presto %>%
    dplyr::group_by(Cluster) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::select(Gene.Name, Cluster) %>%
    dplyr::inner_join(common_genes, by = "Gene.Name") %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Cluster) %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE) %>%
    dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
    dplyr::mutate(mp = replace_na(mp, "")) %>%
    identity()

  mymarkers <-
    test0 %>%
    dplyr::pull(Gene.Name)

  seu$scna <- factor(seu$scna)
  levels(seu$scna)[1] <- "none"

  row_ha <- ComplexHeatmap::rowAnnotation(mp = rev(test0$mp))

  ggplotify::as.ggplot(seu_complex_heatmap(seu, features = mymarkers, group.by = group.by, col_arrangement = col_arrangement, cluster_rows = FALSE, right_annotation = row_ha)) +
    labs(title = sample_id)


  # seu_heatmap <- Seurat::DoHeatmap(seu, features = mymarkers, group.by = "SCT_snn_res.0.4")

  heatmap_file <- glue("results/{sample_id}_{label}heatmap.pdf")

  ggsave(heatmap_file, height = 8, width = 8)

  return(heatmap_file)
}

find_diffex_bw_divergent_clusters <- function(sample_id, tumor_id, seu, mynb, to_SCT_snn_res. = "SCT_snn_res.1", to_clust = c("1", "10"), clone_comparisons) {
  possible_make_cluster_comparison <- possibly(make_cluster_comparison)

  clone_comparisons <- clone_comparisons[[tumor_id]] %>%
    tibble::enframe("comparison", "segment") %>%
    dplyr::mutate(clones = str_extract(comparison, "[0-9]_v_[0-9]")) %>%
    dplyr::mutate(clones = str_split(clones, "_v_")) %>%
    dplyr::rowwise() %>%
    dplyr::filter(all(clones %in% seu$clone_opt)) %>%
    dplyr::select(comparison, segment) %>%
    tibble::deframe() %>%
    identity()

  message(clone_comparisons)

  diffex <- imap(clone_comparisons, make_cluster_comparison, seu, mynb, to_SCT_snn_res., to_clust) %>%
    identity()

  return(diffex)
}

find_diffex_bw_clones_for_each_cluster_integrated <- function(seu_path, kept_samples = c("SRX11133594", "SRX11133593", "SRX11133592"), clone_comparisons = list("2_v_1_16q-" = c("16c", "16b"), "3_v_2_1q+" = c("1b")), location = "cis", scna_of_interest = "1q") {

  seu <- readRDS(seu_path)

  numbat_dir <- "numbat_sridhar"

  mynbs <- glue("output/{numbat_dir}/{kept_samples}_numbat.rds") %>%
    map(readRDS)

  location <- "cis"

  seu <- seu[, !is.na(seu$clone_opt)]

  clone_diff_per_cluster <- function(cluster_for_diffex, seu, group.by) {
    #
    seu0 <- seu[, seu[["clusters"]] == cluster_for_diffex]

    Idents(seu0) <- seu0$clone_opt

    # diffex <- FindAllMarkers(seu0) %>%
    #   dplyr::rename(clone_opt = cluster)

    diffex <- imap(clone_comparisons, make_clone_comparison_integrated, seu0, mynbs, location = location)

    return(diffex)
  }

  myclusters <- sort(unique(seu$clusters)) %>%
    set_names(.)

  possible_clone_diff_per_cluster <- possibly(clone_diff_per_cluster)

  diffex <- map(myclusters, possible_clone_diff_per_cluster, seu) %>%
    compact() %>%
    map(bind_rows, .id = "clone_comparison") %>%
    bind_rows(.id = "cluster") %>%
    # dplyr::arrange(cluster, p_val_adj) %>%
    identity()

  kept_samples_slug <- paste(kept_samples, collapse = "_")

  diffex_path <- glue("{numbat_dir}/{kept_samples_slug}_cluster_clone_comparison_diffex_{location}.csv")
  write_csv(diffex, diffex_path)

  return(diffex_path)
}

drop_bad_cells <- function(seu_path, bad_cell_types = c("RPCs", "Late RPCs", c("Red Blood Cells", "Microglia", "Muller Glia", "RPE", "Horizontal Cells", "Rod Bipolar Cells", "Pericytes", "Bipolar Cells", "Astrocytes", "Endothelial", "Schwann", "Fibroblasts"))) {
  seu <- readRDS(seu_path)

  seu <- seu[, !seu$type %in% bad_cell_types]

  add_hash_metadata(seu = seu, filepath = str_replace(seu_path, "_seu.rds", "_dropped_cells_seu.rds"))

  retainedcells_dimplot <- DimPlot(
    seu,
    group.by = "type"
  ) +
    labs(title = fs::path_file(seu_path))

  return(seu_path)
}