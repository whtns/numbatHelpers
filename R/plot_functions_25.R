# Plot Functions (124)

#' Create a plot visualization
#'
#' @param mymarkers Parameter for mymarkers
#' @param seu_paths File path
#' @param plot_type Parameter for plot type
#' @param group_by Character string (default: "gene_snn_res.0.2")
#' @param extension Character string (default: "_filtered")
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging
# - rownames_round_trip: Avoid unnecessary rownames conversions - keep as column when possible

plot_putative_marker_across_samples <- function(mymarkers, seu_paths, plot_type = FeaturePlot, group_by = "gene_snn_res.0.2", extension = "_filtered") {
  print(mymarkers)

  sample_ids <- str_extract(seu_paths, "SR[RX][0-9]+")

  myplots <- map(seu_paths, plot_markers_in_sample, mymarkers = mymarkers, plot_type = plot_type, group_by = group_by) %>%
    set_names(sample_ids)


  myplots0 <-
    myplots %>%
    transpose() %>%
    imap(~ {
      patchwork::wrap_plots(.x) +
        plot_annotation(
          title = .y
        )
    })

  if (identical(plot_type, VlnPlot)) {
    plot_type_label <- "VlnPlot"
  } else if (identical(plot_type, FeaturePlot)) {
    plot_type_label <- "FeaturePlot"
  }


  plot_paths <- glue("results/numbat_sridhar/gene_plots/{names(myplots0)}_{plot_type_label}_{group_by}_{extension}.pdf")

  map2(plot_paths, myplots0, ~ ggsave(.x, .y, height = 8, width = 14))

  plot_path <- qpdf::pdf_combine(plot_paths, glue("results/numbat_sridhar/gene_plots/{plot_type_label}_{group_by}_{extension}.pdf"))

  fs::file_delete(plot_paths)

  return(plot_path)
}

#' Perform differential expression analysis
#'
#' @param numbat_rds_file File path
#' @param cluster_dictionary Cluster information
#' @param ident.1 Cell identities or groups
#' @param ident.2 Cell identities or groups
#' @return ggplot2 plot object
#' @export
find_diffex_bw_clusters_for_each_clone <- function(numbat_rds_file, cluster_dictionary, ident.1 = "G2M", ident.2 = "cone") {
  #
  sample_id <- str_extract(numbat_rds_file, "SR[RX][0-9]+")


  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[, !is.na(seu$clone_opt)]

  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")

  seu <- AddMetaData(seu, test0)

  
#' Perform clustering analysis
#'
#' @param clone_for_diffex Parameter for clone for diffex
#' @param seu Seurat object
#' @return Function result
#' @export
cluster_diff_per_clone <- function(clone_for_diffex, seu) {
    #
    seu0 <- seu[, seu$clone_opt == clone_for_diffex]

    Idents(seu0) <- seu0$abbreviation

    diffex <- FindAllMarkers(seu0, group.by = "abbreviation")
  }

  myclones <- sort(unique(seu$clone_opt)) %>%
    set_names(.)

  possible_cluster_diff_per_clone <- possibly(cluster_diff_per_clone)

  diffex <- map(myclones, possible_cluster_diff_per_clone, seu)

  diffex0 <- map(diffex, compact) %>%
    map(tibble::rownames_to_column, "symbol") %>%
    dplyr::bind_rows(.id = "clone") %>%
    dplyr::select(-symbol) %>%
    dplyr::rename(symbol = gene) %>%
    dplyr::mutate(sample_id = sample_id) %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::select(symbol, description, everything()) %>%
    dplyr::distinct(symbol, clone, .keep_all = TRUE) %>%
    dplyr::arrange(cluster, clone, p_val_adj) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    identity()

  diffex_path <- glue("results/{numbat_dir}/{sample_id}_clone_diffex.csv")
  write_csv(diffex0, diffex_path)


  
#' Perform differential expression analysis
#'
#' @param diffex_bw_clusters_for_each_clone Cluster information
#' @return Differential expression results
#' @export
compare_diffex_cluster_by_clone <- function(diffex_bw_clusters_for_each_clone) {
    test0 <-
      diffex_bw_clusters_for_each_clone %>%
      dplyr::group_by(clone, cluster) %>%
      dplyr::slice_head() %>%
      dplyr::arrange(cluster, clone)
  }

  diffex1 <- compare_diffex_cluster_by_clone(diffex0)

  cluster_clone <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, clone_opt, abbreviation) %>%
    tidyr::unite(cluster_clone, abbreviation, clone_opt) %>%
    dplyr::select(cell, cluster_clone) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- AddMetaData(seu, cluster_clone)

  diffex_cluster_by_clone_dotplot <-
    DotPlot(seu, features = unique(diffex1$symbol), group.by = "cluster_clone") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    coord_flip() +
    labs(title = sample_id)

  dotplot_path <- glue("results/{numbat_dir}/{sample_id}_cluster_clone_diffex_dotplot.pdf")
  ggsave(dotplot_path, diffex_cluster_by_clone_dotplot, width = 10, height = 12)

  trend_genes <- diffex1 %>%
    slice_head(n = 1)

  diffex_cluster_by_clone_trendplot <-
    plot_gene_clone_trend(seu, trend_genes$symbol) +
    labs(title = sample_id)

  trendplot_path <- glue("results/{numbat_dir}/{sample_id}_gene_trend.pdf")
  ggsave(trendplot_path, diffex_cluster_by_clone_trendplot, width = 10, height = 30)

  return(list("diffex" = diffex_path, "plot" = dotplot_path))
}
#' Perform clustering analysis
#'
#' @param cluster_for_diffex Cluster information
#' @param seu Seurat object
#' @param group.by Parameter for group.by
#' @return Function result
#' @export
clone_diff_per_cluster <- function(cluster_for_diffex, seu, group.by) {
    #
    seu0 <- seu[, seu[["clusters"]] == cluster_for_diffex]

    Idents(seu0) <- seu0$clone_opt

    # diffex <- FindAllMarkers(seu0) %>%
    #   dplyr::rename(clone_opt = cluster)

    diffex <- imap(clone_comparisons, make_clone_comparison_integrated, seu0, mynbs, location = location)

    return(diffex)
  }

plot_fig_09_10 <- function(seu_subset, corresponding_seus, corresponding_clusters_diffex, corresponding_clusters_enrichments, plot_path = "result/dotplot_recurre.pdf", widths = rep(8, 3), heights = rep(12, 3), common_seus = NULL, ...){
	
	common_seus <- common_seus %||% seu_subset
	
	names(corresponding_seus) <- fs::path_file(unlist(corresponding_seus))
	names(corresponding_clusters_diffex) <- names(corresponding_seus)
	names(corresponding_clusters_enrichments) <- names(corresponding_seus)
	
	seu_subset <- corresponding_seus[corresponding_seus %in% seu_subset]
	
	corresponding_clusters_diffex <- corresponding_clusters_diffex[names(seu_subset)]
	
	corresponding_clusters_enrichments <- corresponding_clusters_enrichments[names(seu_subset)]
	
	enrichment_files <- compare_corresponding_enrichments(corresponding_clusters_enrichments, common_seus, plot_path = tempfile(fileext = ".pdf"), width = widths[[1]]*2, height = heights[[1]])
	
	all_dotplot <- corresponding_clusters_diffex |> 
		map_depth(2, "all") |> 
		dotplot_recurrent_genes(...)
	
	cis_dotplot <- corresponding_clusters_diffex |> 
		map_depth(2, "cis") |> 
		dotplot_recurrent_genes(...)
	
	trans_dotplot <- corresponding_clusters_diffex |> 
		map_depth(2, "trans") |> 
		dotplot_recurrent_genes(...)
	
	fig_09_10_panels <- list(
		"all" = all_dotplot,
		"cis" = cis_dotplot,
		"trans" = trans_dotplot
	) |> 
		imap(~{.x + labs(title = .y)}) %>%
		identity()
	
	tmpplots <- 
		list(
		"plot" = fig_09_10_panels, 
		"height" = heights, 
		"width" = widths) |> 
		pmap(~ggsave(tempfile(fileext = ".pdf"), plot = ..1, height = ..2, width = ..3))
	
	qpdf::pdf_combine(c(tmpplots, enrichment_files), plot_path)
	
	return(plot_path)
	
}