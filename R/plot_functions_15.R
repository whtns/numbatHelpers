# Plot Functions (107)

#' Create a plot visualization
#'
#' @param seu_list Parameter for seu list
#' @param cluster_orders Cluster information
#' @param plot_path File path
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - rownames_round_trip: Avoid unnecessary rownames conversions - keep as column when possible
# - map_bind_rows: Use map_dfr() instead of map() %>% bind_rows() for better performance

sample_specific_analyses_of_tumors_with_scna_subclones_after_integration <- function(seu_list, cluster_orders, plot_path = "results/fig_s04.pdf"){
  
  
	
	plot_path <- seu_list |> 
		imap(~make_clone_distribution_figure(.x, cluster_orders, height = 12, width = 20, plot_path = tempfile(pattern = .y, fileext = ".pdf"), heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"))) |> 
		qpdf::pdf_combine(plot_path)
	
	return(plot_path)
}

#' Create a plot visualization
#'
#' @return ggplot2 plot object
#' @export
plot_fig_s09 <- function(){
  
  
	seu <- readRDS("output/seurat/integrated_1q/integrated_seu_1q_complete.rds")
	
	comparisons <- 
		tibble::tribble(
			~ident.1, ~ident.2,
			"g1_6", "g1_3",
			# "g1_0", "g1_3",
			"g1_6", "g1_0",
			# "g2_m_5", "s_g2_4",
			# "s_star_10", "s_1"
		) |> 
		dplyr::mutate(comparison = glue("{ident.1} v. {ident.2}"))
	
	diffexes <- map2(comparisons$ident.1, comparisons$ident.2, ~FindMarkers(seu, group.by = "clusters", ident.1 = .x, .y)) |> 
		set_names(comparisons$comparison)
	
	c6_enrichments <- diffexes |> 
		map(enrichment_analysis, analysis_method = "ora", gene_set = "C6")
	
	c6_enrichment_plots <-
		c6_enrichments |> 
		map(plot_enrichment, analysis_method = "ora") |> 
		imap(~{.x + labs(title = .y)})
	
	h_enrichments <- diffexes |> 
		map(enrichment_analysis, analysis_method = "ora", gene_set = "H")
	
	h_enrichments_plots <- 
		h_enrichments |> 
		map(plot_enrichment, analysis_method = "ora") |> 
		imap(~{.x + labs(title = .y)})
	
	
	enrichment_plots <- c(rbind(c6_enrichment_plots, h_enrichments_plots)) |> 
		map(~{.x + scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "_", " "), width = 15)) })
	
	plot_path <- "results/fig_s09.pdf"
	pdf(plot_path, w = 6, h = 6)
	print(enrichment_plots)
	dev.off()
	
	diffex_path <- 
		diffexes |> 
		map(tibble::rownames_to_column, "symbol") |> 
		writexl::write_xlsx("results/fig_s09.xlsx")
	
	return(list("diffex" = diffex_path, "enrichment" = plot_path))
	
}

#' Create a plot visualization
#'
#' @param plot_path File path
#' @param integrated_seu_path File path
#' @param direction Character string (default: "up")
#' @param p_val_adj_threshold Threshold value for filtering
#' @param recurrence_threshold Threshold value for filtering
#' @param n_genes Gene names or identifiers
#' @param ident.1 Cell identities or groups
#' @param ident.2 Cell identities or groups
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_diffex_genes_on_split_integrated_seu <- function(plot_path = "results/fig_07d.pdf", integrated_seu_path = "output/seurat/integrated_1q/integrated_seu_1q_complete.rds", direction = "up", p_val_adj_threshold = 0.05, recurrence_threshold = 3, n_genes = 10, ident.1 = "w_scna", ident.2 = "wo_scna", ...){

	integrated_seu <- readRDS(integrated_seu_path)
	seu_list <- integrated_seu |> 
		SplitObject(split.by = "batch")
	
	for(batch_name in names(seu_list)){
		DefaultAssay(seu_list[[batch_name]]) <- "gene"
	}
	
	diffex_list <- map(seu_list, ~FindMarkers(.x, group.by = "scna", ident.1 = ident.1, ident.2 = ident.2)) |> 
		map(tibble::rownames_to_column, "symbol") |> 
		map(dplyr::filter, !symbol %in% unlist(Seurat::cc.genes.updated.2019)) |> 
		map(tibble::column_to_rownames, "symbol")
	
	gene_list <- 
		diffex_list |> 
		map(tibble::rownames_to_column, "symbol") |> 
		dplyr::bind_rows(.id = "sample_id") |>
		dplyr::filter(if (direction == "up") {
			avg_log2FC > 0
		} else if (direction == "down") {
			avg_log2FC < 0
		} else {
			avg_log2FC
		}) |> 
		dplyr::filter(p_val_adj <= p_val_adj_threshold) |>
		dplyr::group_by(symbol) %>%
		dplyr::filter(all_same_sign(avg_log2FC)) |>
		dplyr::mutate(abs_mean_FC = abs(mean(avg_log2FC))) %>%
		dplyr::mutate(recurrence = dplyr::n_distinct(sample_id)) %>%
		dplyr::mutate(sign = sign(avg_log2FC)) |> 
		dplyr::arrange(sign, desc(abs_mean_FC)) |>
		dplyr::filter(recurrence >= recurrence_threshold) |> 
		dplyr::distinct(sample_id, .keep_all = TRUE) |> 
		dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
		identity()
	
	gene_list$symbol <- factor(gene_list$symbol, levels = unique(gene_list$symbol))
	
	top_genes <- gene_list |> 
		dplyr::group_by(sign) |> 
		dplyr::distinct(symbol) |> 
		dplyr::slice_head(n = n_genes) |> 
		dplyr::pull(symbol)
	
	diffex_dotplot <- 
	gene_list |> 
		dplyr::filter(symbol %in% top_genes) |> 
		dplyr::mutate(mean_FC = abs_mean_FC*sign) |> 
		dplyr::arrange(mean_FC) |> 
		dplyr::mutate(tumor_id = str_extract(sample_id, "SR[RX][0-9]+")) |> 
		ggplot(aes(y = symbol, x = tumor_id, size = neg_log_p_val_adj, color = avg_log2FC)) + 
		geom_point() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
		scale_color_gradient2()
	
	ggsave(plot_path, ...)

	return(plot_path)
	
}

#' Create a plot visualization
#'
#' @param plot_path File path
#' @param integrated_seu_path File path
#' @param direction Character string (default: "down")
#' @param p_val_adj_threshold Threshold value for filtering
#' @param recurrence_threshold Threshold value for filtering
#' @param n_genes Gene names or identifiers
#' @param ident.1 Cell identities or groups
#' @param ident.2 Cell identities or groups
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_fig_08d <- function(plot_path = "results/fig_08d.pdf", integrated_seu_path = "output/seurat/integrated_16q/integrated_seu_16q_complete.rds", direction = "down", p_val_adj_threshold = 0.1, recurrence_threshold = 2, n_genes = 10, ident.1 = "16q-", ident.2 = ".diploid", ...){
	
	seu_list <- integrated_seu_path |> 
		SplitObject(split.by = "batch")
	
	for(batch_name in names(seu_list)){
		DefaultAssay(seu_list[[batch_name]]) <- "gene"
	}
	
	diffex_list <- map(seu_list, ~FindMarkers(.x, group.by = "scna", ident.1 = ident.1, ident.2 = ident.2)) |> 
		map(tibble::rownames_to_column, "symbol") |> 
		map(dplyr::filter, !symbol %in% unlist(Seurat::cc.genes.updated.2019)) |> 
		map(tibble::column_to_rownames, "symbol")
	
	gene_list <- 
		diffex_list |> 
		map(tibble::rownames_to_column, "symbol") |> 
		dplyr::bind_rows(.id = "sample_id") |>
		dplyr::filter(if (direction == "up") {
			avg_log2FC > 0
		} else if (direction == "down") {
			avg_log2FC < 0
		} else {
			avg_log2FC
		}) |> 
		dplyr::filter(p_val_adj <= p_val_adj_threshold) |>
		dplyr::group_by(symbol) %>%
		dplyr::filter(all_same_sign(avg_log2FC)) |>
		dplyr::mutate(abs_mean_FC = abs(mean(avg_log2FC))) %>%
		dplyr::mutate(recurrence = dplyr::n_distinct(sample_id)) %>%
		dplyr::mutate(sign = sign(avg_log2FC)) |> 
		dplyr::arrange(sign, desc(abs_mean_FC)) |>
		dplyr::filter(recurrence >= recurrence_threshold) |> 
		dplyr::distinct(sample_id, .keep_all = TRUE) |> 
		dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
		identity()
	
	gene_list$symbol <- factor(gene_list$symbol, levels = unique(gene_list$symbol))
	
	top_genes <- gene_list |> 
		dplyr::group_by(sign) |> 
		dplyr::distinct(symbol) |> 
		dplyr::slice_head(n = n_genes) |> 
		dplyr::pull(symbol)
	
	diffex_dotplot <- 
		gene_list |> 
		dplyr::filter(symbol %in% top_genes) |> 
		dplyr::mutate(mean_FC = abs_mean_FC*sign) |> 
		dplyr::arrange(mean_FC) |> 
		dplyr::mutate(tumor_id = str_extract(sample_id, "SR[RX][0-9]+")) |> 
		ggplot(aes(y = symbol, x = tumor_id, size = neg_log_p_val_adj, color = avg_log2FC)) + 
		geom_point() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
		scale_color_gradient2()
	
	ggsave(plot_path, ...)
	
	return(plot_path)
	
}

simplify_gt_col <- function(gt_val, scna_key) {
    #
    if (gt_val != "") {
      gt_vals <-
        gt_val %>%
        str_split(pattern = ",") %>%
        tibble::enframe("name", "seg") %>%
        tidyr::unnest(seg) %>%
        dplyr::left_join(scna_key, by = "seg") %>%
        dplyr::filter(!is.na(scna)) %>%
        dplyr::mutate(scna = paste(scna, collapse = ",")) %>%
        dplyr::pull(scna) %>%
        unique()
      if (length(gt_vals) == 0) gt_vals <- ""
    } else {
      gt_vals <- ""
    }
    return(gt_vals)
  }

  check_tcga_peaks <- function(cohort, arm_ranges = list(
  "1q" = list(chrom = "1", start = 121700000, end = 248956422, direction = "Amp"),
  "2p" = list(chrom = "2", start = 1, end = 93400000, direction = "Amp"),
  "6p" = list(chrom = "6", start = 1, end = 61000000, direction = "Amp"),
  "16q" = list(chrom = "16", start = 33100000, end = 90354753, direction = "Del")
)) {

    # browser()

    tcga_cohort <- 
    cohort |>
    TCGAgistic::tcga_gistic_load(source = "Firehose", cnLevel = "all")


    df <- tcga_cohort@cytoband.summary |>
    separate(Wide_Peak_Limits, into = c("chrom", "range"), sep = ":") %>%
    separate(range, into = c("start", "end"), sep = "-") %>%
    mutate(start = as.integer(start), end = as.integer(end)) |>
    mutate(chrom = str_remove(chrom, "chr"))
    found <- sapply(names(arm_ranges), function(arm) {
      arm_info <- arm_ranges[[arm]]
    #   browser()
      any(
        df$chrom == arm_info$chrom &
        df$qvalues < 0.05 &
        df$start <= arm_info$end &
        df$end >= arm_info$start &
        df$Variant_Classification == arm_info$direction
      )
    })
    data.frame(abbreviation = cohort, t(found))
}

#' Make integrated numbat plots
#' @export
make_integrated_numbat_plots <- function(seu_path, extension = "") {
  #

  sample_id <- str_extract(seu_path, "SRX.*(?=_seu.rds)")

  dir_create(glue("results/numbat_sridhar/{sample_id}"))

  seu <- readRDS(seu_path) %>%
    filter_sample_qc()

  seu <- seu[, !seu$abbreviation %in% c("APOE", "MALAT1")]

  seu@meta.data$abbreviation <-
    seu@meta.data$abbreviation %>%
    str_replace_all("ARL1IP1", "ARL6IP1") %>%
    str_replace_all("PCLAF", "TFF1")


  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  seu <- seu[, !is.na(seu$clone_opt)]

  # seu <- seurat_cluster(seu, resolution = 0.1)
  #
  # plot_markers(seu, metavar = "integrated_snn_res.0.1", marker_method = "presto", return_plotly = FALSE, hide_technical = "all") +
  #   ggplot2::scale_y_discrete(position = "left")

  markerplot <- plot_markers(seu, metavar = "abbreviation", marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 3) +
    ggplot2::scale_y_discrete(position = "left")

  # ggsave(glue("results/numbat_sridhar/{sample_id}/{sample_id}_sample_marker{extension}.pdf"), height = 8, width = 6)

  dimplot <- DimPlot(seu, group.by = c("abbreviation", "Phase")) +
    plot_annotation(title = sample_id)
  # ggsave(glue("results/numbat_sridhar/{sample_id}/{sample_id}_dimplot{extension}.pdf"), width = 8, height = 4)


  ## clone distribution ------------------------------
  distplot <- plot_distribution_of_clones_across_clusters(seu, sample_id)

  # ggsave(glue("results/numbat_sridhar/{sample_id}/{sample_id}_clone_distribution{extension}.pdf"), width = 4, height = 4)

  seu_list <- list("markerplot" = markerplot, "dimplot" = dimplot, "distplot" = distplot)

  return(seu_list)
}

make_clustree_for_clone_comparison <- function(seu, sample_id, clone_set, mylabel = "sample_id", assay = "SCT", resolutions = seq(0.2, 1.0, by = 0.2), fisher_p_val_threshold = 0.1) {
	
	
	
  seu <- seu[, seu$scna %in% clone_set]

  resolutions <-
    glue("{assay}_snn_res.{resolutions}") %>%
    set_names(.)

  seu_meta <- 
  	seu@meta.data %>%
  	dplyr::mutate(across(where(is.character), ~na_if(.x, ""))) %>%
    dplyr::mutate(scna = replace_na(scna, "diploid")) %>%
    identity()

  speckle_proportions <- map(resolutions, ~ seu_meta[, c(.x, "scna")]) %>%
    map(set_names, c("cluster", "scna")) %>%
    map(janitor::tabyl, cluster, scna) %>%
    map(janitor::adorn_percentages) %>%
    map(dplyr::rename, samples = cluster) %>%
    imap(~ dplyr::mutate(.x, clusters = paste0(.y, "C", samples))) %>%
    dplyr::bind_rows() %>%
    janitor::clean_names() %>%
    na.omit() %>%
    identity()

  clustree_plot <- clustree::clustree(seu, assay = assay, show_axis = TRUE)

  clustree_meta <- seu@meta.data[, str_subset(colnames(seu@meta.data), glue("{assay}_snn.res.*"))]

  clustree_graph <- clustree:::build_tree_graph(
    clusterings = clustree_meta,
    prefix = glue("{assay}_snn.res."),
    metadata = clustree_meta,
    node_aes_list = list(colour = list(value = glue("{assay}_snn.res."), aggr = NULL), size = list(
      value = "size", aggr = NULL
    ), alpha = list(value = 1, aggr = NULL)),
    prop_filter = 0.1,
    count_filter = 0
  )

  from_res_col = as.character(glue("from_{assay}_snn.res."))
  to_res_col = as.character(glue("to_{assay}_snn.res."))
  
  daughter_clusters <-
  	clustree_graph %>%
  	tidygraph::activate(edges) %>%
  	data.frame() %>%
  	dplyr::group_by(.data[[from_res_col]], from_clust) %>%
  	dplyr::arrange(.data[[from_res_col]], from_clust, to_clust) %>%
  	dplyr::filter(n_distinct(to_clust) > 1) %>%
  	dplyr::mutate({{from_res_col}} := as.character(.data[[from_res_col]])) %>%
  	dplyr::mutate(from_clust = as.character(from_clust)) %>%
  	split(.[[from_res_col]]) %>%
  	map(~ split(.x, .x[[from_res_col]])) %>%
  	identity()

  for (resolution in names(daughter_clusters)) {
    for (from_clust in names(daughter_clusters[[resolution]])) {
      daughter_clusters[[resolution]][[from_clust]] <- chi_sq_daughter_clusters(seu, daughter_clusters, resolution = resolution, from_clust = from_clust, assay = assay)
    }
  }
  
  

  daughter_clusters <-
    daughter_clusters %>%
    map(dplyr::bind_rows) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(p.value < fisher_p_val_threshold) %>%
    # dplyr::group_by(clone_comparison) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(clone_comparison, p.value) %>%
    dplyr::slice_min(p.value, n = 20, by = clone_comparison) %>%
    dplyr::distinct(from_clust, to_clust, clone_comparison, .data[[from_res_col]], .keep_all = TRUE) %>%
    dplyr::select(all_of(c(
      "to_clust", to_res_col, "from_clust", from_res_col,
      "count", "in_prop", "clone_comparison", "p.value", "method"
    ))) %>%
    dplyr::group_by(from_clust, .data[[from_res_col]], clone_comparison) %>%
    dplyr::mutate(to_clust = paste(to_clust, collapse = "_")) %>%
    dplyr::distinct(from_clust, .data[[from_res_col]], clone_comparison, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from = "clone_comparison", values_from = "p.value") %>%
    identity()

  # clone_comparisons <- str_subset(colnames(daughter_clusters), "_v_") %>%
  #   set_names(.)
  
  clone_comparisons <-
  	paste0(clone_set, collapse= "_v_")

  res_col = glue("{assay}_snn_res.")
  
  # clustree_plot$data <-
  #   dplyr::left_join(clustree_plot$data, speckle_proportions, by = c("node" = "clusters")) %>%
  #   dplyr::left_join(daughter_clusters, by = c(res_col = from_res_col, "cluster" = "from_clust")) %>%
  #   # dplyr::mutate(signif = ifelse(is.na(method), 0, 1)) %>%
  #   identity()

  clustree_plot$data <-
  clustree_plot$data |> 
  	dplyr::left_join(speckle_proportions, by = c("node" = "clusters")) %>%
  	dplyr::left_join(daughter_clusters, by = join_by({{res_col}} == {{from_res_col}}, "cluster" == "from_clust")) |>
  	identity()
  
  clustree_plot$layers[[2]] <- NULL
  
  # browser()
  
  clustree_res <- map(clone_comparisons, plot_clustree_per_comparison, clustree_plot, speckle_proportions, sample_id)

  
  clustree_plots <- map(clustree_res, "plot")

  clustree_plot_path <- glue("results/{mylabel}_{clone_comparisons}_clustree.pdf")

  pdf(clustree_plot_path, width = 8, height = 10)
  print(clustree_plots)
  dev.off()

  clustree_table_path <- glue("results/{mylabel}_{clone_comparisons}_clustree.xlsx")

  clustree_tables <- map(clustree_res, "table") %>%
    purrr::flatten()

  writexl::write_xlsx(clustree_tables, clustree_table_path)

  return(clustree_plot_path)
}
