# Figure Assembly Functions
# Top-level functions that assemble multiple component plots into final paper figures.
# Sources: plot_functions_13, 14, 15, 16, 18, 25, 30

# --- Main figures ---

#' Create a plot visualization
#'
#' @param cluster_orders Cluster information
#' @param plot_path File path
#' @return ggplot2 plot object
#' @export
plot_fig_04_afterall <- function(cluster_orders, plot_path = "results/fig_02afterall.pdf"){
    integrated_1q_seus <- dir_ls("output/seurat/integrated_1q_16q/", regexp = ".*SR[RX][0-9].*_1q_filtered_seu.rds") |> as.character()

    batch_corrected_seu <- "output/seurat/integrated_1q_16q/integrated_seu_1q_afterall.rds"

    names(integrated_1q_seus) <- fs::path_ext_remove(fs::path_file(integrated_1q_seus))

    fig_2a_c <- make_clone_distribution_figure(batch_corrected_seu, cluster_orders,
                                               height = 12, width = 20, plot_path = "results/fig_2a_c.pdf", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"), heatmap_arrangement = c("clusters", "scna", "batch"))

    possibly_plot_clone_cc_plots <- possibly(plot_clone_cc_plots)
    fig_2d_cc <- map(c(batch_corrected_seu, integrated_1q_seus), plot_clone_cc_plots, var_y = "clusters", scna_of_interest = "1q", labeled_values = c("g2_m"))

    fig_2d_cc_table <- map(fig_2d_cc, ~.x$data) |>
        set_names(c("integrated", str_extract(integrated_1q_seus, "SR[RX][0-9]+"))) |>
        write_xlsx(path = "results/fig_2d_cc.xlsx")

    for(plot_i in 2:length(fig_2d_cc)){
        fig_2d_cc[[plot_i]] <- drop_y_axis(fig_2d_cc[[plot_i]])
    }

    fig_2d_cc <-
        fig_2d_cc |>
        wrap_plots() +
        plot_layout(
            nrow = 1,
            guides = "collect",
            axis_titles = "collect") +
        NULL

    fig_2d_cc_path <- ggsave("results/fig_2d_cc_afterall.pdf", fig_2d_cc, width = 12, height = 4)

    possibly_plot_clone_cc_plots <- possibly(plot_clone_cc_plots)

    fig_2d_hypoxia_hsp <- map(integrated_1q_seus, ~possibly_plot_clone_cc_plots(.x, scna_of_interest = "1q", kept_phases = c("hsp", "hypoxia"), labeled_values = c("hsp", "hypoxia")))

    fig_2d_hypoxia_hsp <-
        fig_2d_hypoxia_hsp |>
        wrap_plots() +
        plot_layout(
            nrow = 1,
            guides = "collect",
            axis_titles = "collect") +
        NULL

    fig_2d_hypoxia_hsp_path <- ggsave("results/fig_2d_hypoxia_hsp_afterall.pdf", fig_2d_hypoxia_hsp, width= 14, height = 6)

    qpdf::pdf_combine(
        list(
            fig_2a_c,
            fig_2d_cc_path,
            fig_2d_hypoxia_hsp_path
        ),
        plot_path)

    return(list("plot" = plot_path, "table" = fig_2d_cc_table))
}

#' Create a plot visualization
#'
#' @param cluster_orders Cluster information
#' @param plot_path File path
#' @return ggplot2 plot object
#' @export
plot_fig_03 <- function(cluster_orders, plot_path = "results/fig_03.pdf"){

	integrated_16q_seus <- dir_ls("output/seurat/integrated_16q/", regexp = ".*SR[RX][0-9].*rds") |>
		as.character()

	batch_corrected_seu <- "output/seurat/integrated_16q/integrated_seu_16q_complete.rds"

	names(integrated_16q_seus) <- fs::path_ext_remove(fs::path_file(integrated_16q_seus))

	fig_3a_c <- make_clone_distribution_figure("output/seurat/integrated_16q/integrated_seu_16q_complete.rds", cluster_orders,
																						 height = 12, width = 18, plot_path = "results/fig_3a_c.pdf",
																						 heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"), heatmap_arrangement = c("clusters", "scna", "batch"))

	fig_3d_cc <- map(c(batch_corrected_seu, integrated_16q_seus), plot_clone_cc_plots, scna_of_interest = "16q", var_y = "clusters", labeled_values = c("pm"))

	fig_3d_cc_table <- map(fig_3d_cc, ~.x$data) |>
		set_names(c("integrated", str_extract(integrated_16q_seus, "SR[RX][0-9]+"))) |>
		write_xlsx(path = "results/fig_3d_cc.xlsx")

	for(plot_i in 2:length(fig_3d_cc)){
		fig_3d_cc[[plot_i]] <- drop_y_axis(fig_3d_cc[[plot_i]])
	}

	fig_3d_cc <-
		fig_3d_cc |>
		wrap_plots() +
		plot_layout(
			nrow = 1,
			guides = "collect",
			axis_titles = "collect") +
		NULL

	fig_3d_cc_path <- ggsave("results/fig_3d_cc.pdf", fig_3d_cc, width = 12, height = 6)

	qpdf::pdf_combine(
		list(
			fig_3a_c,
			fig_3d_cc_path
		),
		plot_path)

	return(list("plot" = plot_path, "table" = fig_3d_cc_table))
}

#' Create a plot visualization
#'
#' @param cluster_orders Cluster information
#' @param plot_path File path
#' @return ggplot2 plot object
#' @export
plot_fig_03_afterall <- function(cluster_orders, plot_path = "results/fig_03_afterall.pdf"){

    integrated_16q_seus <- dir_ls("output/seurat/integrated_1q_16q/", regexp = ".*SR[RX][0-9].*_16q_filtered_seu.rds") |>
        as.character()

    batch_corrected_seu <- "output/seurat/integrated_1q_16q/integrated_seu_16q_afterall.rds"

    var_y_levels <- levels(readRDS(batch_corrected_seu)$clusters)

    names(integrated_16q_seus) <- fs::path_ext_remove(fs::path_file(integrated_16q_seus))

    fig_3a_c <- make_clone_distribution_figure(batch_corrected_seu, cluster_orders,
                                               height = 12, width = 18, plot_path = "results/fig_3a_c.pdf",
                                               heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"), heatmap_arrangement = c("clusters", "scna", "batch"))

    fig_3d_cc <- map(c(batch_corrected_seu, integrated_16q_seus), plot_clone_cc_plots, scna_of_interest = "16q", var_y = "clusters", labeled_values = c("pm"), var_y_levels = var_y_levels)

    fig_3d_cc_table <- map(fig_3d_cc, ~.x$data) |>
        set_names(c("integrated", str_extract(integrated_16q_seus, "SR[RX][0-9]+"))) |>
        write_xlsx(path = "results/fig_3d_cc_afterall.xlsx")

    for(plot_i in 2:length(fig_3d_cc)){
        fig_3d_cc[[plot_i]] <- drop_y_axis(fig_3d_cc[[plot_i]])
    }

    fig_3d_cc <-
        fig_3d_cc |>
        wrap_plots() +
        plot_layout(
            nrow = 1,
            guides = "collect",
            axis_titles = "collect") +
        NULL

    fig_3d_cc_path <- ggsave("results/fig_3d_cc_afterall.pdf", fig_3d_cc, width = 12, height = 6)

    qpdf::pdf_combine(
        list(
            fig_3a_c,
            fig_3d_cc_path
        ),
        plot_path)

    return(list("plot" = plot_path, "table" = fig_3d_cc_table))
}

#' Create a plot visualization
#'
#' @param cluster_orders Cluster information
#' @param plot_path File path
#' @return ggplot2 plot object
#' @export
plot_fig_04_07 <- function(cluster_orders, plot_path = "results/fig_04_07.pdf"){

	fig_3a_c <- make_clone_distribution_figure("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", cluster_orders,
																						 height = 12, width = 18, plot_path = tempfile(fileext = ".pdf"),
																						 heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"), heatmap_arrangement = c("clusters", "scna", "batch"))

	qpdf::pdf_combine(
		list(
			fig_3a_c
		),
		plot_path)

	return("plot" = plot_path)
}

#' Create a plot visualization
#'
#' @param seu_paths File path
#' @param integrated_enrichment Parameter for integrated enrichment
#' @param plot_path File path
#' @param integrated_seu_paths File path
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_fig_04_05 <- function(seu_paths, integrated_enrichment, plot_path = "results/fig_04.pdf", integrated_seu_paths, ...){
	panels_a_d <- map(seu_paths, plot_fig_04_05_panels , ...)

	panels_e <-
		integrated_enrichment$enrichment[[1]] |>
		dplyr::filter(!str_detect(ID, "UP")) |>
		dplyr::filter(!str_detect(ID, "DN")) |>
		plot_enrichment(p_val_cutoff = 0.1) +
		facet_grid (.~.sign, scales = "free_y", space = "free_y")

	panel_e_file <- ggsave(tempfile(fileext = ".pdf"), width = 10, height = 3)

	possibly_plot_clone_cc_plots <- possibly(plot_clone_cc_plots)
	panels_f <- map(c(seu_paths, integrated_seu_paths), plot_clone_cc_plots, var_y = "clusters")

	panel_f <-
		panels_f |>
		wrap_plots() +
		plot_layout(
			nrow = 1,
			guides = "collect",
			axis_titles = "collect") +
		NULL

	panel_f_file <- ggsave(tempfile(fileext = ".pdf"), panel_f, width = 14, height = 6)

	paths <- unlist(c(panels_a_d, panel_e_file, panel_f_file))
	paths <- paths[!is.na(paths)]
	qpdf::pdf_combine(paths, plot_path)

	return(plot_path)
}

#' Create a plot visualization
#'
#' @param figure_input Parameter for figure input
#' @param x_var Character string (default: "sample_cluster")
#' @param plot_title Plot title
#' @param p_adj_threshold Threshold value for filtering
#' @param plot_path File path
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_fig_07_08 <- function(figure_input, x_var = "sample_cluster", plot_title = "fig_07", p_adj_threshold = 0.1, plot_path = "results/fig_07.pdf", ...){
	unfiltered_input <-
		unlist(figure_input) |>
		set_names(str_extract(unlist(figure_input), "SR[RX][0-9]+")) |>
		map(read_csv) |>
		purrr::keep(~ nrow(.x) > 0 && "location" %in% names(.x)) |>
		map(~split(.x, .x$cluster)) |>
		purrr::list_transpose() |>
		map(~purrr::compact(.x)) |>
		map(dplyr::bind_rows, .id = "sample_id") |>
		dplyr::bind_rows(.id = "clusters") |>
		(\(df) if ("location" %in% names(df)) dplyr::filter(df, location == "cis") else df)()

	plot_input <-
		unfiltered_input |>
		dplyr::filter(p_val_adj <= p_adj_threshold) |>
		dplyr::group_by(symbol) %>%
		dplyr::filter(all_same_sign(avg_log2FC)) |>
		dplyr::mutate(abs_mean_FC = abs(mean(avg_log2FC))) %>%
		dplyr::mutate(recurrence = dplyr::n_distinct(sample_id)) %>%
		dplyr::arrange(desc(abs_mean_FC)) |>
		dplyr::mutate(sample_cluster = glue("{sample_id}_{cluster}")) |>
		dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
		identity()

	myplot <- ggplot(plot_input, aes(y = symbol, x = .data[[x_var]], size = neg_log_p_val_adj, color = avg_log2FC)) +
		geom_point() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
		scale_color_gradient2() +
		labs(title = plot_title)

	pdf(plot_path, ...)
	print(myplot)
	plot(unfiltered_input$avg_log2FC, unfiltered_input$p_val_adj,
			 pch = 16)
	abline(h = p_adj_threshold, col = "red")
	dev.off()

	return(plot_path)
}

#' Create a plot visualization
#'
#' @param seu_subset Parameter for seu subset
#' @param corresponding_seus Parameter for corresponding seus
#' @param corresponding_clusters_diffex Parameter for corresponding clusters diffex
#' @param corresponding_clusters_enrichments Parameter for corresponding clusters enrichments
#' @param plot_path File path
#' @param widths Parameter for widths
#' @param heights Parameter for heights
#' @param common_seus Parameter for common seus
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
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

# --- Supplementary figures ---

#' Assemble sample-specific supplementary figure (fig_s04)
#'
#' @param seu_list Parameter for seu list
#' @param cluster_orders Cluster information
#' @param plot_path File path
#' @return File path of assembled PDF
#' @export
sample_specific_analyses_of_tumors_with_scna_subclones_after_integration <- function(seu_list, cluster_orders, plot_path = "results/fig_s04.pdf"){

	plot_path <- seu_list |>
		imap(~make_clone_distribution_figure(.x, cluster_orders, height = 12, width = 20, plot_path = tempfile(pattern = .y, fileext = ".pdf"), heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"))) |>
		qpdf::pdf_combine(plot_path)

	return(plot_path)
}

#' Assemble supplementary figure s07 (clone distribution + cc plots)
#'
#' @param seu_list Parameter for seu list
#' @param cluster_orders Cluster information
#' @param plot_path File path
#' @return File path of assembled PDF
#' @export
not_sure_what_this_does <- function(seu_list, cluster_orders, plot_path = "results/fig_s07.pdf"){

	clone_dist_paths <- seu_list |>
		imap(~make_clone_distribution_figure(.x, cluster_orders, height = 12, width = 20, plot_path = tempfile(pattern = .y, fileext = ".pdf"), heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch")))

	possibly_plot_clone_cc_plots <- possibly(plot_clone_cc_plots)
	fig_2d_cc <- map(seu_list, plot_clone_cc_plots, var_y = "clusters", scna_of_interest = "2p", labeled_values = c("g2_m")) |>
		map(~ggsave(tempfile(fileext=".pdf"), .x, w = 5, h = 4))

	plot_path <- qpdf::pdf_combine(c(clone_dist_paths, fig_2d_cc), plot_path)

	return(plot_path)
}

#' Create supplementary figure s06a (SCNA ideograms for all samples)
#'
#' @param plot_path File path
#' @param table_path File path
#' @return List of plot and table paths
#' @export
plot_fig_s06a <- function(plot_path = "results/fig_s06a.pdf", table_path = "results/table_s12.xlsx"){

	nb_paths <- dir_ls("output/numbat_sridhar/", regexp = ".*SR[RX][0-9]+_numbat.rds", recurse = TRUE) |>
		sort()

	ideogram_res <- map(nb_paths, make_rb_scna_ideograms)

	ideogram_tables <-
		ideogram_res |>
		map("table")

	ideogram_tables |>
		set_names(str_extract(names(ideogram_tables), "SR[RX][0-9]+")) |>
		writexl::write_xlsx(path = table_path)

	ideogram_plots <-
		ideogram_res |>
		map("plot") |>
		qpdf::pdf_combine(plot_path)

	return(list(plot_path, table_path))
}

#' Create supplementary figure s08 (1q diffex per cluster)
#'
#' @param plot_path File path
#' @param table_path File path
#' @return List of plot and table paths
#' @export
plot_fig_s08 <- function(plot_path = "results/fig_s08.pdf", table_path = "results/table_s22.csv") {
	segs_1q <- sapply(c("output/numbat_sridhar/SRX10264526_numbat.rds",
											"output/numbat_sridhar/SRX11133594_numbat.rds",
											"output/numbat_sridhar/SRX11133593_numbat.rds",
											"output/numbat_sridhar/SRX11133592_numbat.rds"), pull_scna_segments, chrom = "1") |>
		as("GRangesList") |>
		unlist() |>
		reduce_ranges()

	seu <- readRDS("output/seurat/integrated_1q/integrated_seu_1q_complete.rds")

	myclusters <- seu$clusters |>
		levels() |>
		set_names(identity)

	diffex_1q <- map(myclusters, ~diffex_per_cluster(seu, mycluster = .x, segments = segs_1q)) |>
		dplyr::bind_rows(.id = "cluster")

	diffex_1q <-
		diffex_1q |>
		dplyr::filter(p_val_adj <= 0.1) |>
		dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
		dplyr::filter(avg_log2FC >= 0.5) |>
		dplyr::filter(!cluster %in% c("hsp_8", "hypoxia_2")) |>
		dplyr::arrange(desc(cluster), desc(avg_log2FC)) |>
		dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) |>
		identity()

	plot_input <-
	    diffex_1q |>
	    dplyr::mutate(cluster = factor(cluster, levels = c("g1_3", "g1_0", "g1_6", "s_1", "s_g2_4"))) |>
	    dplyr::arrange(desc(cluster), desc(avg_log2FC)) |>
	    dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) |>
		dplyr::filter(location == "cis")

	plot_input |>
	ggplot(aes(y = symbol, x = cluster, size = neg_log_p_val_adj, color = avg_log2FC)) +
		geom_point() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
		scale_color_gradient2() +
		labs(title = "integrated 1q+ diffex per cluster")

	cis_path = "results/fig_s08a.pdf"
	ggsave(cis_path, w = 4, h = 7) |>
		browseURL()

	all_diffex <-
		diffex_1q |>
		dplyr::filter(!is.na(entrez)) |>
		dplyr::mutate(Gene.Name = symbol) |>
		dplyr::select(-any_of(colnames(annotables::grch38))) |>
		split(~cluster) |>
		map(tibble::column_to_rownames, "Gene.Name") |>
		identity()

	enrichments <-
		map(all_diffex, enrichment_analysis, analysis_method = "ora") |>
		purrr::compact() |>
		map(DOSE::setReadable, org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")

	enrichment_plots <-
		enrichments |>
		map(plot_enrichment, analysis_method = "ora") |>
		imap(~{.x + labs(title = .y)})

	trans_path <- "results/fig_s08b.pdf"
	pdf(trans_path, h = 4, w = 6)
	print(enrichment_plots)
	dev.off()

	qpdf::pdf_combine(c(cis_path, trans_path), plot_path)

	mytable <-
	enrichments |>
		map("result") |>
		dplyr::bind_rows(.id = "cluster") |>
		label_enrichment_by_cis("01", "q") |>
		write_csv(table_path)

	return(list("plot" = plot_path, "table" = table_path))
}

#' Create supplementary figure s09 (1q cluster comparisons: enrichments)
#'
#' @return List of diffex table path and enrichment plot path
#' @export
plot_fig_s09 <- function(){

	seu <- readRDS("output/seurat/integrated_1q/integrated_seu_1q_complete.rds")

	comparisons <-
		tibble::tribble(
			~ident.1, ~ident.2,
			"g1_6", "g1_3",
			"g1_6", "g1_0",
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

#' Create supplementary figure s10 (16q diffex per cluster)
#'
#' @param plot_path File path
#' @param table_path File path
#' @return List of plot and table paths
#' @export
plot_fig_s10 <- function(plot_path = "results/fig_s10.pdf", table_path = "results/table_s23.csv") {
	segs_16q <- sapply(c("output/numbat_sridhar/SRX11133594_numbat.rds",
											"output/numbat_sridhar/SRX11133593_numbat.rds",
											"output/numbat_sridhar/SRX11133592_numbat.rds"), pull_scna_segments, chrom = "16") |>
		as("GRangesList") |>
		unlist() |>
		reduce_ranges()

	seu <- readRDS("output/seurat/integrated_16q/integrated_seu_16q_complete.rds")

	myclusters <- seu$clusters |>
		levels() |>
		set_names(identity)

	diffex_16q <- map(myclusters, ~diffex_per_cluster(seu, "16q-", ".diploid", mycluster = .x, segments = segs_16q)) |>
		dplyr::bind_rows(.id = "cluster")

	diffex_16q <-
		diffex_16q |>
		dplyr::filter(p_val_adj <= 0.1) |>
		dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
		dplyr::filter(abs(avg_log2FC) >= 0.5) |>
		dplyr::arrange(desc(cluster), desc(avg_log2FC)) |>
		dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) |>
		identity()

	diffex_16q |>
		dplyr::filter(location == "cis") |>
		ggplot(aes(y = symbol, x = cluster, size = neg_log_p_val_adj, color = avg_log2FC)) +
		geom_point() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
		scale_color_gradient2() +
		labs(title = "integrated 16q+ diffex per cluster")

	cis_path = "results/fig_s10a.pdf"
	ggsave(cis_path, w = 4, h = 3) |>
		browseURL()

	all_diffex <-
		diffex_16q |>
		dplyr::filter(!is.na(entrez)) |>
		dplyr::mutate(Gene.Name = symbol) |>
		dplyr::select(-any_of(colnames(annotables::grch38))) |>
		split(~cluster) |>
		map(tibble::column_to_rownames, "Gene.Name") |>
		identity()

	enrichments <-
		map(all_diffex, enrichment_analysis, analysis_method = "ora") |>
		purrr::compact() |>
		map(DOSE::setReadable, org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")

	enrichment_plots <-
		enrichments |>
		map(plot_enrichment, analysis_method = "ora") |>
		imap(~{.x + labs(title = .y)})

	trans_path <- "results/fig_s10b.pdf"
	pdf(trans_path, h = 4, w = 6)
	print(enrichment_plots)
	dev.off()

	plot_path <- qpdf::pdf_combine(c(cis_path, trans_path), plot_path)

	mytable <-
		enrichments |>
		map("result") |>
		dplyr::bind_rows(.id = "cluster") |>
		label_enrichment_by_cis("16", "q") |>
		write_csv(table_path)

	return(list("plot" = plot_path, "table" = table_path))
}

#' Create supplementary figure s20 (2p diffex per cluster)
#'
#' @return List of plot and table paths
#' @export
plot_fig_s20 <- function() {
	segs_2p <- sapply(c("output/numbat_sridhar/SRX10264525_numbat.rds",
											 "output/numbat_sridhar/SRX14116944_numbat.rds"), pull_scna_segments, chrom = "2") |>
		as("GRangesList") |>
		unlist() |>
		reduce_ranges()

	seu <- readRDS("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds")

	myclusters <- seu$clusters |>
		levels() |>
		set_names(identity)

	diffex_2p <- map(myclusters, ~diffex_per_cluster(seu, "w_scna", "wo_scna", mycluster = .x, segments = segs_2p)) |>
		dplyr::bind_rows(.id = "cluster")

	janitor::tabyl(seu@meta.data, clusters, scna) |>
		adorn_percentages(denominator = "row") |>
		adorn_pct_formatting() |>
		identity()

	diffex_2p <-
		diffex_2p |>
		dplyr::filter(p_val_adj <= 0.05) |>
		dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
		dplyr::filter(avg_log2FC >= 0.5) |>
		dplyr::arrange(desc(cluster), desc(avg_log2FC)) |>
		dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) |>
		identity()

	diffex_2p |>
		dplyr::group_by(symbol) |>
		dplyr::mutate(recurrence = n_distinct(cluster)) |>
		dplyr::filter(recurrence >=4) |>
		dplyr::filter(location == "cis") |>
		ggplot(aes(y = symbol, x = cluster, size = neg_log_p_val_adj, color = avg_log2FC)) +
		geom_point() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
		scale_color_gradient2() +
		labs(title = "integrated 2p+ diffex per cluster")

	cis_path = "results/fig_s20a.pdf"
	ggsave(cis_path, w = 4, h = 5) |>
		browseURL()

	all_diffex <-
		diffex_2p |>
		dplyr::filter(!is.na(entrez)) |>
		dplyr::mutate(Gene.Name = symbol) |>
		dplyr::select(-any_of(colnames(annotables::grch38))) |>
		split(~cluster) |>
		map(tibble::column_to_rownames, "Gene.Name") |>
		identity()

	enrichments <-
		map(all_diffex, enrichment_analysis, analysis_method = "ora") |>
		purrr::compact() |>
		map(DOSE::setReadable, org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")

	enrichment_plots <-
		enrichments |>
		map(plot_enrichment, analysis_method = "ora") |>
		imap(~{.x + labs(title = .y)})

	trans_path <- "results/fig_s20b.pdf"
	pdf(trans_path, h = 4, w = 6)
	print(enrichment_plots)
	dev.off()

	plot_path <- qpdf::pdf_combine(c(cis_path, trans_path), "results/fig_s20.pdf")

	table_path <- "results/table_s24.csv"
	mytable <-
		enrichments |>
		map("result") |>
		dplyr::bind_rows(.id = "cluster") |>
		label_enrichment_by_cis("02", "p") |>
		write_csv(table_path)

	return(list("plot" = plot_path, "table" = table_path))
}

# --- Multi-panel collage assembly helpers ---

#' Assemble a per-sample figure collage (heatmap + CC scatter + clone tree + distribution)
#'
#' @param seu_path File path to filtered Seurat object RDS
#' @param cluster_order Named list of cluster order data frames
#' @param nb_paths Named character vector of numbat RDS paths
#' @param clone_simplifications Clone simplification mappings
#' @param group.by Seurat metadata column for clustering (default: "SCT_snn_res.0.6")
#' @param assay Seurat assay (default: "SCT")
#' @param height PDF height in inches
#' @param width PDF width in inches
#' @param equalize_scna_clones Downsample clones to equal size
#' @param phase_levels Ordered vector of phase level labels
#' @param kept_phases Subset of phase_levels to retain
#' @param rb_scna_samples Samples with RB SCNA
#' @param large_clone_comparisons Named list of clone comparison specs
#' @param scna_of_interest SCNA arm label (default: "1q")
#' @param min_cells_per_cluster Minimum cells required per cluster
#' @param return_plots Return plot objects instead of saving
#' @param split_columns Metadata column to split heatmap columns
#' @return File path of output PDF
#' @export
plot_figure_collage <- function(seu_path = NULL, cluster_order = NULL, nb_paths = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", height = 10, width = 18, equalize_scna_clones = FALSE, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL, rb_scna_samples, large_clone_comparisons, scna_of_interest = "1q", min_cells_per_cluster = 50, return_plots = FALSE, split_columns = "clusters") {
	kept_phases <- kept_phases %||% phase_levels

	file_id <- fs::path_file(seu_path)

	tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

	message(file_id)
	cluster_order <- cluster_order[[file_id]]

	full_seu <- readRDS(seu_path)

	# subset by retained clones ------------------------------
	clone_comparisons <- names(large_clone_comparisons[[sample_id]])
	clone_comparison <- clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
	retained_clones <- clone_comparison %>%
		str_extract("[0-9]_v_[0-9]") %>%
		str_split("_v_", simplify = TRUE)

	nb_paths <- nb_paths %>%
		set_names(str_extract(., "SR[RX][0-9]+"))

	nb_path <- nb_paths[[tumor_id]]

	plot_paths <- vector(mode = "list", length = length(cluster_order))
	names(plot_paths) <- names(cluster_order)

	file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
	plot_path <- glue("results/{file_slug}_{scna_of_interest}_heatmap_phase_scatter_patchwork.pdf")

	pdf(plot_path, height = height, width = width)

	for (resolution in names(cluster_order)) {
		# start loop ------------------------------

		seu <- full_seu[, full_seu$clone_opt %in% retained_clones]

		if (!is.null(cluster_order)) {
			single_cluster_order <- cluster_order[[resolution]]

			group.by <- unique(single_cluster_order$resolution)
		}

		if (equalize_scna_clones) {
			seu_meta <- seu@meta.data %>%
				tibble::rownames_to_column("cell")

			clones <- table(seu_meta$scna)

			min_clone_num <- clones[which.min(clones)]

			selected_cells <-
				seu_meta %>%
				dplyr::group_by(scna) %>%
				slice_sample(n = min_clone_num) %>%
				pull(cell)

			seu <- seu[, selected_cells]
		}

		if (!is.null(single_cluster_order)) {
			single_cluster_order <-
				single_cluster_order |>
				dplyr::mutate(order = dplyr::row_number()) %>%
				dplyr::filter(!is.na(clusters)) %>%
				dplyr::mutate(clusters = as.character(clusters))

			group.by <- unique(single_cluster_order$resolution)

			seu@meta.data$clusters <- seu@meta.data[[group.by]]

			seu_meta <- seu@meta.data %>%
				tibble::rownames_to_column("cell") %>%
				dplyr::select(-any_of(c("phase_level", "order"))) %>%
				dplyr::left_join(single_cluster_order, by = "clusters") %>%
				dplyr::select(-clusters) %>%
				dplyr::rename(phase_level = phase) %>%
				identity()

			phase_levels <- phase_levels[phase_levels %in% unique(seu_meta$phase_level)]

			seu_meta <-
				seu_meta %>%
				tidyr::unite("clusters", all_of(c("phase_level", group.by)), remove = FALSE) %>%
				dplyr::arrange(phase_level, order) %>%
				dplyr::mutate(clusters = factor(clusters, levels = unique(clusters))) %>%
				tibble::column_to_rownames("cell") %>%
				identity()

			seu@meta.data <- seu_meta[rownames(seu@meta.data), ]

			seu <- seu[, seu$phase_level %in% kept_phases]
			seu <- tryCatch(
				find_all_markers(seu, metavar = "clusters", seurat_assay = "SCT"),
				error = function(e) {
					if (grepl("JoinLayers", conditionMessage(e), fixed = TRUE)) {
						warning("SCT marker JoinLayers failed; using stash_marker_features fallback.")
						seu@misc$markers[["clusters"]] <- seuratTools:::stash_marker_features("clusters", seu, seurat_assay = "SCT")
						return(seu)
					}
					stop(e)
				}
			)

			seu@meta.data$clusters <- forcats::fct_drop(seu@meta.data$clusters)

			heatmap_features <-
				seu@misc$markers[["clusters"]][["presto"]] %>%
				dplyr::filter(Gene.Name %in% VariableFeatures(seu))

			tidy_eval_arrange <- function(.data, ...) {
				.data %>%
					arrange(...)
			}

			single_cluster_order_vec <-
				seu@meta.data %>%
				dplyr::select(clusters, !!group.by) %>%
				dplyr::arrange(clusters, !!sym(group.by)) %>%
				dplyr::select(clusters, !!group.by) |>
				dplyr::distinct(.data[[group.by]], .keep_all = TRUE) |>
				dplyr::mutate(!!group.by := as.character(.data[[group.by]])) |>
				tibble::deframe() |>
				identity()

			heatmap_features[["Cluster"]] <-
				factor(heatmap_features[["Cluster"]], levels = levels(seu_meta$clusters))

			heatmap_features <-
				heatmap_features %>%
				dplyr::group_by(Gene.Name) |>
				dplyr::slice_max(order_by = Average.Log.Fold.Change, n = 1) |>
				dplyr::ungroup() |>
				dplyr::arrange(Cluster, desc(Average.Log.Fold.Change)) |>
				group_by(Cluster) %>%
				slice_head(n = 5) %>%
				dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
				dplyr::ungroup() %>%
				dplyr::distinct(Gene.Name, .keep_all = TRUE) |>
				identity()
		} else {
			heatmap_features <-
				seu@misc$markers[[group.by]][["presto"]]

			single_cluster_order <- levels(seu@meta.data[[group.by]]) %>%
				set_names(.)

			seu@meta.data[[group.by]] <-
				factor(seu@meta.data[[group.by]], levels = single_cluster_order)

			group_by_clusters <- seu@meta.data[[group.by]]

			seu@meta.data$clusters <- names(single_cluster_order[group_by_clusters])

			seu@meta.data$clusters <- factor(seu@meta.data$clusters, levels = unique(setNames(names(single_cluster_order), single_cluster_order)[levels(seu@meta.data[[group.by]])]))

			heatmap_features <-
				heatmap_features %>%
				dplyr::arrange(Cluster) %>%
				group_by(Cluster) %>%
				slice_head(n = 6) %>%
				dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
				dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
				identity()
		}

		large_enough_clusters <-
			seu@meta.data %>%
			dplyr::group_by(clusters) %>%
			dplyr::count() |>
			dplyr::filter(n >= min_cells_per_cluster) %>%
			dplyr::pull(clusters)

		seu <- seu[, seu$clusters %in% large_enough_clusters]

		seu$scna[seu$scna == ""] <- ".diploid"
		seu$scna <- factor(seu$scna)

		giotti_genes <- read_giotti_genes()

		heatmap_features <-
			heatmap_features %>%
			dplyr::ungroup() %>%
			left_join(giotti_genes, by = c("Gene.Name" = "symbol")) %>%
			dplyr::mutate(term = replace_na(term, "")) %>%
			dplyr::distinct(Gene.Name, .keep_all = TRUE)

		row_ha <- ComplexHeatmap::rowAnnotation(term = rev(heatmap_features$term))

		if (!is.null(split_columns)) {
			column_split <- sort(seu@meta.data[[split_columns]])
			column_title <- unique(column_split)
		} else {
			column_split <- split_columns
			column_title <- NULL
		}

		seu_heatmap <- ggplotify::as.ggplot(
			seu_complex_heatmap(seu,
													features = heatmap_features$Gene.Name,
													group.by = c("G2M.Score", "S.Score", "scna", "clusters"),
													col_arrangement = c("clusters", "scna"),
													cluster_rows = FALSE,
													column_split = column_split,
													row_split = rev(heatmap_features$Cluster),
													row_title_rot = 0,
													column_title = column_title,
													column_title_rot = 90
			)
		) +
			labs(title = sample_id) +
			theme()

		labels <- data.frame(clusters = unique(seu[[]][["clusters"]]), label = unique(seu[[]][["clusters"]])) %>%
			identity()

		cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))

		centroid_data <-
			cc_data %>%
			dplyr::group_by(clusters) %>%
			dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
			dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
			dplyr::mutate(centroid = "centroids") %>%
			identity()

		centroid_plot <-
			cc_data %>%
			ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["clusters"]])) +
			geom_point(size = 0.1) +
			theme_light() +
			theme(
				strip.background = element_blank(),
				strip.text.x = element_blank()
			) +
			geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
			guides(fill = "none", color = "none") +
			NULL

		facet_cell_cycle_plot <-
			cc_data %>%
			ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
			geom_point(size = 0.1) +
			geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
			facet_wrap(~ .data[["clusters"]], ncol = 2) +
			theme_light() +
			geom_label(
				data = labels,
				aes(label = label),
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
			NULL

		appender <- function(string) str_wrap(string, width = 40)

		labels <- data.frame(scna = unique(seu$scna), label = str_replace(unique(seu$scna), "^$", "diploid"))

		clone_ratio <- janitor::tabyl(as.character(seu$scna))$percent[[2]]

		comparison_scna <-
			janitor::tabyl(as.character(seu$scna))[2, 1]

		clone_distribution_plot <- plot_distribution_of_clones_across_clusters(
			seu,
			seu_name = glue("{tumor_id} {comparison_scna}"), var_x = "scna", var_y = "clusters", signif = TRUE, plot_type = "clone"
		) +
			coord_flip()

		full_seu@meta.data[[group.by]] <- factor(full_seu@meta.data[[group.by]], levels = single_cluster_order_vec)
		levels(full_seu@meta.data[[group.by]]) <- names(single_cluster_order_vec)
		umap_plots <- make_faded_umap_plots(full_seu, retained_clones, group_by = group.by)

		if (!is.null(nb_path)) {
			clone_tree_plot <-
				plot_clone_tree(seu, tumor_id, nb_path, clone_simplifications, sample_id = sample_id, legend = FALSE, horizontal = FALSE)

			collage_plots <- list(
				"A" = seu_heatmap,
				"B" = facet_cell_cycle_plot,
				"C" = clone_tree_plot,
				"D" = centroid_plot,
				"E" = plot_spacer(),
				"F" = plot_spacer(),
				"G" = plot_spacer(),
				"H" = clone_distribution_plot
			)

			layout <- "
              AAAAAACDEEGG
              AAAAAABBEEGG
              AAAAAABBFFGG
              AAAAAABBFFGG
              AAAAAAHHFFGG
              "

			plot_collage <- wrap_plots(collage_plots) +
				plot_layout(design = layout) +
				plot_annotation(tag_levels = "A") +
				NULL
		} else {
			layout <- "
              AAAAAAAAAABBBBCCCC
              AAAAAAAAAABBBBCCCC
              AAAAAAAAAABBBBDDDD
              AAAAAAAAAABBBBDDDD
              AAAAAAAAAABBBBDDDD
      "

			collage_plots <- list(
				"seu_heatmap" = seu_heatmap,
				"facet_cell_cycle_plot" = facet_cell_cycle_plot,
				"centroid_plot" = centroid_plot,
				"clone_distribution_plot" = clone_distribution_plot
			)

			plot_collage <- wrap_plots(collage_plots) +
				plot_layout(design = layout) +
				plot_annotation(tag_levels = "A") +
				NULL
		}

		print(plot_collage)
		# end loop------------------------------
	}

	dev.off()

	return(plot_path)
}

#' Assemble per-sample figure panels (heatmap + CC scatter + clone tree)
#'
#' @param seu_path File path to filtered Seurat object RDS
#' @param cluster_order Named list of cluster order data frames
#' @param nb_paths Named character vector of numbat RDS paths
#' @param clone_simplifications Clone simplification mappings
#' @param group.by Seurat metadata column for clustering (default: "SCT_snn_res.0.6")
#' @param assay Seurat assay (default: "SCT")
#' @param height PDF height in inches
#' @param width PDF width in inches
#' @param equalize_scna_clones Downsample clones to equal size
#' @param phase_levels Ordered vector of phase level labels
#' @param kept_phases Subset of phase_levels to retain
#' @param rb_scna_samples Samples with RB SCNA
#' @param large_clone_comparisons Named list of clone comparison specs
#' @param scna_of_interest SCNA arm label (default: "1q")
#' @param min_cells_per_cluster Minimum cells required per cluster
#' @param return_plots Return plot objects instead of saving
#' @param split_columns Metadata column to split heatmap columns
#' @param heatmap_groups Metadata columns to show as heatmap annotations
#' @param heatmap_arrangement Metadata columns used to arrange heatmap columns
#' @return File path of output PDF
#' @export
plot_fig_04_05_panels <- function(seu_path = NULL, cluster_order = NULL, nb_paths = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", height = 10, width = 14, equalize_scna_clones = FALSE, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL, rb_scna_samples, large_clone_comparisons, scna_of_interest = "1q", min_cells_per_cluster = 50, return_plots = FALSE, split_columns = "clusters", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"), heatmap_arrangement = c("clusters", "scna", "batch")) {
	kept_phases <- kept_phases %||% phase_levels

	file_id <- fs::path_file(seu_path)

	tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

	message(file_id)
	cluster_order <- cluster_order[[file_id]]

	full_seu <- readRDS(seu_path)

	if(!is.null(nb_paths)){
		nb_paths <- nb_paths %>%
			set_names(str_extract(., "SR[RX][0-9]+"))

		nb_path <- nb_paths[[tumor_id]]
	}

	plot_paths <- vector(mode = "list", length = length(cluster_order))
	names(plot_paths) <- names(cluster_order)

	file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
	plot_path <- glue("results/{file_slug}_{scna_of_interest}_heatmap_phase_scatter_patchwork.pdf")

	pdf(plot_path, height = height, width = width)

	for (resolution in names(cluster_order)) {
		# start loop ------------------------------
		seu <- full_seu

		if (!is.null(cluster_order)) {
			single_cluster_order <- cluster_order[[resolution]]

			group.by <- unique(single_cluster_order$resolution)
		}

		if (equalize_scna_clones) {
			seu_meta <- seu@meta.data %>%
				tibble::rownames_to_column("cell")

			clones <- table(seu_meta$scna)

			min_clone_num <- clones[which.min(clones)]

			selected_cells <-
				seu_meta %>%
				dplyr::group_by(scna) %>%
				slice_sample(n = min_clone_num) %>%
				pull(cell)

			seu <- seu[, selected_cells]
		}

		if (!is.null(single_cluster_order)) {
			single_cluster_order <-
				single_cluster_order |>
				dplyr::mutate(order = dplyr::row_number()) %>%
				dplyr::filter(!is.na(clusters)) %>%
				dplyr::mutate(clusters = as.character(clusters))

			group.by <- unique(single_cluster_order$resolution)

			seu@meta.data$clusters <- seu@meta.data[[group.by]]

			seu_meta <- seu@meta.data %>%
				tibble::rownames_to_column("cell") %>%
				dplyr::select(-any_of(c("phase_level", "order"))) %>%
				dplyr::left_join(single_cluster_order, by = "clusters") %>%
				dplyr::select(-clusters) %>%
				dplyr::rename(phase_level = phase) %>%
				identity()

			phase_levels <- phase_levels[phase_levels %in% unique(seu_meta$phase_level)]

			seu_meta <-
				seu_meta %>%
				tidyr::unite("clusters", all_of(c("phase_level", group.by)), remove = FALSE) %>%
				dplyr::arrange(phase_level, order) %>%
				dplyr::mutate(clusters = factor(clusters, levels = unique(clusters))) %>%
				tibble::column_to_rownames("cell") %>%
				identity()

			seu@meta.data <- seu_meta[rownames(seu@meta.data), ]

			seu <- seu[, seu$phase_level %in% kept_phases]
			seu <- tryCatch(
				find_all_markers(seu, metavar = "clusters", seurat_assay = "SCT"),
				error = function(e) {
					if (grepl("JoinLayers", conditionMessage(e), fixed = TRUE)) {
						warning("SCT marker JoinLayers failed; using stash_marker_features fallback.")
						seu@misc$markers[["clusters"]] <- seuratTools:::stash_marker_features("clusters", seu, seurat_assay = "SCT")
						return(seu)
					}
					stop(e)
				}
			)

			seu@meta.data$clusters <- forcats::fct_drop(seu@meta.data$clusters)

			heatmap_features <-
				seu@misc$markers[["clusters"]][["presto"]] %>%
				dplyr::filter(Gene.Name %in% VariableFeatures(seu))

			tidy_eval_arrange <- function(.data, ...) {
				.data %>%
					arrange(...)
			}

			single_cluster_order_vec <-
				seu@meta.data %>%
				dplyr::select(clusters, !!group.by) %>%
				dplyr::arrange(clusters, !!sym(group.by)) %>%
				dplyr::select(clusters, !!group.by) |>
				dplyr::distinct(.data[[group.by]], .keep_all = TRUE) |>
				dplyr::mutate(!!group.by := as.character(.data[[group.by]])) |>
				tibble::deframe() |>
				identity()

			heatmap_features[["Cluster"]] <-
				factor(heatmap_features[["Cluster"]], levels = levels(seu_meta$clusters))

			heatmap_features <-
				heatmap_features %>%
				dplyr::group_by(Gene.Name) |>
				dplyr::slice_max(order_by = Average.Log.Fold.Change, n = 1) |>
				dplyr::ungroup() |>
				dplyr::arrange(Cluster, desc(Average.Log.Fold.Change)) |>
				group_by(Cluster) %>%
				slice_head(n = 5) %>%
				dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
				dplyr::ungroup() %>%
				dplyr::distinct(Gene.Name, .keep_all = TRUE) |>
				identity()
		} else {
			heatmap_features <-
				seu@misc$markers[[group.by]][["presto"]]

			single_cluster_order <- levels(seu@meta.data[[group.by]]) %>%
				set_names(.)

			seu@meta.data[[group.by]] <-
				factor(seu@meta.data[[group.by]], levels = single_cluster_order)

			group_by_clusters <- seu@meta.data[[group.by]]

			seu@meta.data$clusters <- names(single_cluster_order[group_by_clusters])

			seu@meta.data$clusters <- factor(seu@meta.data$clusters, levels = unique(setNames(names(single_cluster_order), single_cluster_order)[levels(seu@meta.data[[group.by]])]))

			heatmap_features <-
				heatmap_features %>%
				dplyr::arrange(Cluster) %>%
				group_by(Cluster) %>%
				slice_head(n = 6) %>%
				dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
				dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
				identity()
		}

		seu$scna[seu$scna == ""] <- ".diploid"
		seu$scna <- factor(seu$scna)

		giotti_genes <- read_giotti_genes()

		heatmap_features <-
			heatmap_features %>%
			dplyr::ungroup() %>%
			left_join(giotti_genes, by = c("Gene.Name" = "symbol")) %>%
			dplyr::mutate(term = replace_na(term, "")) %>%
			dplyr::distinct(Gene.Name, .keep_all = TRUE)

		row_ha <- ComplexHeatmap::rowAnnotation(term = rev(heatmap_features$term))

		if (!is.null(split_columns)) {
			column_split <- sort(seu@meta.data[[split_columns]])
			column_title <- unique(column_split)
		} else {
			column_split <- split_columns
			column_title <- NULL
		}

		heatmap_arrangement = heatmap_arrangement[heatmap_arrangement %in% colnames(seu@meta.data)]
		heatmap_groups = heatmap_groups[heatmap_groups %in% colnames(seu@meta.data)]

		seu_heatmap <- ggplotify::as.ggplot(
			seu_complex_heatmap(seu,
													features = heatmap_features$Gene.Name,
													group.by = heatmap_groups,
													col_arrangement = heatmap_arrangement,
													cluster_rows = FALSE,
													column_split = column_split,
													row_split = rev(heatmap_features$Cluster),
													row_title_rot = 0,
													column_title = column_title,
													column_title_rot = 90
			)
		) +
			labs(title = sample_id) +
			theme()

		labels <- data.frame(clusters = unique(seu[[]][["clusters"]]), label = unique(seu[[]][["clusters"]])) %>%
			identity()

		cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna")) |>
			dplyr::mutate(cluster_facet = dplyr::case_when(
				str_detect(.data[["clusters"]], "hsp") ~ "stress",
				str_detect(.data[["clusters"]], "hyp") ~ "stress",
				.default = "cc"))

		centroid_data <-
			cc_data %>%
			dplyr::group_by(clusters, cluster_facet) %>%
			dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
			dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
			dplyr::mutate(centroid = "centroids") %>%
			identity()

		centroid_plot <-
			cc_data  |>
			ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["clusters"]])) +
			facet_wrap(~cluster_facet) +
			geom_point(size = 0.1) +
			theme_light() +
			theme(
				strip.background = element_blank(),
				strip.text.x = element_blank()
			) +
			geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
			guides(fill = "none", color = "none") +
			NULL

		facet_cell_cycle_plot <-
			cc_data %>%
			ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
			geom_point(size = 0.1) +
			geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
			facet_wrap(~ .data[["clusters"]], ncol = 2) +
			theme_light() +
			geom_label(
				data = labels,
				aes(label = label),
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
			NULL

		appender <- function(string) str_wrap(string, width = 40)

		labels <- data.frame(scna = unique(seu$scna), label = str_replace(unique(seu$scna), "^$", "diploid"))

		clone_ratio <- janitor::tabyl(as.character(seu$scna))$percent[[2]]

		comparison_scna <-
			janitor::tabyl(as.character(seu$scna))[2, 1]

		clone_distribution_plot <- plot_distribution_of_clones_across_clusters(
			seu,
			seu_name = glue("{tumor_id} {comparison_scna}"), var_x = "scna", var_y = "clusters", signif = TRUE, plot_type = "clone"
		)

		full_seu@meta.data[[group.by]] <- factor(full_seu@meta.data[[group.by]], levels = single_cluster_order_vec)
		levels(full_seu@meta.data[[group.by]]) <- names(single_cluster_order_vec)

		if (!is.null(nb_paths)) {
			clone_tree_plot <-
				plot_clone_tree(seu, tumor_id, nb_path, clone_simplifications, sample_id = sample_id, legend = FALSE, horizontal = FALSE)

			collage_plots <- list(
				"A" = seu_heatmap,
				"B" = clone_tree_plot,
				"C" = centroid_plot,
				"D" = facet_cell_cycle_plot,
				"E" = plot_spacer(),
				"F" = plot_spacer()
			)

			layout <- "
              AAAAAAABCEE
              AAAAAAADDEE
              AAAAAAADDFF
              AAAAAAADDFF
              AAAAAAADDFF
              "

			plot_collage <- wrap_plots(collage_plots) +
				plot_layout(design = layout) +
				plot_annotation(tag_levels = "A") +
				NULL
		} else {
			collage_plots <- list(
				"A" = seu_heatmap,
				"B" = facet_cell_cycle_plot,
				"C" = centroid_plot,
				"D" = clone_distribution_plot,
				"E" = plot_spacer()
			)

			layout <- "
              AAAAAACC
              AAAAAABB
              AAAAAABB
              AAAAAABB
              AAAAAADD
      "

			plot_collage <- wrap_plots(collage_plots) +
				plot_layout(design = layout) +
				plot_annotation(tag_levels = "A") +
				NULL
		}

		print(plot_collage)
		# end loop------------------------------
	}

	dev.off()

	return(plot_path)
}

#' Plot DimPlot + marker dotplot for a sample
#'
#' @param seu_path File path to filtered Seurat object RDS
#' @param cluster_order Named list of cluster order data frames
#' @param phase_levels Ordered vector of phase level labels
#' @return File path of output PDF
#' @export
plot_seu_clusters_and_markers <- function(seu_path, cluster_order, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star")) {
	file_id <- fs::path_file(seu_path)

	tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

	message(file_id)
	cluster_order <- cluster_order[[file_id]]

  seu <- readRDS(seu_path)

  if (!is.null(cluster_order)) {
    group.by <- unique(cluster_order$resolution)

    cluster_order <-
      cluster_order %>%
      dplyr::mutate(order = dplyr::row_number()) %>%
      dplyr::filter(!is.na(clusters)) %>%
      dplyr::mutate(clusters = as.character(clusters))

    seu@meta.data$clusters <- seu@meta.data[[group.by]]

    seu_meta <- seu@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::left_join(cluster_order, by = "clusters") %>%
      dplyr::select(-clusters) %>%
      dplyr::rename(clusters = phase) %>%
      identity()

    phase_levels <- phase_levels[phase_levels %in% unique(seu_meta$clusters)]

    seu_meta <-
      seu_meta %>%
      tidyr::unite("new_clusters", all_of(c("clusters", group.by)), remove = FALSE) %>%
      dplyr::arrange(clusters, new_clusters) %>%
      dplyr::mutate(clusters = factor(new_clusters, levels = unique(new_clusters))) %>%
      tibble::column_to_rownames("cell") %>%
      identity()

    seu@meta.data <- seu_meta[rownames(seu@meta.data), ]
  }

  dplot <- Seurat::DimPlot(seu, group.by = "clusters")

  mplot <- plot_markers(seu, metavar = "clusters", marker_method = "presto", return_plotly = FALSE, hide_technical = "all") +
    ggplot2::scale_y_discrete(position = "left") +
    labs(title = sample_id)

  mypatch <- dplot + mplot

  plot_path <- glue("results/{sample_id}_clusters_and_markers.pdf")

  ggsave(plot_path, mypatch, height = 8, width = 12)

  return(plot_path)
}

#' Create a Venn diagram of recurrent diffex genes
#'
#' @param diffex_list Parameter for diffex list
#' @param mytitle Plot title
#' @param selected_phase Parameter for selected phase
#' @return List with ggplot2 Venn plot and region data
#' @export
plot_venn_w_genes <- function(diffex_list, mytitle = NULL, selected_phase = NULL) {

	if(!is.null(selected_phase)){
		diffex_list <-
			diffex_list |>
			map(~dplyr::filter(.x, str_detect(cluster_comparison, selected_phase)))
	}

	mylist <- diffex_list |>
		purrr::compact() |>
		purrr::list_flatten() |>
		dplyr::bind_rows(.id = "sample_id") |>
		dplyr::group_by(symbol) |>
		dplyr::filter(all_same_sign(avg_log2FC)) |>
		dplyr::mutate(symbol = glue("{symbol}_{sign(avg_log2FC)}")) |>
		identity()

	mylist <- split(mylist, mylist$sample_id) |>
		map(dplyr::pull, symbol)

	myregion <-
		mylist |>
		Venn() |>
		process_data() |>
		venn_region()

	intersecting_genes <-
		myregion |>
		tail(n=1) |>
		dplyr::pull(item)

	venn_plot <- ggVennDiagram(mylist) + scale_fill_gradient(low="grey90",high = "red") +
		labs(title = mytitle) +
		NULL

	return(list("plot" = venn_plot, "genes" = myregion))
}

plot_fig_01 <- function(seu_path, plot_path = "results/fig_01.pdf") {
	seu <- readRDS(seu_path)
	seu$scna <- factor(seu$scna, levels = c("", "16q-", "16q- 1q+"))
	seu$scna_status <- factor(
		ifelse(str_detect(seu$scna, "1q"), "w/ 1q+", "w/o 1q+"),
		levels = c("w/o 1q+", "w/ 1q+")
	)
	seu$clusters <- seu$seurat_clusters
	seu$clusters <- factor(seu$clusters)

	labels <- data.frame(
		clusters = unique(seu[[]][["clusters"]]),
		label    = unique(seu[[]][["clusters"]])
	)

	cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))

	centroid_data <-
		cc_data %>%
		dplyr::group_by(clusters) %>%
		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
		dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
		dplyr::mutate(centroid = "centroids")

	dot_colors <- c("clusters", "Phase", "scna")

	centroid_plots <- list()
	ccplots        <- list()
	for (dot_color in dot_colors) {
		ccplots[[dot_color]] <-
			cc_data %>%
			ggplot(aes(
				x = `S.Score`, y = `G2M.Score`,
				group = .data[["clusters"]], color = .data[[dot_color]]
			)) +
			geom_point(size = 0.1) +
			geom_point(
				data = centroid_data,
				aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]),
				size = 6, alpha = 0.7, shape = 23, colour = "black"
			) +
			facet_wrap(~ .data[["clusters"]], nrow = 2) +
			theme_light() +
			geom_label(
				data        = labels,
				aes(label   = label),
				x           = max(cc_data$S.Score) + 0.05,
				y           = max(cc_data$G2M.Score) - 0.1,
				hjust       = 1, vjust = 1,
				inherit.aes = FALSE
			) +
			theme(strip.background = element_blank(), strip.text.x = element_blank())

		centroid_plots[[dot_color]] <-
			cc_data %>%
			ggplot(aes(
				x = `S.Score`, y = `G2M.Score`,
				group = .data[["clusters"]], color = .data[[dot_color]]
			)) +
			geom_point(size = 0.1) +
			theme_light() +
			theme(strip.background = element_blank(), strip.text.x = element_blank()) +
			geom_point(
				data = centroid_data,
				aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]),
				size = 6, alpha = 0.7, shape = 23, colour = "black"
			) +
			guides(fill = "none", color = "none")
	}

	path_d <- "results/fig_01d.pdf"
	path_e <- "results/fig_01e.pdf"
	path_f <- "results/fig_01f.pdf"

	pdf(path_d, h = 3, w = 4)
	print(centroid_plots)
	dev.off()

	pdf(path_e, h = 4, w = 12)
	print(ccplots)
	dev.off()

	clone_distribution_plot <- plot_distribution_of_clones_across_clusters(
		seu,
		seu_name  = glue::glue("asdf"), var_x = "scna", var_y = "clusters",
		signif    = FALSE, plot_type = "clone"
	)
	ggsave(path_f, plot = clone_distribution_plot, w = 3.5, h = 3)

	qpdf::pdf_combine(c(path_d, path_e, path_f), plot_path)
}

plot_fig_02 <- function(seu_path, numbat_rds_files, large_clone_simplifications, plot_path = NULL) {
	if (is.na(seu_path)) return(NA_character_)

	sample_id <- stringr::str_extract(seu_path, "SR[RX][0-9]+")
	plot_path <- plot_path %||% glue::glue("results/fig_02_{sample_id}.pdf")
	fs::dir_create(fs::path_dir(plot_path))

	# Panel A: phylo heatmap (PDF)
	heatmap_paths <- make_numbat_heatmaps(seu_path, numbat_rds_files, p_min = 0.9, show_segment_names_on_x = TRUE)
	panel_a_path <- heatmap_paths[[1]]

	# Panel B: clone tree (PDF)
	panel_b_path <- save_clone_tree_from_path(
		seu_path, numbat_rds_files, large_clone_simplifications,
		label = "_fig_02_clone_tree", legend = FALSE, horizontal = FALSE
	)

	# Helper: load a PDF page as a ggdraw grob
	pdf_to_grob <- function(pdf_path, trim = FALSE) {
		if (is.null(pdf_path) || is.na(pdf_path)) return(NULL)
		img <- pdftools::pdf_render_page(pdf_path, page = 1, dpi = 150)
		if (trim) img <- magick::image_trim(magick::image_read(img))
		cowplot::ggdraw() + cowplot::draw_image(img)
	}

	grob_a <- pdf_to_grob(panel_a_path)
	grob_b <- pdf_to_grob(panel_b_path, trim = TRUE)

	# Load Seurat for panels C-F
	seu <- readRDS(seu_path)
	seu$clusters <- seu$seurat_clusters
	seu$clusters <- factor(seu$clusters)

	cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))

	centroid_data <- cc_data %>%
		dplyr::group_by(clusters) %>%
		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score), .groups = "drop") %>%
		dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters)))

	cluster_g2m_order <- cc_data %>%
		dplyr::group_by(clusters) %>%
		dplyr::summarise(g2m_prop = mean(Phase == "G2M"), .groups = "drop") %>%
		dplyr::arrange(g2m_prop) %>%
		dplyr::pull(clusters)

	cc_data$clusters <- factor(cc_data$clusters, levels = cluster_g2m_order)
	seu$clusters <- factor(seu$clusters, levels = cluster_g2m_order)
	centroid_data$clusters <- factor(centroid_data$clusters, levels = cluster_g2m_order)

	# Panel C: 3 UMAP plots
	panel_c <- (DimPlot(seu, group.by = "clusters", label = TRUE) + NoLegend()) |
	           DimPlot(seu, group.by = "Phase") |
	           DimPlot(seu, group.by = "scna")

	# Panel D: CC space scatter x3 (clusters, Phase, scna) — same dimensions as panel C (2996x1498)
	make_cc_scatter <- function(color_by, show_fill_legend = TRUE) {
		p <- cc_data %>%
			ggplot(aes(x = S.Score, y = G2M.Score, color = .data[[color_by]])) +
			geom_point(size = 0.1) +
			geom_point(
				data = centroid_data,
				aes(x = mean_x, y = mean_y, fill = clusters),
				size = 6, shape = 23, colour = "black", alpha = 0.7, inherit.aes = FALSE
			) +
			theme_light() +
			labs(title = color_by)
		if (!show_fill_legend) p <- p + guides(fill = "none")
		p
	}
	panel_d <- make_cc_scatter("clusters") |
	           make_cc_scatter("Phase", show_fill_legend = FALSE) |
	           make_cc_scatter("scna", show_fill_legend = FALSE)

	# Panel E: faceted by cluster, colored by scna
	panel_e <- cc_data %>%
		ggplot(aes(x = S.Score, y = G2M.Score, color = scna)) +
		geom_point(size = 0.1) +
		geom_point(
			data = centroid_data,
			aes(x = mean_x, y = mean_y, fill = clusters),
			size = 6, shape = 23, colour = "black", alpha = 0.7, inherit.aes = FALSE
		) +
		facet_wrap(~clusters) +
		geom_label(
			data = data.frame(clusters = unique(cc_data$clusters)),
			aes(label = clusters),
			x = max(cc_data$S.Score),
			y = max(cc_data$G2M.Score),
			hjust = 1, vjust = 1,
			inherit.aes = FALSE
		) +
		theme_light() +
		theme(strip.background = element_blank(), strip.text.x = element_blank()) +
		guides(color = "none", fill = "none")

	# Panel F: clone distribution
	panel_f <- plot_distribution_of_clones_across_clusters(
		seu, seu_name = sample_id, var_x = "scna", var_y = "clusters", plot_type = "clone"
	) +
		theme(legend.position = "none")

	# Assemble rows
	row1 <- cowplot::plot_grid(
		grob_a, grob_b, nrow = 1, rel_widths = c(0.8, 0.2),
		labels = c("A", "B"), label_size = 14
	)
	row4 <- cowplot::plot_grid(
		panel_e, panel_f, nrow = 1, rel_widths = c(0.8, 0.2),
		labels = c("E", "F"), label_size = 14
	)

	page <- cowplot::plot_grid(
		row1, panel_c, panel_d, row4,
		ncol = 1,
		rel_heights = c(1000, 1000, 1000, 1070),
		labels = c("", "C", "D", ""), label_size = 14
	)

	# Page width = panel C width; total height = sum of row heights
	ggsave(plot_path, page, width = 2996, height = 1000 + 1000 + 1000 + 1070, units = "px")
	return(plot_path)
}
