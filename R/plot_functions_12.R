# Plot Functions (104)

#' Perform clustering analysis
#'
#' @param seu Seurat object
#' @return Modified Seurat object
#' @export
# Performance optimizations applied:
# - rownames_round_trip: Avoid unnecessary rownames conversions - keep as column when possible

annotate_cluster_scna_percentage <- function(seu){
	scna_per_cluster <- 
		seu@meta.data[c("clusters", "scna")] |> 
		group_by(clusters, scna) |>
		dplyr::summarize(clones_per_cluster = dplyr::n()) |> 
		mutate(percent_scna = proportions(clones_per_cluster) * 100) %>%
		dplyr::filter(scna == "w_scna") |> 
		dplyr::select(clusters, percent_scna) |> 
		identity()
	
	seu_meta <- seu@meta.data |> 
		tibble::rownames_to_column("cell") |> 
		dplyr::left_join(scna_per_cluster, by = "clusters") |> 
		tibble::column_to_rownames("cell")
	
	seu@meta.data <- seu_meta
	
	return(seu)
}

#' Create a plot visualization
#'
#' @param seu Seurat object
#' @param direction Parameter for direction
#' @param mygenes Gene names or identifiers
#' @param max_genes Gene names or identifiers
#' @param scna_of_interest Character string (default: "02_p")
#' @param in_scna Logical flag (default: TRUE)
#' @param only_scna Logical flag (default: FALSE)
#' @return ggplot2 plot object
#' @export
select_genes_to_plot <- function(seu, direction, mygenes, max_genes = 15, scna_of_interest = "02_p", in_scna = TRUE, only_scna = FALSE) {
	median_expression_values <-
		FetchData(seu, c("clusters", mygenes[[direction]], "percent_scna")) |> 
		tidyr::pivot_longer(-any_of(c("clusters", "percent_scna")), names_to = "symbol", values_to = "expression") |>
		dplyr::group_by(clusters, percent_scna, symbol) |>
		dplyr::summarize(med_expression = median(expression)) |>
		dplyr::filter(str_detect(clusters, "^g1_")) |>
		identity()
	
	if(direction == "-1"){
		first_threshold_genes <- median_expression_values |> 
			dplyr::group_by(symbol) |> 
			dplyr::filter(med_expression == max(med_expression) & percent_scna < 50) |> 
			dplyr::pull(symbol)
		
		second_threshold_genes <- median_expression_values |> 
			dplyr::group_by(symbol) |> 
			dplyr::filter(med_expression == min(med_expression) & percent_scna > 50) |> 
			dplyr::pull(symbol)	
	} else if(direction == "1"){
		first_threshold_genes <- median_expression_values |> 
			dplyr::group_by(symbol) |> 
			dplyr::filter(med_expression == min(med_expression) & percent_scna < 50) |> 
			dplyr::pull(symbol)
		
		second_threshold_genes <- median_expression_values |> 
			dplyr::group_by(symbol) |> 
			dplyr::filter(med_expression == max(med_expression) & percent_scna > 50) |> 
			dplyr::pull(symbol)	
	}
	
	variance_genes <- 
		group_by(median_expression_values, symbol) |> 
		dplyr::summarize(var = var(med_expression)) |> 
		dplyr::filter(var > 0) |> 
		dplyr::arrange(desc(var)) |> 
		dplyr::pull(symbol) |>
		identity()
	
	plotted_genes <- purrr::reduce(list(mygenes[[direction]], first_threshold_genes, second_threshold_genes, variance_genes), base::intersect)
	
	plotted_genes <-
		variance_genes[variance_genes %in% plotted_genes]
	
	plotted_genes <- head(plotted_genes, n = max_genes)
	
	if(only_scna){
		gene_locations <- find_genes_by_arm(mygenes[[direction]]) |> 
			tidyr::unite(seqnames_arm, any_of(c("seqnames", "arm")))
	} else {
		gene_locations <- find_genes_by_arm(plotted_genes) |> 
			tidyr::unite(seqnames_arm, any_of(c("seqnames", "arm")))
	}
	
	if(in_scna){
		plotted_genes <- gene_locations |> 
			dplyr::filter(seqnames_arm == scna_of_interest) |> 
			dplyr::pull(symbol)	
	} else {
		plotted_genes <- gene_locations |> 
			dplyr::filter(seqnames_arm != scna_of_interest) |> 
			dplyr::pull(symbol)	
	}
	
}

#' Perform make table s01 operation
#'
#' @param study_cell_stats Cell identifiers or information
#' @param path File path
#' @return Function result
#' @export
make_table_s01 <- function(study_cell_stats, path = "doc/table_s01.csv"){
	study_cell_stats |> 
		dplyr::group_by(study, sample_id) |> 
		dplyr::summarize(
			mean_umi = mean(nCount_gene),
			mean_genes_detected = mean(nFeature_gene),
			mean_percent_mt = mean(percent.mt),
		) |> 
		write_csv(path) |> 
		identity() 
	
	return(path)
}

#' Create a plot visualization
#'
#' @param cluster_orders Cluster information
#' @param plot_path File path
#' @return ggplot2 plot object
#' @export
plot_integrated_1q_fig <- function(cluster_orders, plot_path = "results/fig_02.pdf"){
	integrated_1q_seus <- dir_ls("output/seurat/integrated_1q/", regexp = ".*SR[RX][0-9].*rds") |> as.character()
	
	batch_corrected_seu <- "output/seurat/integrated_1q/integrated_seu_1q_complete.rds"
	
	names(integrated_1q_seus) <- fs::path_ext_remove(fs::path_file(integrated_1q_seus))
	
	fig_2a_c_trio <- make_clone_distribution_figure("output/seurat/integrated_1q/integrated_seu_1q_trio.rds", cluster_orders,
																						 height = 12, width = 20, plot_path = "results/fig_2a_c_trio.pdf", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"), heatmap_arrangement = c("clusters", "scna", "batch"))
	
	fig_2a_c <- make_clone_distribution_figure(batch_corrected_seu, cluster_orders,
																						 height = 12, width = 20, plot_path = "results/fig_2a_c_fig04.pdf", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"), heatmap_arrangement = c("clusters", "scna", "batch"))
	
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
		
	fig_2d_cc_path <- ggsave("results/fig_2d_cc.pdf", fig_2d_cc, width = 12, height = 4)

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
		
	
	fig_2d_hypoxia_hsp_path <- ggsave("results/fig_2d_hypoxia_hsp.pdf", fig_2d_hypoxia_hsp, width= 14, height = 6)
	
	qpdf::pdf_combine(
		list(
			fig_2a_c_trio, 
			fig_2a_c, 
			fig_2d_cc_path, 
			fig_2d_hypoxia_hsp_path
			),
		plot_path)
	
	return(list("plot" = plot_path, "table" = fig_2d_cc_table))
}

#' Perform drop y axis operation
#'
#' @param myplot Plot object (ggplot2)
#' @return ggplot2 plot object
#' @export
drop_y_axis <- function(myplot){
		myplot[[1]] <- myplot[[1]] + 
			theme(
				axis.title.y = element_blank(),
				axis.text.y = element_blank()
			)
		
		myplot[[2]] <- myplot[[2]] + 
			theme(
				axis.title.y = element_blank(),
				axis.text.y = element_blank()
			)
		
		return(myplot)
	}

	#' Create a plot visualization
#'
#' @param cluster_orders Cluster information
#' @param plot_path File path
#' @return ggplot2 plot object
#' @export
plot_integrated_1q_fig_low_hypoxia <- function(seu_path = "output/seurat/integrated_1q/integrated_seu_1q_complete.rds", cluster_orders, plot_path = tempfile(tmpdir = "results", fileext = ".pdf")){
	
	batch_corrected_seu <- readRDS(seu_path)

	batch_hashes <- batch_corrected_seu@meta.data[c("batch", "batch_hash")] |> 
		distinct() |> 
		tibble::deframe() |> 
		identity()

	batch_cluster_orders = map(batch_hashes, ~ {
		read_cluster_orders_table(hash = .x)
	})
	
	integrated_1q_seus <- map_chr(batch_hashes, ~ {
		read_seu_path(hash = .x) |> as.character()
	}) |> 
		identity()

	names(integrated_1q_seus) <- names(batch_hashes)

	fig_2a_c <- make_clone_distribution_figure(seu_path, cluster_orders, height = 12, width = 20, plot_path = tempfile(tmpdir = "results", fileext = ".pdf"), heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"), heatmap_arrangement = c("clusters", "scna", "batch"))
	
	possibly_plot_clone_cc_plots <- possibly(plot_clone_cc_plots)
	fig_2d_cc <- map(c(seu_path, integrated_1q_seus), plot_clone_cc_plots, var_y = "integrated_snn_res.0.6", scna_of_interest = "1q", labeled_values = c("g2_m"))
	# fig_2d_cc <- map(c(seu_path, integrated_1q_seus), plot_clone_cc_plots, var_y = "clusters", scna_of_interest = "1q", labeled_values = c("g2_m"))
	
	fig_2d_cc_table <- map(fig_2d_cc, ~.x$data) |> 
	set_names(c("integrated", names(integrated_1q_seus)))  |> 
		write_xlsx(path = tempfile(tmpdir = "results", fileext = ".xlsx"))
	
	
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
		
	fig_2d_cc_path <- ggsave(tempfile(tmpdir = "results", fileext = ".pdf"), fig_2d_cc, width = 12, height = 4)

	possibly_plot_clone_cc_plots <- possibly(plot_clone_cc_plots)
	
	# fig_2d_hypoxia_hsp <- map(integrated_1q_seus, ~possibly_plot_clone_cc_plots(.x, scna_of_interest = "1q", kept_phases = c("hsp", "hypoxia"), labeled_values = c("hsp", "hypoxia")))
	
	# fig_2d_hypoxia_hsp <- 
	# 	fig_2d_hypoxia_hsp |> 
	# 	wrap_plots() +
	# 	plot_layout(
	# 		nrow = 1,
	# 		guides = "collect",
	# 		axis_titles = "collect") +
	# 	NULL
		
	# fig_2d_hypoxia_hsp_path <- ggsave(tempfile(tmpdir = "results", fileext = ".pdf"), fig_2d_hypoxia_hsp, width= 14, height = 6)
	
	qpdf::pdf_combine(
		list(
			fig_2a_c, 
			fig_2d_cc_path 
			# fig_2d_hypoxia_hsp_path
			),
		plot_path)
	
	return(list("plot" = plot_path, "table" = fig_2d_cc_table))
}