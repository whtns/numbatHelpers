# Plot Functions (106)

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
	
	# panels_e <- 
	# 	integrated_diffex[[1]]$all |> 
	# 	dplyr::filter(p_val_adj <= 0.05) |> 
	# 	# plot_enrichment(p_val_cutoff = 0.1) |> 
	# 	identity()
	# panel_e_file <- ggsave(tempfile(fileext = ".pdf"), width = 8, height = 6)
	
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

#' Perform clustering analysis
#'
#' @param seu Seurat object
#' @param group.by Character string (default: "integrated_snn_res.0.4")
#' @return Modified Seurat object
#' @export
drop_mt_cluster <- function(seu, group.by = "integrated_snn_res.0.4"){
	correct_mt <- seu@misc$markers[[group.by]]$presto |> 
		group_by(Cluster) |> 
		dplyr::slice_head(n = 5) |> 
		dplyr::summarize(correct_mt = ifelse(any(str_detect(Gene.Name, "^MT-")), FALSE, TRUE)) |> 
		tibble::deframe() |>
		which() |>
		identity()
	
	if(is.factor(seu[[]][[group.by]])){
		seu[[]][[group.by]] <- droplevels(seu[[]][[group.by]])
	}
	
	seu <- seu[,seu@meta.data[[group.by]] %in% names(correct_mt)] |> 
		find_all_markers(metavar = group.by)
	
	return(seu)

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
	# Drop skipped samples (NA_character_) and unwritten files before read_csv() takes
	# an NA as a filename -- see .drop_missing_paths().
	paths <- .drop_missing_paths(figure_input)

	if (length(paths) == 0) {
		message("plot_fig_07_08: no usable input CSVs, returning NULL")
		return(invisible(NULL))
	}

	raw_tables <-
		paths |>
		set_names(str_extract(paths, "SR[RX][0-9]+")) |>
		map(read_csv)

	non_empty <- purrr::keep(raw_tables, ~nrow(.x) > 0)

	if (length(non_empty) == 0) {
		message("plot_fig_07_08: all input CSVs are empty, returning NULL")
		return(invisible(NULL))
	}

	unfiltered_input <-
		non_empty |>
		map(~split(.x, .x$cluster)) |>
		purrr::list_transpose() |>
		map(~purrr::compact(.x)) |>
		map(dplyr::bind_rows, .id = "sample_id") |>
		dplyr::bind_rows(.id = "clusters") |>
		dplyr::filter(location == "cis")
	
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
#' @param diffex_list Parameter for diffex list
#' @param x_var Character string (default: "sample_id")
#' @param recurrence_threshold Threshold value for filtering
#' @param n_genes Gene names or identifiers
#' @param p_val_adj_threshold Threshold value for filtering
#' @return ggplot2 plot object
#' @export
dotplot_recurrent_genes <- function(diffex_list, x_var = "sample_id", recurrence_threshold = 3, n_genes = 50, p_val_adj_threshold = 0.05){
	#
	required_cols <- c("symbol", "avg_log2FC", "p_val_adj")
	valid_diffex <- diffex_list |>
		map(~ dplyr::bind_rows(.x, .id = "comparison")) |>
		purrr::keep(~ is.data.frame(.x) && nrow(.x) > 0 && all(required_cols %in% colnames(.x)))

	if (length(valid_diffex) == 0) {
		return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No valid recurrent diffex input"))
	}

	gene_list <-
		valid_diffex |>
		map(dplyr::filter, p_val_adj <= p_val_adj_threshold) |>
		dplyr::bind_rows(.id = "sample_id") |>
		dplyr::group_by(symbol) %>%
		dplyr::filter(all_same_sign(avg_log2FC)) |>
		dplyr::mutate(abs_mean_FC = abs(mean(avg_log2FC))) %>%
		dplyr::mutate(recurrence = dplyr::n_distinct(sample_id)) %>%
		dplyr::mutate(sign = sign(avg_log2FC)) |> 
		dplyr::arrange(sign, desc(abs_mean_FC)) |>
		dplyr::filter(recurrence >= recurrence_threshold) |> 
		dplyr::distinct(sample_id, .keep_all = TRUE) |> 
		dplyr::mutate(sample_comparison = glue("{sample_id}_{comparison}")) |>
		dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
		identity()

	if (nrow(gene_list) == 0) {
		return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No recurrent genes after filtering"))
	}
	
	gene_list$symbol <- factor(gene_list$symbol, levels = unique(gene_list$symbol))
	
	top_genes <- gene_list |> 
		dplyr::group_by(sign) |> 
		dplyr::slice_head(n = n_genes) |> 
		dplyr::pull(symbol)
	
	# gene_list <- 
		gene_list |> 
		dplyr::filter(symbol %in% top_genes) |> 
		dplyr::mutate(mean_FC = abs_mean_FC*sign) |> 
		dplyr::arrange(mean_FC) |> 
			dplyr::mutate(tumor_id = str_extract(sample_id, "SR[RX][0-9]+")) |> 
		ggplot(aes(y = symbol, x = tumor_id, size = neg_log_p_val_adj, color = avg_log2FC)) + 
		geom_point() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
		scale_color_gradient2()
}

#' Create a plot visualization
#'
#' @param diffex_list Parameter for diffex list
#' @param x_var Character string (default: "sample_id")
#' @param recurrence_threshold Threshold value for filtering
#' @param n_genes Gene names or identifiers
#' @param p_val_adj_threshold Threshold value for filtering
#' @return ggplot2 plot object
#' @export
dotplot_diffex <- function(diffex_list, x_var = "sample_id", recurrence_threshold = 3, n_genes = 50, p_val_adj_threshold = 0.05){
	#
	required_cols <- c("symbol", "avg_log2FC", "p_val_adj")
	valid_diffex <- diffex_list |>
		map(~ dplyr::bind_rows(.x, .id = "comparison")) |>
		purrr::keep(~ is.data.frame(.x) && nrow(.x) > 0 && all(required_cols %in% colnames(.x)))

	if (length(valid_diffex) == 0) {
		return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No valid diffex input"))
	}

	gene_list <-
		valid_diffex |>
		map(dplyr::filter, p_val_adj <= p_val_adj_threshold) |>
		dplyr::bind_rows(.id = "sample_id") |>
		dplyr::group_by(symbol) %>%
		dplyr::filter(all_same_sign(avg_log2FC)) |>
		dplyr::mutate(abs_mean_FC = abs(mean(avg_log2FC))) %>%
		dplyr::mutate(recurrence = dplyr::n_distinct(sample_id)) %>%
		dplyr::mutate(sign = sign(avg_log2FC)) |> 
		dplyr::arrange(sign, desc(abs_mean_FC)) |>
		dplyr::filter(recurrence >= recurrence_threshold) |> 
		dplyr::distinct(sample_id, .keep_all = TRUE) |> 
		dplyr::mutate(sample_comparison = glue("{sample_id}_{comparison}")) |>
		dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
		identity()

	if (nrow(gene_list) == 0) {
		return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No genes after filtering"))
	}
	
	gene_list$symbol <- factor(gene_list$symbol, levels = unique(gene_list$symbol))
	
	top_genes <- gene_list |> 
		dplyr::group_by(sign) |> 
		dplyr::slice_head(n = n_genes) |> 
		dplyr::pull(symbol)
	
	# gene_list <- 
	gene_list |> 
		dplyr::filter(symbol %in% top_genes) |> 
		dplyr::mutate(mean_FC = abs_mean_FC*sign) |> 
		dplyr::arrange(mean_FC) |> 
		dplyr::mutate(tumor_id = str_extract(sample_id, "SR[RX][0-9]+")) |> 
		ggplot(aes(y = symbol, x = tumor_id, size = neg_log_p_val_adj, color = avg_log2FC)) + 
		geom_point() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
		scale_color_gradient2()
}

