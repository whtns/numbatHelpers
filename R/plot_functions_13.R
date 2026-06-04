# Plot Functions (105)

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
#' @param diffex_list Parameter for diffex list
#' @param mytitle Plot title
#' @param selected_phase Parameter for selected phase
#' @return ggplot2 plot object
#' @export
plot_venn_w_genes <- function(diffex_list, mytitle = NULL, selected_phase = NULL) {
	# 
	
	if(!is.null(selected_phase)){
		diffex_list <-
			diffex_list |> 
			map(~dplyr::filter(.x, str_detect(cluster_comparison, selected_phase))) 
	}
	
	mylist <- diffex_list |> 
		purrr::compact() |> 
		# map(~reduce(.x, union)) |> 
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
		# labs(title = "2p g1", subtitle = paste(intersecting_genes[[1]], collapse = "; ")) +
		labs(title = mytitle) +
		NULL
	
	return(list("plot" = venn_plot, "genes" = myregion))
}

