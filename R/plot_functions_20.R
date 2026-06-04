# Plot Functions (112)

#' Perform split vector into overlapping chunks operation
#'
#' @param v Parameter for v
#' @param n Parameter for n
#' @return Function result
#' @export
split_vector_into_overlapping_chunks <- function(v, n) {
	indices <- seq(1, length(v) - n + 1)
	chunks <- lapply(indices, function(i) v[i:(i + n - 1)])
	return(chunks)
}

#' Perform make leng sce operation
#'
#' @param tibble Parameter for tibble
#' @param column_to_rownames Color specification
#' @param dplyr Parameter for dplyr
#' @param mutate Parameter for mutate
#' @param enframe Parameter for enframe
#' @param sample_id Parameter for sample id
#' @param SingleCellExperiment Cell identifiers or information
#' @return Function result
#' @export
make_leng_sce <- function(tibble, column_to_rownames, dplyr, mutate, enframe, sample_id, SingleCellExperiment) {
	test1 <- read_csv("data/GSE64016_H1andFUCCI_normalized_EC.csv.gz") |> 
		tibble::column_to_rownames("...1") |> 
		identity()
	
	sce_coldata <- 
		colnames(test1) |>
		tibble::enframe("rownum", "sample_id") |> 
		dplyr::mutate(name = sample_id) |> 
		dplyr::mutate(phase = str_remove(sample_id, "_.*")) |> 
		tibble::column_to_rownames("name") |> 
		DataFrame() |> 
		identity()
	
	leng_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=test1))
	colData(leng_sce) <- sce_coldata
	
	leng_sce <- leng_sce[,!colData(leng_sce)$phase == "H1"]
	
	leng_sce <- logNormCounts(leng_sce)
}

#' Add annotations to data
#'
#' @param seu Seurat object
#' @return Modified Seurat object
#' @export
annotate_cell_cycle_singleR <- function(seu) {
	
	# # Find genes that are cell cycle-related.
	# library(org.Hs.eg.db)
	# cycle.anno <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", 
	# 																		columns="SYMBOL")[,"SYMBOL"] |> 
	# 	unique()
	
	cycle.anno <- readLines("data/go_0007049_cc_genes.txt")
	
	sce <- Seurat::as.SingleCellExperiment(seu)
	
	test.data <- logcounts(sce)
	
	rownames(test.data) <- rownames(rowData(sce))
	
	library(SingleR)
	
	leng_sce <- readRDS("data/leng_sce.rds")
	
	assignments <- SingleR(test.data, ref=leng_sce, label=leng_sce$phase, 
												 de.method="wilcox", restrict=cycle.anno)
	
	seu$singleR_phase <- assignments$labels
	
	return(seu)
}

#' Perform check mycn mean expression operation
#'
#' @param seu_path File path
#' @return Function result
#' @export
check_MYCN_mean_expression <- function(seu_path){
	seu <- readRDS(seu_path)
	FetchData(seu, "MYCN", assay = "counts") |> 
		dplyr::mutate(counts = cur_data()[[1]]) |> 
		dplyr::summarize(
			mean_expr = mean(counts), 
			sd(counts), 
			count_greater_than_zero = sum(counts > 0),
			total_count = n(),
			percentage_greater_than_zero = (count_greater_than_zero / total_count) * 100
		) |>
		identity()
	
}

#' Create a plot visualization
#'
#' @param seu_path File path
#' @param subtype_markers Parameter for subtype markers
#' @return ggplot2 plot object
#' @export
plot_fig_s25 <- function(seu_path = "output/seurat/integrated_1q/integrated_seu_1q_complete.rds", subtype_markers) {
    seu <- readRDS(seu_path)
    nbin <- floor(length(VariableFeatures(seu)) / 100)
    
    seu <- Seurat::AddModuleScore(seu, features = subtype_markers, name = "subtype", nbin = nbin, ctrl = 100)
    
    plot_input <- seu@meta.data[,c("clusters", "subtype1", "subtype2")] |> 
        tidyr::pivot_longer(-clusters, names_to = "subtype", values_to = "score") |> 
        dplyr::filter(clusters %in% c("g1_3", "g1_0", "g1_6"))
    
    ggpubr::ggviolin(plot_input, x = "clusters", y = "score", fill = "subtype", add = "boxplot")
    
    
    fig_s25 <- "results/fig_s25.pdf"
    
    ggsave(fig_s25, w = 8, h = 4) |> 
        browseURL()
}

make_integrated_collage <- function(seu_path = NULL, cluster_order = NULL, nb_paths = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", height = 10, width = 18, equalize_scna_clones = FALSE, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL, large_clone_comparisons = NULL, scna_of_interest = "1q", min_cells_per_cluster = 50, return_plots = FALSE, split_columns = "clusters", wo_nb = FALSE) {
	kept_phases <- kept_phases %||% phase_levels
	
	file_id <- fs::path_file(seu_path)
	
	tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
	
	message(file_id)
	cluster_order <- cluster_order[[file_id]]
	
	full_seu <- readRDS(seu_path)
	
	# subset by retained clones ------------------------------
	if(any(large_clone_comparisons)){
		clone_comparisons <- names(large_clone_comparisons[[sample_id]])
		clone_comparison <- clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
		retained_clones <- clone_comparison %>%
			str_extract("[0-9]_v_[0-9]") %>%
			str_split("_v_", simplify = TRUE)
		
	} else {
		retained_clones = sort(unique(full_seu$clone_opt)) |> 
			set_names()
	}
	
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
			
			# mysec ------------------------------
			
			heatmap_features <-
				seu@misc$markers[["clusters"]][["presto"]] %>%
				dplyr::filter(Gene.Name %in% VariableFeatures(seu))
			
			tidy_eval_arrange <- function(.data, ...) {
				.data %>%
					arrange(...)
			}
			#
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
		
		# large_enough_clusters <-
		# 	seu@meta.data %>%
		# 	dplyr::group_by(clusters) %>%
		# 	dplyr::count() |>
		# 	dplyr::filter(n >= min_cells_per_cluster) %>%
		# 	dplyr::pull(clusters)
		# 
		# seu <- seu[, seu$clusters %in% large_enough_clusters]
		
		seu$scna[seu$scna == ""] <- ".diploid"
		seu$scna <- factor(seu$scna)
		# levels(seu$scna)[1] <- "none"
		
		giotti_genes <- read_giotti_genes()

		heatmap_features <-
			heatmap_features %>%
			dplyr::ungroup() %>%
			left_join(giotti_genes, by = c("Gene.Name" = "symbol")) %>%
			# select(Gene.Name, term) %>%
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
			# dplyr::rename({{group.by}} := cluster) %>%
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
				# x = Inf,
				# y = -Inf,
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
			# guides(color = "none") +
			NULL
		
		appender <- function(string) str_wrap(string, width = 40)
		
		labels <- data.frame(scna = unique(seu$scna), label = str_replace(unique(seu$scna), "^$", "diploid"))
		
		
		# umap_plots <- DimPlot(full_seu, group.by = c("scna", "clusters"), combine = FALSE) %>%
		# 	# map(~(.x + theme(legend.position = "bottom"))) %>%
		# 	wrap_plots(ncol = 1)
		# full_seu$clusters
		# full_seu[[group.by]] <-
		full_seu@meta.data[[group.by]] <- factor(full_seu@meta.data[[group.by]], levels = single_cluster_order_vec)
		levels(full_seu@meta.data[[group.by]]) <- names(single_cluster_order_vec)
		umap_plots <- make_faded_umap_plots(full_seu, retained_clones, group_by = group.by)
		
		clone_ratio <- janitor::tabyl(as.character(seu$scna))$percent[[2]]
		
		comparison_scna <-
			janitor::tabyl(as.character(seu$scna))[2, 1]
		
		clone_distribution_plot <- plot_distribution_of_clones_across_clusters(
			seu,
			seu_name = glue("{tumor_id} {comparison_scna}"), var_x = "scna", var_y = "clusters", signif = TRUE, plot_type = "clone"
		)
		
		if (!is.null(nb_paths)) {
			
			clone_tree_plot <-
				plot_clone_tree(seu, tumor_id, nb_path, clone_simplifications, sample_id = sample_id, legend = FALSE, horizontal = FALSE)
			
			collage_plots <- list(
				"seu_heatmap" = seu_heatmap,
				"facet_cell_cycle_plot" = facet_cell_cycle_plot,
				"umap_plots" = umap_plots,
				"clone_distribution_plot" = clone_distribution_plot,
				"clone_tree_plot" = clone_tree_plot,
				"centroid_plot" = centroid_plot
			)
			
			layout <- "
              AAAAAAEFCC
              AAAAAABBCC
              AAAAAABBDD
              AAAAAABBDD
              AAAAAABBDD
              "
			
			plot_collage <- wrap_plots(collage_plots) +
				# plot_layout(widths = c(16, 4)) +
				plot_layout(design = layout) +
				plot_annotation(tag_levels = "A") +
				NULL
		} else {
			collage_plots <- list(
				"A" = seu_heatmap,
				"B" = facet_cell_cycle_plot,
				"C" = umap_plots,
				"D" = clone_distribution_plot,
				"E" = plot_spacer(),
				"F" = centroid_plot
			)
			
			layout <- "
              AAAAAAEFCC
              AAAAAABBCC
              AAAAAABBDD
              AAAAAABBDD
              AAAAAABBDD
              "
			
			plot_collage <- wrap_plots(collage_plots) +
				# plot_layout(widths = c(16, 4)) +
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

make_clustrees_for_sample <- function(seu_path, mylabel = "sample_id", assay = "SCT", resolutions = seq(0.2, 1.0, by = 0.2), fisher_p_val_threshold = 0.1) {
  sample_id <- str_extract(seu_path, "SR[RX][0-9]+")
  sample_id <- mylabel

  seu <- readRDS(seu_path)
  
  seu@meta.data <- 
  	seu@meta.data %>%
  	dplyr::mutate(across(where(is.character), ~na_if(.x, ""))) %>%
  	dplyr::mutate(scna = replace_na(scna, ".diploid")) %>%
  	identity()

  scna_counts <-
    seu@meta.data %>%
    dplyr::mutate(scna = factor(scna))

  scna_clones <- levels(scna_counts$scna)

  pairwise_clone_vectors <-
    bind_cols(scna_clones[-length(scna_clones)], scna_clones[-1]) %>%
    t() %>%
    as.data.frame() %>%
    as.list() %>%
    map(as.character) %>%
    identity()

  names(pairwise_clone_vectors) <- map(pairwise_clone_vectors, ~ paste(., collapse = "_v_"))

  possible_make_clustree_for_clone_comparison <- possibly(make_clustree_for_clone_comparison)
  
  
  
  clustree_output <- map(pairwise_clone_vectors, ~ possible_make_clustree_for_clone_comparison(seu, sample_id, .x, mylabel, assay, resolutions, fisher_p_val_threshold))

  return(clustree_output)
}