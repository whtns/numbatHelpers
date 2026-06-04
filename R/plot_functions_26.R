# Plot Functions (125)

#' Perform window operation
#'
#' @param x Parameter for x
#' @param size Size parameter
#' @return Function result
#' @export
window <- function(x, size) {
		lapply(seq_len(length(x) - size + 1), function(i) x[i:(i+size-1)])
	}

#' Perform differential expression analysis
#'
#' @param seu_path File path
#' @param numbat_rds_files File path
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @param cluster_dictionary Cluster information
#' @param location Character string (default: "cis")
#' @param scna_of_interest Parameter for scna of interest
#' @return Differential expression results
#' @export
find_diffex_bw_clones_for_each_cluster <- function(seu_path, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "cis", scna_of_interest = NULL) {
  #
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  numbat_rds_files <- numbat_rds_files %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  mynb <- readRDS(numbat_rds_files[[tumor_id]])

  seu <- readRDS(seu_path)

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  myclusters <- sort(unique(seu@meta.data[["clusters"]])) %>%
    set_names(.)

  clone_diff_per_cluster <- function(cluster_for_diffex, seu, clone_comparisons, location, mynb) {
    seu0 <- seu[, seu[["clusters"]] == cluster_for_diffex]
    Idents(seu0) <- seu0$clone_opt
    diffex <- imap(clone_comparisons, make_clone_comparison, seu0, mynb, location = location) |>
      purrr::compact()
    return(diffex)
  }
  possible_clone_diff_per_cluster <- possibly(clone_diff_per_cluster)


  clone_comparisons <- large_clone_comparisons[[tumor_id]]
  
  if(!is.null(scna_of_interest)){
  	clone_comparisons <-
  		clone_comparisons[str_detect(names(clone_comparisons), scna_of_interest)]	
  	
  	# log cell number per clone per cluster 
  	idents <-
  		names(clone_comparisons) %>%
  		str_extract("[0-9]_v_[0-9]") %>%
  		str_split(pattern = "_v_") %>%
  		unlist()
  	
  	log_path <- glue("results/{path_ext_remove(path_file(seu_path))}_log.csv")
  	
  	clone_per_cluster <- seu@meta.data |> 
  		dplyr::filter(clone_opt %in% idents) |> 
  		janitor::tabyl(scna, clusters) |> 
  		write_csv(log_path)
  }
  # 
  message(sample_id)
  
  diffex <- map(myclusters, possible_clone_diff_per_cluster, seu, clone_comparisons = clone_comparisons, location = location, mynb = mynb) %>%
    compact() %>%
    map(bind_rows, .id = "clone_comparison") %>%
    bind_rows(.id = "cluster") %>%
    # dplyr::arrange(cluster, p_val_adj) %>%
    identity()

  scna_of_interest <- scna_of_interest %||% ""
  
  diffex_path <- glue("results/{sample_id}_cluster_clone_comparison_diffex_{location}_{scna_of_interest}.csv")
  write_csv(diffex, diffex_path)

  return(diffex_path)
}

#' Create a plot visualization
#'
#' @param diffex_path File path
#' @return ggplot2 plot object
#' @export
gse_plot_from_cluster_diffex <- function(diffex_path) {
  sample_id <- str_extract(diffex_path, "SR[RX][0-9]+")

  numbat_dir <- path_split(diffex_path)[[1]][[2]]

  location <- str_extract(diffex_path, "(?<=diffex_).*_segment")

  diffex <-
    diffex_path %>%
    read_csv() %>%
    split(.$clone_comparison) %>%
    map(~ split(.x, .x[["cluster"]])) %>%
    identity()

  
#' Perform clustering analysis
#'
#' @param plot_list Parameter for plot list
#' @param mylabel Parameter for mylabel
#' @return Function result
#' @export
add_cluster_label <- function(plot_list, mylabel) {
    map(plot_list, ~ {
      .x + labs(subtitle = mylabel)
    })
  }

  annotable_cols <- colnames(annotables::grch38)
  annotable_cols <- annotable_cols[!annotable_cols == "symbol"]
  gse_plots <-
    map(diffex, make_gse_plot, sample_id) %>%
    purrr::compact()

  
#' Create a plot visualization
#'
#' @param plot_list Parameter for plot list
#' @return ggplot2 plot object
#' @export
drop_empty_plots <- function(plot_list) {
    #
    plot_content <-
      plot_list %>%
      map(~ {
        dim(.x[["data"]])
      }) %>%
      map_lgl(is.null) %>%
      identity()

    plot_list <- plot_list[!plot_content]
  }

  gse_plots <-
    gse_plots %>%
    map(drop_empty_plots) %>%
    purrr::compact()

  gse_plot_path <- glue("results/{numbat_dir}/{sample_id}_cluster_clone_comparison_diffex_{location}.pdf")

  if (length(gse_plots) > 0) {
    pdf(gse_plot_path)
    gse_plots
    dev.off()
  }

  return(gse_plots)
}

make_clone_distribution_figure <- function(seu_path = NULL, cluster_order = NULL, nb_paths = NULL, clone_simplifications = NULL, group.bys = "SCT_snn_res.0.6", group.by = "SCT_snn_res.0.6", assay = "SCT", height = 10, width = 18, equalize_scna_clones = FALSE, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL, large_clone_comparisons = NULL, scna_of_interest = "1q", min_cells_per_cluster = 50, return_plots = FALSE, split_columns = "clusters", wo_nb = FALSE, plot_path = NULL, heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters"), heatmap_arrangement = c("clusters", "scna"), file_id = NULL) {
	kept_phases <- kept_phases %||% phase_levels
	
	tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
	
	file_id <- file_id %||% fs::path_file(seu_path)
	
	message(sample_id)

    cluster_order <- read_cluster_orders_table(file_id = file_id)[[1]]

	full_seu <- readRDS(seu_path)
	
	    # full_seu <- readRDS(seu_path)
    if ("phase_level" %in% colnames(full_seu@meta.data)){
        full_seu <- full_seu[, !full_seu$phase_level %in% c("other")]
    }
	
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
	
	plot_path <- ifelse(is.null(plot_path), glue("results/{file_slug}_{str_extract(scna_of_interest, '[0-9][a-z]')}_heatmap_phase_scatter_patchwork.pdf"), plot_path)
	
	pdf(plot_path, height = height, width = width)
	on.exit(dev.off(), add = TRUE)

	if(!is.null(cluster_order) & length(cluster_order) > 0){
		iterative <- names(cluster_order)
		} else {
			iterative <- group.bys
		}
	
	
	for(value in iterative){
		# start loop ------------------------------
		
		if(!is.null(cluster_order) & length(cluster_order) > 0) {
			
			seu <- full_seu[, full_seu$clone_opt %in% retained_clones]
			
			single_cluster_order <- cluster_order[[value]]
			
			single_cluster_order <-
				single_cluster_order |>
				dplyr::mutate(order = dplyr::row_number()) %>%
				dplyr::filter(!is.na(clusters)) %>%
				dplyr::mutate(clusters = as.character(clusters))
			
			group.by <- unique(single_cluster_order[["resolution"]])
			
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
				seu %>%
					find_all_markers(metavar = "clusters", seurat_assay = "SCT") %>%
					identity(),
				error = function(e) {
					# seuratTools::find_all_markers uses JoinLayers on the assay and can fail on SCTAssay.
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
			
			cluster_similarity <- find_cluster_pairwise_distance(seu, "clusters")
			
		} else {
				
				group.by = value
			
			seu <- full_seu[, full_seu$clone_opt %in% retained_clones]
			
			cluster_similarity <- find_cluster_pairwise_distance(seu, value)
            # browser()
			
			seu <- find_all_markers(seu, metavar = value, seurat_assay = "SCT")
			
			seu@misc$markers[[value]][["presto"]] <-
				seu@misc$markers[[value]][["presto"]] |> 
				# dplyr::mutate(Cluster = factor(Cluster, levels = levels(seu[[]][[value]]))) |>
				dplyr::mutate(Cluster = factor(Cluster, levels = levels(cluster_similarity$x))) |>
				dplyr::arrange(Cluster)
			
			heatmap_features <-
				seu@misc$markers[[value]][["presto"]]
			
			seu@meta.data$clusters <- factor(seu@meta.data[[value]], levels(cluster_similarity$x))
			
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
			
		}
		
		# cluster similarity ------------------------------
		cluster_similarity_plot <-
			cluster_similarity |> 
			dplyr::mutate(similarity = 1-avg_dist) |> 
			# dplyr::mutate(x = fct_inseq(x), y = fct_inseq(y)) |>
			dplyr::mutate(x = factor(x, levels = rev(levels(x)))) |>
			arrange(x, similarity) |>
			mutate(group_order = forcats::fct_inorder(interaction(x, y))) |> 
			ggplot(mapping = aes(x = x, y = similarity, fill = y, group = group_order)) +
			geom_col(position = position_dodge()) + 
			# geom_text_repel(aes(label = y),
			# 					position = position_dodge(width = 0.9), angle = 0, hjust = 1) +
			coord_flip() + 
			ggsci::scale_fill_igv()
		# end cluster similarity ------------------------------
		
		seu$scna[seu$scna == ""] <- ".diploid"
		seu$scna <- factor(seu$scna)
		# levels(seu$scna)[1] <- "none"
		
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
			# dplyr::rename({{group.by}} := cluster) %>%
			identity()
		
		cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna")) |> 
			dplyr::mutate(cluster_facet = dplyr::case_when(
				str_detect(.data[["clusters"]], "hsp") ~ "stress",
				str_detect(.data[["clusters"]], "hyp") ~ "stress",
				str_detect(.data[["clusters"]], "star") ~ "stress",
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
		
		# asdf
		facet_cell_cycle_plot <-
			cc_data %>%
			ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
			geom_point(size = 0.1, alpha = 0.7) +
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
				strip.text.x = element_blank(),
				legend.title = element_text(size = 16),
				legend.text = element_text(size = 16)
			) +
			guides(
				color = guide_legend(override.aes = list(size = 5))
			) +
			NULL
		
		appender <- function(string) str_wrap(string, width = 40)
		
		labels <- data.frame(scna = unique(seu$scna), label = str_replace(unique(seu$scna), "^$", "diploid"))
		
		single_cluster_order_vec <-
			seu@meta.data %>%
			dplyr::select(clusters, !!group.by) %>%
			dplyr::arrange(clusters, !!sym(group.by)) %>%
			dplyr::select(clusters, !!group.by) |>
			dplyr::distinct(.data[[group.by]], .keep_all = TRUE) |>
			dplyr::mutate(!!group.by := as.character(.data[[group.by]])) |>
			tibble::deframe() |>
			identity()
		
		# full_seu@meta.data[[group.by]] <- factor(full_seu@meta.data[[group.by]], levels = single_cluster_order_vec)
		# levels(full_seu@meta.data[[group.by]]) <- names(single_cluster_order_vec)
		umap_plots <- make_faded_umap_plots(seu, retained_clones, group_by = group.by)
		
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
				"A" = seu_heatmap,
				"B" = facet_cell_cycle_plot,
				"C" = umap_plots,
				"D" = clone_distribution_plot,
				"E" = clone_tree_plot,
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
		} else {
			collage_plots <- list(
				"A" = seu_heatmap,
				"B" = centroid_plot,
				"C" = facet_cell_cycle_plot,
				"D" = clone_distribution_plot,
				"E" = cluster_similarity_plot
			)
# 			layout <- "
#               AAAAAABBDD
#               AAAAAACCDD
#               AAAAAACCEE
#               AAAAAACCEE
#               AAAAAACCEE
#               "
			
			layout <- "
              AAAAAABB
              AAAAAACC
              AAAAAACC
              AAAAAACC
              AAAAAACC
              AAAAAADD
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
	
	return(plot_path)
}