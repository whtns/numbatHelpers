# Plot Functions (129)
#' Perform differential expression analysis
#'
#' @param seu_path File path
#' @param numbat_rds_files File path
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @param location Character string (default: "cis")
#' @return Differential expression results
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O

find_diffex_clones <- function(seu_path, numbat_rds_files, large_clone_comparisons, location = "cis") {
  
  
  #
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  numbat_rds_files <- numbat_rds_files %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  mynb <- readRDS(numbat_rds_files[[tumor_id]])

  seu <- readRDS(seu_path)

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  possible_clone_comparison <- possibly(make_clone_comparison)
  
  clone_comparisons <- large_clone_comparisons[[sample_id]]
  
  actual_clone_comparisons <- split_vector_into_overlapping_chunks(rev(sort(unique(seu$clone_opt))), 2) |> 
  	map(paste, collapse = "_v_")
  
  clone_comparisons <- map(actual_clone_comparisons, \(x) clone_comparisons[str_detect(names(clone_comparisons), x)]) |> 
  	unlist()

  diffex <- imap(clone_comparisons, possible_clone_comparison, seu, mynb, location = location) %>%
    purrr::compact()

  return(diffex)
}

#' Perform differential expression analysis
#'
#' @param seu_path File path
#' @param phase Character string (default: "g1")
#' @param scna_of_interest Character string (default: "2p")
#' @param numbat_rds_files File path
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @param location Character string (default: "cis")
#' @return Differential expression results
#' @export
find_diffex_clones_in_phase <- function(seu_path, phase = "g1", scna_of_interest = "2p", numbat_rds_files, large_clone_comparisons, location = "cis") {
  
  
  #
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  numbat_rds_files <- numbat_rds_files %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  mynb <- numbat_rds_files[[tumor_id]]

  # Some samples in the phase-subset input lack a standalone per-sample filtered
  # seu on disk; skip rather than crash readRDS with "cannot open the connection".
  seu_file <- glue("output/seurat/{sample_id}_filtered_seu.rds")
  if (!fs::file_exists(seu_file)) {
    warning("Missing filtered seu for ", sample_id, " (", seu_file, "); skipping phase diffex.")
    return(NULL)
  }
  seu <- readRDS(seu_file)

  seu <- seu[, seu$phase_level %in% phase]

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  possible_clone_comparison <- possibly(make_clone_comparison)

  large_clone_comparisons <- large_clone_comparisons[[sample_id]]

  large_clone_comparisons <- large_clone_comparisons[str_detect(names(large_clone_comparisons), scna_of_interest)]

  diffex <- imap(large_clone_comparisons, possible_clone_comparison, seu, mynb, location = location) %>%
    purrr::compact()

  #

  diffex <- diffex[str_detect(names(diffex), scna_of_interest)]
  enrichments <- diffex |>
    purrr::list_flatten() |>
    map(enrich_diffex)

  return(diffex)
}

plot_seu_clusters_and_markers <- function(seu_path, cluster_order, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star")) {
	file_id <- fs::path_file(seu_path)
	
	tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
	
	message(file_id)
	cluster_order_list <- cluster_order[[file_id]]
	cluster_order <- if (!is.null(cluster_order_list)) cluster_order_list[["0"]] %||% cluster_order_list[[1]] else NULL

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

  return(
    plot_path
  )
}


plot_plae_celltype_expression <- function(mygenes = c("RXRG", "NRL"), plot_type = "box") {
  #

  celltypes <- c("Retinal Ganglion Cell", "Amacrine Cell", "Horizontal Cell", "RPC", "Early RPC", "Muller Glia", "Bipolar Cell", "Late RPC", "Neurogenic Cell", "B-Cell", "Rod", "Photoreceptor Precursor", "Cone", "RPE", "Microglia", "Red Blood Cell", "Astrocyte", "Rod Bipolar Cell")

  # pseudo_meta <- read_tsv("/dataVolume/storage/scEiad/human_pseudobulk/4000-counts-universe-study_accession-scANVIprojection-15-5-20-50-0.1-CellType-Homo_sapiens.meta.tsv.gz") %>%
  #   tidyr::unite(study_type, study_accession, CellType)

  sub_annotable <-
    annotables::grch38 %>%
    dplyr::filter(symbol %in% mygenes)

  pseudo_counts <-
    "data/plae_pseudobulk_counts.csv" %>%
    read_csv() %>%
    dplyr::inner_join(sub_annotable, by = c("Gene" = "ensgene"), relationship = "many-to-many") %>%
    dplyr::distinct(study, type, symbol, .keep_all = TRUE) %>%
    # dplyr::filter(type %in% celltypes) %>%
    identity()

  if (plot_type == "box") {
    exp_plot <- ggplot(pseudo_counts, aes(
      x = type,
      y = counts,
      color = study
    )) +
      geom_boxplot(color = "black", outlier.shape = NA) +
      ggbeeswarm::geom_quasirandom(groupOnX = TRUE) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_radius(range = c(2, 6)) +
      scale_colour_manual(values = rep(c(pals::alphabet() %>% base::unname()), 20)) +
      theme(legend.position = "bottom") +
      facet_wrap(ncol = 2, scales = "free_y", ~symbol) +
      NULL
  } else if (plot_type == "hmap") {
    exp_plot <-
      pseudo_counts %>%
      group_by(symbol, type) %>%
      dplyr::summarize(median_counts = median(counts)) %>%
      ggplot(aes(
        x = symbol,
        y = type,
        fill = median_counts
      )) +
      geom_tile() +
      # scale_fill_gradient(name = "median_count", trans = "log") +
      NULL
  }

  return(exp_plot)
}

plot_phase_distribution_of_all_samples_by_scna <- function(seu_paths, selected_samples = c(
                                                             "SRX10264520", "SRX10264526", "SRX11133594",
                                                             "SRX11133593", "SRX11133592", "SRX11133585",
                                                             "SRX14116947"
                                                           )) {
  #

  plot_phase_distribution_by_scna <- function(seu_path, seu_name) {
    seu <- readRDS(seu_path)
    seu$Phase <- factor(seu$Phase, levels = c("G1", "S", "G2M"))
    plot_distribution_of_clones_across_clusters(seu, seu_name, var_x = "scna", var_y = "Phase", both_ways = FALSE)
  }

  seu_paths <-
    seu_paths %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  if (!is.null(selected_samples)) {
    seu_paths <- seu_paths[selected_samples]
  }

  seu_plots <-
    seu_paths %>%
    imap(plot_phase_distribution_by_scna)

  test0 <- wrap_plots(seu_plots, ncol = 3) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Proliferation varies directly with scna abundance",
      theme = theme(plot.title = element_text(size = 18))
    )

  plot_path <- "results/relation_between_scna_and_phase.pdf"
  ggsave(plot_path, height = 8, width = 12)

  return(plot_path)
}

plot_mt_v_nUMI <- function(study_cell_stats, plot_path = tempfile(fileext = ".pdf")) {
	# percent_mito_per_read <-
	study_cell_stats %>%
		dplyr::filter(nCount_gene > 1000) %>%
		# dplyr::mutate(mito_per_read = nCount_gene/`percent.mt`) %>%
		ggplot(aes(x = `percent.mt`, y = nCount_gene, group = sample_id)) +
		geom_hex(bins = 70) +
		facet_wrap(~sample_id) +
		geom_vline(xintercept = 5, linetype = "dotted") +
		geom_hline(yintercept = 1e3, linetype = "dotted") +
		# scale_y_discrete(limits=rev, expand = expansion(add = c(0.55, mito_expansion))) +
		labs(title = "percent mito per cell") +
		scale_y_continuous(limits = c(0, 1e5)) +
		xlim(0, 50) +
		ylim(0, 2e5) +
		ggplot2::annotate("rect",
											xmin = 0, xmax = 5, ymin = 1e3, ymax = 2e5,
											alpha = .2, color = "yellow"
		) +
		NULL
	
	ggsave(plot_path, height = 6, width = 8)
}

#' Title
#' @export
plot_markers_by_cell_cycle <- plot_cluster_markers_by_cell_type <- function(seu, checked_cluster_markers) {
  cluster_plots <- map(checked_cluster_markers, ~ VlnPlot(seu, features = .x, group.by = "Phase"))

  cluster_plots <- map2(cluster_plots, names(checked_cluster_markers), ~ (.x + labs(subtitle = .y)))

  return(cluster_plots)
}

plot_gene_clone_trend <- function(seu, mygenes = c("CRABP2", "MEG3")) {
  gene_df <-
    FetchData(seu, vars = c(mygenes, "cluster_clone")) %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(-c("cell", "cluster_clone"), names_to = "gene", values_to = "counts")

  ggplot(gene_df, aes(cluster_clone, counts)) +
    geom_jitter(width = 0.1) +
    facet_wrap(~gene, ncol = 1) +
    NULL
}

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
		
		# seu <- drop_mt_cluster(seu, group.by = group.by)
		
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
			
			seu <- assign_auto_phase_clusters(seu, group.by)

			phase_levels <- phase_levels[phase_levels %in% unique(as.character(seu@meta.data$phase_level))]
			
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
				factor(heatmap_features[["Cluster"]], levels = levels(seu@meta.data$clusters))
			
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
			# dplyr::rename({{group.by}} := cluster) %>%
			identity()
		
		# cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))
		# 
		# centroid_data <-
		# 	cc_data %>%
		# 	dplyr::group_by(clusters) %>%
		# 	dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
		# 	dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
		# 	dplyr::mutate(centroid = "centroids") %>%
		# 	identity()
		# 
		# centroid_plot <-
		# 	cc_data %>%
		# 	ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["clusters"]])) +
		# 	geom_point(size = 0.1) +
		# 	theme_light() +
		# 	theme(
		# 		strip.background = element_blank(),
		# 		strip.text.x = element_blank()
		# 	) +
		# 	geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
		# 	guides(fill = "none", color = "none") +
		# 	NULL
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
		
		# 
		
		
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
		
		#
		
		clone_ratio <- janitor::tabyl(as.character(seu$scna))$percent[[2]]
		
		comparison_scna <-
			janitor::tabyl(as.character(seu$scna))[2, 1]
		
		clone_distribution_plot <- plot_distribution_of_clones_across_clusters(
			seu,
			seu_name = glue("{tumor_id} {comparison_scna}"), var_x = "scna", var_y = "clusters", signif = TRUE, plot_type = "clone"
		)
		
		full_seu@meta.data[[group.by]] <- factor(full_seu@meta.data[[group.by]], levels = single_cluster_order_vec)
		levels(full_seu@meta.data[[group.by]]) <- names(single_cluster_order_vec)
		# umap_plots <- make_faded_umap_plots(full_seu, retained_clones, group_by = group.by)
		
		if (!is.null(nb_paths)) {
			clone_tree_plot <-
				plot_clone_tree(seu, tumor_id, nb_path, clone_simplifications, sample_id = sample_id, legend = FALSE, horizontal = FALSE)
			
			collage_plots <- list(
				"A" = seu_heatmap,
				"B" = clone_tree_plot,
				"C" = centroid_plot,
				"D" = facet_cell_cycle_plot,
				# "E" = umap_plots,
				# "F" = clone_distribution_plot,
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
				# plot_layout(widths = c(16, 4)) +
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