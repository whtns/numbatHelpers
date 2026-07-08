# Plot Functions (111)

#' Perform make table s07 operation
#'
#' @param seu_path File path
#' @param table_path File path
#' @return Function result
#' @export
make_table_s07 <- function(seu_path = "output/seurat/integrated_1q/integrated_seu_1q_complete.rds", table_path = "results/table_s07.csv"){
	seu_integrated <- readRDS(seu_path)
	
	tables_split_by_batch <- janitor::tabyl(seu_integrated@meta.data, clusters, scna, batch) |> 
		adorn_percentages(denominator = "col") |> 
		adorn_pct_formatting() |>
		dplyr::bind_rows(.id = "sample_id") |> 
		identity()
	
	integrated_table <- janitor::tabyl(seu_integrated@meta.data, clusters, scna) |> 
		adorn_percentages(denominator = "col") |> 
		adorn_pct_formatting() |>
		dplyr::mutate(sample_id = "integrated") |> 
		identity()
	
	dplyr::bind_rows(integrated_table, tables_split_by_batch) |> 
		dplyr::arrange(clusters) |> 
		write_csv(table_path) |>
		identity()
	
	return(table_path)
}

#' Perform enrichment analysis
#'
#' @param mytable Parameter for mytable
#' @param mychrom Character string (default: "01")
#' @param myarm Character string (default: "q")
#' @return Enrichment analysis results
#' @export
label_enrichment_by_cis <- function(mytable, mychrom = "01", myarm = "q") {
	test1 <- mytable |> 
		dplyr::mutate(geneID = str_split(geneID, "/")) |> 
		tidyr::unnest(geneID)
	
	genes_1q <-
		find_genes_by_arm(test1$geneID) |> 
		dplyr::filter(seqnames == mychrom, arm == myarm)
	
	test2 <- 
		test1 |> 
		dplyr::filter(test1$geneID %in% genes_1q$symbol) |> 
		dplyr::distinct(cluster, ID, geneID, .keep_all = TRUE) |> 
		dplyr::group_by(cluster, ID) |> 
		summarise(genes_in_cis = paste(geneID, collapse = "/")) |> 
		identity()
	
	test3 <- 
		dplyr::left_join(mytable, test2, by = c("cluster", "ID")) |> 
		dplyr::select(cluster, ID, geneID, genes_in_cis, everything())
}

#' Create a plot visualization
#'
#' @param cohort Character string (default: "READ")
#' @param fullname Character string (default: "READ")
#' @param sig_peaks Parameter for sig peaks
#' @return ggplot2 plot object
#' @export
make_gistic_plot <- function(cohort = "READ", fullname = "READ", sig_peaks){
	
  mycohort = cohort

  peaks <-
  sig_peaks |>
  dplyr::filter(abbreviation == mycohort) |>
  tidyr::pivot_wider(names_from = arm, values_from = is_peak) |>
  unlist() |>
  identity()

	cohort |> 
		TCGAgistic::tcga_gistic_load(source = "Firehose", cnLevel = "all") |> 
		maftools::gisticChromPlot(fdrCutOff = 0.05)
		recordPlot()
	# ggplotify::as.ggplot(gisticChromPlot(cohort))
	title(main = fullname, 
				cex.main = 1,   font.main= 4, col.main= "black"
	)
	
	contig_lens = cumsum(maftools:::getContigLens(build = "hg19"))
	chr_arms = cumsum(get_chr_arms())
	
	boundaries <- par("usr")
	
	rect_col = rgb(red = 0, green = 0, blue = 0, alpha = 0.05)
	
  # browser()
  if(peaks["1q+"] == "TRUE"){
    rect(chr_arms[[1]], boundaries[[3]], chr_arms[[2]], boundaries[[4]], col = rect_col, border = "yellow")
  }

	if(peaks["2p+"] == "TRUE"){
    rect(chr_arms[[2]], boundaries[[3]], chr_arms[[3]], boundaries[[4]], col = rect_col, border = "yellow")
  } 

	if(peaks["6p+"] == "TRUE"){
    rect(chr_arms[[10]], boundaries[[3]], chr_arms[[11]], boundaries[[4]], col = rect_col, border = "yellow")
  }

	if(peaks["16q-"] == "TRUE"){
    rect(chr_arms[[31]], boundaries[[3]], chr_arms[[32]], boundaries[[4]], col = rect_col, border = "yellow")
  } 

}

#' Create a plot visualization
#'
#' @param plot_path File path
#' @return ggplot2 plot object
#' @export
plot_tcga_gistic <- function(plot_path = "results/plot_gistic.pdf") {
	
	tcga_cohorts <- 
		tcga_gistic_available() |> 
		dplyr::filter(CopyNumberLevel == "all") |> 
		dplyr::select(FullName, Cohort) |> 
		dplyr::mutate(FullName = glue("{FullName} ({Cohort})")) |> 
		tibble::deframe() |>
		identity()
	
  possibly_check_tcga_peaks <- possibly(check_tcga_peaks)

	sig_peaks <- 
		tcga_cohorts |> 
		map(possibly_check_tcga_peaks) |>
    dplyr::bind_rows(.id = "cohort") |>
    tidyr::pivot_longer(cols = c("X1q", "X2p", "X6p", "X16q"), names_to = "arm", values_to = "is_peak") |>
    dplyr::mutate(arm = dplyr::case_when(
      arm == "X1q" ~ "1q+",
      arm == "X2p" ~ "2p+",
      arm == "X6p" ~ "6p+",
      arm == "X16q" ~ "16q-",
    )) |>
    identity()

	possibly_gistic_plot <- possibly(make_gistic_plot)
	
	possibly_gistic_plot("BLCA", "BLCA", sig_peaks)
	
	pdf(plot_path, h = 3)
	gistic_plots <- 
		tcga_cohorts |> 
		imap(possibly_gistic_plot, sig_peaks)
	dev.off()
	
	return(plot_path)
}

#' Extract or pull specific data elements
#'
#' @param numbat_rds_files File path
#' @param table_path File path
#' @return Extracted data elements
#' @export
extract_full_segmentation <- function(numbat_rds_files, table_path = "results/table_s11.xlsx") {
	
	segmentation_table <- 
		numbat_rds_files |> 
		set_names(str_extract(numbat_rds_files, "SR[RX][0-9]+")) |> 
		map(retrieve_segmentation) |> 
		# dplyr::bind_rows(.id = "sample_id") |> 
		identity()
	
	writexl::write_xlsx(segmentation_table, table_path)
}

plot_seu_marker_heatmap_by_scna <- function(seu_path = NULL, cluster_order = NULL, nb_paths = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", height = 10, width = 18, equalize_scna_clones = FALSE, phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL, rb_scna_samples, large_clone_comparisons, scna_of_interest = "1q", min_cells_per_cluster = 50, return_plots = FALSE, split_columns = "clusters",  tmp_plot_path = FALSE, hypoxia_expr = NULL, run_hypoxia_clustering = FALSE, cluster_resolutions = seq(0.2, 1, by = 0.2)) {
  kept_phases <- kept_phases %||% phase_levels

  file_id <- fs::path_file(seu_path)
  
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
  
  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
  
  message(file_id)

  cluster_order <- if(is.null(cluster_order[[file_id]])){
    dummy_cluster_order(seu_path, kept_phases = kept_phases)
  } else {
    cluster_order[[file_id]]
  }

  full_seu <- readRDS(seu_path)

  # Optional: subset the full Seurat object by an expression evaluated in meta.data
  if (!is.null(hypoxia_expr)) {
    keep_cells <- tryCatch({
      with(full_seu@meta.data, eval(parse(text = hypoxia_expr)))
    }, error = function(e) {
      warning("Failed to evaluate hypoxia_expr ('", hypoxia_expr, "') on ", seu_path, ": ", e$message)
      NULL
    })
    if (is.logical(keep_cells) && length(keep_cells) == ncol(full_seu)) {
      kept_n <- sum(keep_cells, na.rm = TRUE)
      full_seu <- full_seu[, which(keep_cells)]
      message("Subsetting full_seu by hypoxia_expr '", hypoxia_expr, "' -> kept ", kept_n, " cells")
    } else if (is.numeric(keep_cells)) {
      full_seu <- full_seu[, keep_cells]
      message("Subsetting full_seu by numeric index from hypoxia_expr, kept ", ncol(full_seu), " cells")
    } else {
      warning("hypoxia_expr returned unexpected value; skipping subsetting")
    }
  }

  # Optionally re-run clustering on the hypoxia-subsetted object
  if (!is.null(hypoxia_expr) && isTRUE(run_hypoxia_clustering)) {
    message("Re-running clustering on hypoxia-subsetted object at resolutions: ", paste(cluster_resolutions, collapse = ", "))
    full_seu <- seurat_cluster(full_seu, resolution = cluster_resolutions, seurat_assay = assay)
  }

  if (!sample_id %in% names(large_clone_comparisons)) {
    sample_id <- tumor_id
  }

  # subset by retained clones ------------------------------
  clone_comparisons <- names(large_clone_comparisons[[sample_id]])

  clone_comparison <- clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
  retained_clones <- clone_comparison %>%
    str_extract("[0-9]_v_[0-9]") %>%
    str_split("_v_", simplify = TRUE)

  if(!is.null(nb_paths)){
  	nb_paths <- nb_paths %>%
  		set_names(str_extract(., "SR[RX][0-9]+"))
  	
  	nb_path <- nb_paths[[tumor_id]]	
  }

  plot_paths <- vector(mode = "list", length = length(cluster_order))
  names(plot_paths) <- names(cluster_order)

  file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  plot_path <- if(!tmp_plot_path){
    glue("results/{file_slug}_{scna_of_interest}_heatmap_phase_scatter_patchwork.pdf")
  } else {
    fs::dir_create("tmp")  # tempfile() does not create the tmpdir; ggsave errors "Cannot find directory 'tmp'"
    tempfile(tmpdir = "tmp", fileext = ".pdf")
  }

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
# 
#     large_enough_clusters <-
#       seu@meta.data %>%
#       dplyr::group_by(clusters) %>%
#       dplyr::count() |>
#       dplyr::filter(n >= min_cells_per_cluster) %>%
#       dplyr::pull(clusters)
# 
#     seu <- seu[, seu$clusters %in% large_enough_clusters]

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

    cc_groupby <- c("G2M.Score", "S.Score", "scna", "clusters")
    cc_groupby <- cc_groupby[sapply(cc_groupby, function(col) {
      vals <- seu@meta.data[[col]]
      is.null(vals) || length(unique(vals[!is.na(vals)])) > 1L
    })]

    seu_heatmap <- ggplotify::as.ggplot(
      seu_complex_heatmap(seu,
        features = heatmap_features$Gene.Name,
        group.by = cc_groupby,
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

    #

    clone_ratio <- janitor::tabyl(as.character(seu$scna))$percent[[2]]

    comparison_scna <-
      janitor::tabyl(as.character(seu$scna))[2, 1]

    clone_distribution_plot <- plot_distribution_of_clones_across_clusters(
      seu,
      seu_name = glue("{tumor_id} {comparison_scna}"), var_x = "scna", var_y = "clusters", signif = TRUE, plot_type = "clone"
    )

    # umap_plots <- DimPlot(full_seu, group.by = c("scna", "clusters"), combine = FALSE) %>%
    # 	# map(~(.x + theme(legend.position = "bottom"))) %>%
    # 	wrap_plots(ncol = 1)
    # full_seu$clusters
    # full_seu[[group.by]] <-
    full_seu@meta.data[[group.by]] <- factor(full_seu@meta.data[[group.by]], levels = single_cluster_order_vec)
    levels(full_seu@meta.data[[group.by]]) <- names(single_cluster_order_vec)
    umap_plots <- make_faded_umap_plots(full_seu, retained_clones, group_by = group.by)

    if (!is.null(nb_path)) {
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
      layout <- "
              AAAAABBCC
              AAAAABBCC
              AAAAABBDD
              AAAAABBDD
              AAAAABBDD
      "

      collage_plots <- list(
        "seu_heatmap" = seu_heatmap,
        "facet_cell_cycle_plot" = facet_cell_cycle_plot,
        "centroid_plot" = centroid_plot,
        "clone_distribution_plot" = clone_distribution_plot
      )

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

subset_seu_by_expression <- function(seu_path, hypoxia_expr = NULL,
run_hypoxia_clustering = FALSE, cluster_resolutions = seq(0.2, 1, by = 0.2), assay = "SCT", slug = "hypoxia_low") {
  if (is.na(seu_path)) return(NA_character_)

  seu <- readRDS(seu_path)

    # optionally subset the Seurat object by hypoxia expression string (evaluated in @meta.data)
  if (!is.null(hypoxia_expr)) {
    keep_cells <- tryCatch({
      with(seu@meta.data, eval(parse(text = hypoxia_expr)))
    }, error = function(e) {
      warning("Failed to evaluate hypoxia_expr ('", hypoxia_expr, "') on ", seu_path, ": ", e$message)
      NULL
    })
    if (is.logical(keep_cells) && length(keep_cells) == ncol(seu)) {
      seu <- seu[, which(keep_cells)]
      message("Subsetting seu by hypoxia_expr '", hypoxia_expr, "' -> kept ", sum(keep_cells, na.rm = TRUE), " cells")
    } else if (is.numeric(keep_cells)) {
      seu <- seu[, keep_cells]
      message("Subsetting seu by numeric index from hypoxia_expr, kept ", ncol(seu), " cells")
    } else {
      warning("hypoxia_expr returned unexpected value; skipping subsetting")
    }
  }

  # Optionally re-run clustering on the hypoxia-subsetted object
  if (!is.null(hypoxia_expr) && isTRUE(run_hypoxia_clustering)) {
    n_cells <- ncol(seu)
    if (n_cells < 10) {
      warning("Too few cells (", n_cells, ") after hypoxia subsetting; skipping clustering")
    } else {
      message("Re-running clustering on hypoxia-subsetted object at resolutions: ", paste(cluster_resolutions, collapse = ", "))
      Seurat::DefaultAssay(seu) <- assay
      if (!"data" %in% SeuratObject::Layers(seu[[assay]])) {
        seu <- Seurat::NormalizeData(seu, assay = assay, verbose = FALSE)
      }
      k_param <- min(20L, n_cells - 1L)
      graph_names <- paste0(assay, c("_nn", "_snn"))
      seu <- Seurat::FindNeighbors(seu, dims = 1:30, reduction = "pca", graph.name = graph_names, k.param = k_param)
      for (res in cluster_resolutions) {
        seu <- Seurat::FindClusters(seu, graph.name = paste0(assay, "_snn"), resolution = res)
      }
      message("Re-running UMAP on hypoxia-subsetted object")
      n_neighbors <- min(30L, n_cells - 1L)
      seu <- Seurat::RunUMAP(seu, reduction = "pca", dims = 1:30, n.neighbors = n_neighbors)
      cat("DEBUG: After clustering, computing markers for all resolutions...\n")
      seu <- tryCatch(
        find_all_markers(seu, seurat_assay = assay),
        error = function(e) {
          if (grepl("JoinLayers", conditionMessage(e), fixed = TRUE)) {
            warning("SCT marker JoinLayers failed; using PrepSCTFindMarkers fallback.")
            seu <- PrepSCTFindMarkers(seu)
            seu@misc$markers[["clusters"]] <- FindAllMarkers(seu, assay = assay, verbose = FALSE)
            return(seu)
          }
          stop(e)
        }
      )
      cat("DEBUG: Markers computed. Keys in seu@misc$markers:", paste(names(seu@misc$markers), collapse=", "), "\n")
      # When clustering used SCT, gene_snn_res.0.2 is inherited from the parent and
      # safe_plot_markers in plot_effect_of_filtering looks under that key — compute
      # gene-assay markers for it here. Skip when assay == "gene" because the main
      # find_all_markers call above already covered all gene_snn_res.* columns.
      if (assay != "gene" && "gene_snn_res.0.2" %in% colnames(seu[[]])) {
        seu <- tryCatch(
          find_all_markers(seu, metavar = "gene_snn_res.0.2", seurat_assay = "gene"),
          error = function(e) {
            warning("gene-assay markers for gene_snn_res.0.2 failed: ", conditionMessage(e))
            seu
          }
        )
      }
    }
  }

  new_filepath <- str_replace(seu_path, "_seu.*.rds", paste0("_", slug, "_seu.rds"))
  add_hash_metadata(new_filepath, seu = seu)
  sample_id <- stringr::str_extract(new_filepath, "SR[RX][0-9]+")
  save_cell_barcodes_to_db(new_filepath, sample_id, slug, colnames(seu))
  return(new_filepath)
}

plot_seu_marker_heatmap <- function(seu_path = NULL, cluster_order = NULL, 
nb_paths = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", 
assay = "SCT", label = "_filtered_", height = 10, width = 18, 
equalize_scna_clones = FALSE, 
phase_levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), 
kept_phases = NULL, tmp_plot_path = FALSE, hypoxia_expr = NULL,
run_hypoxia_clustering = FALSE, cluster_resolutions = seq(0.2, 1, by = 0.2)) {
  kept_phases <- kept_phases %||% phase_levels

  if (is.na(seu_path)) return(NA_character_)

  file_id <- fs::path_file(seu_path)

  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
  
  message(file_id)

  cluster_order_list <- if(is.null(cluster_order[[file_id]])){
    dummy_cluster_order(seu_path, kept_phases = kept_phases)
  } else {
    cluster_order[[file_id]]
  }
  cluster_order <- cluster_order_list[["0"]] %||% cluster_order_list[[1]]

  seu <- readRDS(seu_path)

  nb_paths <- nb_paths %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  nb_path <- nb_paths[[tumor_id]]

  if (!is.null(cluster_order)) {
    group.by <- unique(cluster_order$resolution)
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

  if (!is.null(cluster_order)) {
    group.by <- unique(cluster_order$resolution)

    # Auto-generate the cell-cycle phase label per cluster from THIS object's own
    # cells (S.Score/G2M.Score/Phase + marker programs), instead of joining the
    # phase assignments from data/scna_cluster_order.csv, which go stale when
    # clustering changes. Emits g1/s/s_star/g2_m/pm/hsp in the existing
    # phase_levels vocabulary. See auto_phase_level()/assign_auto_phase_clusters().
    seu <- assign_auto_phase_clusters(seu, group.by)

    phase_levels <- phase_levels[phase_levels %in% unique(as.character(seu@meta.data$phase_level))]

    seu <- seu[, seu$phase_level %in% kept_phases]
    # A marker heatmap needs >=2 clusters to find markers between. Degenerate
    # low/high-hypoxia subsets can collapse to <2; skip the whole collage for
    # them (returns NA_character_, same convention as the is.na(seu_path) guard
    # above) instead of proceeding to an inevitable ComplexHeatmap/ggplot failure.
    if (length(unique(seu@meta.data$clusters)) < 2L) {
      warning("Fewer than 2 clusters in ", seu_path,
              "; skipping marker heatmap collage.")
      return(NA_character_)
    }
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

    cluster_order_vec <-
      seu@meta.data %>%
      dplyr::select(clusters, !!group.by) %>%
      dplyr::arrange(clusters, !!sym(group.by)) %>%
      dplyr::pull(!!group.by) %>%
      unique() %>%
      as.character() %>%
      identity()

    heatmap_features[["Cluster"]] <-
      factor(heatmap_features[["Cluster"]], levels = levels(seu@meta.data$clusters))

    heatmap_features <-
      heatmap_features %>%
      dplyr::arrange(Cluster) %>%
      dplyr::group_by(Cluster) %>%
      slice_max(Average.Log.Fold.Change, n = 5) %>%
      identity()
  } else {
    heatmap_features <-
      seu@misc$markers[[group.by]][["presto"]]

    cluster_order <- levels(seu@meta.data[[group.by]]) %>%
      set_names(.)

    seu@meta.data[[group.by]] <-
      factor(seu@meta.data[[group.by]], levels = cluster_order)

    group_by_clusters <- seu@meta.data[[group.by]]

    seu@meta.data$clusters <- names(cluster_order[group_by_clusters])

    seu@meta.data$clusters <- factor(seu@meta.data$clusters, levels = unique(setNames(names(cluster_order), cluster_order)[levels(seu@meta.data[[group.by]])]))

    heatmap_features <-
      heatmap_features %>%
      dplyr::arrange(Cluster) %>%
      group_by(Cluster) %>%
      slice_head(n = 6) %>%
      dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
      dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
      identity()
  }

  seu$scna <- factor(seu$scna)
  levels(seu$scna)[1] <- "none"

  # Display numbat clones (clone_opt) rather than the simplified SCNA labels,
  # which are frequently blank for these subsets. NA / unassigned cells -> "none".
  if ("clone_opt" %in% colnames(seu@meta.data)) {
    .clone <- as.character(seu$clone_opt)
    .clone[is.na(.clone) | .clone == ""] <- "none"
    seu$clone <- factor(.clone)
  } else {
    seu$clone <- factor("none")
  }

  giotti_genes <- read_giotti_genes()

  heatmap_features <-
    heatmap_features %>%
    dplyr::ungroup() %>%
    left_join(giotti_genes, by = c("Gene.Name" = "symbol")) %>%
    # select(Gene.Name, term) %>%
    dplyr::mutate(term = replace_na(term, "")) %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE)

  # No marker features survived the VariableFeatures/logFC filters -> a 0-row
  # heatmap makes ComplexHeatmap set row_order to NULL and as.ggplot() yield a
  # non-plot, which then fails `+ labs()`. Skip the collage for this sample.
  if (is.null(heatmap_features) || nrow(heatmap_features) < 1L) {
    warning("No marker heatmap features for ", seu_path, "; skipping collage.")
    return(NA_character_)
  }

  row_ha <- ComplexHeatmap::rowAnnotation(term = rev(heatmap_features$term))

  cc_groupby <- c("G2M.Score", "S.Score", "clone", "clusters")
  cc_groupby <- cc_groupby[sapply(cc_groupby, function(col) {
    vals <- seu@meta.data[[col]]
    is.null(vals) || length(unique(vals[!is.na(vals)])) > 1L
  })]

  seu_heatmap <- ggplotify::as.ggplot(
    seu_complex_heatmap(seu,
      features = heatmap_features$Gene.Name,
      group.by = cc_groupby,
      col_arrangement = c("clusters", "clone"),
      cluster_rows = FALSE,
      column_split = sort(seu@meta.data$clusters),
      row_split = rev(heatmap_features$Cluster),
      row_title_rot = 0,
      # row_split = sort(seu@meta.data$clusters)
    )
  ) +
    labs(title = sample_id) +
    theme()


  #
  labels <- data.frame(clusters = unique(seu[[]][["clusters"]]), label = unique(seu[[]][["clusters"]])) %>%
    # dplyr::rename({{group.by}} := cluster) %>%
    identity()

  cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "clone"))

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
    NULL


  facet_cell_cycle_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["clone"]])) +
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

  labels <- data.frame(clone = unique(seu$clone), label = as.character(unique(seu$clone)))

  clone_distribution_plot <-
    plot_distribution_of_clones_across_clusters(seu, tumor_id, var_x = "clone", var_y = "clusters")

  umap_plots <- DimPlot(seu, group.by = c("clone", "clusters"), combine = FALSE) %>%
    # map(~(.x + theme(legend.position = "bottom"))) %>%
    wrap_plots(ncol = 1)

  clone_tree_plot <- if (!is.null(nb_path)) {
    plot_clone_tree(seu, tumor_id, nb_path, clone_simplifications, sample_id = sample_id, legend = FALSE, horizontal = TRUE)
  } else {
    NULL
  }

  # Only use the clone-tree layout if we actually got a plot back; otherwise fall
  # through to the no-tree layout rather than feeding a NULL panel to wrap_plots().
  if (!is.null(clone_tree_plot)) {
    collage_plots <- list(
    	"A" = seu_heatmap,
    	"B" = facet_cell_cycle_plot,
    	"C" = umap_plots,
    	"D" = clone_distribution_plot,
    	"E" = clone_tree_plot)

    layout <- "
            AAAAAAAAAAAEEEECCCC
            AAAAAAAAAAABBBBCCCC
            AAAAAAAAAAABBBBDDDD
            AAAAAAAAAAABBBBDDDD
            AAAAAAAAAAABBBBDDDD
            "

    wrap_plots(collage_plots) +
      # plot_layout(widths = c(16, 4)) +
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

    collage_plots <- list(seu_heatmap, facet_cell_cycle_plot, centroid_plot, clone_distribution_plot)

    wrap_plots(collage_plots) +
      # plot_layout(widths = c(16, 4)) +
      plot_layout(design = layout) +
      plot_annotation(tag_levels = "A") +
      NULL
  }

  file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  plot_path <- if(!tmp_plot_path){
    glue("results/{file_slug}_{label}heatmap_phase_scatter_patchwork.pdf")
  } else {
    fs::dir_create("tmp")  # tempfile() does not create the tmpdir; ggsave errors "Cannot find directory 'tmp'"
    tempfile(tmpdir = "tmp", fileext = ".pdf")
  }

  ggsave(plot_path, height = height, width = width)
}


dummy_cluster_order <- function(seu_path, kept_phases = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hsp", "hypoxia", "other", "s_star"), integrated = FALSE) {
  file_id <- fs::path_file(seu_path)
  
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")
  
  sample_id <- str_extract(seu_path, "SR[RX][0-9]+")

  if(integrated){
    resolution = "integrated_snn_res.0.4"
  } else {
    resolution = "SCT_snn_res.0.6"
  }

  seu <- readRDS(seu_path)

  clusters <- seu[[]][[resolution]] %>% unique() %>% as.character()

  cluster_order <-
  tibble::tribble(
    ~resolution, ~file_id, ~sample_id, ~preferred_resolution, ~tumor_id, ~phase, ~clusters,
    resolution, file_id, sample_id, 0, tumor_id, sample(kept_phases, length(clusters), replace = TRUE), clusters
  )  |> 
  tidyr::unnest()

  return(list("0" = cluster_order))
}

numeric_col_fun <- function(myvec, color) {
      circlize::colorRamp2(range(myvec), c("white", color))
    }

#' Plot hypoxia score for a Seurat object
#'
#' Create a simple scatter-style plot of the per-cell hypoxia score stored in
#' the Seurat object's metadata. This function expects a column named
#' `hypoxia_score` in `seu@meta.data`.
#'
#' @param seu A Seurat object containing `hypoxia_score` in `@meta.data`.
#' @param mytitle Optional plot title (default: "").
#' @return A ggplot2 object showing per-cell hypoxia scores.
#' @export
plot_hypoxia_score <- function(seu_path, threshold = 0.5) {
    if (is.na(seu_path)) return(NA_character_)
    seu <- readRDS(seu_path)
    sample_id <- str_extract(seu_path, "SR[RX][0-9]+")
    p <- seu@meta.data |> 
        tibble::rownames_to_column("cell") |> 
        dplyr::arrange(hypoxia_score) |> 
        dplyr::mutate(cell = factor(cell, levels = cell)) |> 
        ggplot(aes(x = cell, y = hypoxia_score)) +
        geom_point() + 
        geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
        labs(title = sample_id, x = "Cells ordered by hypoxia score", y = "Hypoxia score")
    fs::dir_create("tmp")  # tempfile() does not create the tmpdir
    out <- tempfile(tmpdir = "tmp", fileext = ".pdf")
    ggplot2::ggsave(filename = out, plot = p, width = 8, height = 3)
    out
}