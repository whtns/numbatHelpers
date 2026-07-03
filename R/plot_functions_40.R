# Plot Functions (139)
#' Create a plot visualization
#'
#' @param unfiltered_seu_path File path
#' @param filtered_seu_path File path
#' @param group.by Character string (default: "SCT_snn_res.0.6")
#' @return ggplot2 plot object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging
# - rownames_round_trip: Avoid unnecessary rownames conversions - keep as column when possible

plot_effect_of_filtering_old <- function(unfiltered_seu_path, filtered_seu_path, group.by = "SCT_snn_res.0.6") {
  
	#
	
	plot_list <- list()
	
	sample_id <- str_extract(unfiltered_seu_path, "SR[RX][0-9]+")
	
	unfiltered_seu <- unfiltered_seu_path
	
	filtered_seu <- filtered_seu_path
	
	fs::dir_create("results/effect_of_filtering")
	
	# distribution ------------------------------
	dir_create("results/effect_of_filtering/distribution")
	plot_distribution_of_clones_across_clusters(filtered_seu, sample_id, var_x = "scna", var_y = group.by)
	fs::dir_create(glue("results/effect_of_filtering/distribution/filtered/"))
	plot_path <- glue("results/effect_of_filtering/distribution/filtered/{sample_id}_filtered_distribution.pdf")
	plot_list["filtered_distribution"] <- plot_path
	ggsave(plot_path, height = 4, width = 8)
	
	filtered_dist_tables <- table_distribution_of_clones_across_clusters(filtered_seu, sample_id, clusters = group.by)
	
	table_path <- glue("results/effect_of_filtering/distribution/filtered/{sample_id}_filtered_distribution.xlsx")
	plot_list["filtered_distribution_tables"] <- table_path
	writexl::write_xlsx(filtered_dist_tables, table_path)
	
	plot_distribution_of_clones_across_clusters(unfiltered_seu, sample_id, var_y = group.by)
	fs::dir_create(glue("results/effect_of_filtering/distribution/filtered"))
	plot_path <- glue("results/effect_of_filtering/distribution/filtered/{sample_id}_unfiltered_distribution.pdf")
	plot_list["unfiltered_distribution"] <- plot_path
	ggsave(plot_path, height = 4, width = 8)
	
	unfiltered_dist_tables <- table_distribution_of_clones_across_clusters(unfiltered_seu, sample_id, clusters = group.by)
	
	table_path <- glue("results/effect_of_filtering/distribution/filtered/{sample_id}_unfiltered_distribution.xlsx")
	plot_list["unfiltered_distribution_tables"] <- table_path
	writexl::write_xlsx(unfiltered_dist_tables, table_path)
	
	# abbreviation markers ------------------------------
	dir_create("results/effect_of_filtering/abbreviation")
	(plot_markers(filtered_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
			labs(title = "filtered")) +
		(plot_markers(unfiltered_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
		 	labs(title = "unfiltered")) +
		plot_annotation(title = sample_id)
	
	plot_path <- glue("results/effect_of_filtering/abbreviation/{sample_id}_abbreviation_markers.pdf")
	plot_list["abbreviation_markers"] <- plot_path
	ggsave(plot_path, height = 12, width = 15)
	
	filtered_marker_tables <- table_cluster_markers(filtered_seu)
	
	table_path <- glue("results/effect_of_filtering/abbreviation/{sample_id}_filtered_markers.xlsx")
	plot_list["filtered_marker_tables"] <- table_path
	writexl::write_xlsx(filtered_marker_tables, table_path)
	
	unfiltered_marker_tables <- table_cluster_markers(unfiltered_seu) %>%
		purrr::compact()
	
	table_path <- glue("results/effect_of_filtering/abbreviation/{sample_id}_unfiltered_markers.xlsx")
	plot_list["unfiltered_marker_tables"] <- table_path
	writexl::write_xlsx(unfiltered_marker_tables, table_path)
	
	heatmap_features <-
		table_cluster_markers(unfiltered_seu) %>%
		pluck(group.by) %>%
		group_by(Cluster) %>%
		slice_head(n = 10) %>%
		dplyr::pull(Gene.Name) %>%
		identity()
	
	ggplotify::as.ggplot(
		seu_complex_heatmap(unfiltered_seu,
												features = heatmap_features,
												group.by = c(group.by, "Phase", "scna"),
												col_arrangement = c(group.by, "Phase", "scna"),
												cluster_rows = FALSE, use_raster = TRUE
		)
	) +
		labs(title = sample_id)
	
	plot_path <- glue("results/effect_of_filtering/abbreviation/{sample_id}_abbreviation_heatmap.pdf")
	plot_list["abbreviation_heatmap"] <- plot_path
	ggsave(plot_path, height = 8, width = 8)
	
	# nCount_gene umaps ------------------------------
	
	unfiltered_seu$log_nCount_gene <- log1p(unfiltered_seu$nCount_gene)
	filtered_seu$log_nCount_gene <- log1p(filtered_seu$nCount_gene)
	
	(FeaturePlot(unfiltered_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
			labs(title = "unfiltered")) +
		(FeaturePlot(filtered_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
		 	labs(title = "filtered")) +
		plot_annotation(title = sample_id)
	
	fs::dir_create(glue("results/effect_of_filtering/nCount_gene"))
	plot_path <- glue("results/effect_of_filtering/nCount_gene/{sample_id}_nCount_gene_umaps.pdf")
	plot_list["nCount_gene_umaps"] <- plot_path
	ggsave(plot_path, height = 8, width = 10)
	
	# abbreviation umaps ------------------------------
	mycols <- scales::hue_pal()(length(unique(unfiltered_seu@meta.data[["abbreviation"]])))
	
	(DimPlot(unfiltered_seu, group.by = "abbreviation", cols = mycols) +
			labs(title = "unfiltered")) +
		(DimPlot(filtered_seu, group.by = "abbreviation", cols = mycols) +
		 	labs(title = "filtered")) +
		# (DimPlot(regressed_seu, group.by = "gene_snn_res.0.2") +
		#   labs(title = "regressed")) +
		plot_annotation(title = sample_id)
	
	fs::dir_create(glue("results/effect_of_filtering/abbreviation"))
	plot_path <- glue("results/effect_of_filtering/abbreviation/{sample_id}_abbreviation_umaps.pdf")
	plot_list["abbreviation_umaps"] <- plot_path
	ggsave(plot_path, height = 8, width = 10)
	
	# scna umaps ------------------------------
	mycols <- scales::hue_pal()(length(unique(unfiltered_seu@meta.data[["scna"]])))
	
	unfiltered_seu@meta.data$scna <- vec_split_label_line(unfiltered_seu@meta.data$scna, 3)
	filtered_seu@meta.data$scna <- vec_split_label_line(filtered_seu@meta.data$scna, 3)
	
	(DimPlot(unfiltered_seu, group.by = "scna", cols = mycols) +
			labs(title = "unfiltered")) +
		(DimPlot(filtered_seu, group.by = "scna", cols = mycols) +
		 	labs(title = "filtered")) +
		# (DimPlot(regressed_seu, group.by = "gene_snn_res.0.2") +
		#     labs(title = "regressed")) +
		plot_annotation(title = sample_id)
	
	fs::dir_create(glue("results/effect_of_filtering/scna"))
	plot_path <- glue("results/effect_of_filtering/scna/{sample_id}_scna_umaps.pdf")
	plot_list["scna_umaps"] <- plot_path
	ggsave(plot_path, height = 8, width = 10)
	
	return(plot_list)
}

#' Create a plot visualization
#'
#' @param unfiltered_seu_path File path
#' @param filtered_seu_path File path
#' @param group.by Character string (default: "gene_snn_res.0.2")
#' @param cluster_dictionary Cluster information
#' @param plot_path File path
#' @return ggplot2 plot object
#' @export
plot_effect_of_filtering <- function(unfiltered_seu_path, filtered_seu_path = NULL, group.by = "gene_snn_res.0.2", cluster_dictionary,
                                     mito_threshold = c(5, 10),
                                     nCount_threshold = c(500, 1000),
                                     nFeature_threshold = c(500, 1000),
                                     low_hypoxia_seu_path = NULL,
                                     high_hypoxia_seu_path = NULL,
                                     plot_path = NULL) {

  # Helper function: split label lines at commas with word wrap
  split_label_line <- function(label, n_comma_values = 3) {
    if (is.na(label)) return(NA_character_)
    label_vec <- label %>% stringr::str_split_1(",")
    label_groups <- ceiling(seq_along(label_vec) / n_comma_values)
    split_label <- split(label_vec, label_groups) %>%
      purrr::map(paste, collapse = ", ") %>%
      paste(collapse = "\n")
  }

  vec_split_label_line <- Vectorize(split_label_line)

  sample_id <- str_extract(unfiltered_seu_path, "SR[RX][0-9]+")

  if (is.null(plot_path)) {
    fs::dir_create("results/effect_of_filtering")
    plot_path <- glue::glue("results/effect_of_filtering/{sample_id}_effect_of_filtering.pdf")
  }

  unfiltered_seu <- readRDS(unfiltered_seu_path)
  original_filtered_seu <- if (!is.null(filtered_seu_path)) readRDS(filtered_seu_path) else NULL

  # Replace missing/empty scna with differential GT_opt labels
  # For each unique GT_opt value, show only the segments not in its nearest parent clone.
  # e.g. if clone A = "1q", clone B = "1q,2p", clone C = "1q,2p,13q":
  #   A -> "1q",  B -> "[1]+2p",  C -> "[2]+13q"
  differential_gt_labels <- function(gt_vals) {
    unique_gts <- unique(gt_vals[!is.na(gt_vals) & gt_vals != ""])
    if (length(unique_gts) == 0) {
      result <- gt_vals
      result[is.na(gt_vals) | gt_vals == ""] <- "diploid"
      return(result)
    }
    seg_sets <- setNames(
      lapply(unique_gts, function(gt) sort(trimws(strsplit(gt, ",")[[1]]))),
      unique_gts
    )
    sorted_gts <- names(seg_sets)[order(lengths(seg_sets))]
    # For each clone, find the closest proper-subset parent
    parent_of <- setNames(rep(NA_character_, length(sorted_gts)), sorted_gts)
    for (i in seq_along(sorted_gts)) {
      gt <- sorted_gts[i]; segs <- seg_sets[[gt]]
      best_len <- -1L
      for (j in seq_len(i - 1L)) {
        cand <- sorted_gts[j]; cand_segs <- seg_sets[[cand]]
        if (length(cand_segs) > best_len && all(cand_segs %in% segs)) {
          parent_of[gt] <- cand; best_len <- length(cand_segs)
        }
      }
    }
    clone_idx <- setNames(seq_along(sorted_gts), sorted_gts)
    display <- setNames(character(length(sorted_gts)), sorted_gts)
    for (gt in sorted_gts) {
      parent <- parent_of[gt]
      if (is.na(parent)) {
        display[gt] <- paste(seg_sets[[gt]], collapse = ",")
      } else {
        new_segs <- sort(setdiff(seg_sets[[gt]], seg_sets[[parent]]))
        display[gt] <- paste0("[", clone_idx[parent], "]+", paste(new_segs, collapse = ","))
      }
    }
    result <- gt_vals
    for (gt in sorted_gts) result[!is.na(result) & result == gt] <- display[gt]
    result[is.na(gt_vals) | (!is.na(gt_vals) & gt_vals == "")] <- "diploid"
    result
  }

  fill_scna_from_gt_opt <- function(seu) {
    if (!"GT_opt" %in% colnames(seu@meta.data)) return(seu)
    if (!"scna" %in% colnames(seu@meta.data)) {
      needs_gt <- rep(TRUE, nrow(seu@meta.data))
    } else {
      needs_gt <- is.na(seu@meta.data$scna) | seu@meta.data$scna == ""
    }
    if (!any(needs_gt)) return(seu)
    if (!"scna" %in% colnames(seu@meta.data)) seu@meta.data$scna <- NA_character_
    gt_diff <- differential_gt_labels(seu@meta.data$GT_opt[needs_gt])
    seu@meta.data$scna[needs_gt] <- gt_diff
    seu
  }
  unfiltered_seu <- fill_scna_from_gt_opt(unfiltered_seu)
  if (!is.null(original_filtered_seu)) original_filtered_seu <- fill_scna_from_gt_opt(original_filtered_seu)

  plot_list <- list()

  # use filter_reason / filter_keep columns if present, otherwise fall back to cluster_dictionary
  if (all(c("filter_reason", "filter_keep") %in% colnames(unfiltered_seu@meta.data))) {
    filtered_seu <- unfiltered_seu[, unfiltered_seu$filter_keep]
  } else {
    removed_clusters <- cluster_dictionary[[sample_id]] |>
      dplyr::filter(remove == 1) |>
      dplyr::pull(`gene_snn_res.0.2`)
    filtered_seu <- unfiltered_seu[, !unfiltered_seu$gene_snn_res.0.2 %in% removed_clusters]
  }

  # threshold sweep ------------------------------
  meta <- unfiltered_seu@meta.data
  sweep_grid <- expand.grid(
    mito_threshold = mito_threshold,
    nCount_threshold = nCount_threshold,
    nFeature_threshold = nFeature_threshold
  )
  sweep_grid$n_kept <- apply(sweep_grid, 1, function(r) {
    sum(
      meta$percent.mt < r[["mito_threshold"]] &
      meta$nCount_gene > r[["nCount_threshold"]] &
      meta$nFeature_gene > r[["nFeature_threshold"]]
    )
  })
  sweep_grid$pct_kept <- sweep_grid$n_kept / nrow(meta) * 100
  sweep_grid$mito_threshold <- factor(sweep_grid$mito_threshold)
  sweep_grid$nCount_threshold <- factor(sweep_grid$nCount_threshold)
  sweep_grid$nFeature_threshold <- factor(sweep_grid$nFeature_threshold)

  p_sweep <- ggplot(sweep_grid, aes(x = nCount_threshold, y = nFeature_threshold, fill = pct_kept)) +
    geom_tile() +
    geom_text(aes(label = n_kept), size = 2.5) +
    facet_wrap(~ mito_threshold, labeller = label_both) +
    scale_fill_viridis_c(name = "% kept", limits = c(0, 100)) +
    labs(title = glue::glue("{sample_id} — threshold sweep"), x = "nCount_threshold", y = "nFeature_threshold") +
    theme_bw()

  tmp_path <- tempfile(fileext = ".pdf")
  plot_list["threshold_sweep"] <- tmp_path
  ggsave(tmp_path, p_sweep, height = 4, width = 3 * length(mito_threshold))

  # abbreviation markers ------------------------------
  # Skip marker computation for very large objects — presto::wilcoxauc
  # (C++) crashes the process and cannot be caught by tryCatch.
  # Samples with >50k cells are too large to run safely.
  # Instead of using buggy seuratTools::plot_markers, build heatmap directly from presto results
  safe_plot_markers <- function(seu, grp = group.by, max_cells = 50000L,
                                filter_to_clusters = NULL,
                                force_recompute = FALSE, ...) {
    cat("safe_plot_markers: ncol =", ncol(seu), "\n")
    if (ncol(seu) > max_cells) {
      warning("Skipping marker plot — too many cells (", ncol(seu), ")")
      return(NULL)
    }

    # Ensure markers are calculated and in seu@misc$markers
    if (force_recompute ||
        is.null(seu@misc$markers) ||
        is.null(seu@misc$markers[[grp]])) {
      cat("Calculating markers for ", grp, "\n")
      seu <- tryCatch({
        if (exists("find_all_markers", mode = "function")) {
          find_all_markers(seu, grp)
        } else {
          seu
        }
      }, error = function(e) {
        warning("find_all_markers failed: ", conditionMessage(e)); seu
      })
    }

    cat("DEBUG safe_plot_markers: seu@misc$markers is NULL?", is.null(seu@misc$markers), "\n")
    if (!is.null(seu@misc$markers)) cat("DEBUG: seu@misc$markers[[grp]] is NULL?", is.null(seu@misc$markers[[grp]]), "\n")

    # Extract top markers from presto table; fall back to SCT_snn_res.0.2 if gene_snn_res.0.2 is absent
    if (is.null(seu@misc$markers[[grp]]) || is.null(seu@misc$markers[[grp]]$presto)) {
      fallback_grp <- "SCT_snn_res.0.2"
      if (!is.null(seu@misc$markers[[fallback_grp]]$presto)) {
        cat("DEBUG: falling back from", grp, "to", fallback_grp, "\n")
        grp <- fallback_grp
      } else {
        warning("No marker table found for ", grp, " and fallback ", fallback_grp, " also absent")
        return(NULL)
      }
    }

    top_markers_df <- seu@misc$markers[[grp]]$presto %>%
      dplyr::group_by(Cluster) %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(Cluster)

    if (!is.null(filter_to_clusters)) {
      top_markers_df <- dplyr::filter(
        top_markers_df, Cluster %in% filter_to_clusters
      )
    }

    # Drop markers for clusters not present in this Seurat object.
    # filtered_seu inherits misc$markers from the unfiltered object, so removed
    # clusters are still in the presto table. Without this filter, those genes
    # appear on the y-axis but fail the inner_join in marker_boxes (their cluster
    # ID isn't in df2), leaving present-cluster columns with no circles.
    present_clusters <- unique(as.character(seu@meta.data[[grp]]))
    top_markers_df <- dplyr::filter(
      top_markers_df, as.character(Cluster) %in% present_clusters
    )

    cat("DEBUG: top_markers_df nrow =", nrow(top_markers_df), "\n")

    if (nrow(top_markers_df) == 0) {
      warning("No markers found in presto table for ", grp)
      return(NULL)
    }

    # Keep cluster information and sort by cluster, then by gene name
    # This creates a diagonal pattern in the dotplot
    top_markers_df <- top_markers_df %>%
      dplyr::arrange(Cluster, Gene.Name)
    
    # Create feature list with cluster information
    features_to_plot <- unique(top_markers_df$Gene.Name)
    cat("DEBUG: features_to_plot length =", length(features_to_plot), "\n")

    # Filter to genes that exist in the Seurat object, preserving cluster-based order
    valid_features <- features_to_plot[features_to_plot %in% rownames(seu)]
    cat("DEBUG: valid_features length =", length(valid_features), "\n")
    if (length(valid_features) == 0) {
      warning("No valid marker genes found in Seurat object")
      return(NULL)
    }
    
    # Primary cluster per gene (first occurrence after arrange by Cluster)
    gene_primary_cluster <- top_markers_df %>%
      dplyr::distinct(Gene.Name, .keep_all = TRUE) %>%
      { setNames(as.character(.$Cluster), .$Gene.Name) }
    valid_features <- valid_features[
      order(gene_primary_cluster[valid_features], valid_features)
    ]
    cat("DEBUG: valid_features sorted by cluster: ",
        paste(head(valid_features, 10), collapse = ","), "...\n")

    # Create a simple heatmap using DotPlot or DoHeatmap
    tryCatch({
      dp <- DotPlot(seu, features = valid_features, group.by = grp, dot.scale = 8)
      # Manual calculation of percent-expressed and average expression per group
      assay_name <- if ("SCT" %in% names(seu@assays)) "SCT" else Seurat::DefaultAssay(seu)
      cat("DEBUG: assay_name =", assay_name, "\n")
      mat <- tryCatch(as.matrix(Seurat::GetAssayData(seu, assay = assay_name, layer = "data")[valid_features, , drop = FALSE]), error = function(e) {
        cat("ERROR fetching assay data:", conditionMessage(e), "\n"); NULL
      })
      if (is.null(mat)) {
        warning("Failed to fetch assay data for manual dotplot (assay=", assay_name, ")")
        return(NULL)
      }
      cat("DEBUG: mat dim =", dim(mat), "\n")

      groups <- as.character(seu@meta.data[[grp]])
      group_levels <- unique(groups[!is.na(groups)])
      # Naturally sort group levels: try numeric first, fall back to alphanumeric
      tryCatch({
        group_levels <- as.character(sort(as.numeric(group_levels)))
      }, error = function(e) {
        # If not all numeric, use alphanumeric (stringi for natural sort if available)
        if (requireNamespace("stringi", quietly = TRUE)) {
          group_levels <<- stringi::stri_sort(group_levels, numeric = TRUE)
        } else {
          group_levels <<- sort(group_levels)
        }
      })
      cat("DEBUG: group_levels (sorted) =", paste(group_levels, collapse=","), "\n")

      df_list <- list()
      for (g in group_levels) {
        cells <- which(groups == g)
        if (length(cells) == 0) next
        submat <- mat[valid_features, cells, drop = FALSE]
        for (i in seq_along(valid_features)) {
          gene <- valid_features[i]
          vals <- submat[gene, ]
          pct_expr <- mean(vals > 0, na.rm = TRUE) * 100
          avg_exp <- mean(vals, na.rm = TRUE)
          df_list[[length(df_list) + 1]] <- data.frame(id = g, features.plot = gene, pct_expr = pct_expr, avg.exp = avg_exp, stringsAsFactors = FALSE)
        }
      }
      cat("DEBUG: df_list length =", length(df_list), "\n")
      df2 <- do.call(rbind, df_list)
      cat("DEBUG: df2 nrow =", nrow(df2), "\n")
      # keep ordering consistent with original top markers
      df2$features.plot <- factor(df2$features.plot, levels = valid_features)
      df2$id <- factor(df2$id, levels = group_levels)

      # z-score avg expression per gene, matching Seurat's default DotPlot scale
      df2 <- df2 %>%
        dplyr::group_by(features.plot) %>%
        dplyr::mutate(avg.exp.scaled = as.numeric(scale(avg.exp))) %>%
        dplyr::ungroup()

      # one row per (cluster, gene) pair — circles every cluster where the
      # gene is a top marker, not just the primary one
      marker_boxes <- top_markers_df %>%
        dplyr::filter(Gene.Name %in% valid_features) %>%
        dplyr::transmute(
          id            = as.character(Cluster),
          features.plot = Gene.Name
        ) %>%
        dplyr::inner_join(
          dplyr::mutate(
            df2[, c("id", "features.plot", "pct_expr")],
            id            = as.character(id),
            features.plot = as.character(features.plot)
          ),
          by = c("id", "features.plot")
        ) %>%
        dplyr::mutate(
          id            = factor(id,
                                 levels = levels(df2$id)),
          features.plot = factor(features.plot,
                                 levels = levels(df2$features.plot))
        )

      p <- tryCatch({
        cat("DEBUG: building ggplot...\n")
        p_out <- ggplot2::ggplot(df2, ggplot2::aes(x = id, y = features.plot)) +
          ggplot2::geom_point(ggplot2::aes(size = pct_expr, color = avg.exp.scaled)) +
          ggplot2::geom_point(
            data = marker_boxes,
            ggplot2::aes(size = pct_expr),
            shape = 21, fill = NA, color = "black", stroke = 1.5
          ) +
          ggplot2::scale_size_continuous(range = c(1, 8), breaks = c(0, 25, 50, 75, 100), labels = function(x) paste0(x, "%")) +
          ggplot2::scale_color_gradient(low = "lightgrey", high = "blue", na.value = "grey90", name = "avg.exp.scaled") +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = ggplot2::rel(1.0)),
            plot.margin = ggplot2::unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
            panel.spacing.y = ggplot2::unit(0.5, "lines")
          ) +
          ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
          ggplot2::labs(title = paste("Markers for", grp), x = NULL, y = NULL)
        cat("DEBUG: ggplot created successfully\n")
        return(p_out)
      }, error = function(e) {
        cat("ERROR building ggplot:", conditionMessage(e), "\n")
        warning("manual DotPlot construction failed: ", conditionMessage(e)); NULL
      })

      return(p)
    }, error = function(e) {
      warning("DotPlot failed: ", conditionMessage(e))
      return(NULL)
    })
  }

  p_unfiltered <- safe_plot_markers(unfiltered_seu,
    marker_method = "wilcox", return_plotly = FALSE,
    hide_technical = "all", num_markers = 5)
  gc()
  p_filtered <- safe_plot_markers(filtered_seu,
    force_recompute = TRUE,
    marker_method = "wilcox", return_plotly = FALSE,
    hide_technical = "all", num_markers = 5)
  gc()

  # Load hypoxia seus once here — reused for both marker plots and dimplots below
  low_hypoxia_seu <- NULL
  if (!is.null(low_hypoxia_seu_path) && file.exists(low_hypoxia_seu_path)) {
    low_hypoxia_seu <- tryCatch(readRDS(low_hypoxia_seu_path), error = function(e) NULL)
  }
  high_hypoxia_seu <- NULL
  if (!is.null(high_hypoxia_seu_path) && file.exists(high_hypoxia_seu_path)) {
    high_hypoxia_seu <- tryCatch(readRDS(high_hypoxia_seu_path), error = function(e) NULL)
  }

  # include low hypoxia seus if available
  p_low_hypoxia <- NULL
  if (!is.null(low_hypoxia_seu)) {
    p_low_hypoxia <- safe_plot_markers(low_hypoxia_seu, marker_method = "wilcox", return_plotly = FALSE, hide_technical = "all", num_markers = 5)
    gc()
  }

  # include high hypoxia seus if available
  p_high_hypoxia <- NULL
  if (!is.null(high_hypoxia_seu)) {
    p_high_hypoxia <- safe_plot_markers(high_hypoxia_seu, marker_method = "wilcox", return_plotly = FALSE, hide_technical = "all", num_markers = 5)
    gc()
  }

  cat("DEBUG: p_unfiltered is NULL?", is.null(p_unfiltered), "\n")
  cat("DEBUG: p_filtered is NULL?", is.null(p_filtered), "\n")
  cat("DEBUG: p_low_hypoxia is NULL?", is.null(p_low_hypoxia), "\n")
  cat("DEBUG: p_high_hypoxia is NULL?", is.null(p_high_hypoxia), "\n")

  # Collect available plots and combine horizontally.
  # Annotate each marker panel title with its cell count (n).
  n_label <- function(label, seu) {
    if (is.null(seu)) return(label)
    glue::glue("{label} (n = {format(ncol(seu), big.mark = ',')})")
  }
  plot_items <- list()
  if (!is.null(p_unfiltered)) plot_items[["unfiltered"]] <- p_unfiltered + labs(title = n_label("unfiltered", unfiltered_seu))
  if (!is.null(p_filtered)) plot_items[["filtered"]] <- p_filtered + labs(title = n_label("filtered", filtered_seu))
  if (!is.null(p_low_hypoxia)) plot_items[["low_hypoxia"]] <- p_low_hypoxia + labs(title = n_label("low_hypoxia", low_hypoxia_seu))
  if (!is.null(p_high_hypoxia)) plot_items[["high_hypoxia"]] <- p_high_hypoxia + labs(title = n_label("high_hypoxia", high_hypoxia_seu))

  if (length(plot_items) > 0) {
    grid <- patchwork::wrap_plots(plot_items, ncol = length(plot_items)) +
      patchwork::plot_annotation(title = sample_id) +
      patchwork::plot_layout(guides = "collect")

    tmp_path <- tempfile(fileext = ".pdf")
    plot_list["abbreviation_markers"] <- tmp_path
    tryCatch({
      ggsave(tmp_path, grid, height = 10, width = 4 * length(plot_items), units = "in")
      cat("DEBUG: saved abbreviation_markers to", tmp_path, "\n")
    }, error = function(e) {
      warning("Failed to save abbreviation_markers: ", conditionMessage(e))
    })
  } else {
    cat("DEBUG: abbreviation_markers skipped (no plots generated)\n")
  }

  # filter_reason umap ------------------------------
  if ("filter_reason" %in% colnames(unfiltered_seu@meta.data)) {
    reason_levels <- c("clone_opt_na", "qc_fail", "cluster_remove", "malat1", "manual_exclude")
    unfiltered_seu@meta.data$filter_reason <- factor(unfiltered_seu@meta.data$filter_reason, levels = reason_levels)
    reason_cols <- c(scales::hue_pal()(length(reason_levels)), "grey80") |>
      set_names(c(reason_levels, "kept"))

    # sweep dimplots: re-evaluate qc_fail per threshold combo, carry forward other reasons
    sweep_plots <- lapply(seq_len(nrow(sweep_grid)), function(i) {
      r <- sweep_grid[i, ]
      qc_fail_sweep <-
        meta$percent.mt >= as.numeric(as.character(r$mito_threshold)) |
        meta$nCount_gene <= as.numeric(as.character(r$nCount_threshold)) |
        meta$nFeature_gene <= as.numeric(as.character(r$nFeature_threshold))
      sweep_reason <- dplyr::case_when(
        !is.na(unfiltered_seu@meta.data$filter_reason) & unfiltered_seu@meta.data$filter_reason != "qc_fail" ~ as.character(unfiltered_seu@meta.data$filter_reason),
        qc_fail_sweep ~ "qc_fail",
        TRUE ~ NA_character_
      )
      unfiltered_seu@meta.data$sweep_reason_label <- factor(
        tidyr::replace_na(sweep_reason, "kept"),
        levels = c(reason_levels, "kept")
      )
      DimPlot(unfiltered_seu, group.by = "sweep_reason_label") +
        scale_color_manual(values = reason_cols, limits = names(reason_cols), na.value = "grey90") +
        labs(
          title = glue::glue("mito<{r$mito_threshold}"),
          subtitle = glue::glue("nCount>{r$nCount_threshold} nFeature>{r$nFeature_threshold}"),
          colour = "filter reason"
        )
    })

    # all sweep plots on a single page: ncol = nCount thresholds, rows = mito x nFeature combos
    tmp_path <- tempfile(fileext = ".pdf")
    n_sweep_rows <- length(mito_threshold) * length(nFeature_threshold)
    sweep_ok <- tryCatch({
      pdf(tmp_path, height = 3 * n_sweep_rows, width = 4.5 * length(nCount_threshold))
      on.exit(if (!is.null(dev.list())) dev.off(), add = TRUE)
      grid <- patchwork::wrap_plots(sweep_plots, ncol = length(nCount_threshold)) +
        patchwork::plot_annotation(title = sample_id) +
        patchwork::guide_area() +
        patchwork::plot_layout(guides = "collect")
      print(grid)
      dev.off()
      on.exit(NULL)
      TRUE
    }, error = function(e) {
      if (!is.null(dev.list())) dev.off()
      warning("sweep_dimplots failed: ", conditionMessage(e))
      FALSE
    })
    if (sweep_ok) plot_list["sweep_dimplots"] <- tmp_path

    reason_labels <- tidyr::replace_na(as.character(unfiltered_seu@meta.data$filter_reason), "kept")
    unfiltered_seu@meta.data$filter_reason_label <- factor(reason_labels, levels = c(reason_levels, "kept"))

    tmp_path <- tempfile(fileext = ".pdf")
    fr_ok <- tryCatch({
      pdf(tmp_path, height = 3, width = 4.5)
      on.exit(if (!is.null(dev.list())) dev.off(), add = TRUE)
      print(DimPlot(unfiltered_seu, group.by = "filter_reason_label") +
        scale_color_manual(values = reason_cols, limits = names(reason_cols), na.value = "grey90") +
        labs(title = sample_id, colour = "filter reason"))
      dev.off()
      on.exit(NULL)
      TRUE
    }, error = function(e) {
      if (!is.null(dev.list())) dev.off()
      warning("filter_reason_umap failed: ", conditionMessage(e))
      FALSE
    })
    if (fr_ok) plot_list["filter_reason_umap"] <- tmp_path
  }

  # dimplot grids ------------------------------
  scna_cols <- scales::hue_pal()(length(unique(unfiltered_seu@meta.data[["scna"]])))

  unfiltered_seu@meta.data$scna <- vec_split_label_line(unfiltered_seu@meta.data$scna, 3)
  filtered_seu@meta.data$scna   <- vec_split_label_line(filtered_seu@meta.data$scna, 3)

  cluster_labels <- unique(as.character(unfiltered_seu@meta.data[[group.by]]))
  cluster_labels <- stringr::str_sort(cluster_labels, numeric = TRUE)
  group_cols  <- scales::hue_pal()(length(cluster_labels)) |> set_names(cluster_labels)

  n_cols <- if (!is.null(original_filtered_seu)) 3L else 2L

  # row 1: scna
  scna_plots <- list(
    DimPlot(unfiltered_seu, group.by = "scna", cols = scna_cols) + labs(title = "unfiltered"),
    DimPlot(filtered_seu,   group.by = "scna", cols = scna_cols) + labs(title = "filtered")
  )
  if (!is.null(original_filtered_seu)) {
    original_filtered_seu@meta.data$scna <- vec_split_label_line(original_filtered_seu@meta.data$scna, 3)
    scna_plots[[3]] <- DimPlot(original_filtered_seu, group.by = "scna", cols = scna_cols) + labs(title = "re-clustered")
  }

  # row 2: cluster
  cluster_plots <- list(
    DimPlot(unfiltered_seu, group.by = group.by, cols = group_cols) + labs(title = "unfiltered"),
    DimPlot(filtered_seu,   group.by = group.by, cols = group_cols) + labs(title = "filtered")
  )
  if (!is.null(original_filtered_seu)) {
    cluster_plots[[3]] <- DimPlot(original_filtered_seu, group.by = group.by, cols = group_cols) + labs(title = "re-clustered")
  }

  dimplot_grid <- patchwork::wrap_plots(c(scna_plots, cluster_plots), ncol = n_cols) +
    patchwork::plot_annotation(title = sample_id)

  tmp_path <- tempfile(fileext = ".pdf")
  plot_list["abbreviation_umaps"] <- tmp_path
  ggsave(tmp_path, dimplot_grid, height = 6, width = 4.5 * n_cols)
  rm(unfiltered_seu, filtered_seu); gc()

  # hypoxia filtering dimplots ------------------------------
  # Show which cells are dropped between filtered_seu and seus_low_hypoxia.
  # Uses original_filtered_seu (pipeline seu with UMAP) as the "before" object.
  # low_hypoxia_seu already loaded above — no second readRDS needed.
  if (!is.null(low_hypoxia_seu)) {
    {
      source_seu <- if (!is.null(original_filtered_seu)) original_filtered_seu else filtered_seu
      kept_barcodes <- colnames(low_hypoxia_seu)
      source_seu@meta.data$hypoxia_filter_label <- ifelse(
        colnames(source_seu) %in% kept_barcodes, "kept_low_hypoxia", "dropped_high_hypoxia"
      )
      hypoxia_filter_cols <- c(kept_low_hypoxia = "#2166ac", dropped_high_hypoxia = "#d73027")

      p_hypoxia_source <- DimPlot(source_seu, group.by = "hypoxia_filter_label") +
        scale_color_manual(values = hypoxia_filter_cols) +
        labs(title = glue::glue("{sample_id} filtered → hypoxia split"), colour = NULL)

      hypoxia_plot_list <- list(p_hypoxia_source)

      if ("hypoxia_score" %in% colnames(source_seu@meta.data)) {
        p_hypoxia_score <- FeaturePlot(source_seu, features = "hypoxia_score") +
          labs(title = "hypoxia score (filtered seu)")
        hypoxia_plot_list[[2]] <- p_hypoxia_score
      }

      # Gradient FeaturePlot of hypoxia_score (not a categorical scna DimPlot):
      # shows the residual hypoxia signal across the retained low-hypoxia cells.
      p_low_hypoxia_scna <- if ("hypoxia_score" %in% colnames(low_hypoxia_seu@meta.data)) {
        FeaturePlot(low_hypoxia_seu, features = "hypoxia_score", reduction = "umap") +
          labs(title = "low hypoxia seu (hypoxia_score)")
      } else {
        DimPlot(low_hypoxia_seu, group.by = "scna", reduction = "umap") +
          labs(title = "low hypoxia seu (scna)")
      }
      hypoxia_plot_list[[length(hypoxia_plot_list) + 1]] <- p_low_hypoxia_scna

      p_low_hypoxia_clusters <- DimPlot(
        low_hypoxia_seu, group.by = group.by, reduction = "umap"
      ) + labs(title = "low hypoxia seu (clusters)")
      hypoxia_plot_list[[length(hypoxia_plot_list) + 1]] <- p_low_hypoxia_clusters

      p_low_hypoxia_phase <- DimPlot(
        low_hypoxia_seu, group.by = "Phase", reduction = "umap"
      ) + labs(title = "low hypoxia seu (Phase)")
      hypoxia_plot_list[[length(hypoxia_plot_list) + 1]] <- p_low_hypoxia_phase

      hypoxia_grid <- patchwork::wrap_plots(hypoxia_plot_list, ncol = length(hypoxia_plot_list)) +
        patchwork::plot_annotation(title = glue::glue("{sample_id} — hypoxia filtering"))

      tmp_path <- tempfile(fileext = ".pdf")
      hyp_ok <- tryCatch({
        ggsave(tmp_path, hypoxia_grid, height = 4, width = 4.5 * length(hypoxia_plot_list))
        TRUE
      }, error = function(e) {
        warning("hypoxia_filtering_dimplots failed: ", conditionMessage(e))
        FALSE
      })
      if (hyp_ok) plot_list["hypoxia_filtering_dimplots"] <- tmp_path
    }
  }
  rm(low_hypoxia_seu); gc()

  # high hypoxia filtering dimplots ------------------------------
  # high_hypoxia_seu already loaded above — no second readRDS needed.
  if (!is.null(high_hypoxia_seu)) {
    {
      source_seu <- if (!is.null(original_filtered_seu)) original_filtered_seu else filtered_seu
      kept_barcodes <- colnames(high_hypoxia_seu)
      source_seu@meta.data$hypoxia_filter_label <- ifelse(
        colnames(source_seu) %in% kept_barcodes, "kept_high_hypoxia", "dropped_low_hypoxia"
      )
      hypoxia_filter_cols_hh <- c(kept_high_hypoxia = "#d73027", dropped_low_hypoxia = "#2166ac")

      p_hypoxia_source_hh <- DimPlot(source_seu, group.by = "hypoxia_filter_label") +
        scale_color_manual(values = hypoxia_filter_cols_hh) +
        labs(title = glue::glue("{sample_id} filtered → high hypoxia"), colour = NULL)

      hypoxia_plot_list_hh <- list(p_hypoxia_source_hh)

      if ("hypoxia_score" %in% colnames(source_seu@meta.data)) {
        p_hypoxia_score_hh <- FeaturePlot(source_seu, features = "hypoxia_score") +
          labs(title = "hypoxia score (filtered seu)")
        hypoxia_plot_list_hh[[2]] <- p_hypoxia_score_hh
      }

      # Gradient FeaturePlot of hypoxia_score (not a categorical scna DimPlot):
      # shows the residual hypoxia signal across the retained high-hypoxia cells.
      p_high_hypoxia_scna <- if ("hypoxia_score" %in% colnames(high_hypoxia_seu@meta.data)) {
        FeaturePlot(high_hypoxia_seu, features = "hypoxia_score", reduction = "umap") +
          labs(title = "high hypoxia seu (hypoxia_score)")
      } else {
        DimPlot(high_hypoxia_seu, group.by = "scna", reduction = "umap") +
          labs(title = "high hypoxia seu (scna)")
      }
      hypoxia_plot_list_hh[[length(hypoxia_plot_list_hh) + 1]] <- p_high_hypoxia_scna

      p_high_hypoxia_clusters <- DimPlot(
        high_hypoxia_seu, group.by = group.by, reduction = "umap"
      ) + labs(title = "high hypoxia seu (clusters)")
      hypoxia_plot_list_hh[[length(hypoxia_plot_list_hh) + 1]] <- p_high_hypoxia_clusters

      p_high_hypoxia_phase <- DimPlot(
        high_hypoxia_seu, group.by = "Phase", reduction = "umap"
      ) + labs(title = "high hypoxia seu (Phase)")
      hypoxia_plot_list_hh[[length(hypoxia_plot_list_hh) + 1]] <- p_high_hypoxia_phase

      hypoxia_grid_hh <- patchwork::wrap_plots(hypoxia_plot_list_hh, ncol = length(hypoxia_plot_list_hh)) +
        patchwork::plot_annotation(title = glue::glue("{sample_id} — high hypoxia filtering"))

      tmp_path <- tempfile(fileext = ".pdf")
      hyp_ok_hh <- tryCatch({
        ggsave(tmp_path, hypoxia_grid_hh, height = 4, width = 4.5 * length(hypoxia_plot_list_hh))
        TRUE
      }, error = function(e) {
        warning("high_hypoxia_filtering_dimplots failed: ", conditionMessage(e))
        FALSE
      })
      if (hyp_ok_hh) plot_list["high_hypoxia_filtering_dimplots"] <- tmp_path
    }
  }
  rm(high_hypoxia_seu, original_filtered_seu); gc()

  # marker heatmaps at gene_snn_res.0.2 for available inputs ------------------
  maybe_heatmap_file <- function(seu_path, label) {
    if (is.null(seu_path) || !file.exists(seu_path)) return(NULL)
    file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")
    out_path <- glue::glue("results/{file_slug}_{label}heatmap_phase_scatter_patchwork.pdf")
    tryCatch({
      plot_seu_marker_heatmap(seu_path, nb_paths = numbat_rds_files, clone_simplifications = large_clone_simplifications, group.by = "gene_snn_res.0.2", label = label, tmp_plot_path = FALSE)
      if (file.exists(out_path)) out_path else NULL
    }, error = function(e) {
      warning("marker heatmap generation failed for ", seu_path, ": ", conditionMessage(e))
      NULL
    })
  }

  hm_unfiltered <- maybe_heatmap_file(unfiltered_seu_path, "unfiltered_")
  if (!is.null(hm_unfiltered)) plot_list["marker_heatmap_unfiltered"] <- hm_unfiltered

  hm_filtered <- maybe_heatmap_file(filtered_seu_path, "filtered_")
  if (!is.null(hm_filtered)) plot_list["marker_heatmap_filtered"] <- hm_filtered

  hm_low_hypoxia <- maybe_heatmap_file(low_hypoxia_seu_path, "low_hypoxia_")
  if (!is.null(hm_low_hypoxia)) plot_list["marker_heatmap_low_hypoxia"] <- hm_low_hypoxia

  hm_high_hypoxia <- maybe_heatmap_file(high_hypoxia_seu_path, "high_hypoxia_")
  if (!is.null(hm_high_hypoxia)) plot_list["marker_heatmap_high_hypoxia"] <- hm_high_hypoxia

  # Filter to only valid, non-empty PDF files before combining
  valid_pdfs <- Filter(function(f) !is.null(f) && file.exists(f) && file.size(f) > 100L,
                       unlist(plot_list))
  cat("DEBUG: plot_list keys:", names(plot_list), "\n")
  cat("DEBUG: valid_pdfs to combine:", valid_pdfs, "\n")
  if (length(valid_pdfs) == 0) {
    warning("No valid PDF pages produced for ", sample_id)
    return(NULL)
  }
  plot_path <- tryCatch(
    qpdf::pdf_combine(valid_pdfs, plot_path),
    error = function(e) {
      warning("pdf_combine failed: ", conditionMessage(e))
      NULL
    }
  )

  return(plot_path)
}

#' Create a plot visualization
#'
#' @param filtered_seu_path File path
#' @param regressed_seu_path File path
#' @param resolution Parameter for resolution
#' @param filter_dropped_cluster Cluster information
#' @param regress_dropped_cluster Cluster information
#' @param n_features Parameter for n features
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
plot_effect_of_regression <- function(filtered_seu_path, regressed_seu_path, resolution = 0.4, filter_dropped_cluster = NULL, regress_dropped_cluster = NULL, n_features = 2, ...) {
  
  
  
  
  
  
  #

  sample_id <- str_extract(filtered_seu_path, "SR[RX][0-9]+")

  regressed_seu <- readRDS(regressed_seu_path)


  regressed_seu$scna[regressed_seu$scna == ""] <- ".diploid"
  regressed_seu$scna <- factor(regressed_seu$scna)

  regressed_meta <- regressed_seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, starts_with("SCT_snn_res")) %>%
    tibble::column_to_rownames("cell") %>%
    dplyr::rename_with(
      ~ str_replace(.x, "SCT_snn_res.", "regressed."),
      starts_with("SCT_")
    ) %>%
    identity()

  filtered_seu <- readRDS(filtered_seu_path)

  filtered_seu$scna[filtered_seu$scna == ""] <- ".diploid"
  filtered_seu$scna <- factor(filtered_seu$scna)

  filtered_seu <- AddMetaData(filtered_seu, regressed_meta)

  # Drop cells that were not in regressed_seu (AddMetaData gives NA for missing cells)
  cells_with_regressed <- !is.na(filtered_seu@meta.data[[colnames(regressed_meta)[1]]])
  filtered_seu <- filtered_seu[, cells_with_regressed]

  if (inherits(filtered_seu[["gene"]], "Assay5")) {
    filtered_seu <- JoinLayers(filtered_seu, assay = "gene")
  }

  valid_metavars <- colnames(regressed_meta)[
    sapply(colnames(regressed_meta), function(m) {
      dplyr::n_distinct(filtered_seu@meta.data[[m]], na.rm = TRUE) > 1
    })
  ]

  if (length(valid_metavars) == 0) {
    return(NA_character_)
  }

  filtered_seu <- tryCatch({
    find_all_markers(filtered_seu, valid_metavars)
  }, error = function(e) {
    warning("Skipping regression effect plot due to marker calculation failure: ", e$message)
    return(NULL)
  })

  if (is.null(filtered_seu)) {
    return(NA_character_)
  }

  regressed_cluster <- glue("regressed.{resolution}")
  
  # if(!is.null(regress_dropped_cluster)){
  # 	regressed_seu <- regressed_seu[,!regressed_seu@meta.data[[regressed_cluster]] %in% regress_dropped_cluster]
  # }

  filtered_cluster <- glue("SCT_snn_res.{resolution}")
  
  # browser()
  if(!is.null(filter_dropped_cluster)){
  	filtered_seu <- filtered_seu[,!filtered_seu@meta.data[[filtered_cluster]] %in% filter_dropped_cluster]
  }

  regressed_features <-
    filtered_seu@misc$markers[[regressed_cluster]]$presto %>%
    group_by(Cluster) %>%
    slice_head(n = n_features) %>%
    dplyr::select(regressed_cluster = Cluster, Gene.Name) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE) %>%
    dplyr::filter(Gene.Name %in% rownames(Seurat::GetAssayData(filtered_seu, assay = "SCT", layer = "scale.data"))) %>%
    identity()

  filtered_features <-
    filtered_seu@misc$markers[[filtered_cluster]]$presto %>%
    group_by(Cluster) %>%
    slice_head(n = n_features) %>%
    dplyr::select(Cluster, Gene.Name) %>%
    dplyr::mutate(filtered_cluster = Cluster) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE) %>%
    dplyr::filter(Gene.Name %in% rownames(Seurat::GetAssayData(filtered_seu, assay = "SCT", layer = "scale.data"))) %>%
    identity()

  # filtered ------------------------------

  heatmap_features <- filtered_features %>%
    dplyr::left_join(regressed_features, by = "Gene.Name") %>%
    dplyr::mutate(regressed_cluster = replace_na(regressed_cluster, "NA")) %>%
    identity()

  filtered_heatmap <- ggplotify::as.ggplot(
    seu_complex_heatmap(filtered_seu,
      features = heatmap_features$Gene.Name,
      group.by = c("G2M.Score", "S.Score", "scna", regressed_cluster),
      col_arrangement = c(filtered_cluster, "scna"),
      cluster_rows = FALSE,
      column_split = sort(filtered_seu@meta.data[[filtered_cluster]]),
      row_split = rev(heatmap_features$filtered_cluster),
      row_title_rot = 0,
      use_raster = TRUE,
      # right_annotation = row_ha
    )
  ) +
    labs(title = "filtered") +
    theme()

  cc_data <- FetchData(filtered_seu, c(filtered_cluster, regressed_cluster, "G2M.Score", "S.Score", "Phase", "scna"))

  centroid_data <-
    cc_data %>%
    dplyr::group_by(.data[[filtered_cluster]]) %>%
    dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
    # dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
    dplyr::mutate({{ filtered_cluster }} := as.factor(.data[[filtered_cluster]])) %>%
    dplyr::mutate(centroid = "centroids") %>%
    identity()

  filtered_facet_cell_cycle_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[filtered_cluster]], color = .data[["scna"]])) +
    geom_point(size = 0.1) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[filtered_cluster]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
    facet_wrap(~ .data[[filtered_cluster]], ncol = 3) +
    theme_light() +
    # geom_label(data = labels,
    # 					 aes(label = label),
    # 					 # x = Inf,
    # 					 # y = -Inf,
    # 					 x = max(cc_data$S.Score)+0.05,
    # 					 y = max(cc_data$G2M.Score)-0.1,
    # 					 hjust=1,
    # 					 vjust=1,
    # 					 inherit.aes = FALSE) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    labs(title = "filtered") +
    # guides(color = "none") +
    NULL

  # regressed ------------------------------

  heatmap_features <- regressed_features %>%
    dplyr::left_join(filtered_features, by = "Gene.Name") %>%
    dplyr::mutate(filtered_cluster = replace_na(filtered_cluster, "NA")) %>%
    identity()

  regressed_heatmap <- ggplotify::as.ggplot(
    seu_complex_heatmap(filtered_seu,
      features = heatmap_features$Gene.Name,
      group.by = c("G2M.Score", "S.Score", "scna", filtered_cluster),
      col_arrangement = c(regressed_cluster, "scna"),
      cluster_rows = FALSE,
      column_split = sort(filtered_seu@meta.data[[regressed_cluster]]),
      row_split = rev(heatmap_features$regressed_cluster),
      row_title_rot = 0,
      use_raster = TRUE
      # right_annotation = row_ha
    )
  ) +
    labs(title = "regressed") +
    theme()

  cc_data <- FetchData(filtered_seu, c(filtered_cluster, regressed_cluster, "G2M.Score", "S.Score", "Phase", "scna"))

  centroid_data <-
    cc_data %>%
    dplyr::group_by(.data[[regressed_cluster]]) %>%
    dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
    # dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
    dplyr::mutate({{ regressed_cluster }} := as.factor(.data[[regressed_cluster]])) %>%
    dplyr::mutate(centroid = "centroids") %>%
    identity()

  regressed_facet_cell_cycle_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[regressed_cluster]], color = .data[["scna"]])) +
    geom_point(size = 0.1) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[regressed_cluster]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
    facet_wrap(~ .data[[regressed_cluster]], ncol = 3) +
    theme_light() +
    # geom_label(data = labels,
    # 					 aes(label = label),
    # 					 # x = Inf,
    # 					 # y = -Inf,
    # 					 x = max(cc_data$S.Score)+0.05,
    # 					 y = max(cc_data$G2M.Score)-0.1,
    # 					 hjust=1,
    # 					 vjust=1,
    # 					 inherit.aes = FALSE) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    labs(title = "regressed") +
    # guides(color = "none") +
    NULL


  # patchworks ------------------------------
  
  plot_list = list(
  	"A" = filtered_heatmap,
  	"B" = regressed_heatmap,
  	"C" = filtered_facet_cell_cycle_plot,
  	"D" = regressed_facet_cell_cycle_plot,
  	"E" = plot_spacer(),
  	"F" = plot_spacer()
  )
  
  layout <- "
              AAAABBBB
              AAAABBBB
              AAAABBBB
              CCCEDDDF
              CCCEDDDF
              "
  
  mypatch <- wrap_plots(
  	plot_list
  ) +
  	plot_layout(design = layout) +
    plot_annotation(title = sample_id)

  plot_path <- glue("results/{sample_id}_regression_effects.pdf")
  ggsave(plot_path, plot = mypatch, ...)

  # markerplot <- plot_seu_marker_heatmap(filtered_seu_path, nb_path = numbat_rds_file, clone_simplifications = large_clone_simplifications, ...)
  return(plot_path)
}

#' Create a plot visualization
#'
#' @param filtered_seu_path File path
#' @param regressed_seu_path File path
#' @param resolution Character string (default: "0.4")
#' @param group.by Character string (default: "SCT_snn_res.0.6")
#' @return ggplot2 plot object
#' @export
plot_effect_of_regression_old <- function(filtered_seu_path, regressed_seu_path, resolution = "0.4", group.by = "SCT_snn_res.0.6") {
  
  
  
  
  
  
  #

  plot_list <- list()

  sample_id <- str_extract(filtered_seu_path, "SR[RX][0-9]+")

  filtered_seu <- filtered_seu_path

  regressed_seu <- regressed_seu_path

  fs::dir_create("results/effect_of_regression")

  # percent.mt umaps ------------------------------

  (FeaturePlot(filtered_seu, features = "percent.mt") +
    labs(title = "filtered")) +
    (FeaturePlot(regressed_seu, features = "percent.mt") +
      labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/percent_mt"))
  plot_path <- glue("results/effect_of_regression/percent_mt/{sample_id}_percent_mt_umaps.pdf")
  plot_list["percent_mt_umaps"] <- plot_path
  ggsave(plot_path, height = 8, width = 10)

  # distribution ------------------------------

  plot_distribution_of_clones_across_clusters(filtered_seu, sample_id, var_y = group.by)
  fs::dir_create(glue("results/effect_of_regression/distribution/filtered/"))
  plot_path <- glue("results/effect_of_regression/distribution/filtered/{sample_id}_filtered_distribution.pdf")
  plot_list["filtered_distribution"] <- plot_path
  ggsave(plot_path, height = 4, width = 8)

  filtered_dist_tables <- table_distribution_of_clones_across_clusters(filtered_seu, sample_id, clusters = group.by)

  table_path <- glue("results/effect_of_regression/distribution/filtered/{sample_id}_filtered_distribution.xlsx")
  plot_list["filtered_distribution_tables"] <- table_path
  writexl::write_xlsx(filtered_dist_tables, table_path)

  plot_distribution_of_clones_across_clusters(regressed_seu, sample_id, var_y = group.by)
  fs::dir_create(glue("results/effect_of_regression/distribution/regression"))
  plot_path <- glue("results/effect_of_regression/distribution/regression/{sample_id}_regressed_distribution.pdf")
  plot_list["regressed_distribution"] <- plot_path
  ggsave(plot_path, height = 4, width = 8)

  regressed_dist_tables <- table_distribution_of_clones_across_clusters(regressed_seu, sample_id, clusters = group.by)

  table_path <- glue("results/effect_of_regression/distribution/regression/{sample_id}_regressed_distribution.xlsx")
  plot_list["regressed_distribution_tables"] <- table_path
  writexl::write_xlsx(regressed_dist_tables, table_path)

  # abbreviation markers ------------------------------
  (plot_markers(filtered_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
    labs(title = "filtered")) +
    (plot_markers(regressed_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
      labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  plot_path <- glue("results/effect_of_regression/abbreviation/{sample_id}_abbreviation_markers.pdf")
  plot_list["abbreviation_markers"] <- plot_path
  ggsave(plot_path, height = 12, width = 15)

  filtered_marker_tables <- table_cluster_markers(filtered_seu)

  table_path <- glue("results/effect_of_regression/abbreviation/{sample_id}_filtered_markers.xlsx")
  plot_list["filtered_marker_tables"] <- table_path
  writexl::write_xlsx(filtered_marker_tables, table_path)

  regressed_marker_tables <- table_cluster_markers(regressed_seu) %>%
    purrr::compact()

  table_path <- glue("results/effect_of_regression/abbreviation/{sample_id}_regressed_markers.xlsx")
  plot_list["regressed_marker_tables"] <- table_path
  writexl::write_xlsx(regressed_marker_tables, table_path)

  heatmap_features <-
    table_cluster_markers(regressed_seu) %>%
    pluck(group.by) %>%
    group_by(Cluster) %>%
    slice_head(n = 10) %>%
    dplyr::pull(Gene.Name) %>%
    identity()

  ggplotify::as.ggplot(
    seu_complex_heatmap(regressed_seu,
      features = heatmap_features,
      group.by = c(group.by, "Phase", "scna"),
      col_arrangement = c(group.by, "Phase", "scna"),
      cluster_rows = FALSE
    )
  ) +
    labs(title = sample_id)

  plot_path <- glue("results/effect_of_regression/abbreviation/{sample_id}_abbreviation_heatmap.pdf")
  plot_list["abbreviation_heatmap"] <- plot_path
  ggsave(plot_path, height = 8, width = 8)

  # # nCount_gene umaps ------------------------------
  #
  # filtered_seu$log_nCount_gene <- log1p(filtered_seu$nCount_gene)
  # regressed_seu$log_nCount_gene <- log1p(regressed_seu$nCount_gene)
  #
  # (FeaturePlot(filtered_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
  #     labs(title = "filtered")) +
  #   (FeaturePlot(regressed_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
  #      labs(title = "regressed")) +
  #   plot_annotation(title = sample_id)
  #
  # fs::dir_create(glue("results/effect_of_regression/nCount_gene"))
  # plot_path = glue("results/effect_of_regression/nCount_gene/{sample_id}_nCount_gene_umaps.pdf")
  # plot_list["nCount_gene_umaps"] = plot_path
  # ggsave(plot_path, height = 8, width = 10)

  # # cluster umaps ------------------------------
  # mycols = scales::hue_pal()(length(unique(filtered_seu@meta.data[[group.by]])))
  #
  # (DimPlot(filtered_seu, group.by = group.by, cols = mycols) +
  #   labs(title = "filtered")) +
  # (DimPlot(regressed_seu, group.by = group.by, cols = mycols) +
  #   labs(title = "regressed")) +
  #   plot_annotation(title = sample_id)
  #
  # fs::dir_create(glue("results/effect_of_regression/cluster"))
  # plot_path = glue("results/effect_of_regression/cluster/{sample_id}_cluster_umaps.pdf")
  # plot_list["cluster_umaps"] = plot_path
  # ggsave(plot_path, height = 8, width = 10)
  #
  # abbreviation umaps ------------------------------
  mycols <- scales::hue_pal()(length(unique(filtered_seu@meta.data[[group.by]])) + 3)

  (DimPlot(filtered_seu, group.by = group.by, cols = mycols) +
    labs(title = "filtered")) +
    (DimPlot(regressed_seu, group.by = group.by, cols = mycols) +
      labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/abbreviation"))
  plot_path <- glue("results/effect_of_regression/abbreviation/{sample_id}_abbreviation_umaps.pdf")
  plot_list["abbreviation_umaps"] <- plot_path
  ggsave(plot_path, height = 8, width = 10)

  # scna umaps ------------------------------
  mycols <- scales::hue_pal()(length(unique(filtered_seu@meta.data[["scna"]])) + 3)

  # filtered_seu@meta.data$scna <- vec_split_label_line(filtered_seu@meta.data$scna, 3)
  # regressed_seu@meta.data$scna <- vec_split_label_line(regressed_seu@meta.data$scna, 3)

  (DimPlot(filtered_seu, group.by = "scna", cols = mycols) +
    labs(title = "filtered")) +
    (DimPlot(regressed_seu, group.by = "scna", cols = mycols) +
      labs(title = "regressed")) +
    # (DimPlot(regressed_filtered_seu, group.by = group.by) +
    #     labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/scna"))
  plot_path <- glue("results/effect_of_regression/scna/{sample_id}_scna_umaps.pdf")
  plot_list["scna_umaps"] <- plot_path
  ggsave(plot_path, height = 8, width = 10)

  # # scna markers ------------------------------
  # regressed_seu$scna[regressed_seu$scna == ""] <- "none"
  #
  # plot_markers(regressed_seu, "scna", marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10, unique_markers = TRUE) +
  #   labs(title = sample_id)
  #
  # plot_path = glue("results/effect_of_regression/{sample_id}_scna_markers.pdf")
  # plot_list["scna_markers"] = plot_path
  # ggsave(plot_path, height = 8, width = 6)

  # browseURL(glue("results/effect_of_regression/{sample_id}_markers_by_scna.pdf"))
  #
  # #   original_mt_plot <-
  # #     FeaturePlot(regressed_seu, features = "percent.mt") +
  # #     labs(title = "filtered")
  # #
  # #   regressed_mt_plot <- FeaturePlot(regressed_seu0, features = "percent.mt") +
  # #   labs(title = "regressed")
  # #
  # #   wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  # #     plot_annotation(title = sample_id)
  # #
  # # ggsave(glue("results/effect_of_regression/{sample_id}_percent_mt.pdf"), heigh = 6, width = 10)
  # # browseURL(glue("results/effect_of_regression/{sample_id}_percent_mt.pdf"))
  #
  # regressed_seu <- AddModuleScore(regressed_seu, subtype_markers)
  # regressed_seu0 <- AddModuleScore(regressed_seu0, subtype_markers)
  #
  #   original_mt_plot <-
  #     FeaturePlot(regressed_seu, features = "Cluster1") +
  #     labs(title = "filtered")
  #
  #   regressed_mt_plot <- FeaturePlot(regressed_seu0, features = "Cluster1") +
  #   labs(title = "regressed")
  #
  #   wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  #     plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_subtype1.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_subtype1.pdf"))
  #
  # original_mt_plot <-
  #   FeaturePlot(regressed_seu, features = "Cluster2") +
  #   labs(title = "filtered")
  #
  # regressed_mt_plot <- FeaturePlot(regressed_seu0, features = "Cluster2") +
  #   labs(title = "regressed")
  #
  # wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_subtype2.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_subtype2.pdf"))

  # DimPlot(regressed_seu, group.by = group.by) +
  #   labs(title = "filtered") +
  #   DimPlot(regressed_seu0, group.by = group.by) +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_louvain.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_louvain.pdf"))
  #
  # DimPlot(regressed_seu, group.by = "scna") +
  #   labs(title = "filtered") +
  #   DimPlot(regressed_seu0, group.by = "scna") +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_scna.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_scna.pdf"))
  #
  # DimPlot(regressed_seu, group.by = "Phase") +
  #   labs(title = "filtered") +
  #   DimPlot(regressed_seu0, group.by = "Phase") +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_phase.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_phase.pdf"))

  # phase umaps ------------------------------

  (DimPlot(filtered_seu, group.by = "Phase") +
    labs(title = "filtered")) +
    (DimPlot(regressed_seu, group.by = "Phase") +
      labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/phase"))
  plot_path <- glue("results/effect_of_regression/phase/{sample_id}_phase_umaps.pdf")
  plot_list["phase_umaps"] <- plot_path
  ggsave(plot_path, height = 8, width = 10)


  return(plot_list)
}
