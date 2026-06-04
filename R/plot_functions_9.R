# Plot Functions (101)

#' Create a heatmap visualization
#'
#' @param seu_path File path
#' @param child_seu_paths File path
#' @return ggplot2 plot object
#' @export
plot_merged_heatmap <- function(seu_path = "output/seurat/merged_1q_filtered_seu.rds", child_seu_paths = debranched_seus_1q) {
  seu <- check_merged_metadata(seu_path, child_seu_paths)

  # activity ------------------------------

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")

  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  usages <- dir_ls(glue("output/mosaicmpi/{sample_id}/"), glob = "*usage_k*.txt") %>%
    set_names(str_extract_all(., "(?<=k)[0-9]*")) %>%
    map(read_tsv) %>%
    map(tibble::column_to_rownames, "...1") %>%
    imap(~ set_names(.x, paste0(.y, "_", str_pad(colnames(.x), width = 2, pad = "0")))) %>%
    identity()

  usages <-
    usages[as.character(sort(as.numeric(names(usages))))]

  seu0 <- AddMetaData(seu, bind_cols(usages))

  seu0$sample <- seu0$sample_id

  myscna <- str_extract(seu_path, "(?<=merged_).*(?=_filtered)")

  scna_status <- glue("status of {myscna}")

  seu0[[scna_status]] <- str_detect(seu0$scna, myscna)

  phase_levels <- str_remove(seu0$clusters, "_[0-9]*$") |>
    str_replace("g1_stress", "hypoxia") |>
    str_replace("post_mitotic", "pm") |>
    str_replace("s_2", "s_star")

  phase_levels <- factor(phase_levels, levels = c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m", "hypoxia", "hsp", "other", "s_star"))

  seu0$phase_level <- phase_levels

  seu0 <-
    seu0 |>
    ScaleData() %>%
    find_all_markers(metavar = "phase_level") |>
    identity()

  # debug(seu_gene_heatmap)

  cluster_plot <- seu_gene_heatmap(seu0,
    marker_col = "phase_level",
    group.by = c("phase_level", scna_status, "S.Score", "G2M.Score", "sample"),
    col_arrangement = c("phase_level", scna_status),
    column_split = "phase_level",
    # hide_legends = colnames(usages[["6"]]),
    seurat_assay = "gene"
  )

  print(cluster_plot)

  ggsave(glue("results/{str_replace(fs::path_file(seu_path), '.rds', '_heatmap.pdf')}"), w = 8, h = 10)
}

#' Create a plot visualization
#'
#' @param seu Seurat object
#' @param var_y Character string (default: "clusters")
#' @return ggplot2 plot object
#' @export
make_cc_plot <- function(seu, var_y = "clusters") {
	
	seu$scna[seu$scna == ""] <- ".diploid"
	seu$scna <- factor(seu$scna)
	
	cc_data <- FetchData(seu, c(var_y, "G2M.Score", "S.Score", "Phase", "scna", "sample_id"))

	centroid_data <-
		cc_data %>%
		dplyr::group_by(.data[[var_y]]) %>%
		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score), n_cells = dplyr::n()) %>%
		# group_by(.data[[var_y]]) %>%
		mutate(percent = proportions(n_cells) * 100) %>%
		# dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
		dplyr::mutate(centroid = "centroids") %>%
		identity()
		
	split_centroid_data <-
		cc_data %>%
		dplyr::group_by(.data[[var_y]], scna) %>%
		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score), n_cells = dplyr::n()) %>%
		# group_by(.data[[var_y]]) %>%
		group_by(scna) %>%
		mutate(percent = proportions(n_cells) * 100) %>%
		# dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
		dplyr::mutate(centroid = "centroids") %>%
		identity()
	
	centroid_plot <-
		cc_data %>%
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[var_y]], color = .data[["scna"]])) +
		geom_point(size = 0.1) +
		theme_light() +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank()
		) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[var_y]], size = percent), alpha = 0.7, shape = 23, colour = "black") +
		geom_text_repel(data = centroid_data, aes(x = mean_x, y = mean_y, label = .data[[var_y]]), size = 4, color = "black") +
		scale_size_continuous(range = c(0,10), limits = c(1, 100), breaks = c(1,10,30,100)) + 
		guides("fill" = FALSE) +
		NULL
	
	facet_centroid_plot <- 
		centroid_plot + 
		facet_grid(.data[[var_y]]~scna)
	
	cc_plot <-
		cc_data %>%
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[var_y]], color = .data[["scna"]])) +
		geom_point(size = 0.1) +
		theme_light() +
		# theme(
		# 	strip.background = element_blank(),
		# 	strip.text.x = element_blank()
		# ) +
		guides("fill" = FALSE) +
		facet_grid(sample_id~.data[[var_y]]) + 
		geom_point(data = split_centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[var_y]], size = percent), alpha = 0.7, shape = 23, colour = "black") +
		geom_text_repel(data = centroid_data, aes(x = mean_x, y = mean_y, label = .data[[var_y]]), size = 4, color = "black") +
		scale_size_continuous(range = c(0,10), limits = c(1, 100), breaks = c(1,10,30,100)) + 
		NULL
	
	return(list("facet" = centroid_plot, "whole" = cc_plot))
}

#' Generate a volcano plot for differential expression data
#'
#' @param diffex_list Parameter for diffex list
#' @param seu_path File path
#' @param ... Additional arguments passed to other functions
#' @return Data frame
#' @export
plot_corresponding_clusters_diffex_volcanos <- function(diffex_list, seu_path, ...) {
  # asdf
  file_name <- fs::path_ext_remove(fs::path_file(seu_path))

  res <- diffex_list |>
    list_flatten() |>
    map(as_tibble) |>
    imap(~ dplyr::mutate(.x, diffex_comparison = .y)) |>
    map(dplyr::mutate, diffex_comparison = str_replace_all(diffex_comparison, "_", "-")) |>
    map(dplyr::mutate, diffex_comparison = str_replace_all(diffex_comparison, " v. ", "_")) |>
    map(dplyr::mutate, clone_comparison = diffex_comparison) |>
    # map(dplyr::select, -symbol) |>
    map(tibble::column_to_rownames, "symbol") |>
    identity()

  myplots <-
    imap(res, ~ make_volcano_plots(.x, sample_id = .y, mysubtitle = .y, color_by_chrom = FALSE)) %>%
    identity()

  pdf_path <- tempfile(tmpdir = "results", fileext = ".pdf")
  pdf(pdf_path)
  print(myplots)
  dev.off()

  return(pdf_path)
}

#' Create a heatmap visualization
#'
#' @param diffex_list Parameter for diffex list
#' @param seu_path File path
#' @param corresponding_states_dictionary Parameter for corresponding states dictionary
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @param numbat_rds_files File path
#' @param location Character string (default: "all")
#' @return ggplot2 plot object
#' @export
plot_corresponding_clusters_diffex_heatmaps <- function(diffex_list, seu_path, corresponding_states_dictionary, large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all") {

  file_name <- fs::path_ext_remove(fs::path_file(seu_path))

  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  numbat_rds_files <- numbat_rds_files %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  mynb <- readRDS(numbat_rds_files[[tumor_id]])

  seu <- readRDS(seu_path)

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  possible_cluster_comparison <- possibly(make_cluster_comparison)

  corresponding_states <-
    corresponding_states_dictionary |>
    dplyr::mutate(sample = str_remove(file_name, "_filtered_seu.*$")) |>
    dplyr::filter(file_name == fs::path_file(seu_path)) |>
    dplyr::select(w_scna, wo_scna, scna_of_interest) |>
    dplyr::mutate(comparison = glue("{w_scna} v. {wo_scna}")) |>
    identity()

  scna_of_interest <- unique(corresponding_states$scna_of_interest)

  large_clone_comparisons <- large_clone_comparisons[[sample_id]]

  large_clone_comparisons <- large_clone_comparisons[str_detect(names(large_clone_comparisons), scna_of_interest)]

  retained_clones <- names(large_clone_comparisons) %>%
    str_extract("[0-9]_v_[0-9]") %>%
    str_split("_v_", simplify = TRUE)

  seu <- seu[, seu$clone_opt %in% retained_clones]

  test0 <-
    diffex_list |>
    names() |>
    set_names() |>
    tibble::enframe("vs", "comparison") |>
    tidyr::separate(comparison, c("w_scna", "wo_scna"), " v. ") |>
    dplyr::mutate(
      w_scna = str_split(w_scna, "-"),
      wo_scna = str_split(wo_scna, "-")
    ) |>
    # str_split(" v. ") |>
    # flatten() |>
    identity()

  make_little_heatmap <- function(vs, w_scna, wo_scna, ...) {
    # 
    seu$scna_status <-
      seu@meta.data |>
      tibble::rownames_to_column("cell") |>
      dplyr::mutate(scna_status = dplyr::case_when(
        seu$clusters %in% w_scna ~ "w_scna",
        seu$clusters %in% wo_scna ~ "wo_scna"
      )) |>
      dplyr::mutate(scna_status = factor(scna_status, levels = c("wo_scna", "w_scna"))) |>
      dplyr::pull(scna_status)

    seu <-
      seu[, !is.na(seu$scna_status)]


    features <-
      diffex_list[[vs]][["all"]] |>
      dplyr::filter(p_val_adj < 0.05) |>
      dplyr::filter(symbol %in% VariableFeatures(seu)) %>%
      dplyr::mutate(fc_sign = ifelse(avg_log2FC < 0, "-", "+")) |>
      dplyr::group_by(fc_sign) |>
      dplyr::arrange(desc(abs(avg_log2FC))) |>
      dplyr::slice_head(n = 25) |>
      dplyr::pull(symbol) |>
      identity()

    test1 <- ggplotify::as.ggplot(seu_complex_heatmap2(seu, features = features, group.by = c("scna", "clusters"), col_arrangement = "scna")) + 
    labs(title = vs)

    # test1 <- ggplotify::as.ggplot(seu_complex_heatmap2(seu, features = features, group.by = c("scna", "scna_status", "clusters"))) +
    #   labs(title = vs)

    return(test1)
  }

  little_heatmaps <- pmap(test0, make_little_heatmap, diffex_list, seu)

  pdf_path <- tempfile(tmpdir = "results", fileext = ".pdf")
  pdf(pdf_path, w = 8, h = 10.5)
  print(little_heatmaps)
  dev.off()

  return(pdf_path)
}
