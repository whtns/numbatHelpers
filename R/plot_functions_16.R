# Plot Functions (108)

#' Perform make corresponding states dictionary operation
#'
#' @param my_tbl Parameter for my tbl
#' @return Function result
#' @export
make_corresponding_states_dictionary <- function(my_tbl){
	my_tbl <- 
		my_tbl |> 
		dplyr::mutate(file_name = factor(file_name, levels = unique(file_name)))
	
	split(my_tbl, my_tbl$file_name)
}

#' Perform make annohighlight from consensus operation
#'
#' @param ideogramPlot Parameter for ideogramPlot
#' @param ploty Parameter for ploty
#' @param chrom Parameter for chrom
#' @param chromstart Parameter for chromstart
#' @param chromend Parameter for chromend
#' @param fill Parameter for fill
#' @param width Parameter for width
#' @param pg_width Parameter for pg width
#' @return Function result
#' @export
make_annoHighlight_from_consensus <- function(ideogramPlot, ploty, chrom, chromstart, chromend, fill, width, pg_width){
	
	region <- pgParams(chrom = chrom, chromstart = chromstart, chromend = chromend)
	annoHighlight(
		plot = ideogramPlot, params = region,
		fill = fill,
		linecolor = "black",
		y = ploty, height = 0.2, just = c("left", "center"), default.units = "inches"
	)
	
	label_x = unit(0.15, "in") + ((chromstart+(chromend-chromstart)/2)/ideogramPlot$chromend*ideogramPlot$width)
	
	plotText(
		plot = ideogramPlot, params = region,
		label = glue("{round(width/1e6)} Mb"),
		x = label_x,
		y = ploty-0.2,
		fontsize = 8,
		just = c("center", "top"),
		fontcolor = fill
	)
	
}

#' Perform make rb scna ideograms operation
#'
#' @param nb_path File path
#' @param midline_threshold Threshold value for filtering
#' @param suffix Optional suffix inserted before the karyogram filename stem
#' @return ggplot2 plot object
#' @export
#' Generate ideograms for multiple suffixes from a single RDS load
#'
#' @param nb_path Numbat RDS file path
#' @param suffixes Named character vector mapping result names to filename suffixes
#' @param midline_threshold Threshold for midline filtering
#' @param filter_midline Whether to apply midline segment filtering
#' @return Named list: one list(plot=path, table=seg_table) per suffix
#' @export
make_rb_scna_ideograms_multi <- function(nb_path,
                                         suffixes = c(unfiltered = "", filtered = "_subset", low_hypoxia = "_low_hypoxia"),
                                         midline_threshold = 0.4,
                                         filter_midline = FALSE) {
	chrom_lengths <- seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)
	maxChromSize <- max(chrom_lengths)
	hg38_assembly <- assembly("hg38")
	tumor_id <- str_extract(nb_path, "SR[RX][0-9]+")
	dir.create("results/karyograms", showWarnings = FALSE, recursive = TRUE)

	mynb <- readRDS(nb_path)

	if (filter_midline) {
		retained_segs <- mynb$joint_post |>
			dplyr::mutate(at_midline = dplyr::case_when(
				dplyr::between(p_cnv, 0.3, 0.7) ~ 1,
				.default = 0
			)) |>
			group_by(seg) |>
			dplyr::summarise(percent_at_midline = sum(at_midline)/dplyr::n()) |>
			dplyr::filter(percent_at_midline <= midline_threshold) |>
			dplyr::arrange(desc(percent_at_midline)) |>
			dplyr::pull(seg)
	} else {
		retained_segs <- unique(mynb$joint_post$seg)
	}

	segmentation_table <- mynb$joint_post |>
		dplyr::filter(seg %in% retained_segs) |>
		dplyr::filter(p_cnv > 0.9) |>
		dplyr::distinct(CHROM, seg, cnv_state, .keep_all = TRUE) |>
		dplyr::mutate(CHROM = paste0("chr", CHROM)) |>
		dplyr::filter(!is.na(seg)) |>
		dplyr::mutate(fill = dplyr::case_when(cnv_state == "amp" ~ "red",
		                                      cnv_state == "bamp" ~ "pink",
		                                      cnv_state == "del" ~ "blue",
		                                      cnv_state == "loh" ~ "green")) |>
		dplyr::rename(chrom = CHROM, chromstart = seg_start, chromend = seg_end) |>
		dplyr::mutate(width = chromend - chromstart) |>
		dplyr::select(chrom, chromstart, chromend, fill, width, cnv_state) |>
		dplyr::mutate(chrom = factor(chrom, levels = paste0("chr", seq(1, 22)))) |>
		dplyr::arrange(chrom, chromstart)

	pg_width <- 4.25
	pg_height <- 8.85
	primary_suffix <- suffixes[[1]]
	primary_path <- glue("results/karyograms/{tumor_id}{primary_suffix}_karyogram.pdf")

	pdf(primary_path, width = pg_width, height = pg_height)
	pageCreate(width = pg_width, height = pg_height, default.units = "inches",
	           showGuides = FALSE, xgrid = 0, ygrid = 0)
	plotText(label = tumor_id, x = 1, y = 0.25, fontsize = 18)
	plotLegend(legend = c("gain", "balanced gain", "loss", "cnloh"),
	           fill = c("red", "pink", "blue", "green"),
	           border = FALSE, x = 0.1, y = 0.3,
	           width = pg_width - 1, height = 0.7, fontsize = 12,
	           orientation = "h", just = c("left", "top"), default.units = "inches")

	yCoord <- 1
	for (chr in paste0("chr", seq(1, 22))) {
		width <- (4 * chrom_lengths[[chr]]) / maxChromSize
		ideogramPlot <- plotIdeogram(chrom = chr, assembly = hg38_assembly,
		                             orientation = "h", x = 0.15, y = yCoord,
		                             height = 0.2, width = width, just = "left")
		plotText(label = gsub("chr", "", chr), x = 0.05, y = yCoord, fontsize = 10, rot = 90)
		if (chr %in% segmentation_table$chrom) {
			segmentation_table |>
				dplyr::filter(chrom %in% chr) |>
				pmap(~make_annoHighlight_from_consensus(ideogramPlot, yCoord, ..1, ..2, ..3, ..4, ..5, pg_width))
		}
		yCoord <- yCoord + 0.35
	}
	dev.off()

	# Copy the generated PDF to each additional suffix path (content is identical)
	result <- purrr::imap(suffixes, function(suffix, name) {
		dest_path <- glue("results/karyograms/{tumor_id}{suffix}_karyogram.pdf")
		if (dest_path != primary_path) file.copy(primary_path, dest_path, overwrite = TRUE)
		list("plot" = dest_path, "table" = segmentation_table)
	})
	result
}

make_rb_scna_ideograms <- function(nb_path, midline_threshold = 0.4, suffix = "", filter_midline = TRUE) {

	chrom_lengths = seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)
	maxChromSize <- max(chrom_lengths)

	tumor_id <- str_extract(nb_path, "SR[RX][0-9]+")

	dir.create("results/karyograms", showWarnings = FALSE, recursive = TRUE)
	plot_path <- glue("results/karyograms/{tumor_id}{suffix}_karyogram.pdf")

	mynb <- readRDS(nb_path)

	if (filter_midline) {
		retained_segs <- mynb$joint_post |>
			dplyr::mutate(at_midline = dplyr::case_when(
				dplyr::between(p_cnv, 0.3, 0.7) ~ 1,
				.default = 0
			)) |>
			group_by(seg) |>
			dplyr::summarise(percent_at_midline = sum(at_midline)/dplyr::n()) |>
			dplyr::filter(percent_at_midline <= midline_threshold) |>
			dplyr::arrange(desc(percent_at_midline)) |>
			dplyr::pull(seg) |>
			identity()
	} else {
		retained_segs <- unique(mynb$joint_post$seg)
	}

	segmentation_table <- mynb$joint_post |>
		dplyr::filter(seg %in% retained_segs) |>
		dplyr::filter(p_cnv > 0.9) |> 
		dplyr::distinct(CHROM, seg, cnv_state, .keep_all = TRUE) |> 
		dplyr::mutate(CHROM= paste0("chr", CHROM)) |> 
		dplyr::filter(!is.na(seg)) |> 
		# dplyr::rowwise() |> 
		# dplyr::mutate(p_cnv = max(across(any_of(c("p_amp", "p_bamp", "p_del", "p_bel", "p_loh"))))) |> 
		# dplyr::filter(p_cnv > 0.8) |> 
		dplyr::mutate(fill = dplyr::case_when(cnv_state == "amp" ~ "red",
																					cnv_state == "bamp" ~ "pink",
																					cnv_state == "del" ~ "blue",
																					cnv_state == "loh" ~ "green")) |> 
		dplyr::rename(chrom = CHROM, chromstart = seg_start, chromend = seg_end) |> 
		dplyr::mutate(width = chromend - chromstart) |> 
		dplyr::select(chrom, chromstart, chromend, fill, width, cnv_state) |> 
		dplyr::mutate(chrom = factor(chrom, levels = paste0("chr", seq(1,22)))) |> 
		dplyr::arrange(chrom, chromstart)
	
	pg_width = 4.25
	pg_height = 8.85
	hg38_assembly <- assembly("hg38")
	pdf(plot_path, width = pg_width, height = pg_height)

	pageCreate(
		width = pg_width, height = pg_height, default.units = "inches",
		showGuides = FALSE, xgrid = 0, ygrid = 0
	)

	plotText(
		label = tumor_id,
		x = 1, y = 0.25,
		fontsize = 18
	)

	legendPlot <- plotLegend(
		legend = c("gain", "balanced gain", "loss", "cnloh"),
		fill = c("red", "pink", "blue", "green"),
		border = FALSE,
		x = 0.1, y = 0.3, width = pg_width-1, height = 0.7,
		fontsize = 12,
		orientation = "h",
		just = c("left", "top"),
		default.units = "inches"
	)

	yCoord <- 1
	for (chr in c(paste0("chr", seq(1,22)))) {
		width <- (4 * chrom_lengths[[chr]]) / maxChromSize
		ideogramPlot <- plotIdeogram(
			chrom = chr, assembly = hg38_assembly,
			orientation = "h",
			x = 0.15, y = yCoord,
			height = 0.2, width = width,
			just = "left"
		)
		plotText(
			label = gsub("chr", "", chr),
			x = 0.05, y = yCoord, fontsize = 10, rot = 90
		)

		if (chr %in% segmentation_table$chrom) {
			segmentation_table |>
				dplyr::filter(chrom %in% chr) |>
				pmap(~make_annoHighlight_from_consensus(ideogramPlot, yCoord, ..1, ..2, ..3, ..4, ..5, pg_width))
		}

		yCoord <- yCoord + 0.35
	}
	
	axis_length <- (4 * chrom_lengths[["chr1"]]) / maxChromSize
	
	dev.off()
	
	return(list("plot" = plot_path, "table" = segmentation_table))
	
}

#' Perform make rb scna ideograms old operation
#'
#' @param consensus_file File path
#' @param suffix Optional suffix inserted before the karyogram filename stem
#' @return ggplot2 plot object
#' @export
make_rb_scna_ideograms_old <- function(consensus_file, suffix = "") {
	
	tumor_id <- str_extract(consensus_file, "SR[RX][0-9]+")
	
	plot_path <- glue("results/{tumor_id}{suffix}_karyogram.pdf")
	
	test0 <- read_tsv(consensus_file) |> 
		dplyr::mutate(CHROM= paste0("chr", CHROM)) |> 
		dplyr::filter(!is.na(seg)) |> 
		dplyr::mutate(fill = dplyr::case_when(cnv_state == "amp" ~ "red",
																					cnv_state == "del" ~ "blue")) |> 
		dplyr::rename(chrom = CHROM, chromstart = seg_start, chromend = seg_end) |> 
		dplyr::select(chrom, chromstart, chromend, fill)
	
	
	pdf(plot_path)
	
	pageCreate(
		width = 6.25, height = 5.25, default.units = "inches",
		showGuides = FALSE, xgrid = 0, ygrid = 0
	)
	plotText(
		label = tumor_id, 
		x = 0.5, y = 0.5,
		fontsize = 24
	)
	
	mychrom = "chr1"
	ploty = 1.5
	# chr1------------------------------
	ideogramPlot <- plotIdeogram(
		chrom = mychrom, assembly = "hg38",
		orientation = "h",
		x = 0.25, y = ploty, width = 5.75, height = 0.3, just = "left"
	)
	plotText(
		label = gsub("chr", "", mychrom),
		x = 0.25, y = ploty - 0.25, fontsize = 10
	)
	
	if(mychrom %in% test0$chrom){
		print("yes")
		test0 |> 
			dplyr::filter(chrom %in% mychrom) |> 
			pmap(~make_annoHighlight_from_consensus(ideogramPlot, ploty, ..1, ..2, ..3, ..4))
	}
	
	
	# chr2 ------------------------------
	mychrom = "chr2"
	ploty = 2.5
	ideogramPlot <- plotIdeogram(
		chrom = mychrom, assembly = "hg38",
		orientation = "h",
		x = 0.25, y = ploty, width = 5.75, height = 0.3, just = "left"
	)
	plotText(
		label = gsub("chr", "", mychrom),
		x = 0.25, y = ploty - 0.25, fontsize = 10
	)
	
	if(mychrom %in% test0$chrom){
		print("yes")
		test0 |> 
			dplyr::filter(chrom %in% mychrom) |> 
			pmap(~make_annoHighlight_from_consensus(ideogramPlot, ploty, ..1, ..2, ..3, ..4))
	}
	
	# chr6 ------------------------------
	mychrom = "chr6"
	ploty = 3.5
	ideogramPlot <- plotIdeogram(
		chrom = mychrom, assembly = "hg38",
		orientation = "h",
		x = 0.25, y = ploty, width = 5.75, height = 0.3, just = "left"
	)
	plotText(
		label = gsub("chr", "", mychrom),
		x = 0.25, y = ploty - 0.25, fontsize = 10
	)
	
	if(mychrom %in% test0$chrom){
		print("yes")
		test0 |> 
			dplyr::filter(chrom %in% mychrom) |> 
			pmap(~make_annoHighlight_from_consensus(ideogramPlot, ploty, ..1, ..2, ..3, ..4))
	}
	
	# chr16 ------------------------------
	mychrom = "chr16"
	ploty = 4.5
	ideogramPlot <- plotIdeogram(
		chrom = mychrom, assembly = "hg38",
		orientation = "h",
		x = 0.25, y = ploty, width = 5.75, height = 0.3, just = "left"
	)
	plotText(
		label = gsub("chr", "", mychrom),
		x = 0.25, y = ploty - 0.25, fontsize = 10
	)
	
	if(mychrom %in% test0$chrom){
		print("yes")
		test0 |> 
			dplyr::filter(chrom %in% mychrom) |> 
			pmap(~make_annoHighlight_from_consensus(ideogramPlot, ploty, ..1, ..2, ..3, ..4))
	}
	
	dev.off()
	
	return(plot_path)
}

make_volcano_plots <- function(myres, mysubtitle, sample_id, color_by_chrom = TRUE, force_genes = c("MYCN")) {
  diffex_comparison <- str_split(unique(myres[["diffex_comparison"]]), "_", simplify = TRUE)

  right_label <- diffex_comparison[[1]] %||% "right"
  left_label <- diffex_comparison[[2]] %||% "left"

  myres <-
    myres %>%
    dplyr::mutate(chr = case_when(
      chr == "X" ~ "23",
      chr == "Y" ~ "24",
      TRUE ~ as.character(chr)
    )) %>%
    dplyr::mutate(chr = str_pad(chr, side = "left", pad = "0", width = 2)) %>%
    dplyr::mutate(clone_comparison = str_replace_all(clone_comparison, "_", " ")) %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::mutate(rownames = symbol) %>%
    tibble::column_to_rownames("rownames")

  mytitle <- sample_id
  mysubtitle <- mysubtitle

  ref_var <-
    myres$chr %>%
    set_names(.)

  chrs <- str_pad(as.character(1:24), side = "left", pad = "0", width = 2)

  mypal <- scales::hue_pal()(length(chrs))
  names(mypal) <- chrs

  if(!color_by_chrom){
  	custom_cols <- NULL
  } else {
  	custom_cols <- mypal[ref_var]
  }

  FCcutoff <- summary(abs(myres$avg_log2FC))[[5]]

  selected_genes <-
    myres %>%
    dplyr::filter(abs(avg_log2FC) > 0.05, p_val_adj < 0.1) %>%
  	dplyr::slice_head(n = 25) |> 
    dplyr::pull(symbol) |>
  	identity()
  
  selected_genes <- c(force_genes, selected_genes)
  
  myplot <- EnhancedVolcano(myres,
    lab = rownames(myres),
    selectLab = selected_genes,
    labSize = 4,
    pointSize = c(ifelse(myres$symbol %in% force_genes, 4, 1)),
    x = "avg_log2FC",
    y = "p_val_adj",
    FCcutoff = FCcutoff,
    pCutoff = 5e-2,
    colCustom = custom_cols,
    max.overlaps = Inf,
    drawConnectors = TRUE,
    # maxoverlapsConnectors = 15,
    min.segment.length = 0.1,
    colConnectors = "grey",
    raster = TRUE
  ) +
    aes(color = chr) +
    # facet_wrap(~chr) +
    labs(title = mytitle, subtitle = mysubtitle)

  layer_scales(myplot)$y$range$range

  # plot_ymax = max(-log(myplot$data$p_val_adj, base = 10))
  plot_ymax <- max(ggplot_build(myplot)$layout$panel_params[[1]]$y.range)

  myplot <-
    myplot +
    annotation_custom(
      text_grob(
        left_label,
        size = 13,
        color = "red",
        face = "bold"
      ),
      xmin = -Inf,
      xmax = -Inf,
      ymin = plot_ymax,
      ymax = plot_ymax
    ) +
    annotation_custom(
      text_grob(
        right_label,
        size = 13,
        color = "red",
        face = "bold"
      ),
      xmin = Inf,
      xmax = Inf,
      ymin = plot_ymax,
      ymax = plot_ymax
    ) +
    theme(plot.margin = unit(c(1, 3, 1, 1), "lines")) +
    coord_cartesian(clip = "off") +
    NULL


  return(myplot)
}

#' Title
#' @export
make_numbat_heatmaps_old <- function(numbat_rds_file, filter_expressions = NULL, cluster_dictionary, p_min = 0.9, line_width = 0.1, extension = ""){
  #

  sample_id <- str_extract(numbat_rds_file, "SR[RX][0-9]+")

  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}/"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  myannot <- mynb$clone_post[, c("cell")]

  if (!is.null(filter_expressions)) {
    # filter cells
    phylo_heatmap_data <- mynb$clone_post %>%
      dplyr::select(cell, clone_opt) %>%
      dplyr::left_join(mynb$joint_post, by = "cell") %>%
      dplyr::left_join(myannot, by = "cell")

    excluded_cells <- map(filter_expressions[[sample_id]], pull_cells_matching_expression, phylo_heatmap_data) %>%
      unlist()
  }

	keep_cells <- colnames(seu) %in% myannot$cell
	if (exists("excluded_cells")) {
		keep_cells <- keep_cells & !(colnames(seu) %in% excluded_cells)
	}

  clusters_to_remove <-
    cluster_dictionary[[sample_id]] %>%
    dplyr::filter(remove == "1") %>%
    dplyr::pull(`gene_snn_res.0.2`)

	keep_cells <- keep_cells & !(seu$gene_snn_res.0.2 %in% clusters_to_remove)
	seu <- seu[, keep_cells]

  myannot <-
    seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, clone_opt, nCount_gene) %>%
    identity()

  ## numbat ------------------------------
  numbat_heatmap <- safe_plot_numbat(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = p_min, line_width = line_width)[["result"]]
  # numbat_heatmap <- plot_numbat(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = p_min, line_width = line_width)[["result"]]

  scna_variability_plot <- plot_variability_at_SCNA(numbat_heatmap[[3]][["data"]])
  # patchwork::wrap_plots(numbat_heatmap, scna_variability_plot, ncol = 1)
  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_numbat_heatmap{extension}.pdf"), numbat_heatmap)

  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_numbat_probability{extension}.pdf"), scna_variability_plot)

  # ## numbat phylo ------------------------------
  # numbat_heatmap_w_phylo <- safe_plot_numbat_w_phylo(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9)[["result"]]
  #
  # scna_variability_plot_w_phylo <- plot_variability_at_SCNA(numbat_heatmap_w_phylo[[3]][["data"]])
  #
  # patchwork::wrap_plots(numbat_heatmap_w_phylo, scna_variability_plot_w_phylo, ncol = 1)
  # ggsave(glue("results/{numbat_dir}/{sample_id}_numbat_phylo_probability{extension}.pdf"))

  plot_types <- c("numbat_heatmap", "numbat_probability")

  plot_files <- glue("results/{numbat_dir}/{sample_id}/{sample_id}_{plot_types}{extension}.pdf") %>%
    set_names(plot_types)

  return(plot_files)
  # return(numbat_heatmap)
  
}