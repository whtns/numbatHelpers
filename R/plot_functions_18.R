# Plot Functions (110)

#' Extract or pull specific data elements
#'
#' @param nb_path File path
#' @param chrom Character string (default: "1")
#' @return Extracted data elements
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging
# - map_bind_rows: Use map_dfr() instead of map() %>% bind_rows() for better performance

pull_scna_segments <- function(nb_path = "output/numbat_sridhar/SRX11133594_numbat.rds", chrom = "1"){
	mynb <- readRDS(nb_path)
	
	segments <- mynb$clone_post %>%
		dplyr::left_join(mynb$joint_post, by = "cell") %>%
		# dplyr::filter(clone_opt %in% idents) %>%
		# dplyr::filter(seg %in% mysegs) %>%
		dplyr::distinct(CHROM, seg, seg_start, seg_end, cnv_state_map) %>%
		dplyr::mutate(seqnames = CHROM, start = seg_start, end = seg_end) %>%
		dplyr::filter(!cnv_state_map == "neu") %>%
		dplyr::filter(seqnames == chrom) |> 
		plyranges::as_granges() %>%
		identity()
}

#' Perform differential expression analysis
#'
#' @param seu Seurat object
#' @param ident_1 Cell identities or groups
#' @param ident_2 Cell identities or groups
#' @param mycluster Cluster information
#' @param segments Parameter for segments
#' @return Differential expression results
#' @export
diffex_per_cluster <- function(seu, ident_1 = "w_scna", ident_2 = "wo_scna", mycluster = 'g1_1', segments = NULL){
	seu <- seu[,seu$clusters %in% mycluster]
	
	diffex <- FindMarkers(seu, ident.1 = ident_1, ident.2 = ident_2, group = "scna", test.use = "wilcox", assay = "gene", logfc.threshold = 0.1) |> 
		tibble::rownames_to_column("symbol") %>%
		dplyr::left_join(annotables::grch38, by = "symbol") %>%
		dplyr::distinct(ensgene, .keep_all = TRUE) %>%
		dplyr::mutate(seqnames = chr) %>%
		dplyr::filter(!is.na(start), !is.na(end)) %>%
		plyranges::as_granges() |>
		identity()
	
	# cis 
	cis_diffex <- 
		diffex |> 
		plyranges::join_overlap_intersect(segments) %>%
		as_tibble() %>%
		dplyr::select(-any_of(c(
			"CHROM", "seg_start",
			"seg_end", "cnv_state_map", "log2_sign"
		))) %>%
		dplyr::filter(!str_detect(chr, "CHR_")) %>%
		dplyr::distinct(symbol, .keep_all = TRUE)
	
	out_of_segment_ranges <-
		diffex %>%
		plyranges::setdiff_ranges(segments)
	
	trans_diffex <-
		diffex %>%
		plyranges::join_overlap_intersect(out_of_segment_ranges) %>%
		as_tibble() %>%
		dplyr::filter(!str_detect(chr, "CHR_")) %>%
		dplyr::distinct(symbol, .keep_all = TRUE)
	
	diffex <- 
	list("cis" = cis_diffex, "trans" = trans_diffex) |> 
		dplyr::bind_rows(.id = "location")
	
	return(diffex)
	
}

#' Create a plot visualization
#'
#' @return ggplot2 plot object
#' @export
plot_fig_s20 <- function() {
	# 2p  ------------------------------
	
	segs_2p <- sapply(c("output/numbat_sridhar/SRX10264525_numbat.rds",
											 "output/numbat_sridhar/SRX14116944_numbat.rds"), pull_scna_segments, chrom = "2") |> 
		as("GRangesList") |> 
		unlist() |> 
		reduce_ranges()
	
	seu <- readRDS("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds")
	
	myclusters <- seu$clusters |> 
		levels() |> 
		set_names(identity)
	
	diffex_2p <- map(myclusters, ~diffex_per_cluster(seu, "w_scna", "wo_scna", mycluster = .x, segments = segs_2p)) |> 
		dplyr::bind_rows(.id = "cluster")
	
	# percent_2p <- 
		janitor::tabyl(seu@meta.data, clusters, scna) |> 
			adorn_percentages(denominator = "row") |> 
			adorn_pct_formatting() |>
			identity()
	
	diffex_2p <- 
		diffex_2p |> 
		dplyr::filter(p_val_adj <= 0.05) |> 
		dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
		dplyr::filter(avg_log2FC >= 0.5) |>
		# dplyr::filter(!cluster %in% c("hsp_8", "hypoxia_2")) |> 
		dplyr::arrange(desc(cluster), desc(avg_log2FC)) |> 
		dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) |> 
		identity()
	
	diffex_2p |> 
		dplyr::group_by(symbol) |> 
		dplyr::mutate(recurrence = n_distinct(cluster)) |> 
		dplyr::filter(recurrence >=4) |> 
		dplyr::filter(location == "cis") |> 
		ggplot(aes(y = symbol, x = cluster, size = neg_log_p_val_adj, color = avg_log2FC)) + 
		geom_point() + 
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
		scale_color_gradient2() +
		labs(title = "integrated 2p+ diffex per cluster")
	
	cis_path = "results/fig_s20a.pdf"
	ggsave(cis_path, w = 4, h = 5) |> 
		browseURL()
	
	all_diffex <- 
		diffex_2p |> 
		dplyr::filter(!is.na(entrez)) |> 
		dplyr::mutate(Gene.Name = symbol) |> 
		dplyr::select(-any_of(colnames(annotables::grch38))) |> 
		split(~cluster) |> 
		map(tibble::column_to_rownames, "Gene.Name") |> 
		identity()
	
	enrichments <- 
		map(all_diffex, enrichment_analysis, analysis_method = "ora") |> 
		purrr::compact() |> 
		map(DOSE::setReadable, org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
	
	enrichment_plots <- 
		enrichments |> 
		map(plot_enrichment, analysis_method = "ora") |> 
		imap(~{.x + labs(title = .y)})
	
	trans_path <- "results/fig_s20b.pdf"
	pdf(trans_path, h = 4, w = 6)
	print(enrichment_plots)
	dev.off()
	
	plot_path <- qpdf::pdf_combine(c(cis_path, trans_path), "results/fig_s20.pdf")
	
	table_path <- "results/table_s24.csv"
	mytable <- 
		enrichments |> 
		map("result") |> 
		dplyr::bind_rows(.id = "cluster") |> 
		label_enrichment_by_cis("02", "p") |> 
		write_csv(table_path)
	
	return(list("plot" = plot_path, "table" = table_path))
	
	
}

#' Create a plot visualization
#'
#' @param plot_path File path
#' @param table_path File path
#' @return ggplot2 plot object
#' @export
plot_fig_s10 <- function(plot_path = "results/fig_s10.pdf", table_path = "results/table_s23.csv") {
	# 16q  ------------------------------
	
	segs_16q <- sapply(c("output/numbat_sridhar/SRX11133594_numbat.rds",
											"output/numbat_sridhar/SRX11133593_numbat.rds",
											"output/numbat_sridhar/SRX11133592_numbat.rds"), pull_scna_segments, chrom = "16") |> 
		as("GRangesList") |> 
		unlist() |> 
		reduce_ranges()
	
	seu <- readRDS("output/seurat/integrated_16q/integrated_seu_16q_complete.rds")
	
	myclusters <- seu$clusters |> 
		levels() |> 
		set_names(identity)
	
	diffex_16q <- map(myclusters, ~diffex_per_cluster(seu, "16q-", ".diploid", mycluster = .x, segments = segs_16q)) |> 
		dplyr::bind_rows(.id = "cluster")
	
	diffex_16q <- 
		diffex_16q |> 
		dplyr::filter(p_val_adj <= 0.1) |> 
		dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
		dplyr::filter(abs(avg_log2FC) >= 0.5) |>
		# dplyr::filter(!cluster %in% c("hsp_8", "hypoxia_2")) |> 
		dplyr::arrange(desc(cluster), desc(avg_log2FC)) |> 
		dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) |> 
		identity()
	
	diffex_16q |> 
		dplyr::filter(location == "cis") |> 
		ggplot(aes(y = symbol, x = cluster, size = neg_log_p_val_adj, color = avg_log2FC)) + 
		geom_point() + 
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
		scale_color_gradient2() +
		labs(title = "integrated 16q+ diffex per cluster")
	
	cis_path = "results/fig_s10a.pdf"
	ggsave(cis_path, w = 4, h = 3) |> 
		browseURL()
	
	all_diffex <- 
		diffex_16q |> 
		dplyr::filter(!is.na(entrez)) |> 
		dplyr::mutate(Gene.Name = symbol) |> 
		dplyr::select(-any_of(colnames(annotables::grch38))) |> 
		split(~cluster) |> 
		map(tibble::column_to_rownames, "Gene.Name") |> 
		identity()
	
	enrichments <- 
		map(all_diffex, enrichment_analysis, analysis_method = "ora") |> 
		purrr::compact() |> 
		map(DOSE::setReadable, org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
	
	enrichment_plots <- 
		enrichments |> 
		map(plot_enrichment, analysis_method = "ora") |> 
		imap(~{.x + labs(title = .y)})
	
	trans_path <- "results/fig_s10b.pdf"
	pdf(trans_path, h = 4, w = 6)
	print(enrichment_plots)
	dev.off()
	
	plot_path <- qpdf::pdf_combine(c(cis_path, trans_path), plot_path)
	

	mytable <- 
		enrichments |> 
		map("result") |> 
		dplyr::bind_rows(.id = "cluster") |> 
		label_enrichment_by_cis("16", "q") |> 
		write_csv(table_path)
	
	return(list("plot" = plot_path, "table" = table_path))
	
	
}

#' Create a plot visualization
#'
#' @param plot_path File path
#' @param table_path File path
#' @return ggplot2 plot object
#' @export
plot_fig_s08 <- function(plot_path = "results/fig_s08.pdf", table_path = "results/table_s22.csv") {
	# 1q  ------------------------------
	
	segs_1q <- sapply(c("output/numbat_sridhar/SRX10264526_numbat.rds",
											"output/numbat_sridhar/SRX11133594_numbat.rds",
											"output/numbat_sridhar/SRX11133593_numbat.rds",
											"output/numbat_sridhar/SRX11133592_numbat.rds"), pull_scna_segments, chrom = "1") |> 
		as("GRangesList") |> 
		unlist() |> 
		reduce_ranges()
	
	seu <- readRDS("output/seurat/integrated_1q/integrated_seu_1q_complete.rds")
	
	myclusters <- seu$clusters |> 
		levels() |> 
		set_names(identity)
	
	diffex_1q <- map(myclusters, ~diffex_per_cluster(seu, mycluster = .x, segments = segs_1q)) |> 
		dplyr::bind_rows(.id = "cluster")
	
	diffex_1q <- 
		diffex_1q |> 
		dplyr::filter(p_val_adj <= 0.1) |> 
		dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
		dplyr::filter(avg_log2FC >= 0.5) |>
		dplyr::filter(!cluster %in% c("hsp_8", "hypoxia_2")) |> 
		dplyr::arrange(desc(cluster), desc(avg_log2FC)) |> 
		dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) |> 
		identity()
	
	plot_input <- 
	    diffex_1q |> 
	    dplyr::mutate(cluster = factor(cluster, levels = c("g1_3", "g1_0", "g1_6", "s_1", "s_g2_4"))) |> 
	    dplyr::arrange(desc(cluster), desc(avg_log2FC)) |> 
	    dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) |> 
		dplyr::filter(location == "cis")
	
	plot_input |> 
	ggplot(aes(y = symbol, x = cluster, size = neg_log_p_val_adj, color = avg_log2FC)) + 
		geom_point() + 
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
		scale_color_gradient2() +
		labs(title = "integrated 1q+ diffex per cluster")
	
	cis_path = "results/fig_s08a.pdf"
	ggsave(cis_path, w = 4, h = 7) |> 
		browseURL()
	
	all_diffex <- 
		diffex_1q |> 
		dplyr::filter(!is.na(entrez)) |> 
		dplyr::mutate(Gene.Name = symbol) |> 
		dplyr::select(-any_of(colnames(annotables::grch38))) |> 
		split(~cluster) |> 
		map(tibble::column_to_rownames, "Gene.Name") |> 
		identity()
	
	enrichments <- 
		map(all_diffex, enrichment_analysis, analysis_method = "ora") |> 
		purrr::compact() |> 
		map(DOSE::setReadable, org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
	
	enrichment_plots <- 
		enrichments |> 
		map(plot_enrichment, analysis_method = "ora") |> 
		imap(~{.x + labs(title = .y)})
	
	trans_path <- "results/fig_s08b.pdf"
	pdf(trans_path, h = 4, w = 6)
	print(enrichment_plots)
	dev.off()
	
	qpdf::pdf_combine(c(cis_path, trans_path), plot_path)
	
	mytable <- 
	enrichments |> 
		map("result") |> 
		dplyr::bind_rows(.id = "cluster") |> 
		label_enrichment_by_cis("01", "q") |> 
		write_csv(table_path)
	
	return(list("plot" = plot_path, "table" = table_path))
	
		
}

