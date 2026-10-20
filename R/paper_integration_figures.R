# Paper figures on the filtered integrations (integrated_seu_filtered_<scna>).
#
# The older figure functions (plot_integrated_1q_fig, plot_fig_03, plot_fig_s08/s09/
# s10/s20, make_table_s07) read hand-made 2024-25 objects that carry `clusters`
# ("<phase>_<cluster>") and `phase_level` columns. annotate_integrated_states() gives
# the new integrations the same columns from reviewed labels kept in
# batch_hashes.sqlite (file_id "paperint_<scna>"; R/paper_cluster_orders_db.R), and
# the functions below take paths instead of hard-coding them.

#' Label the clusters of a paper integration with cell states
#'
#' Takes the `paperint_<scna>` labels (read_paper_cluster_order()), adds `phase_level`
#' and `clusters` ("<phase>_<cluster>", ordered by phase then stored order), and
#' writes the integrated object and one object per batch to `out_dir`. The same
#' cluster order is registered in the sqlite `cluster_orders` table under each
#' written file name, which is how make_clone_distribution_figure() finds it.
#'
#' @param seu_path Path to the integrated Seurat object
#' @param scna "1q", "2p", "6p" or "16q"
#' @param cluster_order Output of read_paper_cluster_order(scna)
#' @param out_dir Output directory
#' @param sqlite_path Hash database
#' @return Named character vector of written paths, the integrated object first
#' @export
annotate_integrated_states <- function(seu_path, scna,
                                       cluster_order = read_paper_cluster_order(scna, sqlite_path),
                                       out_dir = file.path("output/seurat", paste0("integrated_filtered_", scna)),
                                       sqlite_path = "batch_hashes.sqlite") {
	key <- paste0("paperint_", scna)
	co <- cluster_order |>
		dplyr::mutate(label = paste(phase, clusters, sep = "_"))
	rc <- unique(co$resolution)

	seu <- readRDS(seu_path)
	if (!rc %in% colnames(seu@meta.data)) stop(rc, " not in ", seu_path)
	# the gene assay comes out of integration split into one layer per batch;
	# FindMarkers / presto / GetAssayData need it joined
	seu[["gene"]] <- SeuratObject::JoinLayers(seu[["gene"]])
	# the per-sample Phase calls are incomplete after merging; rescore on the joined
	# gene assay, as src/propose_integrated_cluster_states.R did for the review
	default_assay <- Seurat::DefaultAssay(seu)
	Seurat::DefaultAssay(seu) <- "gene"
	seu <- Seurat::CellCycleScoring(seu, s.features = Seurat::cc.genes.updated.2019$s.genes,
	                                g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)
	Seurat::DefaultAssay(seu) <- default_assay
	cl <- as.character(seu@meta.data[[rc]])
	missing <- setdiff(unique(cl), co$clusters)
	if (length(missing)) stop(key, ": clusters without a state at ", rc, ": ", paste(missing, collapse = ", "))

	seu$phase_level <- factor(as.character(co$phase)[match(cl, co$clusters)], levels = levels(co$phase))
	seu$clusters    <- factor(co$label[match(cl, co$clusters)], levels = co$label)

	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	int_path <- file.path(out_dir, paste0(key, "_filtered_seu.rds"))
	saveRDS(seu, int_path)

	per_paths <- Seurat::SplitObject(seu, split.by = "batch") |>
		purrr::imap_chr(function(x, b) {
			colnames(x) <- stringr::str_remove(colnames(x), "_.*")
			p <- file.path(out_dir, paste0(b, "_", key, "_filtered_seu.rds"))
			saveRDS(x, p)
			p
		})

	paths <- c(integrated = int_path, per_paths)
	ids <- fs::path_file(paths)
	# register the stored labels under each written file name for
	# make_clone_distribution_figure(), which looks them up by file name
	stored <- db_retry(read_cluster_orders_table(sqlite_path, file_id = key))[[key]]
	db_retry(encode_cluster_order_to_hash_table(rlang::set_names(rep(list(stored), length(ids)), ids), sqlite_path))
	paths
}

#' Pick the integrated object / the per-sample objects from annotate_integrated_states() output
#'
#' Selected by file name, since targets may drop the vector's names.
#'
#' @param paths Output of annotate_integrated_states()
#' @return A path; or a character vector named by sample
#' @export
paper_integrated_path <- function(paths) {
	p <- unname(paths[grepl("^paperint_", fs::path_file(paths))])
	stopifnot(length(p) == 1)
	p
}

#' @rdname paper_integrated_path
#' @export
paper_sample_paths <- function(paths) {
	p <- unname(paths[!grepl("^paperint_", fs::path_file(paths))])
	rlang::set_names(p, stringr::str_extract(fs::path_file(p), "SR[RX][0-9]+"))
}

#' Integrated clusters at alternative resolutions (Fig. S10-S12)
#'
#' A clustree of integrated_snn_res.*, then per resolution UMAPs by cluster, next to
#' the reviewed cell states (`clusters`), batch and SCNA status.
#'
#' @param seu_path Annotated integrated object
#' @param scna Label for titles
#' @param resolutions Resolutions to show
#' @param plot_path Output pdf
#' @return plot_path
#' @export
plot_integrated_resolutions <- function(seu_path, scna, resolutions = seq(0.2, 1.0, by = 0.2),
                                        plot_path = paste0("results/fig_", scna, "_alt_resolutions.pdf")) {
	seu <- readRDS(seu_path)
	rescols <- intersect(paste0("integrated_snn_res.", resolutions), colnames(seu@meta.data))
	md <- seu@meta.data[, rescols, drop = FALSE]
	pdf(plot_path, width = 16, height = 11)
	on.exit(dev.off(), add = TRUE)
	# clustree's edge_colourbar guide fails with this ggplot2/ggraph; drop the edge
	# colour legend, and lose only this page if it still fails
	tryCatch(print(clustree::clustree(md, prefix = "integrated_snn_res.") +
	                 guides(edge_colour = "none", edge_alpha = "none") + ggtitle(paste("integrated", scna))),
	         error = function(e) message("clustree failed: ", conditionMessage(e)))
	fixed <- Seurat::DimPlot(seu, group.by = "clusters", label = TRUE, repel = TRUE, raster = TRUE) +
		ggtitle("cell states")
	for (rc in rescols) {
		p <- Seurat::DimPlot(seu, group.by = rc, label = TRUE, raster = TRUE) + Seurat::NoLegend() + ggtitle(rc)
		print(p + fixed +
		      Seurat::DimPlot(seu, group.by = "batch", raster = TRUE) +
		      Seurat::DimPlot(seu, group.by = "scna", raster = TRUE) +
		      patchwork::plot_layout(ncol = 2) +
		      patchwork::plot_annotation(title = paste("integrated", scna, rc)))
	}
	plot_path
}

#' Integration overview figure (Fig. 4, Fig. 5, Fig. S21)
#'
#' Clone distribution figure of the integrated object, then the cell-cycle clone
#' plots of the integrated object and of each sample side by side.
#'
#' @param annotated_paths Output of annotate_integrated_states()
#' @param scna "1q", "2p", "6p" or "16q"
#' @param labeled_values Phase levels to label in the clone cc plots
#' @param plot_path Output pdf
#' @return list(plot, table)
#' @export
plot_paper_integration_fig <- function(annotated_paths, scna, labeled_values = "g2_m",
                                       plot_path = paste0("results/fig_", scna, "_integrated.pdf")) {
	int_path <- paper_integrated_path(annotated_paths)
	sample_paths <- paper_sample_paths(annotated_paths)

	page_a <- make_clone_distribution_figure(int_path, height = 12, width = 20,
		plot_path = tempfile(tmpdir = "results", fileext = ".pdf"),
		heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "batch"),
		heatmap_arrangement = c("clusters", "scna", "batch"))

	cc <- purrr::map(c(int_path, sample_paths), plot_clone_cc_plots, var_y = "clusters",
	                 scna_of_interest = scna, labeled_values = labeled_values)
	names(cc) <- c("integrated", names(sample_paths))
	table_path <- sub("\\.pdf$", "_cc.xlsx", plot_path)
	writexl::write_xlsx(purrr::map(cc, ~ .x$data), table_path)
	for (i in seq_along(cc)[-1]) cc[[i]] <- drop_y_axis(cc[[i]])
	cc_plot <- patchwork::wrap_plots(cc) +
		patchwork::plot_layout(nrow = 1, guides = "collect", axis_titles = "collect")
	page_b <- tempfile(tmpdir = "results", fileext = ".pdf")
	ggsave(page_b, cc_plot, width = 3 + 2.2 * length(cc), height = 5)

	qpdf::pdf_combine(c(page_a, page_b), plot_path)
	unlink(c(page_a, page_b))
	list(plot = plot_path, table = table_path)
}

#' Clone differential expression within each cell-state cluster (Fig. S13, S19, S24)
#'
#' w_scna vs wo_scna within each `clusters` level, cis genes (inside the samples'
#' SCNA segments on the chromosome) as a dot plot, ORA enrichment per cluster.
#' Gains keep genes up in w_scna; losses keep both directions.
#'
#' @param seu_path Annotated integrated object
#' @param nb_paths Numbat RDS paths of the integration samples
#' @param scna "1q", "2p", "6p" or "16q"
#' @param plot_path,table_path Outputs
#' @param p_val_adj_max,min_log2fc Filters
#' @param min_cells Clusters with fewer cells in either group are skipped
#' @param drop_phases Cell states left out of the tests (contaminant/artifact clusters)
#' @param min_recurrence,max_genes The cis dot plot shows genes significant in at least
#'   `min_recurrence` tested clusters, at most `max_genes` of them (most recurrent first,
#'   then by mean |log2FC|); every gene stays in the _genes.csv table
#' @return list(plot, table)
#' @export
plot_paper_clone_diffex <- function(seu_path, nb_paths, scna,
                                    plot_path = paste0("results/fig_", scna, "_clone_diffex_within_clusters.pdf"),
                                    table_path = paste0("results/table_", scna, "_clone_diffex_enrichment.csv"),
                                    p_val_adj_max = 0.1, min_log2fc = 0.5, min_cells = 20,
                                    drop_phases = "other", min_recurrence = 1, max_genes = Inf) {
	chrom <- stringr::str_extract(scna, "^[0-9]+")
	arm   <- stringr::str_extract(scna, "[pq]$")
	gain  <- scna != "16q"

	segs <- lapply(nb_paths, pull_scna_segments, chrom = chrom) |>
		methods::as("GRangesList") |> unlist() |> plyranges::reduce_ranges()

	seu <- readRDS(seu_path)
	Seurat::DefaultAssay(seu) <- "gene"
	n <- table(seu$clusters, seu$scna)
	keep <- rownames(n)[n[, "w_scna"] >= min_cells & n[, "wo_scna"] >= min_cells]
	dropped <- unique(as.character(seu$clusters[seu$phase_level %in% drop_phases]))
	keep <- setdiff(keep, dropped)
	message("clusters tested: ", paste(keep, collapse = ", "))

	diffex <- purrr::map(rlang::set_names(keep),
		~ diffex_per_cluster(seu, "w_scna", "wo_scna", mycluster = .x, segments = segs)) |>
		dplyr::bind_rows(.id = "cluster") |>
		dplyr::filter(p_val_adj <= p_val_adj_max,
		              if (gain) avg_log2FC >= min_log2fc else abs(avg_log2FC) >= min_log2fc) |>
		dplyr::mutate(neg_log_p_val_adj = -log10(p_val_adj),
		              cluster = factor(cluster, levels = keep))
	readr::write_csv(diffex, sub("\\.csv$", "_genes.csv", table_path))

	cis <- diffex |> dplyr::filter(location == "cis")
	shown <- cis |>
		dplyr::group_by(symbol) |>
		dplyr::summarise(recurrence = dplyr::n_distinct(cluster), mean_abs_fc = mean(abs(avg_log2FC))) |>
		dplyr::filter(recurrence >= min_recurrence) |>
		dplyr::arrange(dplyr::desc(recurrence), dplyr::desc(mean_abs_fc)) |>
		utils::head(max_genes)
	message("cis genes: ", dplyr::n_distinct(cis$symbol), "; plotted: ", nrow(shown))
	cis <- cis |>
		dplyr::filter(symbol %in% shown$symbol) |>
		dplyr::mutate(symbol = factor(symbol, levels = rev(shown$symbol)))
	cis_plot <- ggplot(cis, aes(y = symbol, x = cluster, size = neg_log_p_val_adj, color = avg_log2FC)) +
		geom_point() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
		scale_color_gradient2() +
		labs(title = paste0("integrated ", scna, " clone diffex per cluster (cis)"),
		     subtitle = sprintf("cis genes significant in >= %d clusters%s", min_recurrence,
		                        if (is.finite(max_genes)) sprintf(" (top %d)", max_genes) else ""))
	cis_path <- tempfile(tmpdir = "results", fileext = ".pdf")
	ggsave(cis_path, cis_plot, width = 5 + 0.4 * length(keep),
	       height = min(49, 2 + 0.14 * dplyr::n_distinct(cis$symbol)), limitsize = FALSE)

	all_diffex <- diffex |>
		dplyr::filter(!is.na(entrez)) |>
		dplyr::mutate(Gene.Name = symbol) |>
		dplyr::select(-any_of(colnames(annotables::grch38))) |>
		dplyr::group_by(cluster, Gene.Name) |>
		dplyr::slice_max(avg_log2FC, n = 1, with_ties = FALSE) |>
		dplyr::ungroup() |>
		split(~cluster) |>
		purrr::map(tibble::column_to_rownames, "Gene.Name")
	enrichments <- purrr::map(all_diffex, purrr::possibly(enrichment_analysis), analysis_method = "ora") |>
		purrr::compact() |>
		purrr::map(DOSE::setReadable, org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")

	trans_path <- tempfile(tmpdir = "results", fileext = ".pdf")
	pdf(trans_path, height = 4, width = 6)
	enrichments |>
		purrr::map(plot_enrichment, analysis_method = "ora") |>
		purrr::imap(~ .x + labs(title = .y)) |>
		print()
	dev.off()

	qpdf::pdf_combine(c(cis_path, trans_path), plot_path)
	unlink(c(cis_path, trans_path))

	enrichments |>
		purrr::map("result") |>
		dplyr::bind_rows(.id = "cluster") |>
		label_enrichment_by_cis(sprintf("%02d", as.integer(chrom)), arm) |>
		readr::write_csv(table_path)

	list(plot = plot_path, table = table_path)
}

#' Diffex and enrichment between chosen cell-state clusters (Fig. S20)
#'
#' @param seu_path Annotated integrated object
#' @param comparisons Data frame with columns ident.1, ident.2 (`clusters` levels)
#' @param plot_path,xlsx_path Outputs
#' @return list(diffex, enrichment)
#' @export
plot_paper_cluster_diffex_of_interest <- function(seu_path, comparisons,
                                                  plot_path = "results/fig_1q_cluster_diffex_of_interest.pdf",
                                                  xlsx_path = "results/fig_1q_cluster_diffex_of_interest.xlsx") {
	seu <- readRDS(seu_path)
	Seurat::DefaultAssay(seu) <- "gene"
	comparisons <- dplyr::mutate(comparisons, comparison = paste(ident.1, "v.", ident.2))

	diffexes <- purrr::map2(comparisons$ident.1, comparisons$ident.2,
		~ Seurat::FindMarkers(seu, group.by = "clusters", ident.1 = .x, ident.2 = .y)) |>
		rlang::set_names(comparisons$comparison)

	plots <- purrr::map(c("C6", "H"), function(gs) {
		diffexes |>
			purrr::map(enrichment_analysis, analysis_method = "ora", gene_set = gs) |>
			purrr::map(plot_enrichment, analysis_method = "ora") |>
			purrr::imap(~ .x + labs(title = paste(.y, gs)))
	})
	plots <- c(rbind(plots[[1]], plots[[2]])) |>
		purrr::map(~ .x + scale_y_discrete(labels = function(x) stringr::str_wrap(stringr::str_replace_all(x, "_", " "), width = 15)))

	pdf(plot_path, width = 6, height = 6)
	print(plots)
	dev.off()

	diffexes |>
		purrr::map(tibble::rownames_to_column, "symbol") |>
		writexl::write_xlsx(xlsx_path)

	list(diffex = xlsx_path, enrichment = plot_path)
}

#' % of w_scna / wo_scna cells per cell-state cluster (Table S10A-C)
#'
#' make_table_s07() on the annotated integration, leaving out `drop_phases` clusters
#' (contaminant/artifact states) before the percentages are taken.
#'
#' @param seu_path Annotated integrated object
#' @param table_path Output CSV
#' @param drop_phases Cell states left out
#' @return table_path
#' @export
make_paper_clone_per_cluster_table <- function(seu_path, table_path, drop_phases = "other") {
	seu <- readRDS(seu_path)
	md <- seu@meta.data[!seu$phase_level %in% drop_phases, , drop = FALSE]
	md$clusters <- droplevels(md$clusters)
	by_batch <- janitor::tabyl(md, clusters, scna, batch) |>
		janitor::adorn_percentages(denominator = "col") |>
		janitor::adorn_pct_formatting() |>
		dplyr::bind_rows(.id = "sample_id")
	janitor::tabyl(md, clusters, scna) |>
		janitor::adorn_percentages(denominator = "col") |>
		janitor::adorn_pct_formatting() |>
		dplyr::mutate(sample_id = "integrated") |>
		dplyr::bind_rows(by_batch) |>
		dplyr::arrange(clusters) |>
		readr::write_csv(table_path)
	table_path
}
