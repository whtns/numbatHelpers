integration_by_scna_clones <- function(seu_paths, scna_of_interest = "1q", clone_comparisons, ...){
	seus <- seu_paths |> 
		set_names()  |> 
		set_names(str_extract, "SR[RX][0-9]+")
	
	seus <- seus |> 
		map(readRDS)

	seus <- imap(seus,
		subset_seu_by_clones, scna_of_interest = scna_of_interest, clone_comparisons = clone_comparisons, ...)
	
	integrated_seu <- seuratTools::integration_workflow(seus)

	seu_path <- tempfile(pattern = paste0("integrated_", scna_of_interest, "_"), tmpdir = "output/seurat", fileext = "_filtered_seu.rds")

	integrated_seu <- ScaleData(integrated_seu)
	integrated_seu <- seurat_reduce_dimensions(integrated_seu)
	integrated_seu <- seurat_cluster(integrated_seu, resolution = seq(0.2, 2.0, by = 0.2), seurat_assay = "integrated")

	scna_string <- 
		c("1q" = "1q+",
		"2p" = "2p+",
		"6p" = "6p+",
		"16q" = "16q-")[scna_of_interest]

	# integrated_seu$scna <- 
	integrated_seu$scna <- factor(ifelse(str_detect(integrated_seu$scna, pattern = scna_string), "w_scna", "wo_scna"), levels = c("wo_scna", "w_scna"))
	
	seus <- integrated_seu  |> 
	SplitObject(split.by = "batch")

	seus <- map(seus, function(seu) {
		colnames(seu) <- str_remove(colnames(seu), "_.*")
		seu
	})  |> 
	map(~add_hash_metadata(seu = .x))

	add_batch_hash_metadata(seu = integrated_seu, filepath = seu_path)

	return(seu_path)
}

subset_seu_by_clones <- function(seu, sample_id, scna_of_interest = "1q", clone_comparisons, filter_expr = NULL){
		clone_comparisons <- names(clone_comparisons[[sample_id]])
		retained_clones <- clone_comparisons |>
			str_extract("[0-9]_v_[0-9]") |>
			str_split("_v_", simplify = TRUE)

		mode(retained_clones) <- "integer"

		retained_clones <- retained_clones[which.min(rowSums(retained_clones)),]

		seu <- seu[, seu$clone_opt %in% retained_clones]

        # optionally subset the Seurat object by filter expression string (evaluated in @meta.data)
        if (!is.null(filter_expr)) {
            keep_cells <- tryCatch({
            rlang::eval_tidy(rlang::parse_expr(filter_expr), data = seu@meta.data)
            }, error = function(e) {
            warning("Failed to evaluate filter_expr ('", filter_expr, "') on ", sample_id, ": ", e$message)
            NULL
            })
            if (is.logical(keep_cells) && length(keep_cells) == ncol(seu)) {
            seu <- seu[, which(keep_cells)]
            message("Subsetting seu by filter_expr '", filter_expr, "' -> kept ", sum(keep_cells, na.rm = TRUE), " cells")
            } else if (is.numeric(keep_cells)) {
            seu <- seu[, keep_cells]
            message("Subsetting seu by numeric index from filter_expr, kept ", ncol(seu), " cells")
            } else {
            warning("filter_expr returned unexpected value; skipping subsetting")
            }
        }

		return(seu)
	}
