# Plotting and annotation functions (4)

#' Filter data based on specified criteria
#'
#' @param numbat_rds_file File path
#' @param seus Parameter for seus
#' @param cluster_dictionary Cluster information
#' @param large_clone_simplifications Parameter for large clone simplifications
#' @param filter_expressions Parameter for filter expressions
#' @param cells_to_remove Cell identifiers or information
#' @param extension Character string (default: "")
#' @param leiden_cluster_file File path
#' @return Modified Seurat object
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O
# - rownames_round_trip: Avoid unnecessary rownames conversions - keep as column when possible
# - multiple_joins: Combine multiple joins into single join operation where possible

annotate_filter_reason <- function(
  seu,
  sample_id,
  cluster_dictionary,
  cells_to_remove = NULL,
  mito_threshold = 10,
  nCount_threshold = 1000,
  nFeature_threshold = 1000
) {
  if (is.null(cluster_dictionary[[sample_id]])) {
    stop("cluster_dictionary has no entry for ", sample_id, "; cannot annotate filter reasons")
  }

  clusters_to_remove <-
    cluster_dictionary[[sample_id]] %>%
    dplyr::filter(remove == "1") %>%
    dplyr::pull(`gene_snn_res.0.2`)

  cells_to_drop <- character(0)
  if (!is.null(cells_to_remove) && !is.null(cells_to_remove[[sample_id]])) {
    cells_to_drop <- cells_to_remove[[sample_id]][["cell"]]
    if (is.null(cells_to_drop)) {
      cells_to_drop <- character(0)
    }
  }

  meta <- seu@meta.data
  has_clone_col <- "clone_opt" %in% colnames(meta)
  has_cluster_col <- "gene_snn_res.0.2" %in% colnames(meta)
  has_abbrev_col <- "abbreviation" %in% colnames(meta)

  clone_missing <- if (has_clone_col) {
    is.na(meta$clone_opt)
  } else {
    rep(TRUE, nrow(meta))
  }

  qc_fail <-
    meta$percent.mt >= mito_threshold |
    meta$nFeature_gene <= nFeature_threshold |
    meta$nCount_gene <= nCount_threshold

  cluster_remove <- if (has_cluster_col) {
    meta$gene_snn_res.0.2 %in% clusters_to_remove
  } else {
    rep(FALSE, nrow(meta))
  }

  malat1_remove <- if (has_abbrev_col) {
    !is.na(meta$abbreviation) & meta$abbreviation == "MALAT1"
  } else {
    rep(FALSE, nrow(meta))
  }

  manual_remove <- if (length(cells_to_drop) > 0) {
    rownames(meta) %in% cells_to_drop
  } else {
    rep(FALSE, nrow(meta))
  }

  filter_reason <- dplyr::case_when(
    cluster_remove ~ "cluster_remove",
    clone_missing ~ "clone_opt_na",
    qc_fail ~ "qc_fail",
    malat1_remove ~ "malat1",
    manual_remove ~ "manual_exclude",
    TRUE ~ NA_character_
  )

  filter_keep <- is.na(filter_reason)

  reason_meta <- data.frame(
    filter_reason = filter_reason,
    filter_keep = filter_keep,
    stringsAsFactors = FALSE,
    row.names = rownames(meta)
  )

  Seurat::AddMetaData(seu, metadata = reason_meta)
}

filter_cluster_save_seu <- function(numbat_rds_file, seus, cluster_dictionary, large_clone_simplifications, filter_expressions = NULL, cells_to_remove, extension = "", leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv", ...) {
  
  sample_id <- str_extract(numbat_rds_file, "SR[RX][0-9]+")
  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]
  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))
  names(seus) <- str_extract(seus, "SR[RX][0-9]+")
  seu_path <- seus[[sample_id]]
  if (is.null(seu_path) || is.na(seu_path)) {
    message("No unfiltered seu for ", sample_id, "; returning NA")
    return(NA_character_)
  }
  seu <- readRDS(seu_path)
  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))
  mynb <- readRDS(numbat_rds_file)
  if (is.null(mynb[["clone_post"]]) || !all(c("clone_opt", "GT_opt") %in% colnames(mynb[["clone_post"]]))) {
    message("clone_post missing required columns (clone_opt, GT_opt) for ", sample_id, "; returning NA")
    return(NA_character_)
  }
  nb_clone_post <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")
  seu <- Seurat::AddMetaData(seu, nb_clone_post)
  all_cells_meta <- seu@meta.data
  if (!"gene_snn_res.0.2" %in% colnames(seu@meta.data)) {
    message("gene_snn_res.0.2 missing for ", sample_id, "; re-running FindNeighbors + seurat_cluster")
    if (!"gene_nn" %in% names(seu@graphs)) {
      seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
    }
    seu <- seurat_cluster(seu = seu, resolution = seq(0.2, 1.0, by = 0.2), reduction = "pca")
  }
  if (is.null(cluster_dictionary[[sample_id]])) {
    stop("cluster_dictionary has no entry for ", sample_id, "; cannot assign abbreviations")
  }
  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")
  seu <- AddMetaData(seu, test0)
  if (is.null(large_clone_simplifications[[sample_id]])) {
    large_clone_simplifications <- tibble::tibble(scna = character(), seg = character())
  } else {
    large_clone_simplifications <-
      tibble::enframe(large_clone_simplifications[[sample_id]], "scna", "seg") %>%
      tidyr::unnest(seg) %>%
      dplyr::mutate(seg = as.character(seg))
  }

  scna_labels <- vapply(
    nb_clone_post$GT_opt,
    function(gt_opt) {
      wrapped <- wrap_scna_labels(simplify_gt_col(gt_opt, large_clone_simplifications))
      if (length(wrapped) == 0) {
        return(NA_character_)
      }
      as.character(wrapped[[1]])
    },
    FUN.VALUE = character(1)
  )

  scna_metadata <- data.frame(
    scna = scna_labels,
    row.names = rownames(nb_clone_post),
    stringsAsFactors = FALSE
  )

  scna_metadata <- scna_metadata[colnames(seu), , drop = FALSE]

  seu <- Seurat::AddMetaData(seu, scna_metadata)

  if (!all(c("filter_reason", "filter_keep") %in% colnames(seu@meta.data))) {
    stop("Missing filter metadata columns in input Seurat object for ", sample_id,
         ". Expected columns: filter_reason, filter_keep.")
  }

  scna_meta <- seu@meta.data[seu$filter_reason %in% c("qc_fail", "cluster_remove", "malat1", "manual_exclude") | is.na(seu$filter_reason), ]
  qc_meta <- seu@meta.data[seu$filter_reason %in% c("cluster_remove", "malat1", "manual_exclude") | is.na(seu$filter_reason), ]
  seu <- seu[, seu$filter_keep]

  filtered_seu_path <- glue("output/seurat/{sample_id}_filtered{extension}_seu.rds")

  # Cache check: if existing filtered_seu has the same post-filter cells, skip
  # expensive recomputation and only refresh metadata.
  meta_cols_to_refresh <- c("clone_opt", "GT_opt", "scna", "filter_reason", "filter_keep", "abbreviation")
  if (file.exists(filtered_seu_path)) {
    cached_seu <- tryCatch(readRDS(filtered_seu_path), error = function(e) NULL)
    if (!is.null(cached_seu) && identical(sort(colnames(cached_seu)), sort(colnames(seu)))) {
      message("Cell barcodes unchanged for ", sample_id, "; refreshing metadata only (skipping SCTransform/PCA/UMAP)")
      refresh_meta <- seu@meta.data[, base::intersect(meta_cols_to_refresh, colnames(seu@meta.data)), drop = FALSE]
      cached_seu <- Seurat::AddMetaData(cached_seu, refresh_meta)
      Project(cached_seu) <- sample_id
      add_hash_metadata(seu = cached_seu, filepath = filtered_seu_path)
      save_cell_barcodes_to_db(filtered_seu_path, sample_id, "filtered", colnames(cached_seu))
      cell_type_meta <- cached_seu@meta.data
      plot_filtering_timeline(all_cells_meta, scna_meta, qc_meta, cell_type_meta, sample_id)
      ggsave(glue("results/{sample_id}_filtering_timeline_{extension}.pdf"), width = 8, height = 4)
      return(filtered_seu_path)
    }
    message("Cell barcodes changed for ", sample_id, "; running full SCTransform/PCA/UMAP pipeline")
  }

  seu <- SCTransform(seu, assay = "gene", verbose = FALSE)
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- seurat_cluster(
    seu = seu, resolution = seq(0.2, 1.0, by = 0.2),
    reduction = "pca"
  )
  seu <- tryCatch(
    find_all_markers(seu, seurat_assay = "SCT"),
    error = function(e) {
      if (grepl("JoinLayers", conditionMessage(e), fixed = TRUE)) {
        warning("SCT marker JoinLayers failed for ", sample_id, "; using FindAllMarkers fallback.")
        seu <- PrepSCTFindMarkers(seu)
        seu@misc$markers[["clusters"]] <- FindAllMarkers(seu, assay = "SCT", verbose = FALSE)
        return(seu)
      }
      stop(e)
    }
  )
  Project(seu) <- sample_id
  add_hash_metadata(seu = seu, filepath = filtered_seu_path)
  save_cell_barcodes_to_db(filtered_seu_path, sample_id, "filtered", colnames(seu))

  cell_type_meta <- seu@meta.data
  plot_filtering_timeline(all_cells_meta, scna_meta, qc_meta, cell_type_meta, sample_id)
  ggsave(glue("results/{sample_id}_filtering_timeline_{extension}.pdf"), width = 8, height = 4)
  return(filtered_seu_path)
}

#' Filter data based on specified criteria
#'
#' @param numbat_rds_file File path
#' @param cluster_dictionary Cluster information
#' @param large_clone_simplifications Parameter for large clone simplifications
#' @param filter_expressions Parameter for filter expressions
#' @param cells_to_remove Cell identifiers or information
#' @param extension Character string (default: "")
#' @return Modified Seurat object
#' @export
prep_unfiltered_seu <- function(numbat_rds_file, cluster_dictionary, large_clone_simplifications, filter_expressions = NULL, cells_to_remove = NULL, extension = "") {
  options(future.globals.maxSize = Inf)
  sample_id <- str_extract(numbat_rds_file, "SR[RX][0-9]+")
  numbat_dir <- fs::path_split(numbat_rds_file)[[1]][[2]]
  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))
  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))
  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))
  mynb <- readRDS(numbat_rds_file)
  if (is.null(mynb[["clone_post"]]) || !all(c("clone_opt", "GT_opt") %in% colnames(mynb[["clone_post"]]))) {
    message("clone_post missing required columns (clone_opt, GT_opt) for ", sample_id, "; returning NA")
    return(NA_character_)
  }
  nb_clone_post <- mynb[["clone_post"]][, c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")
  seu <- Seurat::AddMetaData(seu, nb_clone_post)
  if (!"gene_snn_res.0.2" %in% colnames(seu@meta.data)) {
    message("gene_snn_res.0.2 missing for ", sample_id, "; re-running FindNeighbors + seurat_cluster")
    if (!"gene_nn" %in% names(seu@graphs)) {
      seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
    }
    seu <- seurat_cluster(seu = seu, resolution = c(0.2, 0.4, 0.6), reduction = "pca")
  }
  if (is.null(cluster_dictionary[[sample_id]])) {
    stop("cluster_dictionary has no entry for ", sample_id, "; cannot assign abbreviations")
  }
  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")
  seu <- AddMetaData(seu, test0)
  if (is.null(large_clone_simplifications[[sample_id]])) {
    large_clone_simplifications <- tibble::tibble(scna = character(), seg = character())
  } else {
    large_clone_simplifications <-
      tibble::enframe(large_clone_simplifications[[sample_id]], "scna", "seg") %>%
      tidyr::unnest(seg) %>%
      dplyr::mutate(seg = as.character(seg))
  }

  scna_labels <- vapply(
    nb_clone_post$GT_opt,
    function(gt_opt) {
      wrapped <- wrap_scna_labels(simplify_gt_col(gt_opt, large_clone_simplifications))
      if (length(wrapped) == 0) {
        return(NA_character_)
      }
      as.character(wrapped[[1]])
    },
    FUN.VALUE = character(1)
  )

  scna_metadata <- data.frame(
    scna = scna_labels,
    row.names = rownames(nb_clone_post),
    stringsAsFactors = FALSE
  )

  scna_metadata <- scna_metadata[colnames(seu), , drop = FALSE]
  seu <- Seurat::AddMetaData(seu, scna_metadata)

  seu <- annotate_filter_reason(
    seu = seu,
    sample_id = sample_id,
    cluster_dictionary = cluster_dictionary
    # cells_to_remove = cells_to_remove
  )

  unfiltered_seu_path <- glue("output/seurat/{sample_id}_unfiltered_seu.rds")

  # Cache check: if existing unfiltered_seu has the same cells, skip expensive
  # recomputation (SCTransform/PCA/UMAP/clustering/markers) and only refresh metadata.
  meta_cols_to_refresh <- c("clone_opt", "GT_opt", "scna", "filter_reason", "filter_keep", "abbreviation")
  if (file.exists(unfiltered_seu_path)) {
    cached_seu <- tryCatch(readRDS(unfiltered_seu_path), error = function(e) NULL)
    if (!is.null(cached_seu) && identical(sort(colnames(cached_seu)), sort(colnames(seu)))) {
      message("Cell barcodes unchanged for ", sample_id, "; refreshing metadata only (skipping SCTransform/PCA/UMAP)")
      refresh_meta <- seu@meta.data[, base::intersect(meta_cols_to_refresh, colnames(seu@meta.data)), drop = FALSE]
      cached_seu <- Seurat::AddMetaData(cached_seu, refresh_meta)
      add_hash_metadata(seu = cached_seu, filepath = unfiltered_seu_path)
      return(unfiltered_seu_path)
    }
    message("Cell barcodes changed for ", sample_id, "; running full SCTransform/PCA/UMAP pipeline")
  }

  seu <- SCTransform(seu, assay = "gene", verbose = FALSE)
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- seurat_cluster(
    seu = seu, resolution = c(0.2, 0.4, 0.6),
    reduction = "pca"
  )
  seu <- tryCatch(
    find_all_markers(seu, seurat_assay = "SCT"),
    error = function(e) {
      if (grepl("JoinLayers", conditionMessage(e), fixed = TRUE)) {
        warning("SCT marker JoinLayers failed for ", sample_id, "; using FindAllMarkers fallback.")
        seu <- PrepSCTFindMarkers(seu)
        seu@misc$markers[["clusters"]] <- FindAllMarkers(seu, assay = "SCT", verbose = FALSE)
        return(seu)
      }
      stop(e)
    }
  )

  add_hash_metadata(seu = seu, filepath = unfiltered_seu_path)
  return(unfiltered_seu_path)
}

#' Filter data based on specified criteria
#'
#' @param filtered_seu_path File path
#' @return Modified Seurat object
#' @export
regress_filtered_seu <- function(filtered_seu_path) {
  sample_id <- str_extract(filtered_seu_path, "SR[RX][0-9]+")
  regressed_seu_path <- str_replace(filtered_seu_path, "_filtered", "_regressed")
  regressed_seu <- readRDS(filtered_seu_path)

  # Cache check: regression doesn't change cell barcodes, so if the cached file
  # exists with the same cells as the input, skip the expensive recomputation.
  if (file.exists(regressed_seu_path)) {
    cached_seu <- tryCatch(readRDS(regressed_seu_path), error = function(e) NULL)
    if (!is.null(cached_seu) && identical(sort(colnames(cached_seu)), sort(colnames(regressed_seu)))) {
      message("Cell barcodes unchanged for ", sample_id, "; skipping SCTransform/PCA/UMAP regression")
      return(regressed_seu_path)
    }
    message("Cell barcodes changed for ", sample_id, "; running full regression pipeline")
  }

  regressed_seu <- PercentageFeatureSet(regressed_seu, pattern = "^MT-", col.name = "percent.mt")
  regressed_seu <- SCTransform(regressed_seu, assay = "gene", vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = FALSE)
  regressed_seu <- RunPCA(regressed_seu, verbose = FALSE)
  regressed_seu <- RunUMAP(regressed_seu, dims = 1:30, verbose = FALSE)
  regressed_seu <- FindNeighbors(regressed_seu, dims = 1:30, verbose = FALSE)
  regressed_seu <- seurat_cluster(
    seu = regressed_seu, resolution = seq(0.2, 1.0, by = 0.2),
    reduction = "pca"
  )
  regressed_seu <- tryCatch(
    find_all_markers(regressed_seu, seurat_assay = "SCT"),
    error = function(e) {
      if (grepl("JoinLayers", conditionMessage(e), fixed = TRUE)) {
        warning("SCT marker JoinLayers failed for ", sample_id, "; using FindAllMarkers fallback.")
        regressed_seu <- PrepSCTFindMarkers(regressed_seu)
        regressed_seu@misc$markers[["clusters"]] <- FindAllMarkers(regressed_seu, assay = "SCT", verbose = FALSE)
        return(regressed_seu)
      }
      stop(e)
    }
  )
  add_hash_metadata(seu = regressed_seu, filepath = regressed_seu_path)
  return(regressed_seu_path)
}


filter_sample_qc <- function(seu, mito_threshold = 10, nCount_threshold = 1000, nFeature_threshold = 1000) {
  #
  seu <-
    seu %>%
    subset(subset = percent.mt < mito_threshold & nFeature_gene > nFeature_threshold & nCount_gene > nCount_threshold) %>%
    identity()
}

#' Extract per-cell filtering metadata from an unfiltered Seurat RDS
#'
#' Reads the unfiltered Seurat object (output of prep_unfiltered_seu, which already
#' contains clone_opt, GT_opt, gene_snn_res.0.2, abbreviation, and scna metadata)
#' and returns a flat tibble of per-cell filtering flags.  No SCTransform, PCA,
#' UMAP, or marker detection is performed.
#'
#' @param seu_path Path to an unfiltered Seurat RDS (element of unfiltered_seus target).
#' @param cluster_dictionary Named list of per-sample tibbles with columns
#'   gene_snn_res.0.2, abbreviation, remove.
#' @param cells_to_remove Named list of per-sample tibbles with a "cell" column
#'   listing manually excluded barcodes.
#' @param large_clone_simplifications Named list of per-sample SCNA simplification
#'   vectors (passed through to identify sample_id only; not used for metadata
#'   derivation since scna is already in the unfiltered object).
#' @return A tibble with one row per cell and columns: cell, sample_id,
#'   nCount_gene, nFeature_gene, percent.mt, clone_opt_is_na, cluster,
#'   abbreviation, cluster_remove_flag, is_malat1, in_manual_exclude, scna.
#' @export
extract_filter_metadata <- function(
  seu_path,
  cluster_dictionary,
  cells_to_remove,
  large_clone_simplifications
) {
  if (is.na(seu_path)) return(NA)
  sample_id <- stringr::str_extract(seu_path, "SR[RX][0-9]+")
  seu <- readRDS(seu_path)

  meta <- seu@meta.data |>
    tibble::rownames_to_column("cell") |>
    tibble::as_tibble() |>
    dplyr::mutate(
      sample_id      = sample_id,
      clone_opt_is_na = is.na(clone_opt),
      cluster        = as.numeric(as.character(gene_snn_res.0.2))
    )

  # cluster_remove_flag from cluster_dictionary
  if (!is.null(cluster_dictionary[[sample_id]])) {
    clusters_to_remove <- cluster_dictionary[[sample_id]] |>
      dplyr::filter(remove == "1") |>
      dplyr::pull(gene_snn_res.0.2)
    meta <- meta |>
      dplyr::mutate(cluster_remove_flag = cluster %in% clusters_to_remove)
  } else {
    meta <- meta |> dplyr::mutate(cluster_remove_flag = FALSE)
  }

  meta <- meta |>
    dplyr::mutate(is_malat1 = !is.na(abbreviation) & abbreviation == "MALAT1")

  manual_cells <- cells_to_remove[[sample_id]][["cell"]]
  if (is.null(manual_cells)) manual_cells <- character(0)
  meta <- meta |>
    dplyr::mutate(in_manual_exclude = cell %in% manual_cells)

  meta |>
    dplyr::select(
      cell, sample_id,
      nCount_gene, nFeature_gene, percent.mt,
      clone_opt_is_na, cluster, abbreviation,
      cluster_remove_flag, is_malat1, in_manual_exclude,
      scna,
      dplyr::any_of("filter_keep")
    )
}

#' Simulate filtering pipeline on per-cell metadata and return staged cell counts
#'
#' Pure dplyr operation on the tibble from extract_filter_metadata.  No Seurat I/O.
#' Stages mirror the four snapshots used by plot_filtering_timeline.
#'
#' @param filter_meta Tibble from extract_filter_metadata (one row per cell).
#' @param mito_threshold Maximum percent.mt to retain (default 10).
#' @param nCount_threshold Minimum nCount_gene to retain (default 1000).
#' @param nFeature_threshold Minimum nFeature_gene to retain (default 1000).
#' @param include_cluster_removal Logical; apply cluster_remove_flag step (default TRUE).
#' @param include_malat1_removal Logical; remove MALAT1 cluster cells (default TRUE).
#' @param include_manual_removal Logical; apply manual exclusion list (default TRUE).
#' @return A tibble with columns: stage, scna, n_cells, pct_remaining.
#'   pct_remaining is the percentage of the pre-filter total at that stage.
#' @export
apply_filter_criteria <- function(
  filter_meta,
  mito_threshold          = 10,
  nCount_threshold        = 1000,
  nFeature_threshold      = 1000,
  include_cluster_removal = TRUE,
  include_malat1_removal  = TRUE,
  include_manual_removal  = TRUE
) {
  n_total <- nrow(filter_meta)
  stage_cells <- list()

  stage_cells[["clone_na"]] <- filter_meta |>
    dplyr::filter(!clone_opt_is_na)

  stage_cells[["qc"]] <- stage_cells[["clone_na"]] |>
    dplyr::filter(
      percent.mt    < mito_threshold,
      nCount_gene   > nCount_threshold,
      nFeature_gene > nFeature_threshold
    )

  stage_cells[["cluster_removal"]] <- stage_cells[["qc"]]
  if (include_cluster_removal) {
    stage_cells[["cluster_removal"]] <- stage_cells[["cluster_removal"]] |>
      dplyr::filter(!cluster_remove_flag)
  }
  if (include_malat1_removal) {
    stage_cells[["cluster_removal"]] <- stage_cells[["cluster_removal"]] |>
      dplyr::filter(!is_malat1)
  }

  stage_cells[["manual"]] <- stage_cells[["cluster_removal"]]
  if (include_manual_removal) {
    stage_cells[["manual"]] <- stage_cells[["manual"]] |>
      dplyr::filter(!in_manual_exclude)
  }

  stage_names <- c("clone_na", "qc", "cluster_removal", "manual")

  purrr::map_dfr(stage_names, function(stage) {
    stage_cells[[stage]] |>
      dplyr::count(scna, name = "n_cells") |>
      dplyr::mutate(
        stage         = stage,
        pct_remaining = round(100 * n_cells / n_total, 1)
      )
  }) |>
    dplyr::mutate(stage = factor(stage, levels = stage_names)) |>
    dplyr::select(stage, scna, n_cells, pct_remaining)
}

#' Sweep QC threshold combinations and plot cell survival landscape
#'
#' Calls apply_filter_criteria over a grid of mito/nCount/nFeature thresholds and
#' returns a faceted ggplot.  No Seurat I/O — operates entirely on the tibble from
#' extract_filter_metadata.
#'
#' @param filter_meta Tibble from extract_filter_metadata (can be combined across
#'   samples with dplyr::bind_rows).
#' @param param_grid Optional data frame with columns mito, nCount, nFeature.
#'   Defaults to expand.grid of c(5,10,15,20) x c(500,1000,2000) x c(500,1000,2000).
#' @return A ggplot2 object.
#' @export
plot_filter_sweep <- function(filter_meta, param_grid = NULL) {
  if (is.null(param_grid)) {
    param_grid <- expand.grid(
      mito     = c(5, 10, 15, 20),
      nCount   = c(500, 1000, 2000),
      nFeature = c(500, 1000, 2000),
      stringsAsFactors = FALSE
    )
  }

  results <- purrr::map_dfr(seq_len(nrow(param_grid)), function(i) {
    mito     <- param_grid[["mito"]][i]
    nCount   <- param_grid[["nCount"]][i]
    nFeature <- param_grid[["nFeature"]][i]
    apply_filter_criteria(
      filter_meta,
      mito_threshold     = mito,
      nCount_threshold   = nCount,
      nFeature_threshold = nFeature
    ) |>
      dplyr::filter(stage == "manual") |>
      dplyr::mutate(mito = mito, nCount = nCount, nFeature = nFeature)
  })

  results_long <- results |>
    tidyr::pivot_longer(
      cols      = c("mito", "nCount", "nFeature"),
      names_to  = "qc_metric",
      values_to = "threshold"
    )

  ggplot2::ggplot(results_long, ggplot2::aes(
    x     = threshold,
    y     = pct_remaining,
    color = scna,
    group = interaction(scna, qc_metric)
  )) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ qc_metric, scales = "free_x") +
    ggplot2::labs(
      x     = "Threshold value",
      y     = "% cells remaining (post all filters)",
      color = "SCNA"
    ) +
    ggplot2::theme_bw()
}

#' Read all hypoxia Seurat RDS files and pull hypoxia_score metadata
#'
#' @param dir Directory to search for hypoxia Seurat files. Defaults to "output/seurat".
#' @param pattern Regex pattern to match hypoxia Seurat files. Defaults to "_seu_hypoxia.rds$".
#' @param recursive Whether to search recursively. Defaults to TRUE.
#' @return A tibble with columns: sample_id, file, cell, hypoxia_score
#' @export
read_all_hypoxia_scores <- function(files) {

    res <- purrr::map_dfr(files, function(f) {
  tryCatch({
    sample_id <- stringr::str_replace(basename(f), "^([^_]+)_.*$", "\\1")
    seu <- readRDS(f)                       # errors here are caught
    md <- seu@meta.data                     # errors here are also caught
    if (!"hypoxia_score" %in% colnames(md)) {
      warning("No hypoxia_score in meta.data for ", f)
      return(tibble::tibble(sample_id = sample_id, file = f, cell = rownames(md), hypoxia_score = NA_real_))
    }
    tibble::tibble(
      sample_id = sample_id,
      file = f,
      cell = rownames(md),
      hypoxia_score = md[["hypoxia_score"]]
    )
  }, error = function(e) {
    warning("Failed processing ", f, ": ", e$message)
    # return an empty tibble row (or whatever shape you prefer)
    tibble::tibble(sample_id = stringr::str_replace(basename(f), "^([^_]+)_.*$", "\\1"),
                   file = f,
                   cell = NA_character_,
                   hypoxia_score = NA_real_)
  })
})

    return(res)
}

assemble_diploid_seu <- function(filtered_seus_paths,
                                 out_path = "output/seurat/diploid_seu.rds",
                                 integrate = TRUE) {
  paths <- unlist(filtered_seus_paths)
  paths <- paths[!is.na(paths)]
  sample_ids <- stringr::str_extract(paths, "SR[RX][0-9]+")

  cone_seus <- purrr::imap(
    purrr::set_names(paths, sample_ids),
    function(path, sample_id) {
      seu <- readRDS(path)
      diploid_mask <- is.na(seu$scna) | seu$scna == ""
      if (!any(diploid_mask)) return(NULL)
      seu <- seu[, diploid_mask]
      for (assay_name in SeuratObject::Assays(seu)) {
        if (inherits(seu[[assay_name]], "Assay5"))
          seu[[assay_name]] <- SeuratObject::JoinLayers(seu[[assay_name]])
      }
      seu$sample_source <- sample_id
      seu <- plot_celltype_predictions(seu, sample_id = sample_id, group.by = "gene_snn_res.0.2")$seu
      cone_mask <- tolower(seu$type) == "cones"
      if (!any(cone_mask)) return(NULL)
      seu[, cone_mask]
    }
  ) |>
    purrr::compact()

  if (integrate) {
    old_assay_version <- getOption("Seurat.object.assay.version")
    options(Seurat.object.assay.version = "v3")
    on.exit(options(Seurat.object.assay.version = old_assay_version), add = TRUE)

    cone_seus <- lapply(cone_seus, function(s) {
      Seurat::DefaultAssay(s) <- "gene"
      if (inherits(s[["gene"]], "Assay5")) {
        counts_mat <- GetAssayData(s[["gene"]], layer = "counts")
        s[["gene"]] <- CreateAssayObject(counts = counts_mat)
      }
      s[["SCT"]] <- NULL
      s
    })

    seu <- seuratTools::integration_workflow(cone_seus, resolution = c(0.2, 0.4), find_markers = FALSE)
  } else {
    seus_list <- as.list(cone_seus)
    seu <- seus_list[[1]]
    for (i in seq_along(seus_list)[-1]) {
      seu <- merge(seu, seus_list[[i]])
      for (assay_name in SeuratObject::Assays(seu)) {
        if (inherits(seu[[assay_name]], "Assay5"))
          seu[[assay_name]] <- SeuratObject::JoinLayers(seu[[assay_name]])
      }
    }
  }

  add_hash_metadata(seu = seu, filepath = out_path)
  out_path
}

merge_hypoxia_with_diploid <- function(hypoxia_seu_path, diploid_seu_path, slug) {
  sample_id <- stringr::str_extract(hypoxia_seu_path, "SR[RX][0-9]+")

  old_assay_version <- getOption("Seurat.object.assay.version")
  options(Seurat.object.assay.version = "v3")
  on.exit(options(Seurat.object.assay.version = old_assay_version), add = TRUE)

  prep_v3 <- function(s) {
    Seurat::DefaultAssay(s) <- "gene"
    if (inherits(s[["gene"]], "Assay5")) {
      counts_mat <- GetAssayData(s[["gene"]], layer = "counts")
      s[["gene"]] <- CreateAssayObject(counts = counts_mat)
    }
    s[["SCT"]] <- NULL
    s[["integrated"]] <- NULL
    s
  }

  hypoxia_seu <- prep_v3(readRDS(hypoxia_seu_path))
  diploid_seu  <- prep_v3(readRDS(diploid_seu_path))

  shared_cells <- intersect(colnames(hypoxia_seu), colnames(diploid_seu))
  if (length(shared_cells) > 0) {
    diploid_seu <- diploid_seu[, !colnames(diploid_seu) %in% shared_cells]
  }

  merged <- list(diploid = diploid_seu)
  merged[[sample_id]] <- hypoxia_seu
  seu <- seuratTools::integration_workflow(merged, resolution = c(0.2, 0.4), find_markers = FALSE)

  out_path <- glue::glue("output/seurat/{sample_id}_{slug}_diploid_merged.rds")
  saveRDS(seu, out_path)
  out_path
}

#' Add cell cycle Phase to a filtered Seurat object and save
#'
#' @param seu_path File path to filtered Seurat RDS
#' @param organism Character string (default: "human")
#' @return Path to updated RDS file
#' @export
add_phase_to_filtered_seu <- function(seu_path, organism = "human") {
  if (is.na(seu_path)) return(NA_character_)
  
  seu <- readRDS(seu_path)
  
  # Add cell cycle phase calculation
  seu <- annotate_cell_cycle_without_1q(seu, organism = organism)
  
  # Save back to disk
  saveRDS(seu, seu_path)
  
  seu_path
}

