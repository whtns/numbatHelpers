
load_and_save_hypoxia_score <- function(seu_path) {
    if (is.na(seu_path)) return(NA_character_)

    # determine target path first and skip work if file already exists

    hypoxia_path <- sub("^([^_]+)_.*$", "\\1_seu_hypoxia.rds", seu_path)

    # if (file.exists(hypoxia_path)) {
    #     message("Hypoxia file exists, skipping: ", hypoxia_path)
    #     return(hypoxia_path)
    # }

    seu <- readRDS(seu_path)

    if (!all(c("G2M.Score", "S.Score") %in% colnames(seu@meta.data))) {
        message("CC scores missing in ", basename(seu_path), "; running CellCycleScoring.")
        seu <- tryCatch(
            Seurat::CellCycleScoring(seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE),
            error = function(e) {
                message("CellCycleScoring failed (", e$message, "); setting CC scores to 0.")
                seu$S.Score <- 0
                seu$G2M.Score <- 0
                seu$Phase <- "G1"
                return(seu)
            }
        )
    }

    seu <- add_hypoxia_score(seu)

    add_hash_metadata(seu = seu, filepath = hypoxia_path)

    return(hypoxia_path)

}

#' Compute and save hypoxia score for a Seurat RDS file
#'
#' Given the path to a Seurat RDS file, compute the hypoxia score by calling
#' `add_hypoxia_score()` and save the modified Seurat object to a new RDS. If
#' the target hypoxia RDS already exists the function returns early and does
#' not re-compute.
#'
#' @param seu_path Path to an existing Seurat RDS file (e.g. "SRXxxxxx_seu.rds").
#' @return Path to the hypoxia RDS file that was written (or existing file path
#'   if skipped).
#' @export

subset_to_low_hypoxia <- function(seu_path, threshold = 0.5, slug="", ...) {
    seu <- readRDS(seu_path)

    seu_low_hypoxia <- seu[,seu$hypoxia_score < threshold] |> 
    # seurat_cluster(resolution = seq(0.2, 1, by = 0.2))  |> 
    identity()

    # seu_low_hypoxia <- assign_phase_clusters(seu_low_hypoxia, file_id, ...)

    low_hypoxia_path <-  stringr::str_replace(seu_path, "_seu_hypoxia.rds", "_seu_low_hypoxia.rds")

    add_hash_metadata(seu = seu_low_hypoxia, filepath = low_hypoxia_path)

    return(low_hypoxia_path)
}

#' Subset Seurat object to low-hypoxia cells and save
#'
#' Read a hypoxia-seurat RDS file, subset to cells with `hypoxia_score < threshold`,
#' run clustering, and save the resulting Seurat object to disk.
#'
#' @param seu_path Path to a Seurat RDS file containing `hypoxia_score`.
#' @param threshold Numeric threshold for low hypoxia (default 0.5).
#' @param slug Optional slug appended to output filenames.
#' @param ... Additional arguments passed to downstream clustering/assignment functions.
#' @return Path to the saved low-hypoxia Seurat RDS file.
#' @export

subset_to_high_hypoxia <- function(seu_path, threshold = 0.5, slug="", ...) {
    seu <- readRDS(seu_path)

    seu_high_hypoxia <- seu[,seu$hypoxia_score > threshold] |> 
    # seurat_cluster(resolution = seq(0.2, 1, by = 0.2))  |> 
    identity()

    # seu_high_hypoxia <- assign_phase_clusters(seu_high_hypoxia, file_id, ...)

    high_hypoxia_path <-  stringr::str_replace(seu_path, "_seu_hypoxia.rds", "_seu_high_hypoxia.rds")

    add_hash_metadata(seu = seu_high_hypoxia, filepath = high_hypoxia_path)

    return(high_hypoxia_path)
}

#' Subset Seurat object to high-hypoxia cells and save
#'
#' Read a hypoxia-seurat RDS file, subset to cells with `hypoxia_score > threshold`,
#' run clustering, and save the resulting Seurat object to disk.
#'
#' @param seu_path Path to a Seurat RDS file containing `hypoxia_score`.
#' @param threshold Numeric threshold for high hypoxia (default 0.5).
#' @param slug Optional slug appended to output filenames.
#' @param ... Additional arguments passed to downstream clustering/assignment functions.
#' @return Path to the saved high-hypoxia Seurat RDS file.
#' @export


add_hypoxia_score <- function(seu) {
    # mt_genes <- str_subset(rownames(seu), "MT-.*")
    mt_genes <- c("MT-CO3", "MT-ND3", "MT-CYB", "MT-ATP6", "MT-CO2")
    
    # hypoxia_genes <-
    #     dplyr::filter(seu@misc$markers$clusters$presto, Cluster == "hypoxia_2") |> 
    #     slice_head(n = 100) |>
    #     dplyr::pull(Gene.Name) |>
    #     identity()
    
    # hypoxia_genes <- c("BNIP3", "GAS5")
    
    msig_db_human <- msigdbr::msigdbr(species = "Homo sapiens")

    hypoxia_genes <- msig_db_human %>%
    dplyr::filter(gs_name == "HALLMARK_HYPOXIA")  |> 
    dplyr::pull(gene_symbol)

    # hypoxia_genes <- c("BNIP3", "GAS5", "EPB41L4A-AS1")

    nbin <- max(24L, floor(length(VariableFeatures(seu)) / 100))

    seu <- tryCatch(
      Seurat::AddModuleScore(seu, features = list("hypoxia" = hypoxia_genes, "MT" = mt_genes), name = "hypoxia", nbin = nbin, ctrl = 100),
      error = function(e) {
        message("AddModuleScore failed with nbin=", nbin, ", retrying with nbin=5: ", e$message)
        tryCatch(
          Seurat::AddModuleScore(seu, features = list("hypoxia" = hypoxia_genes, "MT" = mt_genes), name = "hypoxia", nbin = 5L, ctrl = 5L),
          error = function(e2) {
            message("AddModuleScore failed with nbin=5, retrying with nbin=2: ", e2$message)
            Seurat::AddModuleScore(seu, features = list("hypoxia" = hypoxia_genes, "MT" = mt_genes), name = "hypoxia", nbin = 2L, ctrl = 2L)
          }
        )
      }
    )
    
    seu$hypoxia <- seu$hypoxia1
    seu$MT <- seu$hypoxia2
    
    seu$MT <- seu$MT*-1
    
    seu$hypoxia_score <-
        rowMeans(seu@meta.data[c("hypoxia", "MT")])
    
    seu$hypoxia_score = scales::rescale(seu$hypoxia_score, c(0,1))
    
    return(seu)
}

#' Calculate hypoxia score and add to Seurat metadata
#'
#' Compute a hypoxia module score using MSigDB HALLMARK_HYPOXIA genes and a small
#' mitochondrial gene set, then store a rescaled `hypoxia_score` in `seu@meta.data`.
#'
#' @param seu A Seurat object. Variable features should be set.
#' @return The input Seurat object with `hypoxia`, `MT`, and `hypoxia_score` added to `@meta.data`.
#' @export

plot_hypoxia_gene_heatmap <- function(seu_path, group.by = "gene_snn_res.0.2", n_genes = 50) {
  if (is.na(seu_path)) return(NA_character_)

  sample_id <- stringr::str_extract(seu_path, "SR[RX][0-9]+")
  out_dir <- "results/hypoxia_gene_heatmaps"
  fs::dir_create(out_dir)
  out_path <- glue::glue("{out_dir}/{sample_id}_hypoxia_gene_heatmap.pdf")

  seu <- readRDS(seu_path)
  Seurat::DefaultAssay(seu) <- "gene"

  # Same gene set used in add_hypoxia_score
  msig_db_human <- msigdbr::msigdbr(species = "Homo sapiens")
  hypoxia_genes <- msig_db_human %>%
    dplyr::filter(gs_name == "HALLMARK_HYPOXIA") %>%
    dplyr::pull(gene_symbol)

  available_genes <- intersect(hypoxia_genes, rownames(seu))
  if (length(available_genes) == 0) {
    warning("No HALLMARK_HYPOXIA genes found in: ", seu_path)
    return(NA_character_)
  }

  # Rank by cross-cell variance; keep top n_genes
  if (length(available_genes) > n_genes) {
    gene_mat <- Seurat::GetAssayData(seu, assay = "gene", layer = "data")[available_genes, , drop = FALSE]
    gene_vars <- apply(gene_mat, 1, var)
    available_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(n_genes)]
  }

  # Fall back to a valid grouping column if requested one is absent
  if (!group.by %in% colnames(seu@meta.data)) {
    group.by <- grep("_snn_res\\.", colnames(seu@meta.data), value = TRUE)[1]
  }

  seu <- Seurat::ScaleData(seu, assay = "gene", features = available_genes, verbose = FALSE)

  p <- Seurat::DoHeatmap(seu, features = available_genes, group.by = group.by,
                         assay = "gene", layer = "scale.data") +
    ggplot2::ggtitle(glue::glue("{sample_id}: HALLMARK_HYPOXIA (low-hypoxia cells)")) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                   axis.text.y  = ggplot2::element_text(size = 6))

  ggplot2::ggsave(out_path, p,
                  width  = 14,
                  height = max(6, length(available_genes) * 0.18),
                  limitsize = FALSE)
  out_path
}

subset_to_1q <- function(seu, file_id = NULL, slug="", ...) {
  clone_selection <- tribble(
      ~batch, ~clone_opt,
      "SRX10264526", c(1,2),
      "SRX11133594", c(2,3),
      "SRX11133593", c(2,3),
      "SRX11133592", c(2,3),
      "SRX10264523", c(1,2),
      "SRX14116944", c(1,2),
      "SRX10831287", c(1,2)
  ) |>
      tidyr::unnest(clone_opt)
  
  selected_cells <- seu@meta.data |>
      tibble::rownames_to_column("cell") |>
      dplyr::inner_join(clone_selection, by = c("batch", "clone_opt"))
  
  seu_1q <- seu[,selected_cells$cell]
  
  seu_1q$scna <-
      factor(ifelse(str_detect(seu_1q$GT_opt, pattern = "1[a-z]"), "w_scna", "wo_scna"), levels = c("wo_scna", "w_scna"))
  
  seu_1q <- assign_phase_clusters(seu_1q, file_id, ...)
  
  return(seu_1q)
}

#' Subset a Seurat object to 1q clones
#'
#' Select cells from the Seurat object that belong to predefined 1q clone
#' batches. This function expects `seu@meta.data` to contain `batch` and
#' `clone_opt` fields.
#'
#' @param seu A Seurat object.
#' @param file_id Optional file/sample identifier passed to downstream functions.
#' @param slug Optional slug string.
#' @param ... Additional arguments forwarded to `assign_phase_clusters()`.
#' @return A Seurat object filtered to the selected 1q clone cells.
#' @export

select_1q_clones <- function(seu, slug="", ...) {
    
    seu_1q <- subset_to_1q(seu, ...)
    
    rds_path <-  glue("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6{slug}.rds")
    

    add_hash_metadata(seu = seu_1q, filepath = rds_path)
    
    pdf_path <- make_clone_distribution_figure_debug(seu_1q, rds_path, group.bys = "clusters", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "hypoxia_score", "batch"), ...)
    
    return(pdf_path)
    
}

#' Select 1q clone Seurat objects and save
#'
#' Wrapper that subsets the provided Seurat object to 1q clones, saves the
#' resulting Seurat RDS and creates a clone distribution figure.
#'
#' @param seu A Seurat object.
#' @param slug Optional slug appended to output filenames.
#' @param ... Additional arguments passed to downstream plotting functions.
#' @return Path to the generated PDF clone distribution figure.
#' @export

subset_to_16q <- function(seu, file_id = NULL, slug="", ...) {
    clone_selection <- tribble(
        ~batch, ~clone_opt,
        "SRX11133594", c(1,2),
        "SRX11133593", c(1,2),
        "SRX11133592", c(1,2)
    ) |>
        tidyr::unnest(clone_opt)
    
    selected_cells <- seu@meta.data |>
        tibble::rownames_to_column("cell") |>
        dplyr::inner_join(clone_selection, by = c("batch", "clone_opt"))
    
    seu_16q <- seu[,selected_cells$cell]
    
    seu_16q$scna <-
        factor(ifelse(str_detect(seu_16q$GT_opt, pattern = "16[a-z]"), "w_scna", "wo_scna"), levels = c("wo_scna", "w_scna"))
    
    seu_16q <- assign_phase_clusters(seu_16q, file_id, ...)
    
    return(seu_16q)
}

#' Subset a Seurat object to 16q clones
#'
#' Select cells from the Seurat object that belong to predefined 16q clone
#' batches. This function expects `seu@meta.data` to contain `batch` and
#' `clone_opt` fields.
#'
#' @param seu A Seurat object.
#' @param file_id Optional file/sample identifier.
#' @param slug Optional slug string.
#' @param ... Additional arguments forwarded to `assign_phase_clusters()`.
#' @return A Seurat object filtered to the selected 16q clone cells.
#' @export

select_16q_clones <- function(seu, file_id = NULL, slug="", ...) {
    
    seu_16q <- subset_to_16q(seu, ...)
    
    rds_path <-  glue("output/seurat/integrated_1q_16q/integrated_seu_16q_afterall6{slug}.rds")
    add_hash_metadata(seu = seu_16q, filepath = rds_path)

    pdf_path <- make_clone_distribution_figure_debug(seu_16q, rds_path, group.bys = "clusters", heatmap_groups = c("G2M.Score", "S.Score", "scna", "clusters", "hypoxia_score", "batch"), file_id = file_id, ...)
    
    return(pdf_path)
    
}

#' Select 16q clone Seurat objects and save
#'
#' Wrapper that subsets the provided Seurat object to 16q clones, saves the
#' resulting Seurat RDS and creates a clone distribution figure.
#'
#' @param seu A Seurat object.
#' @param file_id Optional file/sample identifier.
#' @param slug Optional slug appended to output filenames.
#' @param ... Additional arguments passed to downstream plotting functions.
#' @return Path to the generated PDF clone distribution figure.
#' @export

check_integrated_cluster_numbers <- function(seu_path){
    
    seu <- readRDS(seu_path)
    
    n_clusters <- seu@meta.data[c("integrated_snn_res.0.4", "integrated_snn_res.0.6", "integrated_snn_res.0.8")] |> 
        map(dplyr::n_distinct)
    
    return(n_clusters)
}

#' Check number of integrated clusters at several resolutions
#'
#' Read a Seurat RDS and return the number of unique clusters for a set of
#' integrated resolutions (0.4, 0.6, 0.8).
#'
#' @param seu_path Path to a Seurat RDS file.
#' @return A named list (or vector) with distinct cluster counts for each resolution.
#' @export

