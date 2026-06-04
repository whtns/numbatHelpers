# Numbat Functions (1)

#' Extract or pull specific data elements
#'
#' @param clustrees Parameter for clustrees
#' @param divergent_cluster_file File path
#' @return Extracted data elements
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

pull_clustree_tables <- function(clustrees, divergent_cluster_file = "data/clustree_divergent_clusters.csv") {
  clustrees <- map(clustrees, purrr::compact)

  divergent_clusters <- divergent_cluster_file %>%
    read_csv() %>%
    identity()

  clustree_tables <-
    clustrees %>%
    unlist() %>%
    str_replace(".pdf", ".xlsx") %>%
    set_names(str_remove(fs::path_file(.), "_x.*")) %>%
    set_names(str_remove(names(.), "_diploid.*")) %>%
    map(myreadxl) %>%
    map(bind_rows, .id = "clone") %>%
    bind_rows(.id = "sample_id") %>%
    dplyr::mutate(from_SCT_snn_res. = as.double(SCT_snn_res.), from_clust = cluster) %>%
    dplyr::mutate(from_clust = as.integer(str_extract(node, "[0-9]$"))) %>%
    dplyr::select(-to_clust) %>%
    dplyr::inner_join(divergent_clusters, by = c("sample_id", "to_SCT_snn_res.", "from_SCT_snn_res.", "from_clust")) %>%
    dplyr::filter(!is.na(from_clust)) %>%
    dplyr::select(to_clust, everything()) %>%
    dplyr::arrange(sample_id, to_clust) %>%
    dplyr::select(-clone_comparison) %>%
    identity()


  clustree_tables <-
    clustree_tables %>%
    tidyr::pivot_longer(contains("_v_"), names_to = "clone_comparison", values_to = "p.value") %>%
    dplyr::filter(!is.na(p.value)) %>%
    dplyr::filter(!is.na(from_SCT_snn_res.)) %>%
    dplyr::mutate(clone = str_remove(clone, "_[0-9]$")) %>%
    dplyr::select(-clone) %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::select(-starts_with("x")) %>%
    dplyr::select(sample_id, clone_comparison, from_SCT_snn_res., from_clust, to_SCT_snn_res., to_clust, p.value, everything()) %>%
    dplyr::arrange(sample_id, clone_comparison) %>%
    split(.$sample_id) %>%
    identity()

  return(clustree_tables)
}

#' Extract or pull specific data elements
#'
#' @param branch_dictionary_file File path
#' @return Extracted data elements
#' @export
pull_branches <- function(branch_dictionary_file = "data/branch_dictionary.csv") {
  branch_dictionary <- read_csv(branch_dictionary_file) %>%
    dplyr::mutate(branch_members = str_split(branch_members, "_")) %>%
    group_by(sample_id) %>%
    tidyr::nest(branch_termination, branch_members) %>%
    # group_split() %>%
    # dplyr::select(sample_id, branch_members) %>%
    tibble::deframe() %>%
    map(deframe) %>%
    identity()

  # branch_dictionary <- branch_dictionary[map(branch_dictionary, length) > 1]

  return(branch_dictionary)
}

#' Perform debranch seus operation
#'
#' @param filtered_seus Parameter for filtered seus
#' @param branch_dictionary Parameter for branch dictionary
#' @param ... Additional arguments passed to other functions
#' @return Modified Seurat object
#' @export
debranch_seus <- function(filtered_seus, branch_dictionary, ...) {
  filtered_seus <- filtered_seus %>%
    purrr::set_names(str_extract(., "SR[RX][0-9]+"))

  debranched_seus <- map(filtered_seus, split_seu_by_branch, branch_dictionary, ...)

  debranched_seus <- unlist(debranched_seus)

  return(debranched_seus)
}

#' Perform prep seu branch operation
#'
#' @param debranched_seu Parameter for debranched seu
#' @param ... Additional arguments passed to other functions
#' @return Modified Seurat object
#' @export
prep_seu_branch <- function(debranched_seu, ...) {
  DefaultAssay(debranched_seu) <- "gene"

  debranched_seu <-
    debranched_seu %>%
    seurat_preprocess() %>%
    seurat_reduce_dimensions()

  debranched_seu <- FindNeighbors(debranched_seu, dims = 1:30, verbose = FALSE)

  debranched_seu <- seurat_cluster(
    seu = debranched_seu, resolution = seq(0.2, 1.0, by = 0.2),
    reduction = "pca"
  )

  debranched_seu <- find_all_markers(debranched_seu, seurat_assay = "gene")

  DefaultAssay(debranched_seu) <- "SCT"

  debranched_seu <- seurat_cluster(
    seu = debranched_seu, resolution = seq(0.2, 1.0, by = 0.2),
    reduction = "pca"
  )

  debranched_seu <- find_all_markers(debranched_seu, seurat_assay = "SCT")

  debranched_seu <- assign_phase_clusters(debranched_seu, ...)

  return(debranched_seu)
}

#' Perform split seu by branch operation
#'
#' @param seu_path File path
#' @param branches Parameter for branches
#' @param ... Additional arguments passed to other functions
#' @return Modified Seurat object
#' @export
split_seu_by_branch <- function(seu_path, branches, ...) {
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  branches <- branches[[tumor_id]]

  if (length(branches) == 1) {
    seu_paths <- seu_path
  } else {
    seu <- readRDS(seu_path)
    seu_paths <- list()
    for (terminal_clone in names(branches)) {
      branch_path <- str_replace(seu_path, "_filtered_seu.*", glue("_branch_{terminal_clone}_filtered_seu.rds"))
      debranched_seu <- seu[, seu$clone_opt %in% branches[[terminal_clone]]]
      sample_id <- glue("{tumor_id}_branch_{terminal_clone}")

      # Cache check: if branch file exists with the same cells, skip the expensive
      # prep_seu_branch (SCTransform x2 + PCA + UMAP + clustering x2 + markers x2).
      if (file.exists(branch_path)) {
        cached_seu <- tryCatch(readRDS(branch_path), error = function(e) NULL)
        if (!is.null(cached_seu) && identical(sort(colnames(cached_seu)), sort(colnames(debranched_seu)))) {
          message("Cell barcodes unchanged for ", sample_id, "; skipping prep_seu_branch")
          seu_paths[[terminal_clone]] <- branch_path
          next
        }
        message("Cell barcodes changed for ", sample_id, "; running full prep_seu_branch")
      }

      debranched_seu <- prep_seu_branch(debranched_seu, tumor_id = tumor_id, sample_id = sample_id, ...)
      seu_paths[[terminal_clone]] <- branch_path
      add_hash_metadata(seu = debranched_seu, filepath = seu_paths[[terminal_clone]])
    }
  }

  return(seu_paths)
}

