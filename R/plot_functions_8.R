# Plot Functions (100)

#' Perform all same sign operation
#'
#' @param x Parameter for x
#' @return Function result
#' @export
# Performance optimizations applied:
# - rownames_round_trip: Avoid unnecessary rownames conversions - keep as column when possible

all_same_sign <- function(x) {
  OR <- `||`
  AND <- `&&`

  OR(length(x) <= 1L, {
    if (anyNA(x1 <- x[1L])) {
      return(NA)
    }
    if (x1 == 0) {
      AND(
        min(x) == 0,
        max(x) == 0
      )
    } else if (x1 > 0) {
      min(x) > 0
    } else {
      max(x) < 0
    }
  })
}

#' Perform interleave lists operation
#'
#' @param list1 Parameter for list1
#' @param list2 Parameter for list2
#' @return List object
#' @export
interleave_lists <- function(list1, list2) {
  max_length <- max(length(list1), length(list2))
  interleaved <- as.list(c(mapply(c, list1, list2, SIMPLIFY = FALSE)))
  return(c(interleaved, list1[length(list1):max_length], list2[length(list2):max_length]))
}

#' Perform clustering analysis
#'
#' @param seu_path File path
#' @param cluster_orders Cluster information
#' @param resolution_dictionary Parameter for resolution dictionary
#' @return Function result
#' @export
assign_designated_phase_clusters <- function(seu_path, cluster_orders, resolution_dictionary) {
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  file_id <- fs::path_file(seu_path)

  chosen_resolution <- resolution_dictionary[[file_id]]

  assign_phase_cluster_at_resolution(seu_path, cluster_order = cluster_orders, assay = "SCT", resolution = chosen_resolution)
}

#' Create a heatmap visualization
#'
#' @param seu_path File path
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @param scna_of_interest Parameter for scna of interest
#' @param ... Additional arguments passed to other functions
#' @return Data frame
#' @export
plot_seu_gene_heatmap <- function(seu_path, large_clone_comparisons, scna_of_interest, ...) {
  pdf_path <- fs::path("results", str_replace(fs::path_file(seu_path), ".rds", "_heatmap.pdf"))

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  seu <- readRDS(seu_path)

  # subset by retained clones ------------------------------
  clone_comparisons <- names(large_clone_comparisons[[sample_id]])
  clone_comparison <- clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
  retained_clones <- clone_comparison %>%
    str_extract("[0-9]_v_[0-9]") %>%
    str_split("_v_", simplify = TRUE)

  seu <- seu[, seu$clone_opt %in% retained_clones]

  usages <- dir_ls(glue("output/mosaicmpi/{sample_id}/"), glob = "*usage_k*.txt") %>%
    set_names(str_extract_all(., "(?<=k)[0-9]*")) %>%
    map(read_tsv) %>%
    map(tibble::column_to_rownames, "...1") %>%
    imap(~ set_names(.x, paste0(.y, "_", colnames(.x)))) %>%
    identity()

  usages <-
    usages[as.character(sort(as.numeric(names(usages))))]

  colnames(usages[["6"]]) <- paste0("factor.", 1:6)

  seu0 <- AddMetaData(seu, bind_cols(usages))

  seu0$sample <- sample_id

  cluster_plot <- seu_gene_heatmap(seu0,
    marker_col = "clusters",
    group.by = c("clusters", "scna", "S.Score", "G2M.Score", colnames(usages[["6"]])),
    col_arrangement = c("clusters", "scna"),
    column_split = "clusters",
    hide_legends = colnames(usages[["6"]])
  ) + labs(title = sample_id, subtitle = "6")
  print(cluster_plot)

  ggsave(pdf_path, cluster_plot, ...)

  return(pdf_path)
}

#' Merge or join datasets
#'
#' @param merged_seu_path File path
#' @param child_seu_paths File path
#' @return Modified Seurat object
#' @export
check_merged_metadata <- function(merged_seu_path, child_seu_paths) {
  child_seu_meta <- child_seu_paths |>
    map(
      ~ {
        readRDS(.x)@meta.data |>
          dplyr::select(-any_of(c("cell"))) |>
          tibble::rownames_to_column("cell")
      }
    ) |>
    bind_rows(.id = "sample") |>
    dplyr::select(sample, cell, clusters)

  merged_seu <- readRDS(merged_seu_path)

  merged_seu_meta <-
    merged_seu@meta.data |>
    tibble::rownames_to_column("cell") |>
    dplyr::mutate(sample_id = str_remove(sample_id, "_.*")) |>
    dplyr::mutate(truncated_cell = str_remove(cell, "_[0-9]*")) |>
    dplyr::select(-any_of("clusters")) |>
    dplyr::left_join(child_seu_meta, by = c("sample_id" = "sample", "truncated_cell" = "cell")) |>
    dplyr::select(cell, clusters) |>
    tibble::column_to_rownames("cell") |>
    # tibble::deframe() |>
    # dim() |>
    identity()

  merged_seu$clusters <- NULL

  merged_seu <- Seurat::AddMetaData(merged_seu, merged_seu_meta)

  return(merged_seu)
}

