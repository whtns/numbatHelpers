# Plot Functions (149)

#' Load or read data from file
#'
#' @return Loaded data object
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging
# - rownames_round_trip: Avoid unnecessary rownames conversions - keep as column when possible

read_liu_lu_supp_tables <- function() {
  
#' Load or read data from file
#'
#' @param excel_file File path
#' @param skip Parameter for skip
#' @return Loaded data object
#' @export
myreadxl <- function(excel_file, skip = 0) {
    sheets <- readxl::excel_sheets(excel_file) %>% set_names(.)

    map(sheets, ~ readxl::read_excel(excel_file, .x, skip = skip))
  }

  liu_supp_tables <- "data/liu_lu_supp_data/supp_table_4.xlsx" %>%
    myreadxl(skip = 1) %>%
    map(janitor::clean_names) %>%
    identity()
}
#' Perform set final seus operation
#'
#' @param interesting_samples Parameter for interesting samples
#' @return Function result
#' @export
set_final_seus <- function(interesting_samples) {
  final_seus <- fs::dir_ls("output/seurat/", regexp = ".*SR[RX][0-9]+_filtered_seu.rds") %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  final_seus[names(final_seus) %in% interesting_samples]
}

#' Perform run velocity operation
#'
#' @param debranched_seu Parameter for debranched seu
#' @return Function result
#' @export
run_velocity <- function(debranched_seu) {
  adata <- chevreul::run_scvelo(debranched_seu)
  chevreul::plot_scvelo(adata, mode = "dynamical")
}

#' Perform clustering analysis
#'
#' @param seu Seurat object
#' @param file_id File path
#' @param cluster_orders Cluster information
#' @param resolution Character string (default: "0")
#' @param group.by Character string (default: "SCT_snn_res.0.6")
#' @return Modified Seurat object
#' @export
assign_phase_clusters <- function(seu = NULL, file_id = NULL, cluster_orders = NULL, resolution = "0", group.by = "SCT_snn_res.0.6") {
  phase_levels <- c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star")
  
  if(is.null(file_id)){
      seu$clusters <- factor(as.numeric(seu@meta.data[[group.by]]))
      
      levels(seu$clusters) <- paste0("g1_", levels(seu$clusters))
      
      return(seu)
  }
  
  message(file_id)

  cluster_order <-
    cluster_orders[[file_id]][[resolution]]

  if (!is.null(cluster_order)) {
    group.by <- unique(cluster_order[["resolution"]])
  }

  cluster_order <-
    cluster_order %>%
    dplyr::mutate(order = dplyr::row_number()) %>%
    dplyr::filter(!is.na(clusters)) %>%
    dplyr::mutate(clusters = as.character(clusters)) |> 
  	dplyr::select(-any_of(c("resolution", "tumor_id")))

  seu@meta.data$clusters <- seu@meta.data[[group.by]]

  seu_meta <- seu@meta.data %>%
  	dplyr::select(-starts_with("order")) |> 
    tibble::rownames_to_column("cell") %>%
    dplyr::left_join(cluster_order, by = "clusters") %>%
    dplyr::select(-clusters) %>%
    dplyr::select(-any_of(c("phase_level"))) %>%
    dplyr::rename(phase_level = phase) %>%
  	dplyr::select(-ends_with(".x")) |> 
  	dplyr::select(-ends_with(".y")) |> 
    identity()

  phase_levels <- phase_levels[phase_levels %in% unique(seu_meta$phase_level)]

  seu_meta <-
    seu_meta %>%
    tidyr::unite("clusters", all_of(c("phase_level", group.by)), remove = FALSE) %>%
    dplyr::arrange(phase_level, order) %>%
    dplyr::mutate(clusters = factor(clusters, levels = unique(clusters))) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu@meta.data <- seu_meta[rownames(seu@meta.data), ]

  seu <- find_all_markers(seu, metavar = "clusters", seurat_assay = "SCT")

  seu@meta.data$clusters <- forcats::fct_drop(seu@meta.data$clusters)

  return(seu)
}

