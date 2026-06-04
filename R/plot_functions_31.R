# Plot Functions (130)

#' Perform differential expression analysis
#'
#' @param seu_path File path
#' @param corresponding_states_dictionary Parameter for corresponding states dictionary
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @param numbat_rds_files File path
#' @param location Character string (default: "cis")
#' @param scna_of_interest Character string (default: "2p")
#' @param w_scna Parameter for w scna
#' @param wo_scna Parameter for wo scna
#' @return Differential expression results
#' @export
# Performance optimizations applied:
# - repeated_file_reads: Cache file reads to avoid redundant I/O

find_diffex_clones_between_corresponding_states <- function(seu_path, corresponding_states_dictionary, large_clone_comparisons = large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "cis", scna_of_interest = "2p", w_scna = NULL, wo_scna = NULL) {
  
  
  
  
  #
  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  numbat_rds_files <- numbat_rds_files %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  mynb <- numbat_rds_files[[tumor_id]]

  seu <- readRDS(seu_path)

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  large_clone_comparisons <- large_clone_comparisons[[sample_id]]

  large_clone_comparisons <- large_clone_comparisons[str_detect(names(large_clone_comparisons), scna_of_interest)]

  corresponding_states <-
    corresponding_states_dictionary |>
    dplyr::mutate(sample = str_remove(file_name, "_filtered_seu.*$")) |>
    dplyr::filter(sample == sample_id) |>
    dplyr::select(w_scna, wo_scna) |>
    dplyr::mutate(comparison = glue("{w_scna} v. {wo_scna}"))
  identity()

  possible_cluster_comparison <- possibly(make_cluster_comparison)
  diffex <- pmap(corresponding_states, ~ possible_cluster_comparison(large_clone_comparisons, names(large_clone_comparisons), seu, mynb, .x, .y))

  diffex <- diffex[str_detect(names(diffex), scna_of_interest)]

  names(diffex) <-
    corresponding_states$comparison

  return(diffex)
}

#' Perform differential expression analysis
#'
#' @param seu_path File path
#' @param corresponding_states_dictionary Parameter for corresponding states dictionary
#' @param large_clone_comparisons Parameter for large clone comparisons
#' @param numbat_rds_files File path
#' @param location Character string (default: "cis")
#' @return Differential expression results
#' @export
find_diffex_clusters_between_corresponding_states <- function(seu_path, corresponding_states_dictionary, large_clone_comparisons = large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "cis") {
  
  file_name <- fs::path_ext_remove(fs::path_file(seu_path))

  tumor_id <- str_extract(seu_path, "SR[RX][0-9]+")

  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.*")

  numbat_rds_files <- numbat_rds_files %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  seu <- readRDS(seu_path)
  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  corresponding_states <-
  	corresponding_states_dictionary[[1]] |>
  	dplyr::mutate(sample = str_remove(file_name, "_filtered_seu.*$")) |>
  	dplyr::filter(file_name == fs::path_file(seu_path)) |>
  	dplyr::select(w_scna, wo_scna, scna_of_interest) |>
  	dplyr::mutate(comparison = glue("{w_scna} v. {wo_scna}")) |>
  	identity()
  
  scna_of_interest <- unique(corresponding_states$scna_of_interest)
  
  if("integrated" %in% names(seu@assays)){
  	sample_ids <- unique(as.character(seu$batch))
  	mynb <- numbat_rds_files[sample_ids] |> 
  		map(readRDS)
  	
  	large_clone_comparisons <- 
  		large_clone_comparisons[sample_ids] |> 
  		map(~.x[str_detect(names(.x), scna_of_interest)]) |> 
  		purrr::list_flatten()
  	
  } else {
  	mynb <- numbat_rds_files[[tumor_id]]	
  	
  	large_clone_comparisons <- large_clone_comparisons[[sample_id]]
  	
  	large_clone_comparisons <- large_clone_comparisons[str_detect(names(large_clone_comparisons), scna_of_interest)]
  	
  }
  
  possible_cluster_comparison <- possibly(make_cluster_comparison)
  diffex <- pmap(corresponding_states, ~ possible_cluster_comparison(large_clone_comparisons, names(large_clone_comparisons), seu, mynb, .x, .y))

  names(diffex) <-
    corresponding_states$comparison
  
  tables_path <- diffex |>
    purrr::list_flatten() |>
    purrr::compact() |>
    writexl::write_xlsx(glue("results/corresponding_states_diffex_{file_name}.xlsx"))

  return(diffex)
}

#' Title
#' @export
compare_enrichment <- function(diffex, ...) {
  #
  enrich_diffex <- function(df) {
    df |>
      tibble::column_to_rownames("symbol") |>
      dplyr::select(-any_of(colnames(annotables::grch38))) |>
      enrichment_analysis(...)
  }

  enrichments <- diffex |>
    purrr::list_flatten() |>
    map(enrich_diffex)

  return(enrichments)
}

#' Perform differential expression analysis
#'
#' @param seu_path File path
#' @param kept_samples Parameter for kept samples
#' @param clone_comparisons Parameter for clone comparisons
#' @param location Character string (default: "cis")
#' @return Differential expression results
#' @export
find_diffex_clones_integrated <- function(seu_path, kept_samples, clone_comparisons, location = "cis") {
  
  
  
  
  #

  seu <- seu_path

  seu <- seu[, !is.na(seu$clone_opt)]

  mynbs <- glue("output/numbat_sridhar/{kept_samples}_numbat.rds") %>%
    map(readRDS)


  diffex <- imap(clone_comparisons, make_clone_comparison_integrated, seu, mynbs, location = location)


  return(diffex)
}

