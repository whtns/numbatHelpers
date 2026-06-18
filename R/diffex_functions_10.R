# Diffex Functions (10)

#' Load or read data from file
#'
#' @param mp_file File path
#' @return Loaded data object
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

read_mps <- function(mp_file = "/dataVolume/storage/Homo_sapiens/3ca/ITH_hallmarks/MPs_distribution/MP_list.RDS") {
  mps <- readRDS(mp_file)

  mps <- purrr::map(mps, ~ purrr::set_names(.x, janitor::make_clean_names(names(.x))))
}

#' Perform differential expression analysis
#'
#' @param oncoprint_input_by_scna_unfiltered Parameter for oncoprint input by scna unfiltered
#' @return List object
#' @export
tally_num_diffex <- function(oncoprint_input_by_scna_unfiltered) {
  symbol_tally <- purrr::map_depth(oncoprint_input_by_scna_unfiltered, 2, ~ {
    table(.x$symbol)
  })
  sample_tally <- purrr::map_depth(oncoprint_input_by_scna_unfiltered, 2, ~ {
    table(.x$sample_id)
  })
  return(list("symbol" = symbol_tally, "sample" = sample_tally))
}

#' Perform clustering analysis
#'
#' @param diffex_1 Parameter for diffex 1
#' @param diffex_2 Parameter for diffex 2
#' @param new_col_name Color specification
#' @return Function result
#' @export
annotate_cluster_membership <- function(diffex_1, diffex_2, new_col_name) {
  diffex_2 <-
    diffex_2 %>%
    dplyr::ungroup() %>%
    dplyr::select(clone_comparison, symbol) %>%
    dplyr::mutate({{ new_col_name }} := TRUE) %>%
    identity()

  diffex_1 <-
    diffex_1 %>%
    dplyr::ungroup() %>%
    dplyr::left_join(diffex_2, by = c("clone_comparison", "symbol"))

  return(diffex_1)
}

#' Perform tally kooi candidates operation
#'
#' @param cis_diffex_clones Character string (default: "results/diffex_bw_clones_large_in_segment_by_chr.xlsx")
#' @param trans_diffex_clones Character string (default: "results/diffex_bw_clones_trans_by_chr.xlsx")
#' @return List object
#' @export
tally_kooi_candidates <- function(cis_diffex_clones = "results/diffex_bw_clones_large_in_segment_by_chr.xlsx", trans_diffex_clones = "results/diffex_bw_clones_trans_by_chr.xlsx") {
  #

  cis_diffex_clones <-
    cis_diffex_clones %>%
    excel_sheets() %>%
    set_names() %>%
    map(read_excel, path = cis_diffex_clones) %>%
    map(dplyr::filter, !is.na(kooi_region))

  cis_diffex_clones <-
    cis_diffex_clones[map_lgl(cis_diffex_clones, ~ (nrow(.x) > 0))] %>%
    dplyr::bind_rows(.id = "chr")

  trans_diffex_clones <-
    trans_diffex_clones %>%
    excel_sheets() %>%
    set_names() %>%
    map(read_excel, path = trans_diffex_clones) %>%
    map(dplyr::filter, !is.na(kooi_region))

  trans_diffex_clones <-
    trans_diffex_clones[map_lgl(trans_diffex_clones, ~ (nrow(.x) > 0))] %>%
    dplyr::bind_rows(.id = "chr")

  return(list("cis" = cis_diffex_clones, "trans" = trans_diffex_clones))
}

#' Perform retrieve cell stats operation
#'
#' @param seu_path File path
#' @param sqlite_path Path to Seurat metadata SQLite database
#' @return Function result
#' @export
retrieve_cell_stats <- function(seu_path, sqlite_path = "batch_hashes.sqlite") {
    if (file.exists(sqlite_path)) {
      con <- connect_hash_db(sqlite_path)
      on.exit(DBI::dbDisconnect(con), add = TRUE)

      if (DBI::dbExistsTable(con, "cell_qc_values")) {
        stats <- DBI::dbGetQuery(
          con,
          "SELECT cell, nCount_gene, nFeature_gene, percent_mt
           FROM cell_qc_values
           WHERE filepath = ?
           ORDER BY cell",
          params = list(seu_path)
        )

        if (nrow(stats) > 0) {
          stats <- stats |>
            dplyr::rename(`percent.mt` = percent_mt)
          return(stats)
        }
      }
    }

    seu <- readRDS(seu_path)

    stats <- seu@meta.data[c("nCount_gene", "nFeature_gene", "percent.mt")] |>
      tibble::rownames_to_column("cell")

    return(stats)
  }

#' Perform collect study metadata operation
#'
#' @return Function result
#' @export
collect_study_metadata <- function() {
  #

  seus <-
    dir_ls("output/seurat/", regexp = "\\/SR[RX][0-9]+_seu.rds") %>%
    set_names(str_extract(., "SR[RX][0-9]+"))

  # collin ------------------------------

  collin_cell_stats <- seus[c("SRX10031191", "SRX10031192", "SRX10031193", "SRX10031194")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  # field ------------------------------

  field_cell_stats <- seus[c("SRX14116948", "SRX14116947", "SRX14116946", "SRX14116945", "SRX14116944")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  # wu ------------------------------

  wu_cell_stats <- seus[c("SRX10264517", "SRX10264518", "SRX10264519", "SRX10264520", "SRX10264521", "SRX10264522", "SRX10264523", "SRX10264524", "SRX10264525", "SRX10264526")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  # yang ------------------------------

  yang_cell_stats <- seus[c("SRX11133594", "SRX11133593", "SRX11133592", "SRX11133591", "SRX11133590", "SRX11133589", "SRX11133588", "SRX11133587", "SRX11133586", "SRX11133585")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  liu_cell_stats <- seus[c("SRX22868105", "SRX22868104", "SRX22868103", "SRX22868102")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  # combined ------------------------------

  study_cell_stats <- dplyr::bind_rows(list("collin" = collin_cell_stats, "field" = field_cell_stats, "wu" = wu_cell_stats, "yang" = yang_cell_stats, "liu" = liu_cell_stats), .id = "study")

  write_csv(study_cell_stats, "results/study_cell_stats.csv")

  return(study_cell_stats)
}

