## Helper functions for targets

generate_filtering_cell_counts <- function(filtered_seus, seus_low_hypoxia, filter_inspection_metadata, sqlite_path = "batch_hashes.sqlite", out_csv = "results/filtering_cell_counts.csv") {
  filter_meta <- Filter(
    function(x) inherits(x, "data.frame") && nrow(x) > 0,
    filter_inspection_metadata
  )

  if (length(filter_meta) == 0) {
    base <- tibble::tibble(
      sample_id = character(),
      n_unfiltered = integer(),
      n_annotation_filtered = integer()
    )
  } else {
    base <- dplyr::bind_rows(filter_meta) |>
      dplyr::group_by(sample_id) |>
      dplyr::summarise(
        n_unfiltered = dplyr::n(),
        n_annotation_filtered = if ("filter_keep" %in% names(dplyr::pick(dplyr::everything()))) {
          sum(filter_keep, na.rm = TRUE)
        } else {
          sum(
            !clone_opt_is_na &
              percent.mt < 10 & nCount_gene > 1000 &
              nFeature_gene > 1000 &
              !cluster_remove_flag & !is_malat1 & !in_manual_exclude
          )
        },
        .groups = "drop"
      )
  }

  filtered_paths  <- na.omit(unlist(filtered_seus))
  lh_paths        <- na.omit(unlist(seus_low_hypoxia))

  # Fallback: for samples in base not covered by pipeline targets, check DB directly
  covered_filtered <- stringr::str_extract(filtered_paths, "SR[RX][0-9]+")
  missing_filtered <- setdiff(base$sample_id, covered_filtered)
  missing_filtered <- missing_filtered[!is.na(missing_filtered)]
  if (length(missing_filtered) > 0) {
    fb <- file.path("output/seurat", paste0(missing_filtered, "_filtered_seu.rds"))
    filtered_paths <- c(filtered_paths, fb[file.exists(fb)])
  }

  covered_lh <- stringr::str_extract(lh_paths, "SR[RX][0-9]+")
  missing_lh <- setdiff(base$sample_id, covered_lh)
  missing_lh <- missing_lh[!is.na(missing_lh)]
  if (length(missing_lh) > 0) {
    fb <- file.path("output/seurat", paste0(missing_lh, "_hypoxia_low_seu.rds"))
    lh_paths <- c(lh_paths, fb[file.exists(fb)])
  }

  all_paths       <- c(filtered_paths, lh_paths)

  db_rows <- if (length(all_paths) == 0) {
    tibble::tibble(filepath = character(), n_cells = integer())
  } else {
    con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    placeholders <- paste(rep("?", length(all_paths)), collapse = ", ")
    DBI::dbGetQuery(
      con,
      glue::glue(
        "SELECT filepath, n_cells FROM hashes ",
        "WHERE filepath IN ({placeholders})"
      ),
      params = as.list(unname(all_paths))
    )
  }

  extract_counts <- function(paths, col) {
    out <- db_rows[db_rows$filepath %in% paths, ] |>
      dplyr::mutate(
        sample_id = stringr::str_extract(filepath, "SR[RX][0-9]+")
      ) |>
      dplyr::select(sample_id, n_cells)
    names(out)[names(out) == "n_cells"] <- col
    out
  }

  result <- base |>
    dplyr::left_join(
      extract_counts(filtered_paths, "n_pipeline_filtered"),
      by = "sample_id"
    ) |>
    dplyr::left_join(
      extract_counts(lh_paths, "n_low_hypoxia"),
      by = "sample_id"
    ) |>
    dplyr::mutate(
      pct_annotation_filtered =
        round(n_annotation_filtered / n_unfiltered * 100, 1),
      pct_pipeline_filtered =
        round(n_pipeline_filtered / n_unfiltered * 100, 1),
      pct_low_hypoxia =
        round(n_low_hypoxia / n_unfiltered * 100, 1)
    ) |>
    dplyr::arrange(sample_id)

  fs::dir_create(dirname(out_csv))
  readr::write_csv(result, out_csv)
  out_csv
}
