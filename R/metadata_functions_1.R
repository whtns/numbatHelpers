# Parameter and metadata functions (11)

#' Add batch hash metadata column
#'
#' @param filepath Path to save the Seurat object RDS; auto-generated from hash if NULL
#' @param seu Seurat object with a 'batch' metadata column
#' @param sqlite_path Path to the SQLite database file
#' @return Character string: the filepath where the Seurat object was saved
#' @export
add_batch_hash_metadata <- function(filepath = NULL, seu = NULL, sqlite_path = "batch_hashes.sqlite") {
  if (is.null(seu)) {
    seu <- readRDS(filepath)
  }
  # Split by batch
  seu_list <- Seurat::SplitObject(seu, split.by = "batch")

  # Compute hash for each batch
  batch_hashes <- seu_list  |>
  map(colnames)  |>
  map(str_remove, "_.*")  |>
  purrr::map_chr(digest::digest)

  # Create a lookup data frame
  hash_lookup <- data.frame(
    batch = names(batch_hashes),
    batch_hash = base::unname(batch_hashes),
    stringsAsFactors = FALSE
  )

  # Add hash to metadata by matching batch
  seu$batch_hash <- hash_lookup$batch_hash[match(seu$batch, hash_lookup$batch)]
  # add integration hash
  hash <- digest::digest(colnames(seu))
  seu$hash <- hash

  if (is.null(filepath)){
    filepath <- glue::glue("output/seurat/{hash}_seu.rds")
  }

  # Write to sqlite
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA busy_timeout = 30000")
  # Create table if not exists
  DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS hashes (filepath TEXT PRIMARY KEY, hash TEXT)")
  # Insert or replace
  DBI::dbExecute(con, "INSERT OR REPLACE INTO hashes (filepath, hash) VALUES (?, ?)", params = list(filepath, hash))
  saveRDS(seu, filepath)
  return(filepath)
}

add_hash_metadata <- function(filepath = NULL, seu = NULL, sqlite_path = "batch_hashes.sqlite"){
  if (is.null(seu)) {
    seu <- readRDS(filepath)
  }
  hash <- digest::digest(colnames(seu))
  seu$hash <- hash

  if (is.null(filepath)){
    filepath <- glue::glue("output/seurat/{hash}_seu.rds")
  }

  n_cells <- ncol(seu)

  # Write to sqlite
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA busy_timeout = 30000")
  DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS hashes (filepath TEXT PRIMARY KEY, hash TEXT)")
  # Add n_cells column if the table predates this change (idempotent)
  tryCatch(
    DBI::dbExecute(con, "ALTER TABLE hashes ADD COLUMN n_cells INTEGER"),
    error = function(e) invisible(NULL)
  )
  DBI::dbExecute(con,
    "INSERT OR REPLACE INTO hashes (filepath, hash, n_cells) VALUES (?, ?, ?)",
    params = list(filepath, hash, n_cells)
  )
  saveRDS(seu, filepath)
  return(filepath)
}

read_seu_hash <- function(sqlite_path = "batch_hashes.sqlite", filepath = NULL) {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  if (!is.null(filepath)) {
    # Query for specific filepath
    hashes_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes WHERE filepath = ?", params = list(filepath))
  } else {
    # Query all records
    hashes_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes")
  }

  return(hashes_df)
}

make_filepaths_unique_in_hashes_table <- function(){
  con <- DBI::dbConnect(RSQLite::SQLite(), "batch_hashes.sqlite")
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  hashes_df <-
  tbl(con, "hashes")  |>
  dplyr::distinct(hash, .keep_all = TRUE)  |>
  dplyr::collect()

  DBI::dbWriteTable(
      con,
      "hashes", # The name of the table to replace
      hashes_df,
      overwrite = TRUE,
      temporary = FALSE # Set to FALSE for a persistent table
    )

  message("Made filepaths unique in hashes table")
}

read_seu_path <- function(hash = NULL, sqlite_path = "batch_hashes.sqlite") {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  if (!is.null(hash)) {
    # Query for specific hash
    filepath_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes WHERE hash = ?", params = list(hash))
  } else {
    # Query all records
    filepath_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes")
  }

  filepath_df <- dplyr::slice_tail(filepath_df, n = 1)

  return(filepath_df$filepath)
}

read_batch_hashes <- function(filepath = NULL, sqlite_path = "batch_hashes.sqlite") {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  if (!is.null(filepath)) {
    # Query for specific filepath
    hashes_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes WHERE filepath = ?", params = list(filepath))
    seu <- readRDS(hashes_df$filepath)
    return(unique(seu$batch_hash))
  }
}

#' Upsert resolution dictionary into SQLite
#'
#' @param resolution_df Data frame with columns file_id and resolution.
#' @param sqlite_path Path to the SQLite database file.
#' @return Invisibly returns the number of rows written.
#' @export
upsert_resolution_dictionary <- function(resolution_df, sqlite_path = "batch_hashes.sqlite") {
  stopifnot(all(c("file_id", "resolution") %in% colnames(resolution_df)))

  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA busy_timeout = 30000")

  DBI::dbExecute(con, paste0(
    "CREATE TABLE IF NOT EXISTS resolution_dictionary ",
    "(file_id TEXT PRIMARY KEY, resolution TEXT, updated_at TEXT)"
  ))

  for (i in seq_len(nrow(resolution_df))) {
    DBI::dbExecute(
      con,
      paste0(
        "INSERT OR REPLACE INTO resolution_dictionary ",
        "(file_id, resolution, updated_at) VALUES (?, ?, ?)"
      ),
      params = list(
        resolution_df$file_id[[i]],
        as.character(resolution_df$resolution[[i]]),
        format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )
    )
  }

  invisible(nrow(resolution_df))
}

#' Read resolution dictionary from SQLite
#'
#' @param sqlite_path Path to the SQLite database file.
#' @return Named list mapping file_id to resolution.
#' @export
read_resolution_dictionary <- function(sqlite_path = "batch_hashes.sqlite") {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  if (!DBI::dbExistsTable(con, "resolution_dictionary")) {
    return(list())
  }

  resolution_df <- DBI::dbGetQuery(
    con,
    "SELECT file_id, resolution FROM resolution_dictionary ORDER BY file_id"
  )

  if (nrow(resolution_df) == 0) {
    return(list())
  }

  resolution_list <- as.list(resolution_df$resolution)
  names(resolution_list) <- resolution_df$file_id
  resolution_list
}

encode_cluster_order_to_hash_table <- function(
    cluster_orders, sqlite_path = "batch_hashes.sqlite") {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA busy_timeout = 30000")
  # Create table if not exists
  DBI::dbExecute(con, paste0(
    "CREATE TABLE IF NOT EXISTS cluster_orders ",
    "(file_id TEXT PRIMARY KEY, cluster_order TEXT)"
  ))

  for (file_id in names(cluster_orders)) {
    # Convert to JSON
    cluster_order_json <- jsonlite::toJSON(
      cluster_orders[[file_id]], auto_unbox = TRUE
    )
    # Insert or replace
    DBI::dbExecute(con, paste0(
      "INSERT OR REPLACE INTO cluster_orders ",
      "(file_id, cluster_order) VALUES (?, ?)"
    ), params = list(file_id, cluster_order_json))
  }
}

read_cluster_orders_table <- function(
    sqlite_path = "batch_hashes.sqlite",
    file_id = NULL, hash = NULL, as_list = TRUE) {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  if (!DBI::dbExistsTable(con, "cluster_orders")) {
    if (as_list) {
      return(list())
    }
    return(data.frame())
  }

  if (!is.null(file_id)) {
    # Query for specific file_id
    cluster_orders_df <- DBI::dbGetQuery(
      con, "SELECT * FROM cluster_orders WHERE file_id = ?",
      params = list(file_id)
    )
  } else if (!is.null(hash)) {
    # Query for specific hash
    filepaths_df <- DBI::dbGetQuery(
      con, "SELECT * FROM hashes WHERE hash = ?", params = list(hash)
    )
    cluster_orders_df <- DBI::dbGetQuery(
      con, "SELECT * FROM cluster_orders WHERE file_id = ?",
      params = list(fs::path_file(filepaths_df$filepath))
    )
  } else {
    # Query all records
    cluster_orders_df <- DBI::dbGetQuery(con, "SELECT * FROM cluster_orders")
  }

  # Parse JSON if requested
  if (as_list && nrow(cluster_orders_df) > 0) {
    cluster_orders_list <- lapply(cluster_orders_df$cluster_order, function(x) {
      tryCatch({
        jsonlite::fromJSON(x)
      }, error = function(e) {
        warning(
          "Non-JSON format detected for cluster_order entry; skipping. ",
          "Run drop_cluster_orders_table() to reset."
        )
        NULL
      })
    })
    names(cluster_orders_list) <- cluster_orders_df$file_id
    return(cluster_orders_list)
  }

  return(cluster_orders_df)
}

#' Drop and recreate cluster_orders table
#' Use this if you need to clear old non-JSON data
drop_cluster_orders_table <- function(sqlite_path = "batch_hashes.sqlite") {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbExecute(con, "DROP TABLE IF EXISTS cluster_orders")
  message("Dropped cluster_orders table from ", sqlite_path)
}

#' Perform retrieve snakemake params operation
#'
#' @param numbat_rds_file File path
#' @return List object
#' @export
retrieve_snakemake_params <- function(numbat_rds_file) {
  sample_id <- stringr::str_extract(numbat_rds_file, "SR[RX][0-9]+")
  log_file <- fs::path(fs::path_dir(numbat_rds_file), "log.txt")
  log <- readr::read_lines(log_file)[3:26] |>
    stringr::str_split(" = ") |>
    purrr::transpose() |>
    identity()
  params <- log[[2]] |>
    purrr::set_names(log[[1]])
  return(list(sample_id, params))
}

#' Perform retrieve current param operation
#'
#' @param current_params Parameter for current params
#' @param myparam Parameter for myparam
#' @return Named list of param values by sample ID
#' @export
retrieve_current_param <- function(current_params, myparam) {
  sample_ids <- purrr::map(current_params, 1)
  purrr::map(current_params, 2) |>
    purrr::map(myparam) |>
    purrr::set_names(sample_ids)
}

#' Save retained cell barcodes for a Seurat object to the database
#'
#' @param filepath Path to the saved Seurat RDS file (used as primary key)
#' @param sample_id Sample identifier
#' @param seu_type Label for the seu subset type (e.g. "filtered", "hypoxia_low")
#' @param cells Character vector of retained cell barcodes
#' @param sqlite_path Path to the SQLite database
#' @return Invisibly returns filepath
#' @export
save_cell_barcodes_to_db <- function(filepath, sample_id, seu_type, cells, sqlite_path = "batch_hashes.sqlite") {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA busy_timeout = 30000")
  DBI::dbExecute(con, paste0(
    "CREATE TABLE IF NOT EXISTS seu_cells ",
    "(filepath TEXT PRIMARY KEY, sample_id TEXT, seu_type TEXT, cells TEXT)"
  ))
  DBI::dbExecute(con,
    "INSERT OR REPLACE INTO seu_cells (filepath, sample_id, seu_type, cells) VALUES (?, ?, ?, ?)",
    params = list(filepath, sample_id, seu_type, paste(cells, collapse = "\n"))
  )
  invisible(filepath)
}

#' Read retained cell barcodes for a Seurat object from the database
#'
#' @param filepath Path to the Seurat RDS file
#' @param sqlite_path Path to the SQLite database
#' @return Character vector of cell barcodes, or NULL if not found
#' @export
read_cell_barcodes_from_db <- function(filepath, sqlite_path = "batch_hashes.sqlite") {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  if (!DBI::dbExistsTable(con, "seu_cells")) return(NULL)
  result <- DBI::dbGetQuery(con, "SELECT cells FROM seu_cells WHERE filepath = ?", params = list(filepath))
  if (nrow(result) == 0 || is.na(result$cells[[1]])) return(NULL)
  strsplit(result$cells[[1]], "\n")[[1]]
}

#' Compute clone simplifications from a numbat RDS file
#'
#' Derives a named list mapping SCNA labels (e.g. "1q+", "16q-") to primary
#' segment names by reading segs_consensus from the numbat object and
#' classifying each non-neutral segment by chromosomal arm (p = before
#' centromere, q = after) using hg38 centromere midpoints.
#'
#' When multiple segments share the same SCNA label, the alphabetically first
#' segment name is chosen as the primary representative.
#'
#' @param numbat_rds_path Path to a numbat RDS file
#' @return Named list of SCNA_label -> segment_name, or an empty list if no
#'   non-neutral segments are found
#' @export
compute_clone_simplifications <- function(numbat_rds_path) {
  # hg38 centromere midpoints in bp (approximate, from UCSC genome browser)
  centromeres_bp <- c(
    "1" = 123500000L, "2" = 93100000L,  "3" = 92200000L,  "4" = 50700000L,
    "5" = 48300000L,  "6" = 59200000L,  "7" = 59800000L,  "8" = 45000000L,
    "9" = 44400000L,  "10"= 40600000L,  "11"= 52700000L,  "12"= 35600000L,
    "13"= 17000000L,  "14"= 17100000L,  "15"= 18400000L,  "16"= 37300000L,
    "17"= 24900000L,  "18"= 18200000L,  "19"= 25900000L,  "20"= 28200000L,
    "21"= 11900000L,  "22"= 14000000L
  )

  mynb <- readRDS(numbat_rds_path)
  sc <- as.data.frame(mynb$segs_consensus)

  non_neu <- sc[sc$cnv_state != "neu", , drop = FALSE]
  if (nrow(non_neu) == 0) return(list())

  chrom_str <- as.character(non_neu$CHROM)
  midpoint  <- (non_neu$seg_start + non_neu$seg_end) / 2

  arm <- mapply(function(chr, mid) {
    if (!chr %in% names(centromeres_bp)) return("?")
    if (mid < centromeres_bp[[chr]]) "p" else "q"
  }, chrom_str, midpoint, USE.NAMES = FALSE)

  suffix <- dplyr::case_when(
    non_neu$cnv_state %in% c("amp",  "bamp") ~ "+",
    non_neu$cnv_state %in% c("del",  "bdel") ~ "-",
    non_neu$cnv_state %in% c("loh", "cnloh") ~ "cnloh",
    TRUE ~ non_neu$cnv_state
  )

  non_neu$scna_label <- paste0(chrom_str, arm, suffix)
  non_neu$arm <- arm

  # One representative segment per SCNA label: alphabetically first
  result <- non_neu[order(non_neu$scna_label, non_neu$seg), ]
  result <- result[!duplicated(result$scna_label), ]

  out <- as.list(result$seg)
  names(out) <- result$scna_label
  out
}