# Seurat object metadata tracking in SQLite
#
# Schema
# ------
# seurat_objects   - one row per RDS file (identity, size, provenance)
# cell_metadata    - one row per (filepath, column): dtype, n_unique, JSON summary
# cell_qc_values   - one row per (filepath, cell): per-cell QC metrics for fast queries
# cluster_composition - one row per (filepath, cluster): cell counts
# cluster_markers  - one row per (filepath, cluster, rank): top marker genes
# qc_metrics       - one row per (filepath, metric): quantile summary
# hashes           - legacy: filepath → content hash  (kept for compatibility)
# cluster_orders   - legacy: file_id → JSON cluster order (kept for compatibility)

DEFAULT_DB <- "batch_hashes.sqlite"

# ── Schema ─────────────────────────────────────────────────────────────────────

#' Initialise all metadata tables (idempotent)
#'
#' @param sqlite_path Path to the SQLite database file.
#' @return Invisibly returns the database path.
#' @export
init_seu_metadata_db <- function(sqlite_path = DEFAULT_DB) {
  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con))

  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS seurat_objects (
      filepath             TEXT PRIMARY KEY,
      hash                 TEXT,
      n_cells              INTEGER,
      n_features           INTEGER,
      sample_id            TEXT,
      scna_type            TEXT,
      processing_stage     TEXT,
      default_assay        TEXT,
      metadata_columns_json TEXT,
      recorded_at          TEXT
    )")

  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS cell_metadata (
      filepath     TEXT,
      column_name  TEXT,
      dtype        TEXT,
      n_cells      INTEGER,
      n_unique     INTEGER,
      n_na         INTEGER,
      summary_json TEXT,
      PRIMARY KEY (filepath, column_name)
    )")

  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS cell_qc_values (
      filepath      TEXT,
      cell          TEXT,
      sample_id     TEXT,
      nCount_gene   REAL,
      nFeature_gene REAL,
      percent_mt    REAL,
      PRIMARY KEY (filepath, cell)
    )")

  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS cluster_composition (
      filepath TEXT,
      cluster  TEXT,
      n_cells  INTEGER,
      pct_cells REAL,
      PRIMARY KEY (filepath, cluster)
    )")

  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS cluster_markers (
      filepath    TEXT,
      cluster     TEXT,
      marker_rank INTEGER,
      gene_name   TEXT,
      PRIMARY KEY (filepath, cluster, marker_rank)
    )")

  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS qc_metrics (
      filepath TEXT,
      metric   TEXT,
      min      REAL,
      q25      REAL,
      median   REAL,
      mean     REAL,
      q75      REAL,
      max      REAL,
      PRIMARY KEY (filepath, metric)
    )")

  # Legacy tables kept for backward compatibility
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS hashes (
      filepath TEXT PRIMARY KEY,
      hash     TEXT
    )")

  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS cluster_orders (
      file_id       TEXT PRIMARY KEY,
      cluster_order TEXT
    )")

  invisible(sqlite_path)
}

# ── Extraction ─────────────────────────────────────────────────────────────────

#' Extract metadata from a Seurat object and write to SQLite
#'
#' Reads the RDS at \code{filepath} (or uses an in-memory object if \code{seu}
#' is provided), derives the content hash, and upserts six metadata tables.
#' The function is idempotent: re-running on the same file updates existing rows.
#'
#' @param filepath Path to the Seurat RDS file. Used as the primary key.
#' @param seu Optional Seurat object already in memory. If NULL, read from filepath.
#' @param scna_type SCNA type label (e.g. "1q", "2p", "6p", "16q"). Inferred
#'   from \code{filepath} when NULL.
#' @param processing_stage Processing stage label (e.g. "debranched", "integrated").
#'   Inferred from \code{filepath} when NULL.
#' @param sqlite_path Path to the SQLite database.
#' @param qc_cols Character vector of QC numeric columns to summarise.
#' @return Invisibly returns \code{filepath}.
#' @export
extract_seu_metadata <- function(
    filepath,
    seu            = NULL,
    scna_type      = NULL,
    processing_stage = NULL,
    sqlite_path    = DEFAULT_DB,
    qc_cols        = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                       "nCount_SCT", "nFeature_SCT",
                       "log_nCount_gene", "nCount_gene")
) {
  if (is.null(seu)) {
    seu <- readRDS(filepath)
  }

  init_seu_metadata_db(sqlite_path)
  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con))

  md       <- seu@meta.data
  n_cells  <- nrow(md)
  hash     <- digest::digest(colnames(seu))

  # Infer labels from filepath when not supplied
  if (is.null(scna_type)) {
    scna_type <- dplyr::case_when(
      grepl("_1q",  filepath) ~ "1q",
      grepl("_2p",  filepath) ~ "2p",
      grepl("_6p",  filepath) ~ "6p",
      grepl("_16q", filepath) ~ "16q",
      TRUE ~ NA_character_
    )
  }

  if (is.null(processing_stage)) {
    processing_stage <- dplyr::case_when(
      grepl("integrated", filepath) ~ "integrated",
      grepl("branch",     filepath) ~ "debranched",
      grepl("filtered",   filepath) ~ "filtered",
      TRUE ~ "raw"
    )
  }

  sample_id <- stringr::str_extract(filepath, "SR[RX][0-9]+")

  # ── seurat_objects ──────────────────────────────────────────────────────────
  DBI::dbExecute(con, "
    INSERT OR REPLACE INTO seurat_objects
      (filepath, hash, n_cells, n_features, sample_id, scna_type,
       processing_stage, default_assay, metadata_columns_json, recorded_at)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
    params = list(
      filepath,
      hash,
      n_cells,
      nrow(seu),
      sample_id,
      scna_type,
      processing_stage,
      Seurat::DefaultAssay(seu),
      jsonlite::toJSON(colnames(md), auto_unbox = FALSE),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
  )

  # ── hashes (legacy) ─────────────────────────────────────────────────────────
  DBI::dbExecute(con,
    "INSERT OR REPLACE INTO hashes (filepath, hash) VALUES (?, ?)",
    params = list(filepath, hash)
  )

  # ── cell_metadata ───────────────────────────────────────────────────────────
  col_rows <- lapply(colnames(md), function(col) {  # nolint: object_overwrite_linter
    vals  <- md[[col]]
    dtype <- class(vals)[1]
    n_na  <- sum(is.na(vals))

    summary_json <- if (is.numeric(vals)) {
      qs <- quantile(vals, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
      as.character(jsonlite::toJSON(list(
        min    = base::unname(qs[1]),
        q25    = base::unname(qs[2]),
        median = base::unname(qs[3]),
        mean   = mean(vals, na.rm = TRUE),
        q75    = base::unname(qs[4]),
        max    = base::unname(qs[5])
      ), auto_unbox = TRUE))
    } else {
      counts <- sort(table(vals), decreasing = TRUE)
      # cap at 50 most common values to keep JSON small
      as.character(jsonlite::toJSON(as.list(head(counts, 50)), auto_unbox = TRUE))
    }

    list(
      filepath    = filepath,
      column_name = col,
      dtype       = dtype,
      n_cells     = n_cells,
      n_unique    = length(unique(vals)),
      n_na        = n_na,
      summary_json = summary_json
    )
  })

  cell_meta_df <- dplyr::bind_rows(col_rows)
  # upsert row by row via DELETE + INSERT for SQLite compatibility
  for (i in seq_len(nrow(cell_meta_df))) {
    r <- cell_meta_df[i, ]
    DBI::dbExecute(con, "
      INSERT OR REPLACE INTO cell_metadata
        (filepath, column_name, dtype, n_cells, n_unique, n_na, summary_json)
      VALUES (?, ?, ?, ?, ?, ?, ?)",
      params = base::unname(as.list(r))
    )
  }

  # ── cell_qc_values ──────────────────────────────────────────────────────────
  qc_cols_available <- base::intersect(c("nCount_gene", "nFeature_gene", "percent.mt"), colnames(md))
  if (length(qc_cols_available) > 0) {
    qc_df <- md |>
      tibble::rownames_to_column("cell") |>
      dplyr::transmute(
        filepath = filepath,
        cell = cell,
        sample_id = sample_id,
        nCount_gene = if ("nCount_gene" %in% colnames(md)) .data$nCount_gene else NA_real_,
        nFeature_gene = if ("nFeature_gene" %in% colnames(md)) .data$nFeature_gene else NA_real_,
        percent_mt = if ("percent.mt" %in% colnames(md)) .data[["percent.mt"]] else NA_real_
      )

    DBI::dbExecute(con,
      "DELETE FROM cell_qc_values WHERE filepath = ?",
      params = list(filepath)
    )

    if (nrow(qc_df) > 0) {
      DBI::dbWriteTable(con, "cell_qc_values", qc_df, append = TRUE)
    }
  }

  # ── cluster_composition ─────────────────────────────────────────────────────
  cluster_col <- dplyr::case_when(
    "clusters"         %in% colnames(md) ~ "clusters",
    "seurat_clusters"  %in% colnames(md) ~ "seurat_clusters",
    "leiden"           %in% colnames(md) ~ "leiden",
    TRUE ~ NA_character_
  )

  if (!is.na(cluster_col)) {
    counts <- table(md[[cluster_col]])
    comp_df <- data.frame(
      filepath  = filepath,
      cluster   = names(counts),
      n_cells   = as.integer(counts),
      pct_cells = as.numeric(counts) / n_cells * 100,
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(comp_df))) {
      r <- comp_df[i, ]
      DBI::dbExecute(con, "
        INSERT OR REPLACE INTO cluster_composition
          (filepath, cluster, n_cells, pct_cells)
        VALUES (?, ?, ?, ?)",
        params = list(r$filepath, r$cluster, r$n_cells, r$pct_cells)
      )
    }
  }

  # ── cluster_markers ─────────────────────────────────────────────────────────
  markers_tbl <- seu@misc$markers$gene_snn_res.0.2$presto
  if (!is.null(markers_tbl) && nrow(markers_tbl) > 0) {
    gene_col <- dplyr::case_when(
      "Gene.Name" %in% colnames(markers_tbl) ~ "Gene.Name",
      "gene"      %in% colnames(markers_tbl) ~ "gene",
      "feature"   %in% colnames(markers_tbl) ~ "feature",
      TRUE ~ NA_character_
    )

    if (!is.na(gene_col) && "Cluster" %in% colnames(markers_tbl)) {
      markers_df <- markers_tbl |>
        dplyr::mutate(
          cluster = as.character(Cluster),
          gene_name = .data[[gene_col]]
        ) |>
        dplyr::group_by(cluster) |>
        dplyr::slice_head(n = 20) |>
        dplyr::ungroup() |>
        dplyr::group_by(cluster) |>
        dplyr::mutate(marker_rank = dplyr::row_number()) |>
        dplyr::ungroup() |>
        dplyr::transmute(
          filepath = filepath,
          cluster = cluster,
          marker_rank = marker_rank,
          gene_name = as.character(gene_name)
        )

      DBI::dbExecute(con,
        "DELETE FROM cluster_markers WHERE filepath = ?",
        params = list(filepath)
      )

      for (i in seq_len(nrow(markers_df))) {
        r <- markers_df[i, ]
        DBI::dbExecute(con, "
          INSERT OR REPLACE INTO cluster_markers
            (filepath, cluster, marker_rank, gene_name)
          VALUES (?, ?, ?, ?)",
          params = list(r$filepath, r$cluster, r$marker_rank, r$gene_name)
        )
      }
    }
  }

  # ── qc_metrics ──────────────────────────────────────────────────────────────
  present_qc <- base::intersect(qc_cols, colnames(md))
  for (metric in present_qc) {
    vals <- md[[metric]]
    if (!is.numeric(vals)) next
    qs <- quantile(vals, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
    DBI::dbExecute(con, "
      INSERT OR REPLACE INTO qc_metrics
        (filepath, metric, min, q25, median, mean, q75, max)
      VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
      params = list(filepath, metric,
                    base::unname(qs[1]), base::unname(qs[2]), base::unname(qs[3]),
                    mean(vals, na.rm = TRUE),
                    base::unname(qs[4]), base::unname(qs[5]))
    )
  }

  invisible(filepath)
}

# ── Query helpers ───────────────────────────────────────────────────────────────

#' List all tracked Seurat objects
#'
#' @param sqlite_path Path to the SQLite database.
#' @param stage Optional processing stage filter (e.g. "debranched").
#' @param scna Optional SCNA type filter (e.g. "1q").
#' @return A data frame with one row per tracked object.
#' @export
list_seu_objects <- function(sqlite_path = DEFAULT_DB, stage = NULL, scna = NULL) {
  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con))

  q <- "SELECT filepath, hash, n_cells, n_features, sample_id,
               scna_type, processing_stage, default_assay, recorded_at
        FROM seurat_objects WHERE 1=1"
  params <- list()
  if (!is.null(stage)) { q <- paste0(q, " AND processing_stage = ?"); params <- c(params, stage) }
  if (!is.null(scna))  { q <- paste0(q, " AND scna_type = ?");        params <- c(params, scna)  }
  q <- paste0(q, " ORDER BY recorded_at DESC")

  DBI::dbGetQuery(con, q, params = params)
}

#' Summarise metadata columns for one Seurat object
#'
#' @param filepath Path used as primary key in \code{seurat_objects}.
#' @param sqlite_path Path to the SQLite database.
#' @return A data frame with one row per metadata column.
#' @export
get_metadata_summary <- function(filepath, sqlite_path = DEFAULT_DB) {
  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, "
    SELECT column_name, dtype, n_cells, n_unique, n_na, summary_json
    FROM   cell_metadata
    WHERE  filepath = ?
    ORDER  BY column_name",
    params = list(filepath)
  )
}

#' Retrieve per-cell QC values for one Seurat object
#'
#' @param filepath Path used as primary key in \code{seurat_objects}.
#' @param sqlite_path Path to the SQLite database.
#' @return A data frame with cell-level QC columns.
#' @export
get_cell_qc_values <- function(filepath, sqlite_path = DEFAULT_DB) {
  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, "
    SELECT cell, nCount_gene, nFeature_gene, percent_mt
    FROM   cell_qc_values
    WHERE  filepath = ?
    ORDER  BY cell",
    params = list(filepath)
  )
}

#' Retrieve cluster composition for one Seurat object
#'
#' @param filepath Path used as primary key in \code{seurat_objects}.
#' @param sqlite_path Path to the SQLite database.
#' @return A data frame with cluster, n_cells, pct_cells columns.
#' @export
get_cluster_composition <- function(filepath, sqlite_path = DEFAULT_DB) {
  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, "
    SELECT cluster, n_cells, pct_cells
    FROM   cluster_composition
    WHERE  filepath = ?
    ORDER  BY n_cells DESC",
    params = list(filepath)
  )
}

#' Retrieve top marker genes for one Seurat object
#'
#' @param filepath Path used as primary key in \code{seurat_objects}.
#' @param top_n Number of top markers to return per cluster.
#' @param sqlite_path Path to the SQLite database.
#' @return A data frame with cluster, marker_rank, and gene_name columns.
#' @export
get_cluster_markers <- function(filepath, top_n = 5, sqlite_path = DEFAULT_DB) {
  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, "
    SELECT cluster, marker_rank, gene_name
    FROM   cluster_markers
    WHERE  filepath = ?
       AND marker_rank <= ?
    ORDER  BY cluster, marker_rank",
    params = list(filepath, top_n)
  )
}

#' Compare cluster composition across SCNA types or processing stages
#'
#' @param sqlite_path Path to the SQLite database.
#' @param group_by Column from \code{seurat_objects} to group by
#'   (\code{"scna_type"} or \code{"processing_stage"}).
#' @return A data frame with one row per (group, cluster).
#' @export
compare_cluster_composition <- function(sqlite_path = DEFAULT_DB,
                                        group_by = c("scna_type", "processing_stage")) {
  group_by <- match.arg(group_by)
  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, glue::glue("
    SELECT s.{group_by},
           cc.cluster,
           SUM(cc.n_cells)                              AS n_cells,
           SUM(cc.n_cells) * 100.0 / SUM(SUM(cc.n_cells))
             OVER (PARTITION BY s.{group_by})           AS pct_cells
    FROM   cluster_composition cc
    JOIN   seurat_objects       s  ON s.filepath = cc.filepath
    GROUP  BY s.{group_by}, cc.cluster
    ORDER  BY s.{group_by}, n_cells DESC"))
}

#' Retrieve QC metric summaries for one or all objects
#'
#' @param filepath Optional filepath filter.
#' @param metric Optional metric name filter (e.g. "nFeature_RNA").
#' @param sqlite_path Path to the SQLite database.
#' @return A data frame with quantile columns.
#' @export
get_qc_metrics <- function(filepath = NULL, metric = NULL, sqlite_path = DEFAULT_DB) {
  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con))

  q <- "SELECT q.filepath, s.sample_id, s.scna_type, s.processing_stage,
               q.metric, q.min, q.q25, q.median, q.mean, q.q75, q.max
        FROM   qc_metrics q
        JOIN   seurat_objects s ON s.filepath = q.filepath
        WHERE  1=1"
  params <- list()
  if (!is.null(filepath)) { q <- paste0(q, " AND q.filepath = ?"); params <- c(params, filepath) }
  if (!is.null(metric))   { q <- paste0(q, " AND q.metric   = ?"); params <- c(params, metric)   }
  q <- paste0(q, " ORDER BY s.scna_type, s.processing_stage, q.metric")

  DBI::dbGetQuery(con, q, params = params)
}

#' Parse a JSON summary column from get_metadata_summary() into a named list
#'
#' @param summary_json A single JSON string from the \code{summary_json} column.
#' @return A named list or vector.
#' @export
parse_metadata_summary <- function(summary_json) {
  jsonlite::fromJSON(summary_json)
}

#' Bulk-extract metadata from a named character vector of filepaths
#'
#' Designed for use as a \code{tar_target()} command when \code{filepaths} is
#' a targets upstream dependency. Calls \code{extract_seu_metadata()} on each
#' path and returns the filepaths invisibly (so the target stores a small value).
#'
#' @param filepaths Named character vector of Seurat RDS paths.
#' @param sqlite_path Path to the SQLite database.
#' @param ... Additional arguments forwarded to \code{extract_seu_metadata()}.
#' @return Invisibly returns \code{filepaths}.
#' @export
bulk_extract_seu_metadata <- function(filepaths, sqlite_path = DEFAULT_DB, ...) {
  purrr::walk(filepaths, extract_seu_metadata, sqlite_path = sqlite_path, ...)
  invisible(filepaths)
}
