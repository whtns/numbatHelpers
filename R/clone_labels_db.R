#' Subclone labels stored in batch_hashes.sqlite
#'
#' Subclone ("scna") labels are a function of a sample's numbat genotype string
#' (GT_opt, e.g. "1a,16b") and config/large_clone_simplifications.yaml. They used to
#' be written into each Seurat object, but the per-sample targets that build the
#' objects never rebuild on a YAML change (cue depend = FALSE) and the saved objects
#' carry an empty scna column (docs/initial_clone_scna_audit.md). Instead, one
#' target writes the labels to the `clone_labels` table and plotting / integration
#' code joins them onto the cells at load time with [add_scna_labels()].
#'
#' Table `clone_labels`, one row per sample and clone:
#'   sample_id, clone_opt, GT_opt, n_cells, scna, yaml_md5, recorded_at
#' `scna` is the space-joined label ("16q- 1q+"); "" for the diploid clone (GT_opt
#' ""), and the raw GT_opt when none of its segments are in the YAML.
#' @name clone_labels_db
NULL

#' Build the subclone label table from numbat outputs and the clone simplifications
#'
#' @param numbat_rds_files paths to `<sample>_numbat.rds`; clone_post is read from each
#' @param large_clone_simplifications parsed large_clone_simplifications.yaml
#' @return tibble with sample_id, clone_opt, GT_opt, n_cells, scna
#' @export
build_clone_labels <- function(numbat_rds_files, large_clone_simplifications) {
  purrr::map_dfr(numbat_rds_files, function(rds) {
    sample_id <- stringr::str_extract(rds, "SR[RX][0-9]+")
    clone_post <- readRDS(rds)[["clone_post"]]
    if (is.null(clone_post) || !all(c("clone_opt", "GT_opt") %in% colnames(clone_post))) {
      message("no clone_post for ", sample_id, "; skipped")
      return(NULL)
    }
    key <- large_clone_simplifications[[sample_id]]
    lookup <- if (length(key) == 0) {
      tibble::tibble(scna = character(), seg = character())
    } else {
      tibble::enframe(key, "scna", "seg") |>
        tidyr::unnest(seg) |>
        dplyr::mutate(seg = as.character(seg))
    }
    clone_post |>
      dplyr::count(clone_opt, GT_opt, name = "n_cells") |>
      dplyr::mutate(
        sample_id = sample_id,
        clone_opt = as.integer(clone_opt),
        GT_opt = dplyr::coalesce(as.character(GT_opt), ""),
        scna = purrr::map_chr(GT_opt, function(gt) {
          if (gt == "") return("")
          label <- simplify_gt_col(gt, lookup)
          if (length(label) == 0 || label[[1]] == "") return(gt)
          # one label per segment; a label covering several segments appears once
          paste(unique(strsplit(label[[1]], ",")[[1]]), collapse = " ")
        })
      ) |>
      dplyr::select(sample_id, clone_opt, GT_opt, n_cells, scna)
  })
}

#' Replace the stored labels for the samples in `labels`
#'
#' Only the clone_labels_sqlite target should call this (one writer; see the
#' one-tar_make rule in CLAUDE.md).
#'
#' @param labels output of [build_clone_labels()]
#' @param yaml_md5 md5 of the clone simplifications YAML the labels came from
#' @param sqlite_path the hash database
#' @return hash of `labels`, invisibly; a target returns it so dependents rerun only
#'   when the labels change
#' @export
write_clone_labels_db <- function(labels, yaml_md5, sqlite_path = "batch_hashes.sqlite") {
  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  rows <- dplyr::mutate(labels, yaml_md5 = unname(yaml_md5), recorded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  db_retry(DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, paste(
      "CREATE TABLE IF NOT EXISTS clone_labels (",
      "sample_id TEXT, clone_opt INTEGER, GT_opt TEXT, n_cells INTEGER, scna TEXT,",
      "yaml_md5 TEXT, recorded_at TEXT, PRIMARY KEY (sample_id, clone_opt))"
    ))
    for (s in unique(rows$sample_id)) {
      DBI::dbExecute(con, "DELETE FROM clone_labels WHERE sample_id = ?", params = list(s))
    }
    DBI::dbAppendTable(con, "clone_labels", as.data.frame(rows))
  }))
  invisible(rlang::hash(as.data.frame(labels)))
}

#' Read stored subclone labels
#'
#' @param sample_ids optional sample ids to return
#' @param sqlite_path the hash database
#' @return tibble as written by [write_clone_labels_db()]; empty if the table is missing
#' @export
read_clone_labels <- function(sample_ids = NULL, sqlite_path = "batch_hashes.sqlite") {
  con <- connect_hash_db(sqlite_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  if (!DBI::dbExistsTable(con, "clone_labels")) {
    stop("no clone_labels table in ", sqlite_path, "; build the clone_labels_sqlite target")
  }
  out <- tibble::as_tibble(DBI::dbReadTable(con, "clone_labels"))
  if (!is.null(sample_ids)) out <- dplyr::filter(out, sample_id %in% sample_ids)
  out
}

#' Join stored subclone labels onto a Seurat object's cells
#'
#' Matches on GT_opt within the sample, so the join does not depend on clone
#' numbering. Sets `seu$scna` to a factor ordered by clone number. Cells without a
#' clone_opt (not in numbat's clone_post) get NA.
#'
#' @param seu Seurat object with GT_opt in its metadata
#' @param sample_id the object's sample; a vector of the same length as the cells
#'   for merged objects
#' @param diploid_label label for the diploid clone (stored as "")
#' @param sqlite_path the hash database
#' @return `seu` with `scna` replaced
#' @export
add_scna_labels <- function(seu, sample_id, diploid_label = "", sqlite_path = "batch_hashes.sqlite") {
  if (!"GT_opt" %in% colnames(seu@meta.data)) stop("no GT_opt column; cannot label subclones")
  cell_sample <- if (length(sample_id) == 1) rep(sample_id, ncol(seu)) else sample_id
  labels <- read_clone_labels(unique(cell_sample), sqlite_path) |>
    dplyr::distinct(sample_id, GT_opt, .keep_all = TRUE)
  missing_samples <- setdiff(unique(cell_sample), labels$sample_id)
  if (length(missing_samples) > 0) {
    stop("no stored clone labels for ", paste(missing_samples, collapse = ", "))
  }
  # numbat's diploid clone has GT_opt "", which reaches the objects as NA; a cell
  # with a clone_opt but no GT_opt is diploid, a cell without a clone_opt is unknown.
  gt <- as.character(seu$GT_opt)
  if ("clone_opt" %in% colnames(seu@meta.data)) gt[is.na(gt) & !is.na(seu$clone_opt)] <- ""
  cells <- tibble::tibble(sample_id = cell_sample, GT_opt = gt) |>
    dplyr::left_join(labels, by = c("sample_id", "GT_opt"))
  unmatched <- unique(cells$GT_opt[!is.na(cells$GT_opt) & is.na(cells$scna)])
  if (length(unmatched) > 0) {
    warning("GT_opt not in stored clone labels (stale object?): ", paste(unmatched, collapse = "; "))
  }
  display <- ifelse(!is.na(cells$scna) & cells$scna == "", diploid_label, cells$scna)
  lvls <- cells |>
    dplyr::mutate(scna = display) |>
    dplyr::filter(!is.na(scna)) |>
    dplyr::group_by(scna) |>
    dplyr::summarise(first_clone = min(clone_opt), .groups = "drop") |>
    dplyr::arrange(first_clone) |>
    dplyr::pull(scna)
  seu$scna <- factor(display, levels = lvls)
  seu
}
