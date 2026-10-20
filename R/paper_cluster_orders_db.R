# Cell-state labels and Fig. S20 cluster comparisons for the paper integrations,
# kept in batch_hashes.sqlite (2026-10-07) instead of data/scna_cluster_order.csv.
#
# The labels live in the existing cluster_orders table under file_id
# "paperint_<scna>", in the same JSON layout pull_cluster_orders() produces, so
# make_clone_distribution_figure() and read_cluster_orders_table() read them
# unchanged. The cluster_orders_sqlite target re-copies only the CSV's file_ids, so
# these rows are left alone unless a paperint_* row is also added to the CSV (the CSV
# would then win). targets cannot see sqlite edits: the readers below are meant for
# targets with cue = tar_cue("always"), whose downstream rebuilds only when the
# returned rows change.

.paper_phase_levels <- c("pm", "g1", "g1_s", "s", "s_g2", "g2", "g2_m",
                         "hsp", "hypoxia", "other", "s_star")

#' Write the cell-state labels of one paper integration to sqlite
#'
#' @param ... One argument per state, clusters joined by "_" in display order, e.g.
#'   `g1 = "0_7_6", pm = "4", other = "8_10"`. Names must be in the state vocabulary.
#' @param scna "1q", "2p", "6p" or "16q"
#' @param resolution Clustering column, e.g. "integrated_snn_res.0.4"
#' @param sqlite_path Hash database
#' @return The stored cluster order (invisibly)
#' @details `scna` and `resolution` come after `...` so they must be named in full;
#'   otherwise a state such as `s =` would partially match `scna`.
#' @export
set_paper_cluster_order <- function(..., scna, resolution, sqlite_path = "batch_hashes.sqlite") {
	states <- list(...)
	if (is.null(names(states)) || any(names(states) == "")) stop("every state must be named")
	bad <- setdiff(names(states), .paper_phase_levels)
	if (length(bad)) stop("unknown states: ", paste(bad, collapse = ", "))
	clusters <- unlist(strsplit(unlist(states), "_"))
	if (anyDuplicated(clusters)) stop("cluster(s) given twice: ", paste(unique(clusters[duplicated(clusters)]), collapse = ", "))

	key <- paste0("paperint_", scna)
	row <- tibble::tibble(file_id = key, sample_id = key, preferred_resolution = 1,
	                      resolution = resolution, tumor_id = key)
	for (p in .paper_phase_levels) row[[p]] <- if (is.null(states[[p]])) NA_character_ else as.character(states[[p]])

	# go through pull_cluster_orders() so the stored layout matches the CSV-derived rows
	tmp <- tempfile(fileext = ".csv")
	on.exit(unlink(tmp), add = TRUE)
	readr::write_csv(row, tmp, na = "")
	order <- suppressMessages(pull_cluster_orders(tmp))
	db_retry(encode_cluster_order_to_hash_table(order, sqlite_path))
	invisible(order[[key]])
}

#' Read the cell-state labels of one paper integration from sqlite
#'
#' @param scna "1q", "2p", "6p" or "16q"
#' @param sqlite_path Hash database
#' @return A data frame with resolution, phase (factor, state vocabulary order) and
#'   clusters, one row per cluster in display order
#' @export
read_paper_cluster_order <- function(scna, sqlite_path = "batch_hashes.sqlite") {
	key <- paste0("paperint_", scna)
	order <- db_retry(read_cluster_orders_table(sqlite_path, file_id = key))
	if (!length(order) || is.null(order[[key]])) stop("no cluster order for ", key, " in ", sqlite_path)
	order[[key]][[1]] |>
		tibble::as_tibble() |>
		dplyr::filter(!is.na(clusters), clusters != "") |>
		dplyr::mutate(phase = factor(phase, levels = .paper_phase_levels)) |>
		dplyr::arrange(phase) |>
		dplyr::select(resolution, phase, clusters)
}

#' Export the paper cluster orders as one wide row per integration (for review)
#'
#' @param path Output CSV, in the data/scna_cluster_order.csv layout
#' @param scnas Integrations to export
#' @param sqlite_path Hash database
#' @return path
#' @export
export_paper_cluster_orders <- function(path = "results/paper_cluster_orders.csv",
                                        scnas = c("1q", "16q", "2p", "6p"),
                                        sqlite_path = "batch_hashes.sqlite") {
	rows <- lapply(scnas, function(s) {
		co <- tryCatch(read_paper_cluster_order(s, sqlite_path), error = function(e) NULL)
		if (is.null(co)) return(NULL)
		co |>
			dplyr::group_by(resolution, phase) |>
			dplyr::summarise(clusters = paste(clusters, collapse = "_"), .groups = "drop") |>
			tidyr::pivot_wider(names_from = phase, values_from = clusters) |>
			dplyr::mutate(file_id = paste0("paperint_", .env$s), .before = 1)  # .env: a state column is named s
	})
	dplyr::bind_rows(rows) |>
		dplyr::select(file_id, resolution, dplyr::any_of(.paper_phase_levels)) |>
		readr::write_csv(path, na = "")
	path
}

#' Write / read the Fig. S20 cluster comparisons (table paper_cluster_comparisons)
#'
#' @param comparisons Data frame with columns ident.1, ident.2 (`clusters` levels of
#'   the annotated 1q integration, e.g. "g1_7"); replaces the stored rows
#' @param sqlite_path Hash database
#' @export
set_paper_cluster_comparisons <- function(comparisons, sqlite_path = "batch_hashes.sqlite") {
	stopifnot(all(c("ident.1", "ident.2") %in% names(comparisons)))
	con <- connect_hash_db(sqlite_path)
	on.exit(DBI::dbDisconnect(con), add = TRUE)
	db_retry(DBI::dbWriteTable(con, "paper_cluster_comparisons",
	                           as.data.frame(comparisons[c("ident.1", "ident.2")]), overwrite = TRUE))
	invisible(comparisons)
}

#' @rdname set_paper_cluster_comparisons
#' @export
read_paper_cluster_comparisons <- function(sqlite_path = "batch_hashes.sqlite") {
	con <- connect_hash_db(sqlite_path)
	on.exit(DBI::dbDisconnect(con), add = TRUE)
	if (!DBI::dbExistsTable(con, "paper_cluster_comparisons"))
		stop("no paper_cluster_comparisons table in ", sqlite_path)
	tibble::as_tibble(db_retry(DBI::dbReadTable(con, "paper_cluster_comparisons")))
}
