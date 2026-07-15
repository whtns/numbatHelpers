# Shared guard for target outputs that are "a path, or NA if there was nothing to do".
#
# find_diffex_bw_clones_for_each_cluster() / find_diffex_clones() return
# NA_character_ for a sample with nothing to compare -- no numbat run, or no clone
# comparison carrying the requested SCNA (e.g. SRX10831287 has only an 11a+
# comparison, so it has no 1q+). That is a legitimate skip, not a failure.
#
# Every downstream consumer flattens those lists and reads each element as a file,
# so an unfiltered NA becomes a filename and dies with
#   'NA' does not exist in current working directory
# which errors the whole target. This bit three separate consumers
# (plot_fig_07_08, tabulate_diffex_clones, make_volcano_diffex_clones), so the
# filtering lives in one place rather than being rediscovered one crash at a time.

#' Flatten a list of file paths, dropping skips and missing files
#'
#' @param x A path, or list/vector of paths, possibly containing `NA` (a skipped
#'   sample) or paths that were never written.
#' @return A plain character vector of paths that exist on disk (possibly empty).
#' @keywords internal
.drop_missing_paths <- function(x) {
  p <- unlist(x, use.names = FALSE)
  if (is.null(p)) return(character(0))
  p <- as.character(p)
  p[!is.na(p) & nzchar(p) & file.exists(p)]
}

#' Name a per-sample diffex list, erroring clearly on a length mismatch
#'
#' The oncoprint builders name their cis/trans/all diffex lists from a vector of
#' sample ids (`names(x) <- ids`). Base R recycles or errors cryptically when the
#' lengths disagree ("'names' attribute [34] must be the same length as the vector
#' [33]"), and -- worse -- a same-length-but-wrong vector silently mislabels every
#' sample. The list must be named from the SAME target it maps over (the
#' low-hypoxia sample ids), one id per branch, in branch order.
#'
#' @param x A per-branch list (one element per low-hypoxia sample).
#' @param ids Character vector of sample ids, one per element of `x`.
#' @param what Label for the error message.
#' @return `x` with `names(x) <- ids`.
#' @keywords internal
.name_by_sample <- function(x, ids, what = "diffex list") {
  if (length(ids) != length(x)) {
    stop(sprintf(
      paste0("%s has %d elements but %d sample ids were supplied. The naming ",
             "vector must be one id per branch, in branch order -- derive it from ",
             "the same target the diffex list maps over (seus_low_hypoxia), not ",
             "from debranched_ids (which is branch-level and a different length)."),
      what, length(x), length(ids)), call. = FALSE)
  }
  stats::setNames(x, ids)
}
