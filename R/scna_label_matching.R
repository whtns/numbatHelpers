#' Test subclone labels for an arm-level SCNA
#'
#' Subclone labels (the `scna` column, built from
#' `config/large_clone_simplifications.yaml`) are space-joined event names such as
#' "16q- 1q+". Plain substring matching misfires once every segment is labelled:
#' "12p-" contains "2p" and "16pcnloh" contains "6p", while the whole-chromosome
#' gain "2+" does not contain "2p". This matches the chromosome number only when it
#' is not preceded by another digit, and counts a whole-chromosome event of the
#' same sign ("1+", "16-") as the arm event. Focal labels ("6 p+ focal") never match.
#'
#' @param labels character vector of subclone labels
#' @param scna arm event: "1q", "2p", or with a sign, "1q+", "16q-"
#' @return logical vector, FALSE for NA labels
#' @export
scna_label_has <- function(labels, scna) {
  parts <- stringr::str_match(scna, "^([0-9]+)([pq])([+-]?)$")
  if (is.na(parts[1, 1])) {
    return(!is.na(labels) & stringr::str_detect(labels, stringr::fixed(scna)))
  }
  chrom <- parts[1, 2]; arm <- parts[1, 3]; sign <- parts[1, 4]
  pattern <- if (nzchar(sign)) {
    paste0("(?<![0-9])", chrom, arm, "?\\", sign)
  } else {
    paste0("(?<![0-9])", chrom, "(?:", arm, "|(?=[+-]))")
  }
  !is.na(labels) & stringr::str_detect(as.character(labels), pattern)
}
