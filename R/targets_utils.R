# Rename a target in the targets store without remaking it.
#
# Handles simple object targets and dynamic branching targets (pattern +
# branch rows). Updates both the object files in `objects/` and the name,
# parent, and children fields in the meta file.
#
# Args:
#   old_name      : character - current target name
#   new_name      : character - desired target name
#   store         : character - path to the targets store
#   targets_file  : character - path to _targets.R to update (NULL to skip)
#
# Returns invisibly: a character vector of renamed object file paths.
rename_target <- function(old_name, new_name,
                          store = "_targets_r431",
                          targets_file = "_targets.R") {

  meta_path <- file.path(store, "meta", "meta")
  objects_dir <- file.path(store, "objects")

  stopifnot(file.exists(meta_path))

  # ── 0. Collision check ────────────────────────────────────────────────────
  meta_names <- read.table(
    meta_path, sep = "|", header = TRUE,
    quote = "", comment.char = ""
  )$name
  if (new_name %in% meta_names) {
    stop(sprintf(
      "A target named '%s' already exists in the store.", new_name
    ))
  }
  new_obj_exists <- length(list.files(
    objects_dir,
    pattern = paste0("^", new_name, "(_[0-9a-f]+)?$")
  )) > 0
  if (new_obj_exists) {
    stop(sprintf(
      "Object file(s) for '%s' already exist in '%s'.",
      new_name, objects_dir
    ))
  }

  # ── 1. Rename object files ────────────────────────────────────────────────
  # Simple target: objects/old_name
  # Dynamic branches: objects/old_name_<hash>
  old_files <- list.files(
    objects_dir,
    pattern = paste0("^", old_name, "(_[0-9a-f]+)?$"),
    full.names = TRUE
  )

  renamed_files <- character(0)
  for (f in old_files) {
    base <- basename(f)
    new_base <- sub(paste0("^", old_name), new_name, base)
    new_f <- file.path(objects_dir, new_base)
    file.rename(f, new_f)
    renamed_files <- c(renamed_files, new_f)
  }

  # ── 2. Update meta file ───────────────────────────────────────────────────
  lines <- readLines(meta_path)

  # Helper: replace whole-word occurrences of old_name in a field value,
  # where "word" means bounded by | or * (the two delimiters used in meta).
  replace_name_in_field <- function(x) {
    # Replace old_name that appears at start, end, or between delimiters
    gsub(
      paste0("(^|[|*])", old_name, "($|[|*])"),
      paste0("\\1", new_name, "\\2"),
      x,
      perl = TRUE
    )
  }

  updated <- vapply(lines, replace_name_in_field, character(1))

  writeLines(updated, meta_path)

  # ── 3. Update _targets.R ─────────────────────────────────────────────────
  targets_lines_updated <- 0L
  if (!is.null(targets_file) && file.exists(targets_file)) {
    tlines <- readLines(targets_file)
    # Replace whole-word occurrences (bounded by non-identifier characters)
    tupdated <- gsub(
      paste0("(?<![A-Za-z0-9_.])", old_name, "(?![A-Za-z0-9_.])"),
      new_name,
      tlines,
      perl = TRUE
    )
    writeLines(tupdated, targets_file)
    targets_lines_updated <- sum(tupdated != tlines)
  }

  message(sprintf(
    paste0("Renamed '%s' -> '%s': %d object file(s) renamed, ",
           "%d meta line(s) updated, %d _targets.R line(s) updated."),
    old_name, new_name,
    length(renamed_files),
    sum(updated != lines),
    targets_lines_updated
  ))

  invisible(renamed_files)
}
