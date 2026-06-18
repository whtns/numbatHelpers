# Centralized SQLite connection helpers for the hash/metadata database.
#
# Background: batch_hashes.sqlite lives on a networked filesystem (NFS/VAST) and
# is accessed concurrently by many crew workers running on different physical
# nodes. SQLite's WAL journal mode relies on a shared-memory index (the -shm
# file) that does NOT work over networked filesystems, which produced
# intermittent SQLITE_PROTOCOL errors surfacing in R as "locking protocol".
# Those errors caused targets branches (error = "null") to be cached as NULL.
#
# Fixes implemented here:
#   * Use rollback (TRUNCATE) journal mode instead of WAL -> eliminates the
#     SQLITE_PROTOCOL ("locking protocol") root cause on NFS.
#   * Apply a long busy_timeout to EVERY connection (consistently) so SQLite
#     waits on a locked DB rather than erroring immediately (SQLITE_BUSY).
#   * Wrap connection acquisition in bounded exponential-backoff retry to absorb
#     any residual transient lock contention.

#' Retry an expression on transient SQLite locking errors.
#'
#' Retries on SQLITE_BUSY / SQLITE_LOCKED / SQLITE_PROTOCOL and related transient
#' lock messages ("database is locked", "locking protocol", "disk I/O error")
#' with exponential backoff plus jitter. Non-transient errors are re-raised
#' immediately.
#'
#' @param expr Expression to evaluate (lazily; evaluated in the caller's frame).
#' @param max_tries Maximum number of attempts.
#' @param base_sleep Initial backoff in seconds; doubles each retry.
#' @return The value of \code{expr}.
#' @keywords internal
db_retry <- function(expr, max_tries = 8, base_sleep = 0.25) {
  expr <- substitute(expr)
  penv <- parent.frame()
  transient <- paste(
    "locking protocol", "database is locked", "database table is locked",
    "disk I/O error", "SQLITE_BUSY", "SQLITE_LOCKED", "SQLITE_PROTOCOL",
    sep = "|"
  )
  for (attempt in seq_len(max_tries)) {
    res <- tryCatch(
      list(ok = TRUE, value = eval(expr, envir = penv)),
      error = function(e) list(ok = FALSE, error = e)
    )
    if (isTRUE(res$ok)) return(res$value)
    msg <- conditionMessage(res$error)
    if (attempt == max_tries || !grepl(transient, msg, ignore.case = TRUE)) {
      stop(res$error)
    }
    Sys.sleep(base_sleep * (2^(attempt - 1)) + stats::runif(1, 0, base_sleep))
  }
}

#' Open an NFS-safe connection to the hash SQLite database.
#'
#' Connects and applies rollback (TRUNCATE) journal mode and a long busy_timeout
#' on every connection, retrying the acquisition on transient lock errors. Use
#' in place of \code{DBI::dbConnect(RSQLite::SQLite(), path)}; the caller is still
#' responsible for \code{DBI::dbDisconnect()} (typically via \code{on.exit}).
#'
#' @param sqlite_path Path to the SQLite database file.
#' @param busy_timeout_ms Milliseconds SQLite waits on a locked DB before erroring.
#' @param max_tries Maximum connection attempts on transient lock errors.
#' @return A DBIConnection with pragmas applied.
#' @keywords internal
connect_hash_db <- function(sqlite_path = "batch_hashes.sqlite",
                            busy_timeout_ms = 60000, max_tries = 8) {
  db_retry(
    {
      con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
      ok <- FALSE
      on.exit(if (!ok) try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)
      DBI::dbExecute(con, "PRAGMA journal_mode = TRUNCATE")
      DBI::dbExecute(con, sprintf("PRAGMA busy_timeout = %d", as.integer(busy_timeout_ms)))
      ok <- TRUE
      con
    },
    max_tries = max_tries
  )
}
