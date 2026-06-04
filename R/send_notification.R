# Helper to send pipeline completion notifications
#
# Uses blastula if SMTP env vars are provided, otherwise falls back to a
# system `mail` command. Recipient and SMTP credentials are read from
# environment variables and not stored in the repository.

send_pipeline_notification <- function(to = Sys.getenv("PIPELINE_NOTIFY_EMAIL"),
                                       subject = "Pipeline completed",
                                       body = NULL,
                                       from = Sys.getenv("PIPELINE_NOTIFY_FROM", "noreply@localhost")) {
  if (!nzchar(to)) {
    message("No PIPELINE_NOTIFY_EMAIL set; skipping notification")
    return(invisible(FALSE))
  }

  if (is.null(body) || !nzchar(body)) {
    body <- paste0("Targets pipeline completed at ", Sys.time(), " on host ", Sys.info()[["nodename"]], "\nWorking directory: ", normalizePath("."))
  }

  # Try blastula first if installed and SMTP vars are present
  smtp_host <- Sys.getenv("SMTP_HOST")
  smtp_port <- as.integer(Sys.getenv("SMTP_PORT", "587"))
  smtp_user <- Sys.getenv("SMTP_USER")
  smtp_pass <- Sys.getenv("SMTP_PASSWORD")
  use_ssl <- identical(tolower(Sys.getenv("SMTP_USE_SSL", "true")), "true")

  if (requireNamespace("blastula", quietly = TRUE) && nzchar(smtp_host) && nzchar(smtp_user) && nzchar(smtp_pass)) {
    email <- tryCatch({
      blastula::compose_email(body = blastula::md(body))
    }, error = function(e) {
      warning("Failed to compose email with blastula: ", e$message)
      return(NULL)
    })

    if (!is.null(email)) {
      creds <- tryCatch({
        # create credentials on the fly
        blastula::creds(
          host = smtp_host,
          port = smtp_port,
          user = smtp_user,
          pass = smtp_pass,
          use_ssl = use_ssl
        )
      }, error = function(e) {
        warning("Failed to prepare SMTP credentials: ", e$message)
        return(NULL)
      })

      if (!is.null(creds)) {
        tryCatch({
          blastula::smtp_send(email = email, from = from, to = to, subject = subject, creds = creds)
          message("Notification sent to ", to, " via blastula SMTP")
          return(invisible(TRUE))
        }, error = function(e) {
          warning("blastula::smtp_send failed: ", e$message)
        })
      }
    }
  }

  # Fallback to system mail command if available
  mail_cmd <- Sys.which("mail")
  if (nzchar(mail_cmd)) {
    tmp <- tempfile()
    writeLines(body, tmp)
    cmd <- sprintf('%s -s %s %s < %s', shQuote(mail_cmd), shQuote(subject), shQuote(to), shQuote(tmp))
    rc <- system(cmd)
    unlink(tmp)
    if (rc == 0) {
      message("Notification sent to ", to, " via system mail")
      return(invisible(TRUE))
    } else {
      warning("system mail returned non-zero status: ", rc)
    }
  }

  warning("No available method to send notification (blastula SMTP not configured and 'mail' not found)")
  return(invisible(FALSE))
}
