ensure_audit_log <- function(path = AUDIT_LOG_PATH) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) {
    cols <- c(
      "log_id", "timestamp", "variant_id", "pipeline_mode",
      "classification", "criteria_met", "analyst_session", "details"
    )
    empty <- setNames(
      data.frame(matrix(nrow = 0, ncol = length(cols)), stringsAsFactors = FALSE),
      cols
    )
    write.csv(empty, path, row.names = FALSE)
  }
  invisible(path)
}

build_audit_entry <- function(variant_id, pipeline_mode, classification, criteria_met,
                              timestamp, session_id, details = "") {
  data.frame(
    log_id = paste0("LOG-", format(Sys.time(), "%Y%m%d%H%M%S"), "-", sample.int(1e6, 1)),
    timestamp = timestamp,
    variant_id = variant_id,
    pipeline_mode = pipeline_mode,
    classification = classification,
    criteria_met = criteria_met,
    analyst_session = session_id,
    details = details,
    stringsAsFactors = FALSE
  )
}

append_audit_log <- function(entries, path = AUDIT_LOG_PATH) {
  if (is.null(entries) || nrow(entries) == 0) return(invisible(path))
  ensure_audit_log(path)
  existing <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  combined <- if (is.null(existing) || nrow(existing) == 0) entries else rbind(existing, entries)
  write.csv(combined, path, row.names = FALSE)
  invisible(path)
}
