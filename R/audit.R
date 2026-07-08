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

build_audit_entries_from_report <- function(report, session_id = NA_character_) {
  if (is.null(report) || nrow(report) == 0) return(NULL)
  do.call(rbind, lapply(seq_len(nrow(report)), function(i) {
    build_audit_entry(
      variant_id = report$variant_id[i],
      pipeline_mode = report$pipeline_mode[i],
      classification = report$classification[i],
      criteria_met = report$criteria_met[i],
      timestamp = report$classified_at[i],
      session_id = session_id,
      details = report$evidence_summary[i] %||% ""
    )
  }))
}

#' Append audit rows without re-reading the full log file (O(1) per batch).
append_audit_log <- function(entries, path = AUDIT_LOG_PATH) {
  if (is.null(entries) || nrow(entries) == 0) return(invisible(path))
  ensure_audit_log(path)
  utils::write.table(
    entries, file = path, sep = ",",
    row.names = FALSE, col.names = FALSE, append = TRUE, quote = TRUE
  )
  invisible(path)
}
