render_vcf_preview_table <- function(preview_df) {
  DT::datatable(
    preview_df,
    options = list(pageLength = 10, scrollX = TRUE, dom = "tip"),
    rownames = FALSE,
    class = "cell-border stripe hover compact"
  )
}

render_validation_table <- function(checks_df) {
  if (is.null(checks_df) || nrow(checks_df) == 0) {
    return(DT::datatable(
      data.frame(Message = "No validation results."),
      rownames = FALSE, options = list(dom = "t")
    ))
  }

  DT::datatable(
    checks_df,
    options = list(pageLength = 15, scrollX = TRUE, dom = "tip"),
    rownames = FALSE,
    class = "cell-border stripe hover compact"
  ) |>
    DT::formatStyle(
      "status",
      backgroundColor = DT::styleEqual(
        c("PASS", "WARN", "FAIL"),
        c("#d4edda", "#fff3cd", "#f8d7da")
      ),
      fontWeight = "bold"
    )
}

readiness_indicator_ui <- function(validation, label = "Analysis readiness") {
  if (is.null(validation)) {
    return(div(
      class = "readiness-badge readiness-neutral",
      tags$span(class = "readiness-dot readiness-dot-grey"),
      tags$span(class = "readiness-text", icon("hourglass-half"), " Upload files to check readiness")
    ))
  }

  ready <- isTRUE(validation$can_analyze)
  badge_class <- if (ready) "readiness-ready" else "readiness-not-ready"
  dot_class <- if (ready) "readiness-dot-green" else "readiness-dot-red"
  icon_name <- if (ready) "check-circle" else "times-circle"
  status_text <- if (ready) {
    "Ready for analysis — click Run Analysis"
  } else {
    "Not ready for analysis"
  }

  missing_block <- NULL
  missing_items <- validation$missing_items %||% character()
  if (!ready && length(missing_items) > 0) {
    missing_block <- tags$div(
      class = "readiness-missing mt-2",
      tags$strong("Missing:"),
      tags$ul(class = "mb-0 ps-3", lapply(missing_items, tags$li))
    )
  } else if (!ready) {
    missing_block <- tags$div(
      class = "readiness-missing mt-2",
      tags$strong(validation$summary %||% "Fix validation failures before analysis.")
    )
  }

  div(
    class = paste("readiness-badge", badge_class),
    div(
      class = "d-flex align-items-center gap-2",
      tags$span(class = paste("readiness-dot", dot_class)),
      div(
        div(class = "readiness-label", label),
        div(class = "readiness-text", icon(icon_name), status_text)
      )
    ),
    missing_block
  )
}

validation_summary_ui <- function(validation, file_label = "VCF") {
  if (is.null(validation)) {
    return(div(class = "alert alert-secondary", icon("clipboard-check"), " Upload a file to run requirement checks."))
  }
  readiness_indicator_ui(validation, label = paste(file_label, "requirement check"))
}

render_results_table <- function(report_df) {
  display_cols <- intersect(REPORT_COLUMNS, names(report_df))
  dat <- report_df[, display_cols, drop = FALSE]

  DT::datatable(
    dat,
    options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip"),
    rownames = FALSE,
    class = "cell-border stripe hover compact"
  ) |>
    DT::formatStyle(
      "classification",
      backgroundColor = DT::styleEqual(
        c("Pathogenic", "Likely Pathogenic", "VUS", "Benign", "Likely Benign"),
        c("#dc3545", "#fd7e14", "#ffc107", "#28a745", "#6cbf6c")
      ),
      color = DT::styleEqual(
        c("Pathogenic", "Likely Pathogenic", "VUS", "Benign", "Likely Benign"),
        c("#ffffff", "#ffffff", "#212529", "#ffffff", "#ffffff")
      ),
      fontWeight = "bold"
    ) |>
    DT::formatRound(columns = c("gnomad_af", "revel_score"), digits = 4)
}
