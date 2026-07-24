#' Audit Log tab UI.
#' @noRd
audit_nav_panel <- function() {
  bslib::nav_panel(
    title = "Audit Log",
    value = "audit",
    icon = icon("clipboard-list"),
    div(
      class = "container-fluid py-3",
      bslib::card(
        bslib::card_header("Analysis Provenance Log"),
        bslib::card_body(
          p("Every classification verdict is timestamped and persisted to ",
            tags$code("logs/analysis_log.csv"), " for thesis-grade data provenance."),
          actionButton("refresh_audit", "Refresh Log", class = "btn-secondary mb-3"),
          DT::DTOutput("audit_table")
        )
      )
    )
  )
}
