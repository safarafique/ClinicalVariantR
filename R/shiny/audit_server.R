#' Audit log tab and shared load_audit helper.
register_audit_server <- function(ctx) {
  input <- ctx$input
  output <- ctx$output
  audit_data <- ctx$audit_data

  ctx$load_audit <- function() {
    path <- AUDIT_LOG_PATH
    audit_data(if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else data.frame())
  }

  observe({ ctx$load_audit() })
  observeEvent(input$refresh_audit, { ctx$load_audit() })

  output$audit_table <- DT::renderDT({
    DT::datatable(
      audit_data(),
      options = list(pageLength = 20, scrollX = TRUE, order = list(list(1, "desc"))),
      rownames = FALSE
    )
  })

  invisible(ctx)
}
