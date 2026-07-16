#' Evidence Explorer tab UI.
explorer_nav_panel <- function() {
  bslib::nav_panel(
    title = "Evidence Explorer",
    value = "evidence",
    icon = icon("microscope"),
    div(
      class = "container-fluid py-3",
      bslib::card(
        bslib::card_header("Interactive ACMG Evidence Report"),
        bslib::card_body(
          p("Transparent criterion-level explanations with observed values, thresholds, strength, and confidence."),
          fluidRow(
            column(3, selectInput("explorer_source", "Results source",
              choices = c("Group A" = "a", "Group B" = "b", "Group C — Gene panel" = "c"))),
            column(3, selectInput("explorer_pdf_scope", "Export scope",
              choices = c(
                "All variants" = "all",
                "Current search filter" = "filtered",
                "Summary only (1 page)" = "summary"
              ),
              selected = "all")),
            column(4, uiOutput("explorer_status")),
            column(2,
              div(class = "text-end", style = "padding-top: 26px;",
                downloadButton(
                  "download_explorer_pdf",
                  "Save PDF",
                  class = "btn-outline-primary btn-sm mb-1",
                  icon = icon("file-pdf")
                ),
                tags$br(),
                downloadButton(
                  "download_explorer_csv",
                  "Save CSV",
                  class = "btn-outline-secondary btn-sm",
                  icon = icon("file-csv")
                )
              )
            )
          ),
          DT::DTOutput("explorer_results"),
          hr(),
          uiOutput("explorer_selected_variant"),
          DT::DTOutput("explorer_evidence")
        )
      )
    )
  )
}
