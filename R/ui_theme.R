#' Shiny UI theme and reusable layout components.
#' @noRd
clinical_theme <- bslib::bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#2c3e50",
  success = "#18bc9c",
  warning = "#f39c12",
  danger = "#e74c3c",
  base_font = bslib::font_google("Source Sans 3"),
  heading_font = bslib::font_google("Source Sans 3")
)

landing_card <- function(id, title, description, color) {
  bslib::card(
    class = "pipeline-card h-100",
    bslib::card_header(class = paste0("bg-", color, " text-white"), title),
    bslib::card_body(
      p(description),
      actionButton(
        inputId = paste0("select_", id),
        label = paste("Enter", title),
        class = paste0("btn btn-", color, " w-100")
      )
    )
  )
}
