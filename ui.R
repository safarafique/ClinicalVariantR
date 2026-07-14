# ACMGamp — thin UI entry (panels in R/shiny/ui/*_ui.R)
ui <- tagList(
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
  bslib::page_navbar(
  title = tagList(icon("dna"), APP_TITLE),
  theme = clinical_theme,
  id = "main_nav",
    home_nav_panel(),
    group_a_nav_panel(),
    group_b_nav_panel(),
    group_c_nav_panel(),
    explorer_nav_panel(),
    audit_nav_panel(),
  bslib::nav_spacer(),
  bslib::nav_item(tags$span(class = "navbar-text me-2", textOutput("auth_label", inline = TRUE))),
  bslib::nav_item(tags$span(class = "navbar-text me-3", textOutput("session_label", inline = TRUE)))
  )
)
