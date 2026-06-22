# ACMGamp — ACMG/AMP Variant Classification Shiny Application
# Launch with: shiny::runApp("path/to/ACMGamp")

source("global.R", local = TRUE)
source("ui.R", local = TRUE)
source("server.R", local = TRUE)

shinyApp(ui = ui, server = server)
