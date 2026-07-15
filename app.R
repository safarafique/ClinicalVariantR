# ClinicalVariantR: ACMG/AMP Variant Classification Shiny Application
# Launch with: shiny::runApp("path/to/ClinicalVariantR")

source("global.R", local = TRUE)
source("ui.R", local = TRUE)
source("server.R", local = TRUE)

# Ensure 1-hour request timeout is active for this process (large VCF / idle use).
options(shiny.http.timeout = if (exists("SESSION_IDLE_TIMEOUT_SEC")) SESSION_IDLE_TIMEOUT_SEC else 3600L)

shinyApp(ui = ui, server = server)


