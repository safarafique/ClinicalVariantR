#' Home / landing tab UI.
home_nav_panel <- function() {
  bslib::nav_panel(
    title = "Home",
    value = "home",
    icon = icon("home"),
    div(
      class = "container-fluid py-4",
      div(
        class = "text-center mb-4",
        h2("ClinicalVariantR - Variant Prediction Platform"),
        p(
          class = "lead text-muted",
          "ACMG/AMP evidence-based variant classification with streaming VCF analysis, ",
          "per-variant curator workflow, and structured prediction reports."
        ),
        tags$span(class = "badge bg-primary me-2", "Prediction mode"),
        tags$small(class = "text-muted", paste("Version", APP_VERSION)),
        uiOutput("reference_readiness_home")
      ),
      div(
        class = "row justify-content-center mb-4",
        div(class = "col-lg-10", expert_csv_handoff_ui())
      ),
      fluidRow(
        column(4, landing_card("group_a", "Group A - Clinical Prediction",
          "VCF + clinical + pedigree. All 28 criteria evaluated; per-variant curator overrides. Best for sign-out workflows.", "primary")),
        column(4, landing_card("group_b", "Group B - Automated Prediction",
          "VCF-only batch prediction with 18 automated criteria and full evidence per variant.", "success")),
        column(4, landing_card("group_c", "Group C - Gene Panel Prediction",
          "Same automated engine scoped to your gene panel - ideal for targeted NGS panels.", "warning"))
      )
    )
  )
}
