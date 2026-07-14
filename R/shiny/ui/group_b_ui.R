#' Group B — automated prediction tab UI.
group_b_nav_panel <- function() {
  bslib::nav_panel(
    title = "Group B — Automated Prediction",
    value = "group_b",
    icon = icon("bolt"),
    div(
      class = "container-fluid py-3",
      fluidRow(
        column(
          width = 4,
          bslib::card(
            bslib::card_header("Required Input"),
            bslib::card_body(
              fileInput("vcf_b", "VCF File (max 1 GB)", accept = c(".vcf", ".vcf.gz", ".bcf")),
              selectInput("profile_b", "Disease profile", choices = profile_choices(), selected = DEFAULT_PROFILE_ID),
              disease_profile_help_ui(),
              uiOutput("readiness_b"),
              div(
                class = "run-action-panel run-action-ready mt-2",
                actionButton(
                  "run_b",
                  "Run Analysis",
                  class = "btn-success btn-lg w-100",
                  icon = icon("play")
                ),
                tags$p(class = "text-muted small mt-2 mb-0", "Click to start automated ACMG classification.")
              ),
              uiOutput("run_b_state"),
              hr(),
              h5("Complete VCF Analysis"),
              checkboxInput("complete_vcf_b", "Analyze entire VCF (no row limit)", TRUE),
              checkboxInput("pass_only_b", "Load passing-filter rows only (PASS or .)", FALSE),
              numericInput("min_qual_b", "Minimum QUAL", value = 0, min = 0, step = 1),
              checkboxInput("use_bcftools_b", "Use bcftools (Ubuntu/WSL — faster)", bcftools_available()),
              checkboxInput("skip_audit_b", "Skip audit log (faster analysis)", FALSE),
              numericInput("chunk_size_b", "Chunk size (variants per batch)", value = 10000, min = 1000, step = 1000),
              helpText(textOutput("engine_status_b", inline = TRUE)),
              p(class = "text-muted small",
                "Runs ACMGamp Pro automated criteria from VEP CSQ or SnpEff ANN fields."),
              p(class = "text-muted small",
                tags$strong("Prediction v2.7:"),
                " sensitivity rules upgrade strong VUS to Likely Pathogenic when PVS1/PM/PP evidence supports it.",
                " Use the classification summary table after run for Pathogenic / Likely Pathogenic counts.")
            )
          )
        ),
        column(
          width = 8,
          bslib::card(
            bslib::card_header("VCF Requirement Check"),
            bslib::card_body(
              uiOutput("validation_status_b"),
              DT::DTOutput("validation_b")
            )
          ),
          bslib::card(
            class = "mt-3",
            bslib::card_header("VCF Preview"),
            bslib::card_body(
              uiOutput("preview_status_b"),
              DT::DTOutput("preview_b")
            )
          ),
          bslib::card(
            class = "mt-3",
            bslib::card_header("Analysis Results"),
            bslib::card_body(
              uiOutput("results_placeholder_b"),
              uiOutput("status_b"),
              uiOutput("expert_handoff_b"),
              h5(class = "mt-2", "Classification summary"),
              p(class = "text-muted small", "Click a row to view variants in that ACMG category."),
              uiOutput("category_hint_b"),
              DT::DTOutput("classification_summary_b"),
              hr(),
              h5("Variant details"),
              DT::DTOutput("results_b"),
              br(),
              bslib::accordion(
                id = "evidence_accordion_b",
                open = FALSE,
                bslib::accordion_panel(
                  "Variant evidence detail",
                  p(class = "text-muted small", "Click a row in the results table to inspect criterion-level evidence."),
                  uiOutput("selected_variant_b"),
                  DT::DTOutput("evidence_detail_b")
                ),
                bslib::accordion_panel(
                  "Expert review checklist",
                  expert_review_checklist_static_ui(show_curation = FALSE),
                  uiOutput("expert_review_stats_b"),
                  div(
                    class = "d-flex flex-wrap gap-2 mt-2",
                    downloadButton(
                      "download_b_lp_plus",
                      "Export LP+ only (CSV)",
                      class = "btn-outline-danger btn-sm",
                      icon = icon("file-csv")
                    ),
                    downloadButton(
                      "download_b_worklist",
                      "Export expert worklist (CSV)",
                      class = "btn-outline-warning btn-sm",
                      icon = icon("clipboard-check")
                    )
                  ),
                  p(
                    class = "text-muted small mt-2 mb-0",
                    "Manual criteria (PS2, PP4, etc.) must be reviewed off-platform or in Group A."
                  )
                ),
                bslib::accordion_panel(
                  "Reproducibility metadata",
                  uiOutput("repro_b")
                )
              ),
              br(),
              downloadButton("download_b", "Download full prediction report (CSV)",
                             class = "btn-outline-success", icon = icon("download"))
            )
          )
        )
      )
    )
  )
}
