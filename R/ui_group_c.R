#' Group C - gene panel prediction tab UI.
group_c_nav_panel <- function() {
  bslib::nav_panel(
    title = "Group C - Gene Panel Prediction",
    value = "group_c",
    icon = icon("filter"),
    div(
      class = "container-fluid py-3",
      fluidRow(
        column(
          width = 4,
          bslib::card(
            bslib::card_header("Gene-targeted analysis"),
            bslib::card_body(
              fileInput("vcf_c", "Patient VCF (max 1 GB)", accept = c(".vcf", ".vcf.gz", ".bcf")),
              textAreaInput(
                "genes_c",
                "Gene symbol(s)",
                placeholder = "ABL1, BCR, MYC, ENO1\n(HGNC symbols; BCR-ABL expands to ABL1 + BCR)",
                rows = 4,
                resize = "vertical"
              ),
              helpText("Matches any VEP/SnpEff annotation for these genes across the full VCF. Preview shows first 50 rows only; run analysis for the complete gene panel."),
              selectInput("profile_c", "Disease profile", choices = profile_choices(), selected = DEFAULT_PROFILE_ID),
              disease_profile_help_ui(),
              uiOutput("readiness_c"),
              hr(),
              h5("Analysis options"),
              checkboxInput("complete_vcf_c", "Analyze entire VCF (no row limit)", TRUE),
              checkboxInput("pass_only_c", "Load passing-filter rows only (PASS or .)", FALSE),
              numericInput("min_qual_c", "Minimum QUAL", value = 0, min = 0, step = 1),
              checkboxInput("use_bcftools_c", "Use bcftools (Ubuntu/WSL - faster)", bcftools_available()),
              checkboxInput("skip_audit_c", "Skip audit log (faster analysis)", FALSE),
              numericInput("chunk_size_c", "Chunk size (variants per batch)", value = 10000, min = 1000, step = 1000),
              helpText(textOutput("engine_status_c", inline = TRUE)),
              p(class = "text-muted small",
                "Runs the same automated ACMG criteria as Group B, limited to your gene list."),
              uiOutput("run_c_ui")
            )
          )
        ),
        column(
          width = 8,
          bslib::card(
            bslib::card_header("VCF Requirement Check"),
            bslib::card_body(
              uiOutput("validation_status_c"),
              DT::DTOutput("validation_c")
            )
          ),
          bslib::card(
            class = "mt-3",
            bslib::card_header("VCF Preview"),
            bslib::card_body(
              uiOutput("preview_status_c"),
              uiOutput("gene_preview_hint_c"),
              DT::DTOutput("preview_c")
            )
          ),
          bslib::card(
            class = "mt-3",
            bslib::card_header("Gene panel results"),
            bslib::card_body(
              uiOutput("results_placeholder_c"),
              uiOutput("status_c"),
              uiOutput("expert_handoff_c"),
              h5(class = "mt-2", "Classification summary"),
              p(class = "text-muted small", "Click a row to view variants in that ACMG category."),
              uiOutput("category_hint_c"),
              DT::DTOutput("classification_summary_c"),
              hr(),
              h5("Variant details"),
              DT::DTOutput("results_c"),
              br(),
              bslib::accordion(
                id = "evidence_accordion_c",
                open = FALSE,
                bslib::accordion_panel(
                  "Variant evidence detail",
                  p(class = "text-muted small", "Click a row in the results table to inspect criterion-level evidence."),
                  uiOutput("selected_variant_c"),
                  DT::DTOutput("evidence_detail_c")
                ),
                bslib::accordion_panel(
                  "Expert review checklist",
                  expert_review_checklist_static_ui(show_curation = FALSE),
                  uiOutput("expert_review_stats_c"),
                  div(
                    class = "d-flex flex-wrap gap-2 mt-2",
                    downloadButton(
                      "download_c_lp_plus",
                      "Export LP+ only (CSV)",
                      class = "btn-outline-danger btn-sm",
                      icon = icon("file-csv")
                    ),
                    downloadButton(
                      "download_c_worklist",
                      "Export expert worklist (CSV)",
                      class = "btn-outline-warning btn-sm",
                      icon = icon("clipboard-check")
                    )
                  ),
                  p(
                    class = "text-muted small mt-2 mb-0",
                    "Worklist = Pathogenic/Likely Pathogenic plus VUS with >=2 pathogenic evidence; excludes PM2-only VUS."
                  )
                ),
                bslib::accordion_panel(
                  "Reproducibility metadata",
                  uiOutput("repro_c")
                )
              ),
              br(),
              downloadButton("download_c", "Download full gene panel report (CSV)",
                             class = "btn-outline-warning", icon = icon("download"))
            )
          )
        )
      )
    )
  )
}
