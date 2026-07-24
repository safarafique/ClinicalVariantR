#' Group A - clinical prediction tab UI.
#' @noRd
group_a_nav_panel <- function() {
  bslib::nav_panel(
    title = "Group A - Clinical Prediction",
    value = "group_a",
    icon = icon("layer-group"),
    div(
      class = "container-fluid py-3",
      fluidRow(
        column(
          width = 4,
          bslib::card(
            bslib::card_header("Required Inputs"),
            bslib::card_body(
              fileInput("vcf_a", "VCF File (max 1 GB)", accept = c(".vcf", ".vcf.gz", ".bcf")),
              fileInput("clinical_a", "Clinical Logs (CSV)", accept = ".csv"),
              fileInput("pedigree_a", "Pedigree Data (CSV)", accept = ".csv"),
              selectInput("profile_a", "Disease profile", choices = profile_choices(), selected = DEFAULT_PROFILE_ID),
              disease_profile_help_ui(),
              uiOutput("readiness_a"),
              hr(),
              h5("Complete VCF Analysis"),
              checkboxInput("complete_vcf_a", "Analyze entire VCF (no row limit)", TRUE),
              checkboxInput("pass_only_a", "Load passing-filter rows only (PASS or .)", FALSE),
              numericInput("min_qual_a", "Minimum QUAL", value = 0, min = 0, step = 1),
              checkboxInput("use_bcftools_a", "Use bcftools (Ubuntu/WSL - faster)", bcftools_available()),
              checkboxInput("skip_audit_a", "Skip audit log (faster analysis)", TRUE),
              numericInput("chunk_size_a", "Chunk size (variants per batch)", value = 10000, min = 1000, step = 1000),
              helpText(textOutput("engine_status_a", inline = TRUE)),
              p(class = "text-muted small",
                "Prediction mode uses stricter PVS1 and ClinVar rules. ",
                "After analysis, select a variant and apply curator evidence in the evidence panel."),
              uiOutput("run_a_ui")
            )
          )
        ),
        column(
          width = 8,
          bslib::card(
            bslib::card_header("VCF Requirement Check"),
            bslib::card_body(
              uiOutput("validation_status_a"),
              DT::DTOutput("validation_a")
            )
          ),
          bslib::card(
            class = "mt-3",
            bslib::card_header("VCF Preview"),
            bslib::card_body(
              uiOutput("preview_status_a"),
              DT::DTOutput("preview_a")
            )
          ),
          bslib::card(
            class = "mt-3",
            bslib::card_header("Classification Results"),
            bslib::card_body(
              uiOutput("results_placeholder_a"),
              uiOutput("status_a"),
              uiOutput("expert_handoff_a"),
              h5(class = "mt-2", "Classification summary"),
              p(class = "text-muted small", "Click a row to view variants in that ACMG category."),
              uiOutput("category_hint_a"),
              DT::DTOutput("classification_summary_a"),
              hr(),
              h5("Variant details"),
              DT::DTOutput("results_a"),
              br(),
              bslib::accordion(
                id = "evidence_accordion_a",
                open = FALSE,
                bslib::accordion_panel(
                  "Variant evidence detail",
                  p(class = "text-muted small", "Select a variant row, review 28 criteria, then apply per-variant curator evidence below."),
                  uiOutput("selected_variant_a"),
                  div(
                    class = "border rounded p-3 mb-3 bg-light",
                    h6("Per-variant curator evidence"),
                    p(class = "text-muted small", "Select a variant row above, set criteria, then apply to reclassify that variant only."),
                    checkboxInput("cur_ps3", "PS3 - Functional studies support damaging effect", FALSE),
                    checkboxInput("cur_pp4", "PP4 - Phenotype matches gene/disease", FALSE),
                    checkboxInput("cur_ps4", "PS4 - Case-control enrichment", FALSE),
                    checkboxInput("cur_ps2", "PS2 - De novo (confirmed)", FALSE),
                    checkboxInput("cur_pm6", "PM6 - De novo (assumed)", FALSE),
                    checkboxInput("cur_pp1", "PP1 - Co-segregation", FALSE),
                    checkboxInput("cur_pp2", "PP2 - Missense mechanism", FALSE),
                    actionButton("apply_curation_a", "Apply curation & reclassify", class = "btn-primary btn-sm"),
                    uiOutput("curation_status_a")
                  ),
                  DT::DTOutput("evidence_detail_a")
                ),
                bslib::accordion_panel(
                  "Expert review checklist",
                  expert_review_checklist_static_ui(show_curation = TRUE),
                  uiOutput("expert_review_stats_a"),
                  div(
                    class = "d-flex flex-wrap gap-2 mt-2",
                    downloadButton(
                      "download_a_lp_plus",
                      "Export LP+ only (CSV)",
                      class = "btn-outline-danger btn-sm",
                      icon = icon("file-csv")
                    ),
                    downloadButton(
                      "download_a_worklist",
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
                  uiOutput("repro_a")
                )
              ),
              br(),
              downloadButton("download_a", "Download full prediction report (CSV)",
                             class = "btn-outline-primary", icon = icon("download"))
            )
          )
        )
      )
    )
  )
}
