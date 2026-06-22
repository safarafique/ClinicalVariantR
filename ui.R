# CML Clinical Variant Interpreter — User Interface

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

ui <- tagList(
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
  bslib::page_navbar(
  title = tagList(icon("dna"), APP_TITLE),
  theme = clinical_theme,
  id = "main_nav",

  bslib::nav_panel(
    title = "Home",
    value = "home",
    icon = icon("home"),
    div(
      class = "container-fluid py-4",
      div(
        class = "text-center mb-4",
        h2("CML Clinical Variant Interpretation Platform"),
        p(
          class = "lead text-muted",
          "ACMG/AMP-compliant classification for Chronic Myeloid Leukemia research. ",
          "Select a pipeline based on available clinical data."
        ),
        tags$small(class = "text-muted", paste("Version", APP_VERSION))
      ),
      fluidRow(
        column(6, landing_card("group_a", "Group A — Full Pipeline", 
          "Upload VCF, clinical logs, and pedigree data for 28-criteria ACMG/AMP classification with manual curation fields.", "primary")),
        column(6, landing_card("group_b", "Group B — Rapid Pipeline",
          "Upload VCF only for 18 automated ACMG criteria and pathogenicity filtering.", "success"))
      )
    )
  ),

  bslib::nav_panel(
    title = "Group A — Full Pipeline",
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
              uiOutput("readiness_a"),
              hr(),
              h5("Complete VCF Analysis"),
              checkboxInput("complete_vcf_a", "Analyze entire VCF (no row limit)", TRUE),
              checkboxInput("pass_only_a", "Load FILTER=PASS rows only", FALSE),
              numericInput("min_qual_a", "Minimum QUAL", value = 0, min = 0, step = 1),
              checkboxInput("use_bcftools_a", "Use bcftools (Ubuntu/WSL — faster)", bcftools_available()),
              numericInput("chunk_size_a", "Chunk size (variants per batch)", value = 10000, min = 1000, step = 1000),
              helpText(textOutput("engine_status_a", inline = TRUE)),
              hr(),
              h5("Manual ACMG Criteria (Curation)"),
              p(class = "text-muted small",
                "Automated criteria (BA1, PM2, PP3, etc.) run automatically. ",
                "Toggle manual criteria after clinical review."),
              checkboxInput("manual_ps3", "PS3 — Functional studies support damaging effect", FALSE),
              checkboxInput("manual_pp4", "PP4 — Patient phenotype matches gene/disease", FALSE),
              checkboxInput("manual_ps4", "PS4 — Case-control enrichment", FALSE),
              checkboxInput("manual_ps2", "PS2 — De novo (confirmed)", FALSE),
              checkboxInput("manual_pm6", "PM6 — De novo (assumed without confirmation)", FALSE),
              checkboxInput("manual_pp1", "PP1 — Co-segregation with disease", FALSE),
              checkboxInput("manual_pp2", "PP2 — Missense is common disease mechanism", FALSE),
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
              DT::DTOutput("results_a"),
              br(),
              downloadButton("download_a", "Download Final_Clinical_Report.csv",
                             class = "btn-outline-primary", icon = icon("download"))
            )
          )
        )
      )
    )
  ),

  bslib::nav_panel(
    title = "Group B — Rapid Pipeline",
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
              uiOutput("readiness_b"),
              hr(),
              h5("Complete VCF Analysis"),
              checkboxInput("complete_vcf_b", "Analyze entire VCF (no row limit)", TRUE),
              checkboxInput("pass_only_b", "Load FILTER=PASS rows only", FALSE),
              numericInput("min_qual_b", "Minimum QUAL", value = 0, min = 0, step = 1),
              checkboxInput("use_bcftools_b", "Use bcftools (Ubuntu/WSL — faster)", bcftools_available()),
              numericInput("chunk_size_b", "Chunk size (variants per batch)", value = 10000, min = 1000, step = 1000),
              helpText(textOutput("engine_status_b", inline = TRUE)),
              p(class = "text-muted small",
                "Runs 18 automated ACMG criteria using gnomAD v4.1, ClinVar, and REVEL placeholders."),
              checkboxInput("filter_pathogenic_b", "Show only Pathogenic / Likely Pathogenic", TRUE),
              uiOutput("run_b_ui")
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
              DT::DTOutput("results_b"),
              br(),
              downloadButton("download_b", "Download Final_Clinical_Report.csv",
                             class = "btn-outline-success", icon = icon("download"))
            )
          )
        )
      )
    )
  ),

  bslib::nav_panel(
    title = "Audit Log",
    value = "audit",
    icon = icon("clipboard-list"),
    div(
      class = "container-fluid py-3",
      bslib::card(
        bslib::card_header("Analysis Provenance Log"),
        bslib::card_body(
          p("Every classification verdict is timestamped and persisted to ",
            tags$code("logs/analysis_log.csv"), " for thesis-grade data provenance."),
          actionButton("refresh_audit", "Refresh Log", class = "btn-secondary mb-3"),
          DT::DTOutput("audit_table")
        )
      )
    )
  ),

  bslib::nav_spacer(),
  bslib::nav_item(tags$span(class = "navbar-text me-3", textOutput("session_label", inline = TRUE)))
  )
)
