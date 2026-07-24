# Global configuration and module loading for ClinicalVariantR
#
# Runtime dependencies are declared in DESCRIPTION Imports. After installing
# ClinicalVariantR with dependencies = TRUE, users only need:
#   library(ClinicalVariantR)
#   shiny::runApp(ClinicalVariantR())
# Do not install packages from app code (forbidden for Bioconductor packages).

.clinicalvariantr_require_deps <- function() {
  deps <- c("shiny", "bslib", "DT", "data.table", "readr", "jsonlite")
  missing <- deps[!vapply(deps, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0L) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "), ".\n",
      "Reinstall ClinicalVariantR with dependencies = TRUE ",
      "(BiocManager from Bioconductor, or remotes from GitHub/local clone).",
      call. = FALSE
    )
  }
  # Attach only shiny + bslib. Attaching DT/jsonlite with shiny masks
  # renderDataTable, dataTableOutput, and validate. Other Imports are used
  # via DT:: / jsonlite:: / data.table:: / readr:: in application code.
  suppressPackageStartupMessages({
    library(shiny, quietly = TRUE, warn.conflicts = FALSE)
    library(bslib, quietly = TRUE, warn.conflicts = FALSE)
  })
  invisible(TRUE)
}

.clinicalvariantr_require_deps()

# VariantAnnotation is an Import; confirm availability (fallback parser if absent).
.va_available <- requireNamespace("VariantAnnotation", quietly = TRUE)
if (!.va_available) {
  message(
    "VariantAnnotation not available. VCF parsing will use a lightweight fallback. ",
    "Reinstall ClinicalVariantR with dependencies = TRUE to restore full VCF support."
  )
}

APP_TITLE <- "ClinicalVariantR"
APP_VERSION <- "0.99.3"
ACMG_PRO_ENGINE <- "ClinicalVariantR-Prediction-v0.99.3"
ACMG_GUIDELINE_VERSION <- "ACMG/AMP 2015 + ClinGen refinements"

# Maximum upload size for VCF and companion files (1 GiB)
MAX_UPLOAD_SIZE_BYTES <- 1024 * 1024^2
# Keep idle sessions alive for at least 1 hour (local + hosted Shiny)
SESSION_IDLE_TIMEOUT_SEC <- as.integer(Sys.getenv("CLINICALVARIANTR_IDLE_TIMEOUT_SEC", unset = "3600"))
options(
  shiny.maxRequestSize = MAX_UPLOAD_SIZE_BYTES,
  shiny.http.timeout = SESSION_IDLE_TIMEOUT_SEC
)

CONFIG_PATHS <- list(
  thresholds = file.path("config", "acmg_thresholds.csv"),
  criteria   = file.path("config", "acmg_criteria.csv"),
  profiles   = file.path("config", "disease_profiles.csv")
)

DEFAULT_PROFILE_ID <- "general_germline"

ACMG_CLASSIFICATIONS <- c(
  "Pathogenic", "Likely Pathogenic", "VUS", "Likely Benign", "Benign"
)

# Reference data placeholders (gnomAD v4.1, ClinVar, REVEL)
# After `Rscript scripts/install_reference_data.R`, point these at production TSVs.
REFERENCE_PATHS <- list(
  gnomad_v41 = file.path("data", "reference", "gnomad_v41_placeholder.tsv"),
  clinvar    = file.path("data", "reference", "clinvar_placeholder.tsv"),
  revel      = file.path("data", "reference", "revel_placeholder.tsv")
)

AUDIT_LOG_PATH <- file.path("logs", "analysis_log.csv")
REPORT_COLUMNS <- c(
  "variant_id", "chrom", "pos", "ref", "alt", "gene", "consequence",
  "annotation_source", "genome_build_hint",
  "gnomad_af", "revel_score", "cadd_score", "spliceai_max", "alphamissense_score",
  "prediction_scores", "clinvar_classification",
  "criteria_met", "criteria_rationale", "evidence_summary",
  "classification", "confidence_score", "confidence_label",
  "evidence_strength", "pathogenic_evidence_count", "benign_evidence_count",
  "prediction_limitations",
  "disease_profile",
  "pipeline_mode", "classified_at", "analyst_session", "engine",
  "app_version", "acmg_guideline_version", "input_vcf_checksum"
)

EVIDENCE_DETAIL_COLUMNS <- c(
  "variant_id", "gene", "classification", "confidence_score", "confidence_label",
  "criteria_met", "evidence_summary", "prediction_scores", "gnomad_af",
  "disease_profile", "annotation_source", "genome_build_hint"
)

PDF_EXPORT_COLUMNS <- c(
  "variant_id", "gene", "classification", "confidence_score", "confidence_label",
  "criteria_met", "gnomad_af", "disease_profile", "evidence_summary"
)

VARIANT_DETAIL_COLUMNS <- c(
  "variant_id", "chrom", "pos", "gene", "consequence", "classification",
  "evidence_strength", "criteria_met", "confidence_score", "confidence_label",
  "prediction_scores", "prediction_limitations",
  "gnomad_af", "revel_score", "cadd_score", "clinvar_classification"
)

# Source application modules (explicit order for dependencies)
module_files <- c(
  "R/aaa_constants.R",
  "R/classify_variant.R",
  "R/audit.R",
  "R/reference_data.R",
  "R/variant_key.R",
  "R/vcf_unified_parser.R",
  "R/rule_config.R",
  "R/prediction_config.R",
  "R/evidence_report.R",
  "R/clinvar_pathogenic_criteria.R",
  "R/hpo_omim_phenotype.R",
  "R/trio_genotypes.R",
  "R/autopvs1.R",
  "R/ps4_case_control.R",
  "R/clinical_context_criteria.R",
  "R/acmg_vcf_criteria.R",
  "R/reproducibility.R",
  "R/acmg_engine.R",
  "R/variant_rescore.R",
  "R/vcf_stream.R",
  "R/parse_inputs.R",
  "R/vcf_validate.R",
  "R/acmg_pipeline.R",
  "R/benchmark.R",
  "R/intervar_compare.R",
  "R/gene_filter.R",
  "R/expert_review_export.R",
  "R/ui_helpers.R",
  "R/auth_storage.R"
)
invisible(lapply(module_files, source, local = FALSE))

shiny_server_modules <- c(
  "R/shiny/context.R",
  "R/shiny/shared_server.R",
  "R/shiny/upload_server.R",
  "R/shiny/analysis_server.R",
  "R/shiny/audit_server.R",
  "R/shiny/group_a_server.R",
  "R/shiny/group_b_server.R",
  "R/shiny/group_c_server.R",
  "R/shiny/results_server.R",
  "R/shiny/explorer_server.R",
  "R/shiny/auth_server.R"
)
invisible(lapply(shiny_server_modules, source, local = FALSE))

shiny_ui_modules <- c(
  "R/shiny/ui/theme.R",
  "R/shiny/ui/home_ui.R",
  "R/shiny/ui/group_a_ui.R",
  "R/shiny/ui/group_b_ui.R",
  "R/shiny/ui/group_c_ui.R",
  "R/shiny/ui/explorer_ui.R",
  "R/shiny/ui/audit_ui.R"
)
invisible(lapply(shiny_ui_modules, source, local = FALSE))
