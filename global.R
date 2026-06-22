# Global configuration and module loading for ACMGamp

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(data.table)
  library(readr)
})

# VariantAnnotation is optional at startup; loaded when VCF is parsed
.va_available <- requireNamespace("VariantAnnotation", quietly = TRUE)
if (!.va_available) {
  message(
    "VariantAnnotation not installed. VCF parsing will use a lightweight fallback. ",
    "Install with: BiocManager::install('VariantAnnotation')"
  )
}

APP_TITLE <- "ACMGamp"
APP_VERSION <- "1.0.0"

# Maximum upload size for VCF and companion files (1 GiB)
MAX_UPLOAD_SIZE_BYTES <- 1024 * 1024^2
options(shiny.maxRequestSize = MAX_UPLOAD_SIZE_BYTES)

# Reference data placeholders (gnomAD v4.1, ClinVar, REVEL)
REFERENCE_PATHS <- list(
  gnomad_v41 = file.path("data", "reference", "gnomad_v41_placeholder.tsv"),
  clinvar    = file.path("data", "reference", "clinvar_placeholder.tsv"),
  revel      = file.path("data", "reference", "revel_placeholder.tsv")
)

AUDIT_LOG_PATH <- file.path("logs", "analysis_log.csv")
REPORT_COLUMNS <- c(
  "variant_id", "chrom", "pos", "ref", "alt", "gene", "consequence",
  "gnomad_af", "revel_score", "clinvar_classification",
  "criteria_met", "criteria_strength", "classification",
  "pipeline_mode", "classified_at", "analyst_session"
)

# Source application modules (explicit order for dependencies)
module_files <- c(
  "R/classify_variant.R",
  "R/audit.R",
  "R/reference_data.R",
  "R/vcf_stream.R",
  "R/parse_inputs.R",
  "R/vcf_validate.R",
  "R/acmg_pipeline.R",
  "R/ui_helpers.R"
)
invisible(lapply(module_files, source, local = FALSE))
