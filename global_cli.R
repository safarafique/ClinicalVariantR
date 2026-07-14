# CLI bootstrap for ACMGamp scripts (no Shiny / DT / bslib).
# Usage in scripts: source("global_cli.R")

bootstrap_user_r_library <- function() {
  user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
  if (!nzchar(user_lib)) {
    rver <- paste(R.version$major, R.version$minor, sep = ".")
    user_lib <- file.path(
      Sys.getenv("HOME"),
      "R",
      paste0(R.version$platform, "-library"),
      rver
    )
  }
  if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
  invisible(user_lib)
}
bootstrap_user_r_library()

.va_available <- requireNamespace("VariantAnnotation", quietly = TRUE)
if (!.va_available) {
  message(
    "VariantAnnotation not installed — using lightweight VCF parser. ",
    "Optional: BiocManager::install('VariantAnnotation')"
  )
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  user_lib <- bootstrap_user_r_library()
  stop(
    "Package 'data.table' is required for CLI analysis.\n",
    "Install in WSL (one-time):\n",
    "  cd /mnt/e/ACGM/cml_variant_interpreter\n",
    "  Rscript scripts/install_r_cli_deps.R\n",
    "Or via apt (faster, if available):\n",
    "  sudo apt install -y r-cran-data.table r-cran-readr\n",
    if (nzchar(user_lib)) paste0("(Expected user library: ", user_lib, ")\n") else "",
    call. = FALSE
  )
}

APP_TITLE <- "ACMGamp"
APP_VERSION <- "2.7.0"
ACMG_PRO_ENGINE <- "ACMGamp-Prediction-v2.7.0"
ACMG_GUIDELINE_VERSION <- "ACMG/AMP 2015 + ClinGen refinements"

CONFIG_PATHS <- list(
  thresholds = file.path("config", "acmg_thresholds.csv"),
  criteria   = file.path("config", "acmg_criteria.csv"),
  profiles   = file.path("config", "disease_profiles.csv")
)

DEFAULT_PROFILE_ID <- "general_germline"

ACMG_CLASSIFICATIONS <- c(
  "Pathogenic", "Likely Pathogenic", "VUS", "Likely Benign", "Benign"
)

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

cli_module_files <- c(
  "R/aaa_constants.R",
  "R/classify_variant.R",
  "R/audit.R",
  "R/reference_data.R",
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
  "R/benchmark.R",
  "R/intervar_compare.R",
  "R/gene_filter.R",
  "R/expert_review_export.R"
)
invisible(lapply(cli_module_files, source, local = FALSE))
