APP_TITLE <- "ClinicalVariantR"
APP_VERSION <- "0.99.3"
ACMG_PRO_ENGINE <- "ClinicalVariantR-Prediction-v0.99.3"
ACMG_GUIDELINE_VERSION <- "ACMG/AMP 2015 + ClinGen refinements"

ACMG_CLASSIFICATIONS <- c(
  "Pathogenic", "Likely Pathogenic", "VUS", "Likely Benign", "Benign"
)

REFERENCE_PATHS <- list(
  gnomad_v41 = file.path("data", "reference", "gnomad_v41_placeholder.tsv"),
  clinvar = file.path("data", "reference", "clinvar_placeholder.tsv"),
  revel = file.path("data", "reference", "revel_placeholder.tsv")
)

utils::globalVariables(c(
  "ACMG_CLASSIFICATIONS",
  "ACMG_GUIDELINE_VERSION",
  "ACMG_PRO_ENGINE",
  "APP_TITLE",
  "APP_VERSION",
  "AUTH_ENABLED",
  "CONFIG_PATHS",
  "DEFAULT_PROFILE_ID",
  "PREDICTION_MODE",
  "PREDICTION_SETTINGS",
  "REFERENCE_PATHS"
))

FULL_ACMG_CRITERIA <- c(
  "PVS1", "PS1", "PS2", "PS3", "PS4",
  "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
  "PP1", "PP2", "PP3", "PP4", "PP5",
  "BA1", "BS1", "BS2", "BS3", "BS4",
  "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"
)

AUTOMATED_ACMG_CRITERIA <- c(
  "PVS1", "PS1", "PS4",
  "PM1", "PM2", "PM4", "PM5",
  "PP2", "PP3", "PP5",
  "BA1", "BS1", "BS2",
  "BP1", "BP3", "BP4", "BP6", "BP7"
)

CONTEXT_ASSISTED_CRITERIA <- c("PS2", "PP1", "PP4")
MANUAL_ONLY_CRITERIA <- c("PS3", "PM3", "PM6", "BS3", "BS4", "BP2", "BP5")

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
