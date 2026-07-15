APP_TITLE <- "ClinicalVariantR"
APP_VERSION <- "2.7.0"
ACMG_PRO_ENGINE <- "ClinicalVariantR-Prediction-v2.7.0"
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
