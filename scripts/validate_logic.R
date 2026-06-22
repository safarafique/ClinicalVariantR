# Validation script — run before deploying on patient data
# Usage: source("scripts/validate_logic.R")

source("global.R")

cat("=== CML Variant Interpreter — Logic Validation ===\n\n")

# 1. Self-test with dummy variant row (thesis validation command)
cat("1. Full 28-criteria engine test:\n")
test_result <- test_acmg_logic_engine()
print(test_result)
stopifnot("Logic engine must classify high-impact variant as Pathogenic/LP" = test_result$pass)

# 2. End-to-end pipeline with sample VCF
cat("\n2. Sample VCF pipeline test (Group B):\n")
sample_vcf <- file.path("data", "samples", "example_variants.vcf")
variants <- parse_vcf_upload(sample_vcf)
refs <- load_reference_data()
annotated <- annotate_variants(variants, refs)
pipeline_result <- run_pipeline(annotated, mode = "rapid", session_id = "VALIDATION")
print(pipeline_result$report[, c("variant_id", "classification", "criteria_met")])

cat("\n3. Audit log written to:", AUDIT_LOG_PATH, "\n")
cat("\nAll validation checks passed.\n")
