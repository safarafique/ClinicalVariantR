# Verify Group A exports all 28 ACMG criteria when clinical + pedigree files are supplied.
# Usage (WSL): cd /mnt/e/ACGM/cml_variant_interpreter && Rscript scripts/verify_group_a_28.R

source("global_cli.R")
source("R/parse_inputs.R")

sample_vcf <- file.path("data", "samples", "example_variants.vcf")
clinical <- file.path("data", "samples", "example_clinical_logs.csv")
pedigree <- file.path("data", "samples", "example_pedigree.csv")

stopifnot(file.exists(sample_vcf), file.exists(clinical), file.exists(pedigree))

variants <- parse_vcf_upload(sample_vcf)
refs <- tryCatch(load_reference_data(), error = function(e) NULL)
if (!is.null(refs)) variants <- annotate_variants(variants, refs)

clinical_ctx <- parse_clinical_logs(clinical)
pedigree_ctx <- parse_pedigree(pedigree)

manual <- list(
  PS3_functional = TRUE,
  PM6_de_novo = FALSE,
  PS2_de_novo = FALSE,
  PS4_case_control = FALSE,
  PP1_segregation = FALSE,
  PP4_phenotype = FALSE,
  PP2_missense_mechanism = FALSE
)

scored <- score_variants_table(
  variants,
  manual_inputs = manual,
  clinical_context = clinical_ctx,
  pedigree_context = pedigree_ctx,
  evidence_scope = "full"
)

row <- scored[1, , drop = FALSE]
evidence <- parse_evidence_json(row$evidence_json[1])

cat("Group A 28-criteria verification\n")
cat("Variants scored:", nrow(scored), "\n")
cat("Evidence rows:", nrow(evidence), "\n")
cat("Expected criteria:", length(FULL_ACMG_CRITERIA), "\n")

missing <- setdiff(FULL_ACMG_CRITERIA, evidence$criterion)
extra <- setdiff(evidence$criterion, FULL_ACMG_CRITERIA)

if (length(missing) > 0L) {
  stop("Missing criteria in evidence output: ", paste(missing, collapse = ", "))
}
if (length(extra) > 0L) {
  stop("Unexpected criteria in evidence output: ", paste(extra, collapse = ", "))
}
if (nrow(evidence) != length(FULL_ACMG_CRITERIA)) {
  stop("Expected ", length(FULL_ACMG_CRITERIA), " evidence rows, got ", nrow(evidence))
}

ps3_row <- evidence[evidence$criterion == "PS3", , drop = FALSE]
if (!isTRUE(ps3_row$triggered[1])) {
  stop("Manual PS3 should be triggered in test run.")
}

cat("Criteria met:", row$criteria_met[1], "\n")
cat("Classification:", row$classification[1], "\n")
cat("All 28 criteria present in Group A evidence output.\n")
