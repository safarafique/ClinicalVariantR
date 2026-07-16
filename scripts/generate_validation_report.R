#!/usr/bin/env Rscript
# Formal validation report: score expert-classified benchmark variants and emit metrics.
# Usage (from ClinicalVariantR/):
#   Rscript scripts/generate_validation_report.R
#   Rscript scripts/generate_validation_report.R --vcf ../testig/clinicalvariantr_benchmark/clinicalvariantr_group_b_benchmark.vcf

args <- commandArgs(trailingOnly = TRUE)
vcf_path <- if (length(args) >= 2L && args[[1L]] == "--vcf") {
  args[[2L]]
} else {
  file.path("..", "testig", "clinicalvariantr_benchmark", "clinicalvariantr_group_b_benchmark.vcf")
}
validation_tsv <- file.path("data", "validation", "clinicalvariantr_validation_set.tsv")
output_dir <- file.path("results", "validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source("global_cli.R")

if (!file.exists(vcf_path)) {
  stop("Benchmark VCF not found: ", vcf_path, call. = FALSE)
}
if (!file.exists(validation_tsv)) {
  stop("Validation set not found: ", validation_tsv, call. = FALSE)
}

ref <- utils::read.delim(validation_tsv, stringsAsFactors = FALSE, comment.char = "")
ref$variant_key <- variant_key_from_parts(ref$chr, ref$pos, ref$ref, ref$alt)
ref$ref_class <- normalize_acmg_class(ref$acmg_classification_base)
ref$ref_tier <- collapse_to_tier(ref$ref_class)

variants <- parse_testing_vcf(vcf_path, pass_only = TRUE)
if (nrow(variants) == 0L) stop("No variants parsed from VCF")

scored <- score_variants_table(variants, profile_id = DEFAULT_PROFILE_ID, evidence_scope = "automated")
scored$variant_key <- variant_key_from_parts(scored$chrom, scored$pos, scored$ref, scored$alt)
scored$pred_class <- normalize_acmg_class(scored$classification)
scored$pred_tier <- collapse_to_tier(scored$pred_class)

cmp <- merge(
  scored[, c("variant_key", "gene", "consequence", "pred_class", "pred_tier",
             "classification", "criteria_met", "confidence_score", "evidence_strength")],
  ref[, c("variant_key", "benchmark_id", "benchmark_category", "ref_class", "ref_tier",
          "acmg_criteria_base", "benchmark_note", "expert_source")],
  by = "variant_key",
  all.y = TRUE
)

cmp$exact_match <- cmp$pred_class == cmp$ref_class
cmp$tier_match <- cmp$pred_tier == cmp$ref_tier
cmp$missed <- is.na(cmp$classification)

n_compared <- sum(!cmp$missed)
exact_acc <- mean(cmp$exact_match, na.rm = TRUE)
tier_acc <- mean(cmp$tier_match, na.rm = TRUE)

path_ref <- cmp$ref_tier == "pathogenic" & !cmp$missed
ben_ref <- cmp$ref_tier == "benign" & !cmp$missed
path_pred <- cmp$pred_tier == "pathogenic" & !cmp$missed

sensitivity <- if (any(path_ref)) mean(path_pred[path_ref], na.rm = TRUE) else NA_real_
specificity <- if (any(ben_ref)) mean(!path_pred[ben_ref], na.rm = TRUE) else NA_real_

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
csv_out <- file.path(output_dir, paste0("validation_comparison_", stamp, ".csv"))
md_out <- file.path(output_dir, paste0("validation_report_", stamp, ".md"))

utils::write.csv(cmp, csv_out, row.names = FALSE)

md_lines <- c(
  "# ClinicalVariantR Validation Report",
  "",
  sprintf("- **Engine:** %s", ACMG_PRO_ENGINE),
  sprintf("- **App version:** %s", APP_VERSION),
  sprintf("- **Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")),
  sprintf("- **Benchmark VCF:** `%s`", normalizePath(vcf_path, winslash = "/", mustWork = FALSE)),
  sprintf("- **Expert validation set:** `%s` (%d variants)", validation_tsv, nrow(ref)),
  "",
  "## Summary metrics",
  "",
  sprintf("| Metric | Value |"),
  sprintf("|--------|-------|"),
  sprintf("| Variants in expert set | %d |", nrow(ref)),
  sprintf("| Variants scored from VCF | %d |", n_compared),
  sprintf("| Exact class accuracy | %.1f%% |", 100 * exact_acc),
  sprintf("| Tier accuracy (P/LP vs B/LB vs VUS) | %.1f%% |", 100 * tier_acc),
  sprintf("| Pathogenic sensitivity | %s |", ifelse(is.na(sensitivity), "N/A", sprintf("%.1f%%", 100 * sensitivity))),
  sprintf("| Benign specificity | %s |", ifelse(is.na(specificity), "N/A", sprintf("%.1f%%", 100 * specificity))),
  "",
  "## Mismatches",
  ""
)

mismatches <- cmp[!isTRUE(cmp$exact_match) & !cmp$missed, , drop = FALSE]
if (nrow(mismatches) == 0L) {
  md_lines <- c(md_lines, "None — all scored variants match expert classification.")
} else {
  for (i in seq_len(nrow(mismatches))) {
    r <- mismatches[i, , drop = FALSE]
    md_lines <- c(
      md_lines,
      sprintf(
        "- **%s** (%s): expert `%s` vs predicted `%s` — %s",
        r$benchmark_id, r$gene, r$ref_class, r$pred_class, r$benchmark_note
      )
    )
  }
}

md_lines <- c(
  md_lines,
  "",
  "## Notes",
  "",
  paste(
    "This report uses the curated expert validation set in `data/validation/`.",
    "Expand toward 500+ variants by adding ClinGen VCEP or laboratory sign-out cases",
    "with the same TSV schema before production deployment."
  ),
  "",
  sprintf("Detailed comparison CSV: `%s`", csv_out)
)

writeLines(md_lines, md_out)

cat("Validation report written:\n")
cat(" ", md_out, "\n")
cat(" ", csv_out, "\n")
cat(sprintf(
  "Exact accuracy: %.1f%% | Tier accuracy: %.1f%% | n=%d\n",
  100 * exact_acc, 100 * tier_acc, n_compared
))
