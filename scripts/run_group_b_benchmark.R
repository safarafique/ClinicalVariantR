#!/usr/bin/env Rscript
# Score the Group B efficiency benchmark and compare to expected labels.
#
# Usage:
#   Rscript scripts/generate_group_b_benchmark.R   # once, if files missing
#   Rscript scripts/run_group_b_benchmark.R

script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
setwd(project_root)
source("global_cli.R")

bench_dir <- normalizePath(file.path(project_root, "..", "testig", "acmgamp_benchmark"), mustWork = FALSE)
vcf_path <- file.path(bench_dir, "acmgamp_group_b_benchmark.vcf")
tsv_path <- file.path(bench_dir, "acmgamp_group_b_benchmark.acmg.tsv")

if (!file.exists(vcf_path) || !file.exists(tsv_path)) {
  message("Group B benchmark files missing — generating...")
  source(file.path(project_root, "scripts", "generate_group_b_benchmark.R"))
}

out_dir <- file.path(project_root, "..", "results", "acmgamp_benchmark")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== ACMGamp Group B efficiency benchmark ===\n")
cat("Engine:", ACMG_PRO_ENGINE, " (", APP_VERSION, ")\n", sep = "")
cat("VCF:", vcf_path, "\n\n")

res <- benchmark_one_sample(vcf_path, tsv_path, profile_id = "general_germline")

if (!is.null(res$error)) {
  cat("ERROR:", res$error, "\n")
  if (!is.null(res$n_vcf)) cat("VCF variants:", res$n_vcf, " Reference:", res$n_ref, "\n")
  quit(status = 1)
}

cmp <- res$comparison
cmp$match_symbol <- ifelse(cmp$exact_match, "OK", "MISMATCH")
cmp$tier_symbol <- ifelse(cmp$tier_match, "OK", "TIER_MISMATCH")

print(cmp[, c("variant_key", "gene", "pred_class", "ref_class", "criteria_met", "match_symbol")])

cat("\n--- Summary ---\n")
cat("Compared:", res$metrics$n_compared, "/", res$metrics$n_reference_variants, "\n")
cat("Exact accuracy:", round(100 * res$metrics$exact_accuracy, 1), "%\n")
cat("Tier accuracy:", round(100 * res$metrics$tier_accuracy, 1), "%\n")
if (!is.null(res$metrics$sensitivity_pathogenic)) {
  cat("Pathogenic sensitivity:", round(100 * res$metrics$sensitivity_pathogenic, 1), "%\n")
}
if (!is.null(res$metrics$specificity_benign)) {
  cat("Benign specificity (no false pathogenic):", round(100 * res$metrics$specificity_benign, 1), "%\n")
}
cat("Score time:", round(res$metrics$score_seconds, 3), "s\n")
cat("Throughput:", round(res$metrics$variants_per_second, 1), " variants/s\n")

out_csv <- file.path(out_dir, "group_b_benchmark_comparison.csv")
write.csv(cmp, out_csv, row.names = FALSE)
cat("\nWrote:", out_csv, "\n")

mismatches <- cmp[!cmp$exact_match, , drop = FALSE]
if (nrow(mismatches) > 0L) {
  cat("\nMismatches:\n")
  apply(mismatches, 1L, function(r) {
    cat(sprintf("  %s | expected %s | got %s | %s\n",
                r[["variant_key"]], r[["ref_class"]], r[["pred_class"]], r[["criteria_met"]]))
  })
  quit(status = 1)
}

cat("\nAll Group B benchmark variants matched expected classifications.\n")
