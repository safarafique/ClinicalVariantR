#!/usr/bin/env Rscript
# Phase 1 parity benchmark: ACMGamp vs expert labels and optional InterVar-style reference.
#
# Usage (from cml_variant_interpreter/):
#   Rscript scripts/benchmark_vs_intervar.R
#   Rscript scripts/benchmark_vs_intervar.R --vcf ../testig/acmgamp_benchmark/acmgamp_group_b_benchmark.vcf
#   Rscript scripts/benchmark_vs_intervar.R --reference ../testig/testig/sample041.acmg.tsv --min-tier 0.95

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (length(hit) == 0L || hit[[1L]] >= length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

parse_num_arg <- function(flag, default) {
  raw <- parse_arg(flag, NA_character_)
  if (is.na(raw)) return(default)
  suppressWarnings(as.numeric(raw))
}

vcf_path <- parse_arg(
  "--vcf",
  file.path("..", "testig", "acmgamp_benchmark", "acmgamp_group_b_benchmark.vcf")
)
reference_tsv <- parse_arg("--reference", NA_character_)
validation_tsv <- file.path("data", "validation", "acmgamp_validation_set.tsv")
min_tier <- parse_num_arg("--min-tier", 0.95)
output_dir <- file.path("results", "intervar_compare")
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
             "classification", "criteria_met", "confidence_score", "evidence_strength",
             "max_population_af", "popmax_af", "af_known")],
  ref[, c("variant_key", "benchmark_id", "benchmark_category", "ref_class", "ref_tier",
          "acmg_criteria_base", "benchmark_note", "expert_source")],
  by = "variant_key",
  all.y = TRUE
)

cmp$exact_match <- cmp$pred_class == cmp$ref_class
cmp$tier_match <- cmp$pred_tier == cmp$ref_tier
cmp$missed <- is.na(cmp$classification)

exact_acc <- mean(cmp$exact_match, na.rm = TRUE)
tier_acc <- mean(cmp$tier_match, na.rm = TRUE)

path_ref <- cmp$ref_tier == "pathogenic" & !cmp$missed
ben_ref <- cmp$ref_tier == "benign" & !cmp$missed
path_pred <- cmp$pred_tier == "pathogenic" & !cmp$missed

sensitivity <- if (any(path_ref)) mean(path_pred[path_ref], na.rm = TRUE) else NA_real_
specificity <- if (any(ben_ref)) mean(!path_pred[ben_ref], na.rm = TRUE) else NA_real_

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
csv_out <- file.path(output_dir, paste0("phase1_benchmark_", stamp, ".csv"))
summary_out <- file.path(output_dir, paste0("phase1_benchmark_summary_", stamp, ".csv"))

utils::write.csv(cmp, csv_out, row.names = FALSE)

intervar_tier_acc <- NA_real_
intervar_exact_acc <- NA_real_
if (!is.na(reference_tsv) && nzchar(reference_tsv) && file.exists(reference_tsv)) {
  iv <- load_intervar_reference_tsv(reference_tsv)
  acmg_df <- data.frame(
    variant_key = scored$variant_key,
    gene = scored$gene,
    classification = scored$pred_class,
    criteria = scored$criteria_met,
    source = "ACMGamp",
    stringsAsFactors = FALSE
  )
  acmg_df$tier <- collapse_to_tier(acmg_df$classification)
  iv_cmp <- compare_two_classification_tables(acmg_df, iv, left_name = "ACMGamp", right_name = "Reference")
  intervar_tier_acc <- iv_cmp$metrics$tier_accuracy / 100
  intervar_exact_acc <- iv_cmp$metrics$exact_accuracy / 100
  utils::write.csv(
    iv_cmp$comparison,
    file.path(output_dir, paste0("intervar_style_compare_", stamp, ".csv")),
    row.names = FALSE
  )
}

summary <- data.frame(
  engine = ACMG_PRO_ENGINE,
  app_version = APP_VERSION,
  benchmark_vcf = normalizePath(vcf_path, winslash = "/", mustWork = FALSE),
  n_variants = nrow(ref),
  exact_accuracy = exact_acc,
  tier_accuracy = tier_acc,
  pathogenic_sensitivity = sensitivity,
  benign_specificity = specificity,
  intervar_reference = ifelse(is.na(reference_tsv), "", reference_tsv),
  intervar_tier_accuracy = intervar_tier_acc,
  intervar_exact_accuracy = intervar_exact_acc,
  min_tier_required = min_tier,
  pass = tier_acc >= min_tier,
  stringsAsFactors = FALSE
)
utils::write.csv(summary, summary_out, row.names = FALSE)

cat("=== ACMGamp Phase 1 InterVar parity benchmark ===\n")
cat("Engine:", ACMG_PRO_ENGINE, "\n")
cat(sprintf("Exact accuracy: %.1f%% (%d/%d)\n", 100 * exact_acc, sum(cmp$exact_match), nrow(cmp)))
cat(sprintf("Tier accuracy:  %.1f%% (%d/%d)\n", 100 * tier_acc, sum(cmp$tier_match), nrow(cmp)))
cat(sprintf("Pathogenic sensitivity: %s\n", ifelse(is.na(sensitivity), "n/a", sprintf("%.1f%%", 100 * sensitivity))))
cat(sprintf("Benign specificity:     %s\n", ifelse(is.na(specificity), "n/a", sprintf("%.1f%%", 100 * specificity))))
if (!is.na(intervar_tier_acc)) {
  cat(sprintf("InterVar-ref tier agreement: %.1f%%\n", 100 * intervar_tier_acc))
}
cat("\nMismatches:\n")
bad <- cmp[!cmp$tier_match | cmp$missed, ]
if (nrow(bad) == 0L) {
  cat("  (none)\n")
} else {
  for (i in seq_len(nrow(bad))) {
    cat(sprintf(
      "  %s: expert=%s (%s) vs ACMGamp=%s (%s) [%s]\n",
      bad$benchmark_id[i],
      bad$ref_class[i], bad$ref_tier[i],
      ifelse(bad$missed[i], "MISSING", bad$pred_class[i]),
      ifelse(bad$missed[i], "n/a", bad$pred_tier[i]),
      bad$criteria_met[i] %||% ""
    ))
  }
}
cat("\nWrote:", csv_out, "\n")
cat("Wrote:", summary_out, "\n")

if (!isTRUE(summary$pass[[1L]])) {
  cat(sprintf("\nFAIL: tier accuracy %.1f%% < required %.1f%%\n", 100 * tier_acc, 100 * min_tier))
  quit(status = 1L)
}

cat(sprintf("\nPASS: tier accuracy %.1f%% >= %.1f%%\n", 100 * tier_acc, 100 * min_tier))
quit(status = 0L)
