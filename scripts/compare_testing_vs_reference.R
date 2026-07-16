#!/usr/bin/env Rscript
# Score testing VCFs with ClinicalVariantR and compare to bundled .acmg.tsv (InterVar-style reference).
#
# Usage:
#   Rscript scripts/compare_testing_vs_reference.R [test_folder] [output_dir]

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
setwd(project_root)
source("global_cli.R")

test_folder <- if (length(args) >= 1L) normalizePath(args[[1L]], winslash = "/", mustWork = FALSE) else {
  normalizePath(file.path(project_root, "..", "testig", "testig"), winslash = "/", mustWork = FALSE)
}
output_dir <- if (length(args) >= 2L) normalizePath(args[[2L]], winslash = "/", mustWork = FALSE) else {
  file.path(project_root, "..", "results", "intervar_compare")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

vcf_files <- sort(list.files(test_folder, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE, ignore.case = TRUE))
if (length(vcf_files) == 0L) stop("No VCF files in: ", test_folder)

refs <- load_reference_data()
summary_rows <- list()

for (vcf_path in vcf_files) {
  sample_id <- sub("\\.haplotypecaller.*$", "", tools::file_path_sans_ext(basename(vcf_path)))
  tsv_candidates <- list.files(test_folder, pattern = paste0("^", sample_id, ".*\\.acmg\\.tsv$"), full.names = TRUE)
  if (length(tsv_candidates) == 0L) {
    message("Skip (no reference TSV): ", basename(vcf_path))
    next
  }
  tsv_path <- tsv_candidates[[1L]]
  out_prefix <- file.path(output_dir, sample_id)

  message("Analyzing: ", basename(vcf_path))
  result <- analyze_complete_vcf(
    vcf_path = vcf_path,
    mode = "rapid",
    output_csv = paste0(out_prefix, ".clinicalvariantr.csv"),
    pass_only = FALSE,
    use_bcftools = FALSE,
    refs = refs,
    profile_id = DEFAULT_PROFILE_ID
  )

  acmg_df <- load_clinicalvariantr_report_csv(result$output_csv)
  ref_df <- load_intervar_reference_tsv(tsv_path)
  cmp <- compare_two_classification_tables(acmg_df, ref_df, left_name = "ClinicalVariantR", right_name = "Reference")
  paths <- write_comparison_outputs(cmp, out_prefix)

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    sample_id = sample_id,
    n_clinicalvariantr = nrow(acmg_df),
    n_reference = nrow(ref_df),
    n_overlap = cmp$metrics$n_overlap,
    exact_accuracy_pct = cmp$metrics$exact_accuracy,
    tier_accuracy_pct = cmp$metrics$tier_accuracy,
    clinicalvariantr_csv = result$output_csv,
    comparison_csv = paths$comparison_csv,
    stringsAsFactors = FALSE
  )
}

if (length(summary_rows) == 0L) stop("No samples compared.")

summary_df <- do.call(rbind, summary_rows)
summary_path <- file.path(output_dir, "intervar_style_comparison_summary.csv")
utils::write.csv(summary_df, summary_path, row.names = FALSE)

cat("\n=== ClinicalVariantR vs reference .acmg.tsv (InterVar-style) ===\n")
print(summary_df)
cat("\nWrote: ", summary_path, "\n", sep = "")
