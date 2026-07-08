#!/usr/bin/env Rscript
# Compare existing ACMGamp CSV reports to reference .acmg.tsv files (base R only).
# No Shiny, no VCF re-analysis — use when reports were generated on Windows.
#
# Usage:
#   Rscript scripts/compare_reference_only.R [test_folder] [acmgamp_csv_dir] [output_dir]

args <- commandArgs(trailingOnly = TRUE)
script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
setwd(project_root)
source("R/intervar_compare.R")

test_folder <- if (length(args) >= 1L) normalizePath(args[[1L]], winslash = "/", mustWork = FALSE) else {
  normalizePath(file.path(project_root, "..", "testig", "testig"), winslash = "/", mustWork = FALSE)
}
acmg_dir <- if (length(args) >= 2L) normalizePath(args[[2L]], winslash = "/", mustWork = FALSE) else {
  file.path(project_root, "..", "results", "intervar_compare")
}
output_dir <- if (length(args) >= 3L) normalizePath(args[[3L]], winslash = "/", mustWork = FALSE) else acmg_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tsv_files <- sort(list.files(test_folder, pattern = "\\.acmg\\.tsv$", full.names = TRUE))
if (length(tsv_files) == 0L) stop("No .acmg.tsv reference files in: ", test_folder)

summary_rows <- list()
for (tsv_path in tsv_files) {
  sample_id <- sub("\\.acmg\\.tsv$", "", basename(tsv_path))
  acmg_candidates <- list.files(acmg_dir, pattern = paste0("^", sample_id, ".*\\.acmgamp\\.csv$"), full.names = TRUE)
  if (length(acmg_candidates) == 0L) {
    acmg_candidates <- list.files(acmg_dir, pattern = paste0(sample_id, ".*\\.csv$"), full.names = TRUE)
  }
  if (length(acmg_candidates) == 0L) {
    message("Skip (no ACMGamp CSV): ", sample_id)
    next
  }
  acmg_path <- acmg_candidates[[1L]]
  out_prefix <- file.path(output_dir, sample_id)

  acmg_df <- load_acmgamp_report_csv(acmg_path)
  ref_df <- load_intervar_reference_tsv(tsv_path)
  cmp <- compare_two_classification_tables(acmg_df, ref_df, left_name = "ACMGamp", right_name = "Reference")
  paths <- write_comparison_outputs(cmp, out_prefix)

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    sample_id = sample_id,
    n_acmgamp = nrow(acmg_df),
    n_reference = nrow(ref_df),
    n_overlap = cmp$metrics$n_overlap,
    exact_accuracy_pct = cmp$metrics$exact_accuracy,
    tier_accuracy_pct = cmp$metrics$tier_accuracy,
    acmgamp_csv = acmg_path,
    comparison_csv = paths$comparison_csv,
    stringsAsFactors = FALSE
  )
}

if (length(summary_rows) == 0L) {
  stop(
    "No comparisons run. Place ACMGamp CSV files in: ", acmg_dir,
    "\nOr run full pipeline after: Rscript scripts/install_r_cli_deps.R"
  )
}

summary_df <- do.call(rbind, summary_rows)
summary_path <- file.path(output_dir, "intervar_style_comparison_summary.csv")
utils::write.csv(summary_df, summary_path, row.names = FALSE)

cat("\n=== ACMGamp vs reference .acmg.tsv ===\n")
print(summary_df)
cat("\nWrote: ", summary_path, "\n", sep = "")
