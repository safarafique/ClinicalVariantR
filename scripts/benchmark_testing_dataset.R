#!/usr/bin/env Rscript
# Benchmark ACMGamp efficiency against testing dataset with reference .acmg.tsv files.
#
# Usage:
#   Rscript scripts/benchmark_testing_dataset.R [test_folder] [output_dir] [profile_id]
#
# Example:
#   Rscript scripts/benchmark_testing_dataset.R ../testig/testig ../results/benchmark general_germline

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
  file.path(project_root, "..", "results", "benchmark")
}
profile_id <- if (length(args) >= 3L) args[[3L]] else DEFAULT_PROFILE_ID

if (!dir.exists(test_folder)) stop("Test folder not found: ", test_folder)

cat("=== ACMGamp Benchmark ===\n")
cat("Test folder: ", test_folder, "\n", sep = "")
cat("Profile:     ", profile_id, "\n", sep = "")
cat("Output:      ", output_dir, "\n\n", sep = "")

out <- benchmark_testing_folder(test_folder, profile_id = profile_id, output_dir = output_dir)

print(out$summary)

cat("\n--- Interpretation ---\n")
cat("exact_accuracy: 5-class label match (Pathogenic, Likely Pathogenic, VUS, ...)\n")
cat("tier_accuracy:  collapsed match (pathogenic / benign / vus)\n")
cat("variants_per_second: throughput (higher = faster)\n")
cat("\nWrote: ", file.path(output_dir, "benchmark_summary.csv"), "\n", sep = "")
cat("Done.\n")
