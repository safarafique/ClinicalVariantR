#!/usr/bin/env Rscript
# Standalone: compare ACMGamp CSV vs InterVar output or InterVar-style reference TSV.
# Copy this file + intervar_compare.R to any machine with R 4.x.
#
# Official InterVar (Python, not R): https://github.com/WGLab/InterVar
# Web version (exonic SNVs):         https://wintervar.wglab.org/
#
# Usage:
#   Rscript compare_acmgamp_intervar.R --acmgamp report.csv --reference sample.acmg.tsv
#   Rscript compare_acmgamp_intervar.R --acmgamp report.csv --intervar myanno.hg38_multianno.txt.intervar
#   Rscript compare_acmgamp_intervar.R --acmgamp report.csv --reference ref.tsv --out results/compare_042
#
# Optional: run ACMGamp first from project root:
#   Rscript scripts/benchmark_testing_dataset.R ../testig/testig ../results/benchmark

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(acmgamp = NA_character_, reference = NA_character_, intervar = NA_character_, out = NA_character_)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
  if (key %in% c("--acmgamp", "-a") && i < length(args)) {
      out$acmgamp <- args[[i + 1L]]; i <- i + 2L
    } else if (key %in% c("--reference", "-r") && i < length(args)) {
      out$reference <- args[[i + 1L]]; i <- i + 2L
    } else if (key %in% c("--intervar", "-i") && i < length(args)) {
      out$intervar <- args[[i + 1L]]; i <- i + 2L
    } else if (key %in% c("--out", "-o") && i < length(args)) {
      out$out <- args[[i + 1L]]; i <- i + 2L
    } else if (key %in% c("--help", "-h")) {
      cat(
        "Compare ACMGamp vs InterVar or reference TSV\n\n",
        "Required: --acmgamp FILE\n",
        "One of:   --reference FILE  (your .acmg.tsv InterVar-style truth)\n",
        "          --intervar FILE   (InterVar *.intervar output)\n",
        "Optional: --out PREFIX      (default: comparison_<timestamp>)\n\n",
        "InterVar GitHub: https://github.com/WGLab/InterVar\n",
        sep = ""
      )
      quit(status = 0)
    } else {
      stop("Unknown argument: ", key, " (use --help)")
    }
  }
  out
}

opts <- parse_args(args)
if (is.na(opts$acmgamp)) stop("Provide --acmgamp path to ACMGamp CSV report.")
if (is.na(opts$reference) && is.na(opts$intervar)) {
  stop("Provide --reference (.acmg.tsv) or --intervar (.intervar output).")
}

script_dir <- tryCatch({
  dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]), winslash = "/"))
}, error = function(e) getwd())

compare_r <- file.path(script_dir, "..", "R", "intervar_compare.R")
if (!file.exists(compare_r)) compare_r <- file.path(script_dir, "intervar_compare.R")
if (!file.exists(compare_r)) stop("Cannot find intervar_compare.R next to this script.")
source(compare_r, local = FALSE)

acmg_df <- load_acmgamp_report_csv(normalizePath(opts$acmgamp, winslash = "/", mustWork = TRUE))
if (!is.na(opts$intervar)) {
  ref_df <- load_intervar_output(normalizePath(opts$intervar, winslash = "/", mustWork = TRUE))
  ref_name <- "InterVar"
} else {
  ref_df <- load_intervar_reference_tsv(normalizePath(opts$reference, winslash = "/", mustWork = TRUE))
  ref_name <- "Reference"
}

result <- compare_two_classification_tables(acmg_df, ref_df, left_name = "ACMGamp", right_name = ref_name)
prefix <- if (!is.na(opts$out)) {
  normalizePath(opts$out, winslash = "/", mustWork = FALSE)
} else {
  file.path(getwd(), paste0("comparison_", format(Sys.time(), "%Y%m%d_%H%M%S")))
}

paths <- write_comparison_outputs(result, prefix)

cat("=== ACMGamp vs ", ref_name, " ===\n", sep = "")
cat("ACMGamp variants:  ", nrow(acmg_df), "\n", sep = "")
cat("Reference variants:", nrow(ref_df), "\n", sep = "")
cat("Overlap:           ", result$metrics$n_overlap, "\n", sep = "")
cat("Exact 5-class acc: ", result$metrics$exact_accuracy, "%\n", sep = "")
cat("Tier accuracy:     ", result$metrics$tier_accuracy, "%\n\n", sep = "")

cat("ACMGamp category counts:\n")
print(classification_count_table(acmg_df))
cat("\n", ref_name, " category counts:\n", sep = "")
print(classification_count_table(ref_df))
cat("\nPer-class agreement (reference class):\n")
print(result$summary_by_class)

cat("\nWrote:\n")
cat(" ", paths$metrics_txt, "\n", sep = "")
if (!is.na(paths$comparison_csv)) cat(" ", paths$comparison_csv, "\n", sep = "")
if (!is.na(paths$summary_csv)) cat(" ", paths$summary_csv, "\n", sep = "")
