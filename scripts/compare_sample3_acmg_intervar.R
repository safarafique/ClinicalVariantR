#!/usr/bin/env Rscript
# Score Sample 3 pilot VCF with ClinicalVariantR and compare to InterVar when available.
#
# Usage:
#   Rscript scripts/compare_sample3_acmg_intervar.R
#   Rscript scripts/compare_sample3_acmg_intervar.R --intervar path/to/Sample3.pilot.hg38_multianno.txt.intervar
#   Rscript scripts/compare_sample3_acmg_intervar.R --pilot-vcf path/to/pilot.vcf --out-dir results/intervar_compare

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    pilot_vcf = NA_character_,
    intervar = NA_character_,
    out_dir = NA_character_,
    profile = "hematologic_predisposition"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--pilot-vcf", "-v") && i < length(args)) {
      out$pilot_vcf <- args[[i + 1L]]; i <- i + 2L
    } else if (key %in% c("--intervar", "-i") && i < length(args)) {
      out$intervar <- args[[i + 1L]]; i <- i + 2L
    } else if (key %in% c("--out-dir", "-o") && i < length(args)) {
      out$out_dir <- args[[i + 1L]]; i <- i + 2L
    } else if (key == "--profile" && i < length(args)) {
      out$profile <- args[[i + 1L]]; i <- i + 2L
    } else if (key %in% c("--help", "-h")) {
      cat(
        "Compare ClinicalVariantR vs InterVar on Sample 3 pilot VCF\n\n",
        "Options:\n",
        "  --pilot-vcf FILE   Pilot VCF (default: ../results/intervar_compare/Sample3.pilot.vcf)\n",
        "  --intervar FILE    InterVar *.intervar output (auto-detected if present)\n",
        "  --out-dir DIR      Output directory (default: ../results/intervar_compare)\n",
        "  --profile ID       Disease profile (default: hematologic_predisposition)\n",
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

initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", initial_options[grep("^--file=", initial_options)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
data_root <- normalizePath(file.path(project_root, ".."), winslash = "/", mustWork = FALSE)

out_dir <- if (!is.na(opts$out_dir)) {
  normalizePath(opts$out_dir, winslash = "/", mustWork = FALSE)
} else {
  file.path(data_root, "results", "intervar_compare")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pilot_vcf <- if (!is.na(opts$pilot_vcf)) {
  normalizePath(opts$pilot_vcf, winslash = "/", mustWork = TRUE)
} else {
  file.path(out_dir, "Sample3.pilot.vcf")
}
if (!file.exists(pilot_vcf)) {
  stop(
    "Pilot VCF not found: ", pilot_vcf, "\n",
    "Create it with: bash scripts/subset_sample3_vcf.sh"
  )
}

intervar_path <- if (!is.na(opts$intervar)) {
  normalizePath(opts$intervar, winslash = "/", mustWork = FALSE)
} else {
  file.path(out_dir, "Sample3.pilot.hg38_multianno.txt.intervar")
}

setwd(project_root)
source("global.R")

parse_pilot_vcf <- function(path) {
  con <- file(path, "r")
  on.exit(close(con), add = TRUE)
  rows <- list()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (grepl("^#", line)) next
    parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 8L) next
    alt <- strsplit(parts[[5L]], ",")[[1L]][1L]
    rows[[length(rows) + 1L]] <- parse_variant_from_vcf_fields(
      chrom = parts[[1L]], pos = parts[[2L]], ref = parts[[4L]], alt = alt,
      qual = suppressWarnings(as.numeric(parts[[6L]])),
      filter = parts[[7L]], info = parts[[8L]]
    )
  }
  if (length(rows) == 0L) stop("No variants in pilot VCF: ", path)
  do.call(rbind, rows)
}

message("Pilot VCF: ", pilot_vcf)
variants <- parse_pilot_vcf(pilot_vcf)
message("Pilot variants: ", nrow(variants))

lof_panel <- file.path(project_root, "data", "gene_panels", "lof_disease_mechanism_genes.csv")
scored <- score_variants_table(
  variants,
  lof_panel_path = lof_panel,
  profile_id = opts$profile
)

acmg_pilot_path <- file.path(out_dir, "Sample3.pilot.clinicalvariantr.csv")
write.csv(scored, acmg_pilot_path, row.names = FALSE)
message("Wrote ClinicalVariantR pilot results: ", acmg_pilot_path)

cat("\nClinicalVariantR pilot classification counts:\n")
print(table(scored$classification, useNA = "ifany"))

if (!file.exists(intervar_path)) {
  cat("\n")
  cat("=== InterVar output not found ===\n")
  cat("Expected: ", intervar_path, "\n", sep = "")
  cat("\nInterVar requires ANNOVAR. After installing ANNOVAR in WSL, run:\n")
  cat("  bash scripts/run_intervar_sample3.sh\n")
  cat("\nThen re-run this script to compare:\n")
  cat("  Rscript scripts/compare_sample3_acmg_intervar.R --intervar ", intervar_path, "\n", sep = "")
  quit(status = 0)
}

source(file.path(project_root, "R", "intervar_compare.R"), local = FALSE)

acmg_df <- load_clinicalvariantr_report_csv(acmg_pilot_path)
ref_df <- load_intervar_output(intervar_path)
result <- compare_two_classification_tables(acmg_df, ref_df, left_name = "ClinicalVariantR", right_name = "InterVar")
prefix <- file.path(out_dir, "Sample3.pilot_compare")
paths <- write_comparison_outputs(result, prefix)

cat("\n=== ClinicalVariantR vs InterVar (Sample 3 pilot) ===\n")
cat("ClinicalVariantR variants:  ", nrow(acmg_df), "\n", sep = "")
cat("InterVar variants: ", nrow(ref_df), "\n", sep = "")
cat("Overlap:           ", result$metrics$n_overlap, "\n", sep = "")
cat("Exact 5-class acc: ", result$metrics$exact_accuracy, "%\n", sep = "")
cat("Tier accuracy:     ", result$metrics$tier_accuracy, "%\n\n", sep = "")

cat("ClinicalVariantR category counts:\n")
print(classification_count_table(acmg_df))
cat("\nInterVar category counts:\n")
print(classification_count_table(ref_df))
cat("\nPer-class agreement (InterVar class):\n")
print(result$summary_by_class)

mismatches <- result$comparison[
  result$comparison$classification_left != result$comparison$classification_right,
  ,
  drop = FALSE
]
if (nrow(mismatches) > 0L) {
  mismatch_path <- file.path(out_dir, "Sample3.pilot_mismatches.csv")
  utils::write.csv(mismatches, mismatch_path, row.names = FALSE)
  cat("\nMismatches written: ", mismatch_path, " (", nrow(mismatches), " variants)\n", sep = "")
}

cat("\nWrote:\n")
cat(" ", paths$metrics_txt, "\n", sep = "")
if (!is.na(paths$comparison_csv)) cat(" ", paths$comparison_csv, "\n", sep = "")
if (!is.na(paths$summary_csv)) cat(" ", paths$summary_csv, "\n", sep = "")
