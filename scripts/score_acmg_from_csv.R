#!/usr/bin/env Rscript
# Re-score ACMG criteria from an existing parsed variant CSV (no VCF re-parse).

args <- commandArgs(trailingOnly = TRUE)
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", initial_options[grep("^--file=", initial_options)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
input_csv <- if (length(args) >= 1L) args[[1L]] else file.path(project_root, "..", "results", "HMC-1.parsed_variants.csv")
out_csv <- if (length(args) >= 2L) args[[2L]] else sub("\\.csv$", ".acmg_detailed.csv", input_csv)
profile_id <- if (length(args) >= 3L) args[[3L]] else "hematologic_predisposition"

setwd(project_root)
source("global.R")

if (!file.exists(input_csv)) stop("Input CSV not found: ", input_csv)

variants <- read.csv(input_csv, stringsAsFactors = FALSE)
lof_panel <- file.path(project_root, "data", "gene_panels", "lof_disease_mechanism_genes.csv")

message("Scoring: ", nrow(variants), " variants from ", input_csv)
detailed <- score_variants_table(variants, lof_panel_path = lof_panel, profile_id = profile_id)
report <- acmg_pro_to_report(detailed, mode = "rapid", session_id = "CSV-CLI")
write.csv(detailed, out_csv, row.names = FALSE)
write.csv(report, sub("\\.csv$", ".evidence_report.csv", out_csv), row.names = FALSE)
message("Wrote: ", out_csv)

summary <- acmg_pro_criteria_summary(detailed)
summary_path <- sub("\\.csv$", ".acmg_criteria_summary.csv", out_csv)
write.csv(summary, summary_path, row.names = FALSE)
print(summary)
