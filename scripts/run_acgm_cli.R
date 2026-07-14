#!/usr/bin/env Rscript
# Run ACMGamp on one VCF (CLI, no Shiny).
#
# Usage:
#   Rscript scripts/run_acgm_cli.R <vcf> <output_csv> [profile_id] [gene1,gene2,...]
#
# Examples:
#   Rscript scripts/run_acgm_cli.R ../HMC-1.final.vcf ../results/HMC-1.acmgamp.csv hematologic_predisposition
#   Rscript scripts/run_acgm_cli.R sample.vcf out.csv hematologic_predisposition ABL1,BCR,RUNX1,GATA2

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript scripts/run_acgm_cli.R <vcf> <output_csv> [profile_id] [genes_csv]\n",
    call. = FALSE
  )
}

script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
setwd(project_root)
source("global_cli.R")

vcf_path <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
output_csv <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
profile_id <- if (length(args) >= 3L && nzchar(args[[3L]])) args[[3L]] else DEFAULT_PROFILE_ID
gene_filter <- if (length(args) >= 4L && nzchar(args[[4L]])) {
  trimws(unlist(strsplit(args[[4L]], ",", fixed = TRUE)))
} else {
  character()
}

if (!file.exists(vcf_path)) stop("VCF not found: ", vcf_path)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

sample_id <- tools::file_path_sans_ext(basename(vcf_path))
session_id <- paste0(sample_id, "-CLI")

message("ACMGamp CLI: ", basename(vcf_path))
message("Profile: ", profile_id)
if (length(gene_filter) > 0L) message("Gene panel: ", paste(gene_filter, collapse = ", "))

refs <- tryCatch(load_reference_data(), error = function(e) {
  message("Reference load note: ", conditionMessage(e))
  NULL
})

result <- analyze_complete_vcf(
  vcf_path = vcf_path,
  mode = "rapid",
  output_csv = output_csv,
  pass_only = FALSE,
  use_bcftools = FALSE,
  refs = refs,
  profile_id = profile_id,
  session_id = session_id,
  gene_filter = gene_filter
)

message("Rows classified: ", result$rows_classified)
if (length(result$classification_counts) > 0L) {
  for (nm in names(result$classification_counts)) {
    message("  ", nm, ": ", result$classification_counts[[nm]])
  }
}
message("Wrote: ", output_csv)
