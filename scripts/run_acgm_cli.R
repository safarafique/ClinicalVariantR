#!/usr/bin/env Rscript
# Run ClinicalVariantR Group B (rapid) on one VCF — CLI, no Shiny.
#
# Usage:
#   Rscript scripts/run_acgm_cli.R <vcf> <output_csv> [profile_id] [genes_csv]
#     [--pass-only] [--min-qual N] [--skip-audit] [--chunk-size N]
#     [--bcftools|--no-bcftools] [--no-count]
#
# Examples (from ClinicalVariantR/):
#   Rscript scripts/run_acgm_cli.R \
#     BC_sample5.CRCh.38.REVEL.prepared.CADD.prepared.vcf \
#     out/BC_sample5.GroupB.QUAL50.csv \
#     hematologic_predisposition \
#     --pass-only --min-qual 50 --skip-audit --chunk-size 50000 --no-count

raw <- commandArgs(trailingOnly = TRUE)
if (length(raw) < 2L) {
  stop(
    "Usage: Rscript scripts/run_acgm_cli.R <vcf> <output_csv> [profile_id] [genes_csv]\n",
    "  [--pass-only] [--min-qual N] [--skip-audit] [--chunk-size N]\n",
    "  [--bcftools|--no-bcftools] [--no-count]\n",
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

# ---- parse positional + flags ------------------------------------------------
pass_only <- FALSE
min_qual <- 0
write_audit <- TRUE
chunk_size <- 10000L
use_bcftools <- TRUE
no_count <- FALSE
positionals <- character()

i <- 1L
while (i <= length(raw)) {
  a <- raw[[i]]
  if (identical(a, "--pass-only")) {
    pass_only <- TRUE
  } else if (identical(a, "--skip-audit")) {
    write_audit <- FALSE
  } else if (identical(a, "--no-count")) {
    no_count <- TRUE
  } else if (identical(a, "--bcftools")) {
    use_bcftools <- TRUE
  } else if (identical(a, "--no-bcftools")) {
    use_bcftools <- FALSE
  } else if (identical(a, "--min-qual")) {
    i <- i + 1L
    if (i > length(raw)) stop("--min-qual requires a number", call. = FALSE)
    min_qual <- as.numeric(raw[[i]])
  } else if (identical(a, "--chunk-size")) {
    i <- i + 1L
    if (i > length(raw)) stop("--chunk-size requires a number", call. = FALSE)
    chunk_size <- as.integer(raw[[i]])
  } else if (startsWith(a, "--")) {
    stop("Unknown flag: ", a, call. = FALSE)
  } else {
    positionals <- c(positionals, a)
  }
  i <- i + 1L
}

if (length(positionals) < 2L) {
  stop("Need <vcf> and <output_csv> positional arguments.", call. = FALSE)
}

vcf_path <- normalizePath(positionals[[1L]], winslash = "/", mustWork = FALSE)
output_csv <- normalizePath(positionals[[2L]], winslash = "/", mustWork = FALSE)
profile_id <- if (length(positionals) >= 3L && nzchar(positionals[[3L]])) {
  positionals[[3L]]
} else {
  DEFAULT_PROFILE_ID
}
gene_filter <- if (length(positionals) >= 4L && nzchar(positionals[[4L]])) {
  trimws(unlist(strsplit(positionals[[4L]], ",", fixed = TRUE)))
} else {
  character()
}

if (!file.exists(vcf_path)) stop("VCF not found: ", vcf_path, call. = FALSE)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

# Prefer bcftools automatically when available (WSL/Linux). Windows often has none.
if (isTRUE(use_bcftools) && exists("bcftools_available") && is.function(bcftools_available)) {
  use_bcftools <- isTRUE(bcftools_available())
}

sample_id <- tools::file_path_sans_ext(basename(vcf_path))
session_id <- paste0(sample_id, "-CLI")

message("ClinicalVariantR CLI: ", basename(vcf_path))
message("Profile: ", profile_id)
message(
  "Options: pass_only=", pass_only,
  " min_qual=", min_qual,
  " chunk_size=", chunk_size,
  " write_audit=", write_audit,
  " bcftools=", use_bcftools,
  " no_count=", no_count
)
if (length(gene_filter) > 0L) message("Gene panel: ", paste(gene_filter, collapse = ", "))

refs <- tryCatch(load_reference_data(), error = function(e) {
  message("Reference load note: ", conditionMessage(e))
  NULL
})

t0 <- proc.time()[["elapsed"]]
result <- analyze_complete_vcf(
  vcf_path = vcf_path,
  mode = "rapid",
  output_csv = output_csv,
  pass_only = pass_only,
  min_qual = min_qual,
  chunk_size = chunk_size,
  use_bcftools = use_bcftools,
  refs = refs,
  profile_id = profile_id,
  session_id = session_id,
  gene_filter = gene_filter,
  write_audit = write_audit
)
elapsed <- proc.time()[["elapsed"]] - t0

message("Rows classified: ", result$rows_classified %||% NA)
if (!isTRUE(no_count) && length(result$classification_counts) > 0L) {
  for (nm in names(result$classification_counts)) {
    message("  ", nm, ": ", result$classification_counts[[nm]])
  }
}
message(sprintf("Elapsed: %.1f sec", elapsed))
message("Wrote: ", output_csv)
