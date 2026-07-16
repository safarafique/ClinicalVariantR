#!/usr/bin/env Rscript
# Install or validate production reference snapshots for ClinicalVariantR prediction mode.
#
# Usage (from ClinicalVariantR/):
#   Rscript scripts/install_reference_data.R --check
#   Rscript scripts/install_reference_data.R --gnomad /path/gnomad_subset.tsv --clinvar /path/clinvar.tsv --revel /path/revel.tsv
#   Rscript scripts/install_reference_data.R --copy-placeholders
#
# Expected TSV schemas (header row required):
#   gnomAD: chrom, pos, ref, alt, AF, popmax_AF
#   ClinVar: chrom, pos, ref, alt, clinical_significance, review_status
#   REVEL: chrom, pos, ref, alt, REVEL

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (length(hit) == 0L || hit[[1L]] >= length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

ref_dir <- file.path("data", "reference")
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

targets <- list(
  gnomad_v41 = file.path(ref_dir, "gnomad_v41.tsv"),
  clinvar = file.path(ref_dir, "clinvar.tsv"),
  revel = file.path(ref_dir, "revel.tsv")
)

placeholders <- list(
  gnomad_v41 = file.path(ref_dir, "gnomad_v41_placeholder.tsv"),
  clinvar = file.path(ref_dir, "clinvar_placeholder.tsv"),
  revel = file.path(ref_dir, "revel_placeholder.tsv")
)

required_cols <- list(
  gnomad_v41 = c("chrom", "pos", "ref", "alt", "AF"),
  clinvar = c("chrom", "pos", "ref", "alt", "clinical_significance"),
  revel = c("chrom", "pos", "ref", "alt", "REVEL")
)

validate_tsv <- function(path, cols) {
  if (!file.exists(path)) {
    return(list(ok = FALSE, rows = 0L, message = paste("Missing file:", path)))
  }
  hdr <- readLines(path, n = 1L, warn = FALSE)
  if (length(hdr) == 0L || !nzchar(hdr)) {
    return(list(ok = FALSE, rows = 0L, message = "Empty file"))
  }
  names <- strsplit(hdr, "\t", fixed = TRUE)[[1L]]
  missing <- setdiff(cols, names)
  if (length(missing) > 0L) {
    return(list(ok = FALSE, rows = 0L, message = paste("Missing columns:", paste(missing, collapse = ", "))))
  }
  n <- length(readLines(path, warn = FALSE)) - 1L
  list(ok = TRUE, rows = max(0L, n), message = "OK")
}

copy_if_given <- function(src, dest, name) {
  if (is.na(src) || !nzchar(src)) return(invisible(FALSE))
  if (!file.exists(src)) stop("Source not found for ", name, ": ", src, call. = FALSE)
  file.copy(src, dest, overwrite = TRUE)
  cat("Installed", name, "->", dest, "\n")
  invisible(TRUE)
}

if ("--copy-placeholders" %in% args) {
  for (nm in names(targets)) {
    if (file.exists(placeholders[[nm]])) {
      file.copy(placeholders[[nm]], targets[[nm]], overwrite = TRUE)
      cat("Copied placeholder", nm, "->", targets[[nm]], "\n")
    }
  }
}

copy_if_given(parse_arg("--gnomad"), targets$gnomad_v41, "gnomad")
copy_if_given(parse_arg("--clinvar"), targets$clinvar, "clinvar")
copy_if_given(parse_arg("--revel"), targets$revel, "revel")

cat("\nReference validation:\n")
all_ok <- TRUE
for (nm in names(targets)) {
  path <- if (file.exists(targets[[nm]])) targets[[nm]] else placeholders[[nm]]
  chk <- validate_tsv(path, required_cols[[nm]])
  status <- if (chk$ok) "OK" else "FAIL"
  cat(sprintf("  [%s] %s (%d rows) — %s\n", status, path, chk$rows, chk$message))
  if (!chk$ok) all_ok <- FALSE
}

cat("\nTo use production files, set REFERENCE_PATHS in global.R / global_cli.R to:\n")
cat(sprintf('  gnomad_v41 = "%s"\n', targets$gnomad_v41))
cat(sprintf('  clinvar    = "%s"\n', targets$clinvar))
cat(sprintf('  revel      = "%s"\n', targets$revel))

cat("\nSuggested sources:\n")
cat("  gnomAD v4 AF export: https://gnomad.broadinstitute.org/downloads\n")
cat("  ClinVar variant_summary: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz\n")
cat("  REVEL/dbNSFP: https://sites.google.com/site/revelgenomics/downloads\n")

if ("--check" %in% args) {
  quit(status = if (all_ok) 0L else 1L)
}

invisible(all_ok)
