#!/usr/bin/env Rscript
# Export explainable ACMG evidence reports from parsed variant CSV.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript scripts/export_evidence_report.R <parsed_or_scored.csv> [profile_id] [out_dir]")
}

input_csv <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
profile_id <- if (length(args) >= 2L) args[[2L]] else DEFAULT_PROFILE_ID
out_dir <- if (length(args) >= 3L) {
  normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
} else {
  dirname(input_csv)
}

script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
setwd(project_root)
source("global.R")

variants <- read.csv(input_csv, stringsAsFactors = FALSE)
scored <- score_variants_table(variants, profile_id = profile_id)
report <- acmg_pro_to_report(scored, mode = "rapid", session_id = "CLI")

base <- tools::file_path_sans_ext(basename(input_csv))
report_path <- file.path(out_dir, paste0(base, ".evidence_report.csv"))
evidence_path <- file.path(out_dir, paste0(base, ".evidence_detail.csv"))

write.csv(report, report_path, row.names = FALSE)

evidence_rows <- lapply(seq_len(nrow(scored)), function(i) {
  tbl <- parse_evidence_json(scored$evidence_json[i])
  if (nrow(tbl) == 0) return(NULL)
  cbind(variant_id = scored$variant_id[i], tbl)
})
evidence_rows <- Filter(Negate(is.null), evidence_rows)
if (length(evidence_rows) > 0) {
  write.csv(do.call(rbind, evidence_rows), evidence_path, row.names = FALSE)
}

metadata <- build_run_metadata(vcf_path = input_csv, profile_id = profile_id, mode = "rapid", session_id = "CLI")
write_run_metadata_json(metadata, file.path(out_dir, paste0(base, ".metadata.json")))

cat("Wrote:\n", report_path, "\n", sep = "")
if (length(evidence_rows) > 0) cat(evidence_path, "\n", sep = "")
cat("Done.\n")
