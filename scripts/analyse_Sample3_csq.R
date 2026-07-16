#!/usr/bin/env Rscript
# Parse Sample 3 VEP CSQ VCF (GRCh38), filter rare coding variants, run ClinicalVariantR Pro.

args <- commandArgs(trailingOnly = TRUE)
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", initial_options[grep("^--file=", initial_options)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
} else {
  normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
}

acgm_root <- project_root
data_root <- normalizePath(file.path(project_root, ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(project_root, "R"))) {
  acgm_root <- file.path(project_root, "ClinicalVariantR")
  data_root <- project_root
}

vcf_path <- if (length(args) >= 2L) {
  normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
} else {
  file.path(data_root, "Sample3.haplotypecaller.filtered_VEP.ann (2).vcf")
}
out_dir <- if (length(args) >= 3L) {
  normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
} else {
  file.path(data_root, "results")
}

if (!file.exists(vcf_path)) {
  stop("VCF not found: ", vcf_path)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

setwd(acgm_root)
source("global.R")

AF_CUTOFF <- 0.01
CODING_PATTERN <- "missense|synonymous|stop_gained|frameshift|splice|inframe|start_lost|stop_lost"

parse_vcf_csq <- function(path, pass_only = TRUE) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)

  header_cols <- NULL
  rows <- list()

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (grepl("^#CHROM\t", line)) {
      header_cols <- strsplit(sub("^#", "", line), "\t", fixed = TRUE)[[1L]]
      next
    }
    if (grepl("^#", line)) next

    parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 8L) next
    if (isTRUE(pass_only) && !vcf_filter_is_pass(parts[[7L]])) next

    info <- parts[[8L]]
    if (!grepl("(^|;)CSQ=", info, perl = TRUE)) next
    if (!grepl(CODING_PATTERN, info, ignore.case = TRUE)) next

    row <- parse_variant_from_vcf_fields(
      chrom = parts[[1L]], pos = parts[[2L]], ref = parts[[4L]],
      alt = strsplit(parts[[5L]], ",")[[1L]][1L],
      qual = scalar_num(parts[[6L]]),
      filter = parts[[7L]],
      info = info
    )
    rows[[length(rows) + 1L]] <- row
  }

  if (length(rows) == 0L) stop("No coding CSQ variants parsed from: ", path)
  do.call(rbind, rows)
}

message("Parsing: ", vcf_path)
variants <- parse_vcf_csq(vcf_path, pass_only = TRUE)
message("Parsed coding variants: ", nrow(variants))

parsed_path <- file.path(out_dir, "Sample3.parsed_variants.csv")
write.csv(variants, parsed_path, row.names = FALSE)
message("Wrote: ", parsed_path)

lof_panel <- file.path(acgm_root, "data", "gene_panels", "lof_disease_mechanism_genes.csv")
detailed <- score_variants_table(variants, lof_panel_path = lof_panel, profile_id = "hematologic_predisposition")
evidence_report <- acmg_pro_to_report(detailed, mode = "rapid", session_id = "Sample3-CLI")

detailed_path <- file.path(out_dir, "Sample3.acmg_detailed.csv")
write.csv(detailed, detailed_path, row.names = FALSE)
message("Wrote: ", detailed_path)

evidence_report_path <- file.path(out_dir, "Sample3.evidence_report.csv")
write.csv(evidence_report, evidence_report_path, row.names = FALSE)
message("Wrote: ", evidence_report_path)

metadata <- build_run_metadata(vcf_path = vcf_path, profile_id = "hematologic_predisposition", mode = "rapid", session_id = "Sample3-CLI")
write_run_metadata_json(metadata, file.path(out_dir, "Sample3.metadata.json"))
message("Wrote: ", file.path(out_dir, "Sample3.metadata.json"))

rare_detailed <- detailed[is.na(detailed$max_population_af) | detailed$max_population_af <= AF_CUTOFF, , drop = FALSE]
message("Rare variants (AF <= ", AF_CUTOFF, " or missing): ", nrow(rare_detailed))

rare_path <- file.path(out_dir, "Sample3.rare_coding_variants.csv")
write.csv(rare_detailed, rare_path, row.names = FALSE)
message("Wrote: ", rare_path)

acmg_path <- file.path(out_dir, "Sample3.acmg_results.csv")
write.csv(rare_detailed, acmg_path, row.names = FALSE)
message("Wrote: ", acmg_path)

priority <- rare_detailed[order(
  match(rare_detailed$classification, c("Pathogenic", "Likely Pathogenic", "VUS", "Likely Benign", "Benign")),
  rare_detailed$max_population_af
), , drop = FALSE]

top <- priority[priority$classification %in% c("Pathogenic", "Likely Pathogenic", "VUS"), , drop = FALSE]
if (nrow(top) > 200L) top <- top[seq_len(200L), , drop = FALSE]

top_path <- file.path(out_dir, "Sample3.top_candidates.csv")
write.csv(top, top_path, row.names = FALSE)
message("Wrote: ", top_path)

criteria_summary <- acmg_pro_criteria_summary(rare_detailed)
summary_path <- file.path(out_dir, "Sample3.acmg_criteria_summary.csv")
write.csv(criteria_summary, summary_path, row.names = FALSE)
message("Wrote: ", summary_path)

cat("\nSummary\n")
cat("-------\n")
cat("Input variants:     ", nrow(variants), "\n", sep = "")
cat("Rare variants:      ", nrow(rare_detailed), "\n", sep = "")
cat("ACMG classified:    ", nrow(rare_detailed), "\n", sep = "")
print(table(rare_detailed$classification, useNA = "ifany"))
cat("\nCriteria triggered (rare set):\n")
print(criteria_summary)
cat("\nDone.\n")
