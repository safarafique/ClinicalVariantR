#!/usr/bin/env Rscript
# Validate ClinicalVariantR Publishable MVP on HMC-1 ANN and Sample 3 CSQ fixtures.

script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
data_root <- normalizePath(file.path(project_root, ".."), winslash = "/", mustWork = FALSE)
setwd(project_root)
source("global.R")

state <- new.env(parent = emptyenv())
state$pass_n <- 0L
state$fail_n <- 0L

check <- function(name, cond, msg = "") {
  if (isTRUE(cond)) {
    state$pass_n <- state$pass_n + 1L
    cat("[PASS] ", name, "\n", sep = "")
  } else {
    state$fail_n <- state$fail_n + 1L
    cat("[FAIL] ", name, if (nzchar(msg)) paste0(": ", msg) else "", "\n", sep = "")
  }
}

cat("=== ClinicalVariantR Publishable MVP Validation ===\n\n")

cfg <- load_rule_config("general_germline")
check("Rule config loads", length(cfg$thresholds) > 0)
check("Criteria metadata loads", length(cfg$criteria_meta) > 0)
check("Disease profiles load", nrow(cfg$profiles_table) >= 1)

test_row <- data.frame(
  variant_id = "17:41234470:A:G", chrom = "17", pos = 41234470L,
  ref = "A", alt = "G", gene = "BRCA1", consequence = "missense_variant",
  af_1000g = 0.000003, population_af = 0.000003, gnomad_af = 0.000003,
  revel_score = 0.82, cadd = 25, spliceai_max = 0.02,
  alphamissense_score = 0.7, sift = "deleterious", polyphen = "probably_damaging",
  clinvar_classification = "Pathogenic", annotation_source = "CSQ",
  genome_build_hint = "GRCh38", stringsAsFactors = FALSE
)
scored <- score_variants_table(test_row, profile_id = "general_germline")
check("PM2 transparent rationale", grepl("Threshold", scored$PM2_rationale[1]))
check("Confidence score present", !is.na(scored$confidence_score[1]))
check("Evidence JSON present", nzchar(scored$evidence_json[1]))
check("Prediction summary present", nzchar(scored$prediction_scores[1]))

hmc1_vcf <- file.path(data_root, "results", "HMC-1.coding.vcf.gz")
if (!file.exists(hmc1_vcf)) hmc1_vcf <- file.path(data_root, "HMC-1.final.vcf")
if (file.exists(hmc1_vcf)) {
  val <- validate_vcf(hmc1_vcf, mode = "rapid", sample_rows = 20L)
  check("HMC-1 validation ready", isTRUE(val$can_analyze), val$summary)
  line <- NULL
  con <- open_vcf_connection(hmc1_vcf)
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0) break
    if (!grepl("^#", line)) break
  }
  if (!is.null(line)) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    row <- parse_variant_from_vcf_fields(parts[1], parts[2], parts[4], strsplit(parts[5], ",")[[1]][1], info = parts[8])
    check("HMC-1 parser uses ANN", identical(row$annotation_source, "ANN"))
    s <- score_variants_table(row, profile_id = "hematologic_predisposition")
    check("HMC-1 scoring works", !is.na(s$classification[1]))
  }
} else {
  cat("[SKIP] HMC-1 VCF not found\n")
}

sample3_vcf <- file.path(data_root, "Sample3.haplotypecaller.filtered_VEP.ann (2).vcf")
if (file.exists(sample3_vcf)) {
  con_hdr <- file(sample3_vcf, "r")
  header_line <- NULL
  first_variant <- NULL
  repeat {
    l <- readLines(con_hdr, n = 1L, warn = FALSE)
    if (length(l) == 0) break
    if (grepl("^#CHROM\t", l)) header_line <- l
    if (!grepl("^#", l)) { first_variant <- l; break }
  }
  close(con_hdr)
  has_csq <- (!is.null(first_variant) && grepl("CSQ=", first_variant, fixed = TRUE)) ||
    (!is.null(header_line) && grepl("CSQ", header_line, fixed = TRUE))
  check("Sample 3 has CSQ annotation", has_csq)
  con3 <- file(sample3_vcf, "r")
  on.exit(close(con3), add = TRUE)
  line3 <- NULL
  repeat {
    line3 <- readLines(con3, n = 1L, warn = FALSE)
    if (length(line3) == 0) break
    if (!grepl("^#", line3)) break
  }
  if (!is.null(line3)) {
    parts3 <- strsplit(line3, "\t", fixed = TRUE)[[1]]
    row3 <- parse_variant_from_vcf_fields(parts3[1], parts3[2], parts3[4], strsplit(parts3[5], ",")[[1]][1], info = parts3[8])
    check("Sample 3 parser uses CSQ", identical(row3$annotation_source, "CSQ"))
    s3 <- score_variants_table(row3, profile_id = "hematologic_predisposition")
    check("Sample 3 scoring works", !is.na(s3$classification[1]))
    check("Sample 3 has SIFT/PolyPhen fields", nzchar(s3$sift[1]) || nzchar(s3$polyphen[1]))
  }
} else {
  cat("[SKIP] Sample 3 VCF not found\n")
}

cat("\nSummary: ", state$pass_n, " passed, ", state$fail_n, " failed\n", sep = "")
if (state$fail_n > 0) quit(status = 1)
