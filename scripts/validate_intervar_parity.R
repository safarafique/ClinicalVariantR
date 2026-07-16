#!/usr/bin/env Rscript
# Quick validation of InterVar-parity improvements (transcript priority, PS1/PM5/PM1, PVS1).

script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
setwd(project_root)
source("global_cli.R")

test_folder <- normalizePath(Sys.getenv("ACMG_TEST_FOLDER", file.path(project_root, "..", "testig", "testig")),
                             winslash = "/", mustWork = FALSE)

cat("=== ClinicalVariantR InterVar-parity validation ===\n\n")

# 1) TCF4 stop_gained transcript priority (sample 041 chr18:55631301)
vcf_041 <- file.path(test_folder, "260324100041.haplotypecaller_VEP.ann.split.vcf")
if (file.exists(vcf_041)) {
  row <- parse_variant_from_vcf_fields("chr18", 55631301, "G", "A", info = {
    con <- file(vcf_041, "r")
    on.exit(close(con), add = TRUE)
    info <- NA_character_
    repeat {
      line <- readLines(con, n = 1L, warn = FALSE)
      if (length(line) == 0L) break
      if (grepl("\t55631301\t", line, fixed = TRUE)) {
        info <- strsplit(line, "\t", fixed = TRUE)[[1L]][[8L]]
        break
      }
    }
    info
  })
  cat("[TCF4 transcript] gene=", row$gene, " consequence=", row$consequence,
      " protein_coding=", row$is_protein_coding, "\n", sep = "")
  scored <- score_variants_table(row, profile_id = "general_germline")
  cat("[TCF4 classification] ", scored$classification[[1L]],
      " criteria=", scored$criteria_met[[1L]], "\n\n", sep = "")
}

# 2) MECP2 missense PM1/PM5 (sample 041 chrX:154030906)
if (file.exists(vcf_041)) {
  row2 <- parse_variant_from_vcf_fields("chrX", 154030906, "T", "C", info = {
    con <- file(vcf_041, "r")
    on.exit(close(con), add = TRUE)
    info <- NA_character_
    repeat {
      line <- readLines(con, n = 1L, warn = FALSE)
      if (length(line) == 0L) break
      if (grepl("\t154030906\t", line, fixed = TRUE)) {
        info <- strsplit(line, "\t", fixed = TRUE)[[1L]][[8L]]
        break
      }
    }
    info
  })
  cat("[MECP2 transcript] gene=", row2$gene, " hgvs_p=", row2$hgvs_p,
      " aa=", row2$amino_acids, " pos=", row2$protein_position, "\n", sep = "")
  scored2 <- score_variants_table(row2, profile_id = "general_germline")
  cat("[MECP2 classification] ", scored2$classification[[1L]],
      " criteria=", scored2$criteria_met[[1L]], "\n\n", sep = "")
}

# 3) Full benchmark if test folder present
if (dir.exists(test_folder)) {
  out_dir <- file.path(project_root, "..", "results", "intervar_compare_v2")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  res <- benchmark_testing_folder(test_folder, profile_id = "general_germline", output_dir = out_dir)
  print(res$summary)
  cat("\nWrote: ", file.path(out_dir, "benchmark_summary.csv"), "\n", sep = "")
} else {
  cat("Test folder not found: ", test_folder, "\n", sep = "")
}

cat("\nDone.\n")
