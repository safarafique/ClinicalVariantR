#!/usr/bin/env Rscript
script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
setwd(project_root)
source("global_cli.R")

vcf_path <- normalizePath(
  file.path(project_root, "..", "testig", "testig", "260324100042.haplotypecaller_VEP.ann.split.vcf"),
  winslash = "/", mustWork = TRUE
)

cat("VCF:", vcf_path, "\n")
refs <- tryCatch(load_reference_data(), error = function(e) NULL)

result <- tryCatch(
  analyze_complete_vcf(
    vcf_path = vcf_path,
    mode = "rapid",
    pass_only = FALSE,
    chunk_size = 100L,
    use_bcftools = FALSE,
    refs = refs,
    session_id = "repro",
    write_audit = FALSE
  ),
  error = function(e) {
    cat("\nERROR:", conditionMessage(e), "\n")
    traceback()
    quit(status = 1)
  }
)

cat("OK: classified", result$rows_classified, "variants\n")
