#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) stop("need data.table")
})
source("global.R")
vcf <- "../testig/clinicalvariantr_benchmark/clinicalvariantr_group_b_benchmark.vcf"
if (!file.exists(vcf)) vcf <- "testig/clinicalvariantr_benchmark/clinicalvariantr_group_b_benchmark.vcf"
stopifnot(file.exists(vcf))
p <- preview_vcf(vcf)
v <- validate_vcf(vcf, mode = "rapid")
cat("preview_rows:", p$preview_rows, "display:", p$total_display, "\n")
cat("valid:", v$valid, "can_analyze:", v$can_analyze, "\n")
cat("OK\n")
