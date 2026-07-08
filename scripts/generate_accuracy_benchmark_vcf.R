#!/usr/bin/env Rscript
# Build a compact accuracy-check VCF + reference TSV for ACMGamp.
#
# Usage:
#   Rscript scripts/generate_accuracy_benchmark_vcf.R
#
# Output:
#   testig/acmgamp_benchmark/acmgamp_accuracy_benchmark.vcf
#   testig/acmgamp_benchmark/acmgamp_accuracy_benchmark.acmg.tsv
#   testig/acmgamp_benchmark/acmgamp_accuracy_benchmark.key.csv

script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
testig_dir <- normalizePath(file.path(project_root, "..", "testig", "testig"), mustWork = TRUE)
out_dir <- normalizePath(file.path(project_root, "..", "testig", "acmgamp_benchmark"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_col <- "BENCHMARK_SAMPLE"

extract_vcf_line <- function(vcf_path, pattern) {
  con <- file(vcf_path, "r")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (grepl(pattern, line, perl = TRUE)) return(line)
  }
  stop("No line matching: ", pattern, " in ", vcf_path)
}

extract_csq_header <- function(vcf_path) {
  con <- file(vcf_path, "r")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (grepl("^##INFO=<ID=CSQ", line, fixed = TRUE)) return(line)
    if (grepl("^#CHROM\t", line)) break
  }
  stop("CSQ header not found in ", vcf_path)
}

normalize_sample_column <- function(line, sample_col) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
  if (length(parts) < 10L) {
    stop("Expected FORMAT + genotype columns in: ", substr(line, 1L, 80L))
  }
  parts[[9L]] <- "GT:AD:DP:GQ:PL"
  parts[[10L]] <- "0/1:10,2:12:54:54,0,403"
  if (length(parts) > 10L) parts <- parts[1L:10L]
  parts[[10L]] <- sub("^[^:]+$", "0/1:10,2:12:54:54,0,403", parts[[10L]])
  paste(parts, collapse = "\t")
}

vcf_041 <- file.path(testig_dir, "260324100041.haplotypecaller_VEP.ann.split.vcf")
vcf_039 <- file.path(testig_dir, "260324100039.haplotypecaller_VEP.ann.split.vcf")
vcf_042 <- file.path(testig_dir, "260324100042.haplotypecaller_VEP.ann.split.vcf")

variants <- list(
  list(
    file = vcf_041, pattern = "\t55631301\t",
    id = "TCF4_stop_gained", gene = "TCF4", note = "LoF on canonical transcript; must NOT call VUS"
  ),
  list(
    file = vcf_041, pattern = "\t154030906\t",
    id = "MECP2_missense", gene = "MECP2", note = "PM1 domain + PM5 + PS4; must be Likely Pathogenic+"
  ),
  list(
    file = vcf_039, pattern = "\t126544289\t",
    id = "ALDH7A1_common", gene = "ALDH7A1", note = "Common population variant; must be Benign (BA1)"
  ),
  list(
    file = vcf_039, pattern = "\t51687181\t",
    id = "SCN8A_synonymous", gene = "SCN8A", note = "Synonymous + common; must be Benign"
  ),
  list(
    file = vcf_039, pattern = "\t62694680\t",
    id = "BSCL2_rare_missense", gene = "BSCL2", note = "Rare missense PM2 only; expect VUS"
  ),
  list(
    file = vcf_042, pattern = "\t193614765\t",
    id = "OPA1_synonymous", gene = "OPA1", note = "Synonymous rare; expect Likely Benign or Benign"
  ),
  list(
    file = vcf_039, pattern = "\t73223597\t",
    id = "PSEN1_UTR_indel", gene = "PSEN1", note = "Rare UTR indel; expect VUS"
  ),
  list(
    file = vcf_039, pattern = "\t132892556\t",
    id = "TSC1_common", gene = "TSC1", note = "Common 3'UTR; must be Benign"
  )
)

header_lines <- c(
  "##fileformat=VCFv4.2",
  "##FILTER=<ID=PASS,Description=\"All filters passed\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
  "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
  "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
  "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">",
  extract_csq_header(vcf_041),
  paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", sample_col)
)

data_lines <- vapply(variants, function(v) {
  line <- extract_vcf_line(v$file, v$pattern)
  normalize_sample_column(line, sample_col)
}, FUN.VALUE = character(1L))

vcf_out <- file.path(out_dir, "acmgamp_accuracy_benchmark.vcf")
writeLines(c(header_lines, data_lines), vcf_out, useBytes = TRUE)

# Expected ACMGamp v2.3 classifications (ground truth for this benchmark).
truth <- data.frame(
  chr = c("chr18", "chrX", "chr5", "chr12", "chr11", "chr3", "chr14", "chr9"),
  pos = c(55631301L, 154030906L, 126544289L, 51687181L, 62694680L, 193614765L, 73223597L, 132892556L),
  ref = c("G", "T", "G", "C", "G", "A", "TA", "C"),
  alt = c("A", "C", "A", "T", "T", "G", "T", "T"),
  gene_symbol_base = c("TCF4", "MECP2", "ALDH7A1", "SCN8A", "BSCL2", "OPA1", "PSEN1", "TSC1"),
  effect_base = c(
    "stop_gained", "missense_variant", "3_prime_UTR_variant", "synonymous_variant",
    "missense_variant", "synonymous_variant", "3_prime_UTR_variant", "3_prime_UTR_variant"
  ),
  acmg_classification_base = c(
    "Likely Pathogenic", "Likely Pathogenic", "Benign", "Benign",
    "VUS", "Likely Benign", "VUS", "Benign"
  ),
  acmg_criteria_base = c(
    "PVS1,PM2", "PM1,PM2,PM5,PS4", "BA1,BP4,BP6", "BA1,BP4,BP6,BP7",
    "PM2", "PM2,BP4,BP7", "PM2", "BA1,BP4,BP6"
  ),
  benchmark_id = vapply(variants, `[[`, character(1L), "id"),
  benchmark_note = vapply(variants, `[[`, character(1L), "note"),
  stringsAsFactors = FALSE
)

tsv_out <- file.path(out_dir, "acmgamp_accuracy_benchmark.acmg.tsv")
write.table(truth, tsv_out, sep = "\t", row.names = FALSE, quote = FALSE)

key_out <- file.path(out_dir, "acmgamp_accuracy_benchmark.key.csv")
write.csv(truth[, c("benchmark_id", "chr", "pos", "ref", "alt", "gene_symbol_base",
                    "acmg_classification_base", "acmg_criteria_base", "benchmark_note")],
          key_out, row.names = FALSE)

cat("Wrote benchmark VCF:", vcf_out, "\n")
cat("Wrote reference TSV:", tsv_out, "\n")
cat("Wrote key:", key_out, "\n")
cat("Variants:", nrow(truth), "\n")
