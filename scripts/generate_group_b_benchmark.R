#!/usr/bin/env Rscript
# Build Group B efficiency benchmark VCF + reference labels (ACMGamp v2.3.2).
#
# Usage:
#   Rscript scripts/generate_group_b_benchmark.R
#
# Output (testig/acmgamp_benchmark/):
#   acmgamp_group_b_benchmark.vcf
#   acmgamp_group_b_benchmark.acmg.tsv
#   acmgamp_group_b_benchmark.key.csv

script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
testig_dir <- normalizePath(file.path(project_root, "..", "testig", "testig"), mustWork = TRUE)
out_dir <- normalizePath(file.path(project_root, "..", "testig", "acmgamp_benchmark"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_col <- "GROUP_B_BENCHMARK"

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

normalize_sample_column <- function(line) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
  if (length(parts) < 10L) {
    stop("Expected FORMAT + genotype columns in: ", substr(line, 1L, 80L))
  }
  parts[[9L]] <- "GT:AD:DP:GQ:PL"
  parts[[10L]] <- "0/1:10,2:12:54:54,0,403"
  if (length(parts) > 10L) parts <- parts[1L:10L]
  paste(parts, collapse = "\t")
}

vcf_039 <- file.path(testig_dir, "260324100039.haplotypecaller_VEP.ann.split.vcf")
vcf_041 <- file.path(testig_dir, "260324100041.haplotypecaller_VEP.ann.split.vcf")
vcf_042 <- file.path(testig_dir, "260324100042.haplotypecaller_VEP.ann.split.vcf")
vcf_043 <- file.path(testig_dir, "260324100043.clean.vcf")

variants <- list(
  # --- Pathogenic anchors (must detect) ---
  list(file = vcf_041, pattern = "\t55631301\t", id = "TCF4_stop_gained", gene = "TCF4",
       effect = "stop_gained", class = "Likely Pathogenic", criteria = "PVS1,PM2",
       category = "pathogenic_anchor",
       note = "LoF on canonical TCF4; false negative if VUS or Benign"),
  list(file = vcf_041, pattern = "\t154030906\t", id = "MECP2_missense", gene = "MECP2",
       effect = "missense_variant", class = "Pathogenic", criteria = "PS4,PM1,PM2,PM5,PP3",
       category = "pathogenic_anchor",
       note = "Rett hotspot missense; false negative if LP- or below"),
  list(file = vcf_043, pattern = "\t62813296\t", id = "EGR2_stop_gained", gene = "EGR2",
       effect = "stop_gained", class = "Likely Pathogenic", criteria = "PVS1,PM2",
       category = "pathogenic_anchor",
       note = "Rare stop_gained; false negative if VUS"),

  # --- Benign traps (must NOT over-call pathogenic) ---
  list(file = vcf_039, pattern = "\t126544289\t", id = "ALDH7A1_common", gene = "ALDH7A1",
       effect = "3_prime_UTR_variant", class = "Benign", criteria = "BA1,BP4,BP6",
       category = "benign_trap",
       note = "AF ~100%; false positive if Pathogenic"),
  list(file = vcf_039, pattern = "\t51687181\t", id = "SCN8A_synonymous", gene = "SCN8A",
       effect = "synonymous_variant", class = "Benign", criteria = "BA1,BP6,BP7",
       category = "benign_trap",
       note = "Common synonymous; false positive if Pathogenic"),
  list(file = vcf_039, pattern = "\t132892556\t", id = "TSC1_common", gene = "TSC1",
       effect = "3_prime_UTR_variant", class = "Benign", criteria = "BA1,BP6",
       category = "benign_trap",
       note = "Common TSC1 UTR; false positive if Pathogenic"),
  list(file = vcf_041, pattern = "\t100990864\t", id = "TWNK_common_splice", gene = "TWNK",
       effect = "splice_region_variant", class = "Benign", criteria = "BA1,BP6",
       category = "benign_trap",
       note = "AF >50% splice region; false positive if Pathogenic"),
  list(file = vcf_041, pattern = "\t108292760\t", id = "ATM_synonymous", gene = "ATM",
       effect = "synonymous_variant", class = "Benign", criteria = "BA1,BP6,BP7",
       category = "benign_trap",
       note = "Borderline AF ~5%; must stay Benign"),
  list(file = vcf_041, pattern = "\t23539813\t", id = "NPC1_splice_synonymous", gene = "NPC1",
       effect = "synonymous_variant", class = "Benign", criteria = "BA1,BP7",
       category = "splice_trap",
       note = "High AF splice-region synonymous; ClinVar conflict — must stay Benign"),
  list(file = vcf_041, pattern = "\t226881976\t", id = "PSEN2_synonymous", gene = "PSEN2",
       effect = "synonymous_variant", class = "Benign", criteria = "BA1,BP6,BP7",
       category = "benign_trap",
       note = "Common Alzheimer-gene SNP; false positive if Pathogenic"),

  # --- Likely benign ---
  list(file = vcf_042, pattern = "\t193614765\t", id = "OPA1_synonymous_rare", gene = "OPA1",
       effect = "synonymous_variant", class = "Likely Benign", criteria = "PM2,BP4,BP7",
       category = "likely_benign",
       note = "Rare synonymous with benign predictors"),
  list(file = vcf_041, pattern = "\t55224052\t", id = "TCF4_UTR_BS2", gene = "TCF4",
       effect = "3_prime_UTR_variant", class = "Likely Benign", criteria = "BS2,BP6",
       category = "likely_benign",
       note = "Observed in population with ClinVar benign support"),

  # --- VUS (conservative; should not be Pathogenic) ---
  list(file = vcf_039, pattern = "\t62694680\t", id = "BSCL2_rare_missense", gene = "BSCL2",
       effect = "missense_variant", class = "VUS", criteria = "PM2",
       category = "vus_pm2",
       note = "Rare missense PM2 only"),
  list(file = vcf_039, pattern = "\t73223597\t", id = "PSEN1_UTR_indel", gene = "PSEN1",
       effect = "3_prime_UTR_variant", class = "VUS", criteria = "PM2",
       category = "vus_pm2",
       note = "Rare UTR indel PM2 only"),
  list(file = vcf_041, pattern = "\t110155271\t", id = "COL4A1_rare_intron", gene = "COL4A1",
       effect = "intron_variant", class = "VUS", criteria = "PM2",
       category = "vus_pm2",
       note = "Rare intronic PM2 only"),
  list(file = vcf_041, pattern = "\t66508138\t", id = "TK2_rare_UTR", gene = "TK2",
       effect = "3_prime_UTR_variant", class = "VUS", criteria = "PM2",
       category = "vus_pm2",
       note = "Rare UTR insertion PM2 only"),
  list(file = vcf_041, pattern = "\t42914099\t", id = "G6PC1_BS1", gene = "G6PC1",
       effect = "3_prime_UTR_variant", class = "VUS", criteria = "BS1",
       category = "vus_mixed",
       note = "BS1 only — elevated but not BA1"),
  list(file = vcf_041, pattern = "\t55631298\t", id = "TCF4_missense_conflicting", gene = "TCF4",
       effect = "missense_variant", class = "VUS", criteria = "PM1,PM2,PP3,BP4",
       category = "overcall_trap",
       note = "Conflicting in silico on TCF4 missense; must NOT call Pathogenic"),
  list(file = vcf_041, pattern = "\t32156037\t", id = "SPAST_BS2", gene = "SPAST",
       effect = "3_prime_UTR_variant", class = "VUS", criteria = "BS2",
       category = "vus_mixed",
       note = "BS2 only in population"),
  list(file = vcf_041, pattern = "\t50623471\t", id = "ARSA_rare_UTR", gene = "ARSA",
       effect = "3_prime_UTR_variant", class = "VUS", criteria = "PM2",
       category = "vus_pm2",
       note = "Rare UTR SNV PM2 only")
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
  normalize_sample_column(extract_vcf_line(v$file, v$pattern))
}, FUN.VALUE = character(1L))

vcf_out <- file.path(out_dir, "acmgamp_group_b_benchmark.vcf")
writeLines(c(header_lines, data_lines), vcf_out, useBytes = TRUE)

truth <- do.call(rbind, lapply(variants, function(v) {
  parts <- strsplit(extract_vcf_line(v$file, v$pattern), "\t", fixed = TRUE)[[1L]]
  data.frame(
    chr = parts[[1L]],
    pos = as.integer(parts[[2L]]),
    ref = parts[[4L]],
    alt = strsplit(parts[[5L]], ",")[[1L]][1L],
    gene_symbol_base = v$gene,
    effect_base = v$effect,
    acmg_classification_base = v$class,
    acmg_criteria_base = v$criteria,
    benchmark_id = v$id,
    benchmark_category = v$category,
    benchmark_note = v$note,
    stringsAsFactors = FALSE
  )
}))

tsv_out <- file.path(out_dir, "acmgamp_group_b_benchmark.acmg.tsv")
write.table(truth, tsv_out, sep = "\t", row.names = FALSE, quote = FALSE)

key_out <- file.path(out_dir, "acmgamp_group_b_benchmark.key.csv")
write.csv(
  truth[, c("benchmark_id", "benchmark_category", "chr", "pos", "ref", "alt",
            "gene_symbol_base", "acmg_classification_base", "acmg_criteria_base", "benchmark_note")],
  key_out, row.names = FALSE
)

cat("Wrote Group B benchmark VCF:", vcf_out, "\n")
cat("Wrote reference TSV:", tsv_out, "\n")
cat("Wrote key:", key_out, "\n")
cat("Variants:", nrow(truth), "\n")
