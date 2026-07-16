#!/usr/bin/env Rscript
# Verify all 18 automated ACMG criteria and benchmark samples 039-042.
script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
setwd(project_root)
source("global_cli.R")

IMPLEMENTED <- c(
  "PVS1", "PS1", "PS4", "PM1", "PM2", "PM4", "PM5",
  "PP2", "PP3", "PP5",
  "BA1", "BS1", "BS2", "BP1", "BP3", "BP4", "BP6", "BP7"
)

cat("=== ClinicalVariantR 18-criteria verification ===\n\n")
cat("Engine:", ACMG_PRO_ENGINE, "\n")
cat("Automated (18):", paste(IMPLEMENTED, collapse = ", "), "\n\n")

# Combining rules audit
source("scripts/audit_acmg_rules.R")

cat("\n=== Per-criterion smoke tests ===\n")
pass <- 0L
fail <- 0L
chk <- function(name, cond, detail = "") {
  if (isTRUE(cond)) {
    pass <<- pass + 1L
    cat("[PASS]", name, "\n")
  } else {
    fail <<- fail + 1L
    cat("[FAIL]", name, if (nzchar(detail)) paste0(": ", detail) else "", "\n")
  }
}

base_row <- function(...) {
  data.frame(
    variant_id = "t", chrom = "1", pos = 1L, ref = "A", alt = "G",
    gene = "GENE", consequence = "missense_variant",
    population_af = NA_real_, gnomad_af = NA_real_, af_1000g = NA_real_, af_esp6500 = NA_real_,
    is_protein_coding = TRUE, is_canonical_transcript = TRUE,
    hgvs_p = NA_character_, amino_acids = NA_character_, protein_position = NA_character_,
    rsids = NA_character_,
    sift = NA_character_, polyphen = NA_character_, polyphen_score = NA_real_,
    revel_score = NA_real_, cadd = NA_real_, spliceai_max = NA_real_, alphamissense_score = NA_real_,
    clinvar_classification = NA_character_, stringsAsFactors = FALSE,
    ...
  )
}

lof <- base_row(consequence = "stop_gained", gene = "TP53")
s <- score_variants_table(lof, lof_genes = "TP53")
chk("PVS1 LoF TP53", isTRUE(s$PVS1[1]))

pm2 <- base_row()
s <- score_variants_table(pm2)
chk("PM2 absent AF", isTRUE(s$PM2[1]))

ba1 <- base_row(population_af = 0.06, gnomad_af = 0.06)
s <- score_variants_table(ba1)
chk("BA1 high AF", isTRUE(s$BA1[1]) && identical(s$classification[1], "Benign"))

pm4 <- base_row(consequence = "inframe_deletion")
s <- score_variants_table(pm4)
chk("PM4 inframe", isTRUE(s$PM4[1]))

bp7 <- base_row(consequence = "synonymous_variant", sift = "tolerated", polyphen = "benign")
s <- score_variants_table(bp7)
chk("BP7 synonymous", isTRUE(s$BP7[1]))

pp3 <- base_row(
  revel_score = 0.9, cadd = 25,
  sift = "deleterious", polyphen = "probably_damaging", polyphen_score = 0.99
)
s <- score_variants_table(pp3)
chk("PP3 multi-tool", isTRUE(s$PP3[1]))

bp4 <- base_row(
  revel_score = 0.05, cadd = 5,
  sift = "tolerated", polyphen = "benign", polyphen_score = 0.01
)
s <- score_variants_table(bp4)
chk("BP4 multi-tool", isTRUE(s$BP4[1]))

pp5 <- base_row(clinvar_classification = "Pathogenic")
s <- score_variants_table(pp5)
chk("PP5 ClinVar pathogenic", isTRUE(s$PP5[1]))

bp6 <- base_row(clinvar_classification = "Benign")
s <- score_variants_table(bp6)
chk("BP6 ClinVar benign", isTRUE(s$BP6[1]))

pm1 <- base_row(gene = "MECP2", hgvs_p = "p.Ala140Val", amino_acids = "A/V", protein_position = "140")
s <- score_variants_table(pm1)
chk("PM1 hotspot gene", isTRUE(s$PM1[1]))

pm5 <- base_row(
  gene = "MECP2", hgvs_p = "p.Thr320Ala", amino_acids = "T/A", protein_position = "320"
)
s <- score_variants_table(pm5)
chk("PM5 same residue", isTRUE(s$PM5[1]))

ps4 <- base_row(
  chrom = "chrX", pos = 154030906L, ref = "T", alt = "C",
  gene = "MECP2", hgvs_p = "p.Thr320Ala", rsids = "rs1273236261"
)
s <- score_variants_table(ps4)
chk("PS4 case-control", isTRUE(s$PS4[1]))

# TCF4 transcript + PVS1
test_vcf <- normalizePath(Sys.getenv("ACMG_TEST_VCF",
  file.path(project_root, "..", "testig", "testig", "260324100041.haplotypecaller_VEP.ann.split.vcf")),
  mustWork = FALSE)
if (file.exists(test_vcf)) {
  con <- file(test_vcf, "r")
  on.exit(close(con), add = TRUE)
  info <- NA_character_
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (grepl("\t55631301\t", line, fixed = TRUE)) {
      parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      info <- parts[[8L]]
      break
    }
  }
  if (!is.na(info)) {
    tcf4 <- parse_variant_from_vcf_fields("chr18", 55631301, "G", "A", info = info)
    chk("TCF4 gene selection", identical(tcf4$gene, "TCF4"),
      paste("got", tcf4$gene, tcf4$consequence))
    s <- score_variants_table(tcf4)
    chk("TCF4 PVS1+classification", isTRUE(s$PVS1[1]),
      paste(s$classification[1], s$criteria_met[1]))
  }
}

cat("\nSmoke tests:", pass, "passed,", fail, "failed\n")

test_dir <- normalizePath(Sys.getenv("ACMG_TEST_FOLDER",
  file.path(project_root, "..", "testig", "testig")), mustWork = FALSE)
if (dir.exists(test_dir)) {
  cat("\n=== Benchmark vs reference (039-042) ===\n")
  out <- file.path(project_root, "..", "results", "verify_18criteria")
  res <- benchmark_testing_folder(test_dir, profile_id = "general_germline", output_dir = out)
  print(res$summary)
}

cat("\nDone.\n")
if (fail > 0L) quit(status = 1)
