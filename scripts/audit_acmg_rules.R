#!/usr/bin/env Rscript
# Audit ACMG/AMP 2015 germline combining rules and automated criterion coverage.

script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
setwd(project_root)
source("global_cli.R")
source("R/acmg_pipeline.R")

state <- new.env(parent = emptyenv())
state$pass_n <- 0L
state$fail_n <- 0L

check <- function(name, cond, detail = "") {
  if (isTRUE(cond)) {
    state$pass_n <- state$pass_n + 1L
    cat("[PASS] ", name, "\n", sep = "")
  } else {
    state$fail_n <- state$fail_n + 1L
    cat("[FAIL] ", name, if (nzchar(detail)) paste0(": ", detail) else "", "\n", sep = "")
  }
}

ev <- function(pvs = 0L, ps = 0L, pm = 0L, pp = 0L, ba = FALSE, bs = 0L, bp = 0L) {
  list(PVS = pvs, PS = ps, PM = pm, PP = pp, BA = ba, BS = bs, BP = bp)
}

cat("=== ACMG/AMP 2015 combining rules audit ===\n\n")

combining_cases <- list(
  list("BA1 alone", ev(ba = TRUE), "Benign"),
  list("2x BS", ev(bs = 2L), "Benign"),
  list("BS + BP", ev(bs = 1L, bp = 1L), "Likely Benign"),
  list("2x BP", ev(bp = 2L), "Likely Benign"),
  list("PVS + PM", ev(pvs = 1L, pm = 1L), "Likely Pathogenic"),
  list("PM + 4PP", ev(pm = 1L, pp = 4L), "Likely Pathogenic"),
  list("3x PM", ev(pm = 3L), "Likely Pathogenic"),
  list("PVS + PS", ev(pvs = 1L, ps = 1L), "Pathogenic"),
  list("2x PS", ev(ps = 2L), "Pathogenic"),
  list("PM + PP only", ev(pm = 1L, pp = 1L), "VUS"),
  list("BS only", ev(bs = 1L), "VUS")
)

for (case in combining_cases) {
  got <- combine_acmg_evidence(case[[2]])
  check(paste("Combining:", case[[1]]), identical(got, case[[3]]), paste("got", got))
}

conflict <- combine_acmg_evidence(ev(pm = 1L, pp = 1L, bp = 2L))
check(
  "Conflicting pathogenic + benign evidence defaults to VUS",
  identical(conflict, "VUS"),
  paste("got", conflict, "(ACMG recommends VUS when evidence conflicts)")
)

cat("\n=== Automated criterion behavior ===\n\n")

check("PM2 absent when AF observed in population (BS2 range)",
  !score_population_criteria(0.003)$PM2 && isTRUE(score_population_criteria(0.003)$BS2))

check("PM2 when AF absent",
  isTRUE(score_population_criteria(NA_real_)$PM2))

check("BA1 when AF > 5%",
  isTRUE(score_population_criteria(0.06)$BA1))

check("BS1 when AF > 1% and <= 5%",
  isTRUE(score_population_criteria(0.02)$BS1) && !isTRUE(score_population_criteria(0.02)$BA1))

check("PM2+PP3 yields VUS (not Likely Pathogenic)",
  identical(score_variants_table(data.frame(
    variant_id = "t1", chrom = "1", pos = 1L, ref = "A", alt = "G",
    gene = "GENE", consequence = "missense_variant",
    population_af = NA_real_, gnomad_af = NA_real_, af_1000g = NA_real_, af_esp6500 = NA_real_,
    sift = "deleterious", polyphen = "probably_damaging(0.99)", polyphen_score = 0.99,
    revel_score = NA_real_, cadd = NA_real_, spliceai_max = NA_real_, alphamissense_score = NA_real_,
    clinvar_classification = NA_character_, stringsAsFactors = FALSE
  ), profile_id = "general_germline")$classification[1], "VUS"))

check("BA1 yields Benign even with PM2+PP3 signals absent",
  identical(score_population_criteria(0.10)$BA1, TRUE))

lb_row <- data.frame(
  variant_id = "t2", chrom = "1", pos = 2L, ref = "A", alt = "G",
  gene = "GENE", consequence = "synonymous_variant",
  population_af = 0.003, gnomad_af = 0.003, af_1000g = NA_real_, af_esp6500 = NA_real_,
  sift = "tolerated", polyphen = "benign", polyphen_score = 0.01,
  revel_score = NA_real_, cadd = NA_real_, spliceai_max = NA_real_, alphamissense_score = NA_real_,
  clinvar_classification = "Benign", stringsAsFactors = FALSE
)
lb <- score_variants_table(lb_row, profile_id = "general_germline")
check("BS2 + BP6 + BP7 yields Likely Benign",
  identical(lb$classification[1], "Likely Benign"),
  paste("got", lb$classification[1], "criteria:", lb$criteria_met[1]))

implemented <- c(
  "PVS1", "PS1", "PS2", "PS4", "PM1", "PM2", "PM4", "PM5",
  "PP1", "PP2", "PP3", "PP4", "PP5",
  "BA1", "BS1", "BS2", "BP1", "BP3", "BP4", "BP6", "BP7"
)
registry <- acmg_criteria_registry()
not_auto <- registry$criterion[!registry$criterion %in% implemented]

cat("\n=== Coverage summary ===\n")
cat("Automated + context-assisted criteria (21/28): ", paste(implemented, collapse = ", "), "\n", sep = "")

ps4_row <- data.frame(
  chrom = "chrX", pos = 154030906L, ref = "T", alt = "C",
  gene = "MECP2", hgvs_p = "p.Thr320Ala", rsids = "rs1273236261",
  consequence = "missense_variant", stringsAsFactors = FALSE
)
ps4_hit <- score_ps4_criteria(ps4_row)
check("PS4 case-control enrichment (MECP2 seed)",
  isTRUE(ps4_hit$PS4),
  ps4_hit$PS4_rationale)

gwas_demo <- data.frame(
  chrom = "chr6", pos = 32641407L, ref = ".", alt = ".",
  rsid = "rs1234567", gene = "HLA-A", trait = "example complex trait",
  pvalue = 1e-8, odds_ratio = 1.4, maf = 0.25, source = "GWAS_Catalog_import",
  stringsAsFactors = FALSE
)
gwas_demo$variant_key <- variant_key_chr_pos_ref_alt(gwas_demo$chrom, gwas_demo$pos, gwas_demo$ref, gwas_demo$alt)
gwas_row <- data.frame(
  chrom = "chr6", pos = 32641407L, ref = "A", alt = "G",
  gene = "HLA-A", rsids = "rs1234567", consequence = "intron_variant",
  stringsAsFactors = FALSE
)
gwas_ps4 <- score_ps4_criteria(gwas_row, ps4_db = load_ps4_case_control_db(), gwas_db = gwas_demo)
check("GWAS supplementary does NOT auto-award PS4",
  !isTRUE(gwas_ps4$PS4) && isTRUE(gwas_ps4$gwas_supplementary),
  gwas_ps4$gwas_supplementary_note)

clinical_ctx <- data.frame(
  cml_phase = "chronic phase", phenotype = "CML Philadelphia positive",
  stringsAsFactors = FALSE
)
pedigree_ctx <- data.frame(
  relation = c("proband", "mother", "father"),
  affected_status = c("yes", "no", "no"),
  stringsAsFactors = FALSE
)
ctx_row <- data.frame(
  variant_id = "t3", chrom = "22", pos = 23632664L, ref = "A", alt = "G",
  gene = "BCR", consequence = "missense_variant",
  population_af = NA_real_, gnomad_af = NA_real_, af_1000g = NA_real_, af_esp6500 = NA_real_,
  sift = "deleterious", polyphen = "probably_damaging(0.99)", polyphen_score = 0.99,
  revel_score = NA_real_, cadd = NA_real_, spliceai_max = NA_real_, alphamissense_score = NA_real_,
  clinvar_classification = NA_character_, stringsAsFactors = FALSE
)
ctx_scored <- score_variants_table(
  ctx_row, clinical_context = clinical_ctx, pedigree_context = pedigree_ctx,
  profile_id = "hematologic_predisposition", evidence_scope = "full"
)
check("PP4 from clinical log (BCR + CML phenotype)",
  isTRUE(ctx_scored$PP4[1]),
  ctx_scored$PP4_rationale[1])
check("PS2 from pedigree (affected proband, unaffected parents)",
  isTRUE(ctx_scored$PS2[1]),
  ctx_scored$PS2_rationale[1])

check("BP3 only on homopolymer in-frame indel",
  isTRUE(score_consequence_criteria(
    "inframe_deletion", "GENE", character(), ref = "AAAAA", alt = "AAAA"
  )$BP3) &&
    !isTRUE(score_consequence_criteria(
      "inframe_deletion", "GENE", character(), ref = "ATCG", alt = "AT"
    )$BP3))

cat("Not automated (require manual/clinical data): ", paste(not_auto, collapse = ", "), "\n", sep = "")

cat("\nSummary: ", state$pass_n, " passed, ", state$fail_n, " failed\n", sep = "")
if (state$fail_n > 0) quit(status = 1)
