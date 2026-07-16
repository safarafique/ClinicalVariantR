#!/usr/bin/env Rscript
# Verify Group B/C automated criteria coverage and Group C gene filter parity.
#
# Usage (WSL):
#   cd /mnt/e/ACGM/ClinicalVariantR && Rscript scripts/verify_group_b_c.R

source("global_cli.R")

cat("=== Group B / Group C automated criteria audit ===\n\n")
cat("Engine:", ACMG_PRO_ENGINE, "\n")
cat("Automated (18):", paste(AUTOMATED_ACMG_CRITERIA, collapse = ", "), "\n")
cat("Not in B/C (manual/context):", paste(c(CONTEXT_ASSISTED_CRITERIA, MANUAL_ONLY_CRITERIA), collapse = ", "), "\n\n")

state <- new.env(parent = emptyenv())
state$pass <- 0L
state$fail <- 0L
chk <- function(name, cond, detail = "") {
  if (isTRUE(cond)) {
    state$pass <- state$pass + 1L
    cat("[PASS]", name, "\n")
  } else {
    state$fail <- state$fail + 1L
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

# Rapid mode exports all 18 evaluated criteria in evidence JSON
s <- score_variants_table(base_row(revel_score = 0.9, cadd = 25, sift = "deleterious",
                                   polyphen = "probably_damaging", polyphen_score = 0.99),
                         evidence_scope = "automated")
ev <- parse_evidence_json(s$evidence_json[1])
chk("Group B evidence JSON has 18 rows", nrow(ev) == length(AUTOMATED_ACMG_CRITERIA),
    paste("got", nrow(ev)))
missing <- setdiff(AUTOMATED_ACMG_CRITERIA, ev$criterion)
chk("All automated criteria present in evidence", length(missing) == 0L,
    paste(missing, collapse = ", "))

# Group C uses same engine after gene filter
myc_row <- base_row(gene = "MYC", all_genes = "MYC;BCR", consequence = "missense_variant")
filtered <- filter_variants_by_genes(myc_row, "MYC")
chk("Group C gene filter keeps MYC variant", nrow(filtered) == 1L)
s_c <- score_variants_table(filtered, evidence_scope = "automated")
ev_c <- parse_evidence_json(s_c$evidence_json[1])
chk("Group C scoring matches B evidence scope", nrow(ev_c) == length(AUTOMATED_ACMG_CRITERIA))

# Clinical context must not apply in default rapid scoring (no clinical passed)
ctx <- score_variants_table(
  base_row(gene = "BCR"),
  clinical_context = data.frame(phenotype = "CML", cml_phase = "chronic", stringsAsFactors = FALSE),
  pedigree_context = data.frame(relation = "proband", affected_status = "affected", stringsAsFactors = FALSE),
  evidence_scope = "automated"
)
chk("Rapid mode does not auto-trigger PP4 from clinical when scope=automated",
    !isTRUE(ctx$PP4[1]))
chk("Rapid mode does not auto-trigger PS2 from pedigree when scope=automated",
    !isTRUE(ctx$PS2[1]))

# Benchmark if available
bench_vcf <- normalizePath(file.path("..", "testig", "clinicalvariantr_benchmark", "clinicalvariantr_group_b_benchmark.vcf"),
                           mustWork = FALSE)
bench_tsv <- normalizePath(file.path("..", "testig", "clinicalvariantr_benchmark", "clinicalvariantr_group_b_benchmark.acmg.tsv"),
                           mustWork = FALSE)
if (file.exists(bench_vcf) && file.exists(bench_tsv)) {
  cat("\n=== Group B benchmark (20 variants) ===\n")
  res <- benchmark_one_sample(bench_vcf, bench_tsv, profile_id = "general_germline")
  if (!is.null(res$error)) {
    cat("Benchmark error:", res$error, "\n")
    state$fail <- state$fail + 1L
  } else {
    cat("Exact accuracy:", round(100 * res$metrics$exact_accuracy, 1), "%\n")
    cat("Tier accuracy:", round(100 * res$metrics$tier_accuracy, 1), "%\n")
    mism <- res$comparison[!res$comparison$exact_match, , drop = FALSE]
    if (nrow(mism) > 0L) {
      apply(mism, 1L, function(r) {
        cat("  MISMATCH:", r[["variant_key"]], "| expected", r[["ref_class"]],
            "| got", r[["pred_class"]], "\n")
      })
      state$fail <- state$fail + 1L
    } else {
      cat("All 20 benchmark variants matched.\n")
      state$pass <- state$pass + 1L
    }
  }
}

cat("\nSummary:", state$pass, "passed,", state$fail, "failed\n")
if (state$fail > 0L) quit(status = 1)
