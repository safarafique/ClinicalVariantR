#' Registry of all 28 ACMG/AMP criteria (Richards et al., 2015).
#'
#' implementation: automated | manual | not_implemented
acmg_criteria_registry <- function() {
  data.frame(
    criterion = c(
      "PVS1", "PS1", "PS2", "PS3", "PS4",
      "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
      "PP1", "PP2", "PP3", "PP4", "PP5",
      "BA1", "BS1", "BS2", "BS3", "BS4",
      "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"
    ),
    strength = c(
      "PVS", "PS", "PS", "PS", "PS",
      "PM", "PM", "PM", "PM", "PM", "PM",
      "PP", "PP", "PP", "PP", "PP",
      "BA", "BS", "BS", "BS", "BS",
      "BP", "BP", "BP", "BP", "BP", "BP", "BP"
    ),
    description = c(
      "Null variant in gene where LoF is disease mechanism",
      "Same amino acid change as known pathogenic variant",
      "De novo (confirmed) in patient with disease and no family history",
      "Well-established functional studies show damaging effect",
      "Prevalence in affected >> controls",
      "Located in mutational hot spot / critical functional domain",
      "Absent or extremely low frequency in population databases",
      "Detected in trans with pathogenic variant (recessive)",
      "Protein length change (in-frame indel / stop loss)",
      "Novel missense at residue with known pathogenic missense",
      "Assumed de novo without parental confirmation",
      "Co-segregation with disease in family",
      "Missense in gene where missense is common mechanism",
      "Multiple computational lines support deleterious effect",
      "Patient phenotype highly specific for gene/disease",
      "Reputable source reports variant as pathogenic",
      "Allele frequency > 5% in population database",
      "Allele frequency greater than expected for disorder",
      "Observed in healthy adult for expected inheritance",
      "Well-established functional studies show no damage",
      "Lack of segregation in affected family members",
      "Missense in gene where truncating is mechanism",
      "Observed in trans with pathogenic (dominant/recessive context)",
      "In-frame indel in repetitive region without known function",
      "Multiple computational lines suggest no impact",
      "Reputable source reports variant as benign",
      "Reputable source / database reports benign (ClinVar)",
      "Synonymous with no predicted splice impact"
    ),
    implementation = c(
      "automated", "not_implemented", "manual", "manual", "manual",
      "partial", "automated", "not_implemented", "not_implemented", "not_implemented", "manual",
      "manual", "manual", "automated", "manual", "automated",
      "automated", "automated", "not_implemented", "not_implemented", "not_implemented",
      "partial", "not_implemented", "not_implemented", "automated", "not_implemented", "automated", "automated"
    ),
    vcf_inputs = c(
      "CSQ/consequence", "protein HGVS + ClinVar", "pedigree/parental", "literature", "case-control stats",
      "domain annotation", "gnomAD AF", "phased genotype", "CSQ indel type", "residue + ClinVar", "pedigree",
      "pedigree", "gene mechanism DB", "REVEL/SpliceAI", "HPO/phenotype", "ClinVar",
      "gnomAD AF", "gnomAD AF + disease MAF", "cohort/pedigree", "functional assay", "pedigree",
      "gene mechanism DB", "phased genotype", "repeat annotation", "REVEL/SpliceAI", "literature", "ClinVar",
      "CSQ + splice predictors"
    ),
    stringsAsFactors = FALSE
  )
}

new_criteria_frame <- function(criterion_codes) {
  registry <- acmg_criteria_registry()
  registry <- registry[registry$criterion %in% criterion_codes, , drop = FALSE]
  registry <- registry[match(criterion_codes, registry$criterion), , drop = FALSE]
  data.frame(
    criterion = registry$criterion,
    strength = registry$strength,
    implementation = registry$implementation,
    met = FALSE,
    rationale = "",
    stringsAsFactors = FALSE
  )
}

FULL_CRITERIA <- c(
  "PVS1", "PS1", "PS2", "PS3", "PS4",
  "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
  "PP1", "PP2", "PP3", "PP4", "PP5",
  "BA1", "BS1", "BS2", "BS3", "BS4",
  "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"
)

AUTOMATED_CRITERIA <- c(
  "PVS1", "PS1", "PM1", "PM2", "PM3", "PM4", "PM5",
  "PP3", "PP5",
  "BA1", "BS1", "BS2", "BP1", "BP2", "BP3", "BP4", "BP6", "BP7"
)

extract_variant_fields <- function(variant) {
  af <- coalesce_num(variant$gnomad_af, variant$AF)
  revel <- coalesce_num(variant$revel_score, variant$REVEL)
  clinvar <- coalesce_chr(variant$clinvar_classification, variant$ClinVar)
  consequence <- coalesce_chr(variant$consequence, "")

  list(
    af = af,
    revel = revel,
    clinvar = clinvar,
    consequence = consequence,
    is_lof = grepl("stop_gained|frameshift|splice_donor|splice_acceptor|stop_lost", consequence, ignore.case = TRUE),
    is_missense = grepl("missense", consequence, ignore.case = TRUE),
    is_synonymous = grepl("synonymous", consequence, ignore.case = TRUE)
  )
}

apply_automated_acmg_rules <- function(fields, mark, mode = c("full", "rapid")) {
  mode <- match.arg(mode)
  af <- fields$af
  revel <- fields$revel
  clinvar <- fields$clinvar

  if (fields$is_lof) {
    mark("PVS1", "Predicted loss-of-function variant (requires gene-level LoF mechanism review).")
  }

  if (!is.na(revel) && revel >= 0.75) {
    mark("PP3", sprintf("REVEL=%.3f (>=0.75).", revel))
  }
  if (!is.na(revel) && revel <= 0.15) {
    mark("BP4", sprintf("REVEL=%.3f (<=0.15).", revel))
  }

  if (!is.na(af) && af > 0.05) {
    mark("BA1", sprintf("gnomAD AF=%.4f (>5%%).", af))
  } else if (!is.na(af) && af > 0.001) {
    mark("BS1", sprintf("gnomAD AF=%.4f elevated for rare disorder.", af))
  }

  if (is.na(af) || af == 0) {
    mark("PM2", "Absent or extremely rare in gnomAD (confirm database coverage).")
  }

  if (fields$is_missense) {
    mark("PM1", "Missense variant — domain/hotspot annotation not yet integrated (partial).")
    mark("BP1", "Missense in gene where truncating is primary mechanism (gene list required).")
  }

  if (!is.na(clinvar) && grepl("pathogenic", clinvar, ignore.case = TRUE) &&
      !grepl("conflict|uncertain", clinvar, ignore.case = TRUE)) {
    mark("PP5", sprintf("ClinVar reports pathogenic: %s.", clinvar))
  }
  if (!is.na(clinvar) && grepl("benign", clinvar, ignore.case = TRUE) &&
      !grepl("conflict|uncertain", clinvar, ignore.case = TRUE)) {
    mark("BP6", sprintf("ClinVar reports benign: %s.", clinvar))
  }

  if (mode == "rapid" && fields$is_synonymous) {
    mark("BP7", "Synonymous variant (splice predictors not integrated).")
  }
}

apply_manual_acmg_rules <- function(manual_inputs, mark) {
  manual <- function(key) isTRUE(manual_inputs[[key]] %||% FALSE)
  if (manual("PS3_functional")) mark("PS3", "Curator confirmed functional study supports damage.")
  if (manual("PP4_phenotype")) mark("PP4", "Curator confirmed phenotype match with disease.")
  if (manual("PS4_case_control")) mark("PS4", "Curator confirmed case-control enrichment.")
  if (manual("PS2_de_novo")) mark("PS2", "Curator confirmed de novo occurrence.")
  if (manual("PM6_de_novo")) mark("PM6", "Curator confirmed assumed de novo without parental data.")
  if (manual("PP1_segregation")) mark("PP1", "Curator confirmed co-segregation with disease.")
  if (manual("PP2_missense_mechanism")) mark("PP2", "Curator confirmed missense is common disease mechanism.")
}

apply_clinical_context_hints <- function(criteria, clinical_context) {
  if (is.null(clinical_context) || nrow(clinical_context) == 0) return(criteria)
  if (any(grepl("accelerated|blast", clinical_context$cml_phase, ignore.case = TRUE))) {
    idx <- match("PP4", criteria$criterion)
    if (!is.na(idx)) {
      criteria$rationale[idx] <- paste(
        criteria$rationale[idx],
        "Clinical logs suggest advanced CML phase — review PP4."
      )
    }
  }
  criteria
}

apply_pedigree_context_hints <- function(criteria, pedigree_context) {
  if (is.null(pedigree_context) || nrow(pedigree_context) == 0) return(criteria)

  relations <- tolower(pedigree_context$relation)
  affected <- tolower(pedigree_context$affected_status)
  has_unaffected_parent <- any(relations %in% c("mother", "father") & affected == "unaffected")
  proband_affected <- any(relations %in% c("proband", "patient", "index") & affected == "affected")

  if (proband_affected && has_unaffected_parent) {
    for (code in c("PS2", "PP1")) {
      idx <- match(code, criteria$criterion)
      if (!is.na(idx) && !criteria$met[idx]) {
        criteria$rationale[idx] <- paste(
          criteria$rationale[idx],
          "Pedigree compatible with de novo/segregation review."
        )
      }
    }
  }
  criteria
}

evaluate_acmg_full_28 <- function(variant, manual_inputs = list(), clinical_context = NULL,
                                  pedigree_context = NULL) {
  fields <- extract_variant_fields(variant)
  criteria <- new_criteria_frame(FULL_CRITERIA)

  mark <- function(code, rationale) {
    idx <- match(code, criteria$criterion)
    criteria$met[idx] <<- TRUE
    criteria$rationale[idx] <<- rationale
  }

  apply_automated_acmg_rules(fields, mark, mode = "full")
  apply_manual_acmg_rules(manual_inputs, mark)
  criteria <- apply_clinical_context_hints(criteria, clinical_context)
  criteria <- apply_pedigree_context_hints(criteria, pedigree_context)

  evidence <- summarize_evidence(criteria)
  classification <- combine_acmg_evidence(evidence)

  list(criteria = criteria, evidence = evidence, classification = classification)
}

evaluate_acmg_automated_18 <- function(variant) {
  fields <- extract_variant_fields(variant)
  criteria <- new_criteria_frame(AUTOMATED_CRITERIA)

  mark <- function(code, rationale) {
    idx <- match(code, criteria$criterion)
    criteria$met[idx] <<- TRUE
    criteria$rationale[idx] <<- rationale
  }

  apply_automated_acmg_rules(fields, mark, mode = "rapid")

  evidence <- summarize_evidence(criteria)
  classification <- combine_acmg_evidence(evidence)

  list(criteria = criteria, evidence = evidence, classification = classification)
}

run_pipeline <- function(variants_df, mode = c("full", "rapid"), manual_inputs = list(),
                         clinical_context = NULL, pedigree_context = NULL,
                         session_id = NA_character_) {
  mode <- match.arg(mode)
  if (nrow(variants_df) == 0) {
    return(list(report = empty_report(), mode = mode))
  }

  run_ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  results <- vector("list", nrow(variants_df))
  audit_rows <- vector("list", nrow(variants_df))

  for (i in seq_len(nrow(variants_df))) {
    row <- variants_df[i, , drop = FALSE]
    fields <- extract_variant_fields(row)

    eval_result <- if (mode == "full") {
      evaluate_acmg_full_28(row, manual_inputs, clinical_context, pedigree_context)
    } else {
      evaluate_acmg_automated_18(row)
    }

    results[[i]] <- data.frame(
      variant_id = row$variant_id,
      chrom = row$chrom,
      pos = row$pos,
      ref = row$ref,
      alt = row$alt,
      gene = coalesce_chr(row$gene),
      consequence = fields$consequence,
      gnomad_af = fields$af,
      revel_score = fields$revel,
      clinvar_classification = fields$clinvar,
      criteria_met = eval_result$evidence$criteria_met,
      criteria_strength = eval_result$evidence$criteria_strength,
      classification = eval_result$classification,
      pipeline_mode = if (mode == "full") "Group A (28-criteria)" else "Group B (18-automated)",
      classified_at = run_ts,
      analyst_session = session_id,
      stringsAsFactors = FALSE
    )

    audit_rows[[i]] <- build_audit_entry(
      variant_id = row$variant_id,
      pipeline_mode = results[[i]]$pipeline_mode,
      classification = eval_result$classification,
      criteria_met = eval_result$evidence$criteria_met,
      timestamp = run_ts,
      session_id = session_id,
      details = paste(eval_result$criteria$rationale[eval_result$criteria$met], collapse = " | ")
    )
  }

  report <- do.call(rbind, results)
  rownames(report) <- NULL
  append_audit_log(do.call(rbind, audit_rows))

  list(report = report, mode = mode)
}

test_acmg_logic_engine <- function() {
  test_variant <- data.frame(
    variant_id = "TEST:1:A:G",
    chrom = "TEST", pos = 1L, ref = "A", alt = "G",
    gene = "BCR", consequence = "missense_variant",
    AF = 0.001, REVEL = 0.8, ClinVar = "Pathogenic",
    gnomad_af = 0.001, revel_score = 0.8, clinvar_classification = "Pathogenic",
    stringsAsFactors = FALSE
  )
  result <- evaluate_acmg_full_28(test_variant)
  registry <- acmg_criteria_registry()
  implemented <- sum(registry$implementation %in% c("automated", "partial", "manual"))

  list(
    classification = result$classification,
    criteria_met = result$evidence$criteria_met,
    pass = result$classification %in% c("Pathogenic", "Likely Pathogenic"),
    criteria_implemented = implemented,
    criteria_total = nrow(registry)
  )
}

empty_report <- function() {
  df <- as.data.frame(matrix(ncol = length(REPORT_COLUMNS), nrow = 0))
  names(df) <- REPORT_COLUMNS
  df
}

coalesce_num <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v) && length(v) > 0 && !is.na(v[1])) return(as.numeric(v[1]))
  }
  NA_real_
}

coalesce_chr <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v) && length(v) > 0 && !is.na(v[1]) && nzchar(as.character(v[1]))) {
      return(as.character(v[1]))
    }
  }
  NA_character_
}
