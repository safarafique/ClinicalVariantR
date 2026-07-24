#' Structured ACMG evidence tables and confidence scoring.
#' @noRd

EVIDENCE_CRITERIA <- FULL_ACMG_CRITERIA

build_criterion_evidence <- function(
    criterion,
    triggered,
    observed_value,
    threshold_text,
    reason,
    strength = "Unknown",
    description = "",
    review_required = FALSE) {

  status <- if (isTRUE(triggered)) {
    "Triggered"
  } else if (isTRUE(review_required)) {
    "Review required"
  } else {
    "Not triggered"
  }

  data.frame(
    criterion = criterion,
    triggered = isTRUE(triggered),
    status = status,
    observed_value = as.character(observed_value %||% ""),
    threshold = threshold_text,
    reason = reason,
    strength = strength,
    description = description,
    stringsAsFactors = FALSE
  )
}

order_evidence_table <- function(evidence_df) {
  if (is.null(evidence_df) || nrow(evidence_df) == 0L) return(evidence_df)
  evidence_df$criterion <- factor(evidence_df$criterion, levels = FULL_ACMG_CRITERIA)
  evidence_df <- evidence_df[order(evidence_df$criterion), , drop = FALSE]
  evidence_df$criterion <- as.character(evidence_df$criterion)
  rownames(evidence_df) <- NULL
  evidence_df
}

build_variant_evidence_table <- function(scores, thresholds, criteria_meta = list(),
                                         evidence_scope = c("triggered", "automated", "full")) {
  evidence_scope <- match.arg(evidence_scope)
  include_all <- evidence_scope %in% c("automated", "full")
  af <- scores$max_population_af %||% NA_real_
  af_display <- if (is.na(af)) "absent/NA" else sprintf("%.8f", af)

  rows <- list(
    build_criterion_evidence(
      "BA1", scores$BA1, af_display,
      sprintf("> %.4f (%.0f%%)", thresholds$ba1_af, thresholds$ba1_af * 100),
      scores$BA1_rationale %||% "",
      get_criterion_strength("BA1", criteria_meta),
      get_criterion_description("BA1", criteria_meta)
    ),
    build_criterion_evidence(
      "BS1", scores$BS1, af_display,
      sprintf("> %.4f (%.0f%%)", thresholds$bs1_af, thresholds$bs1_af * 100),
      scores$BS1_rationale %||% "",
      get_criterion_strength("BS1", criteria_meta),
      get_criterion_description("BS1", criteria_meta)
    ),
    build_criterion_evidence(
      "BS2", scores$BS2, af_display,
      sprintf(">= %.4f and < %.4f", thresholds$bs2_af %||% 0.001, thresholds$bs1_af),
      scores$BS2_rationale %||% "",
      get_criterion_strength("BS2", criteria_meta),
      get_criterion_description("BS2", criteria_meta)
    ),
    build_criterion_evidence(
      "BS3", scores$BS3, "",
      "Well-established functional studies show no damaging effect",
      scores$BS3_rationale %||% "",
      get_criterion_strength("BS3", criteria_meta),
      get_criterion_description("BS3", criteria_meta),
      review_required = include_all && !isTRUE(scores$BS3)
    ),
    build_criterion_evidence(
      "BS4", scores$BS4, "",
      "Lack of segregation in affected family members",
      scores$BS4_rationale %||% "",
      get_criterion_strength("BS4", criteria_meta),
      get_criterion_description("BS4", criteria_meta),
      review_required = include_all && !isTRUE(scores$BS4) && nzchar(scores$BS4_rationale %||% "")
    ),
    build_criterion_evidence(
      "PM2", scores$PM2, af_display,
      sprintf("< %.4f (%.0f%%)", thresholds$pm2_af, thresholds$pm2_af * 100),
      scores$PM2_rationale %||% "",
      get_criterion_strength("PM2", criteria_meta),
      get_criterion_description("PM2", criteria_meta)
    ),
    build_criterion_evidence(
      "PM3", scores$PM3, "",
      "Detected in trans with pathogenic variant (recessive/compound het)",
      scores$PM3_rationale %||% "",
      get_criterion_strength("PM3", criteria_meta),
      get_criterion_description("PM3", criteria_meta),
      review_required = include_all && !isTRUE(scores$PM3)
    ),
    build_criterion_evidence(
      "PVS1", scores$PVS1, scores$PVS1_rationale %||% "",
      "LoF in LoF-mechanism gene panel or canonical protein-coding transcript",
      scores$PVS1_rationale %||% "",
      get_criterion_strength("PVS1", criteria_meta),
      get_criterion_description("PVS1", criteria_meta)
    ),
    build_criterion_evidence(
      "PS1", scores$PS1, scores$PS1_rationale %||% "",
      "Same amino acid change as established pathogenic variant",
      scores$PS1_rationale %||% "",
      get_criterion_strength("PS1", criteria_meta),
      get_criterion_description("PS1", criteria_meta)
    ),
    build_criterion_evidence(
      "PS2", scores$PS2, "",
      "De novo occurrence in patient with disease and no family history",
      scores$PS2_rationale %||% "",
      get_criterion_strength("PS2", criteria_meta),
      get_criterion_description("PS2", criteria_meta)
    ),
    build_criterion_evidence(
      "PS3", scores$PS3, "",
      "Well-established functional studies show damaging effect",
      scores$PS3_rationale %||% "",
      get_criterion_strength("PS3", criteria_meta),
      get_criterion_description("PS3", criteria_meta),
      review_required = include_all && !isTRUE(scores$PS3) && !nzchar(scores$PS3_rationale %||% "")
    ),
    build_criterion_evidence(
      "PS4", scores$PS4,
      if (isTRUE(scores$PS4)) scores$PS4_rationale %||% "" else scores$gwas_supplementary_note %||% "",
      "Variant-specific case-control enrichment (GWAS is supplementary review only)",
      if (isTRUE(scores$PS4)) {
        scores$PS4_rationale %||% ""
      } else if (isTRUE(scores$gwas_supplementary)) {
        paste0(scores$gwas_supplementary_note %||% "", " [Supplementary - PS4 not auto-applied]")
      } else {
        scores$PS4_rationale %||% ""
      },
      get_criterion_strength("PS4", criteria_meta),
      get_criterion_description("PS4", criteria_meta)
    ),
    build_criterion_evidence(
      "PM1", scores$PM1, "",
      "Missense in curated PM1 hotspot/critical-domain gene panel",
      scores$PM1_rationale %||% "",
      get_criterion_strength("PM1", criteria_meta),
      get_criterion_description("PM1", criteria_meta)
    ),
    build_criterion_evidence(
      "PM4", scores$PM4, "",
      "In-frame indel or stop-loss consequence",
      scores$PM4_rationale %||% "",
      get_criterion_strength("PM4", criteria_meta),
      get_criterion_description("PM4", criteria_meta)
    ),
    build_criterion_evidence(
      "PM5", scores$PM5, scores$PM5_rationale %||% "",
      "Missense at residue with different established pathogenic change",
      scores$PM5_rationale %||% "",
      get_criterion_strength("PM5", criteria_meta),
      get_criterion_description("PM5", criteria_meta)
    ),
    build_criterion_evidence(
      "PM6", scores$PM6, "",
      "Assumed de novo without parental confirmation",
      scores$PM6_rationale %||% "",
      get_criterion_strength("PM6", criteria_meta),
      get_criterion_description("PM6", criteria_meta),
      review_required = include_all && !isTRUE(scores$PM6) && nzchar(scores$PM6_rationale %||% "")
    ),
    build_criterion_evidence(
      "PP1", scores$PP1, "",
      "Co-segregation with disease in multiple affected family members",
      scores$PP1_rationale %||% "",
      get_criterion_strength("PP1", criteria_meta),
      get_criterion_description("PP1", criteria_meta)
    ),
    build_criterion_evidence(
      "PP2", scores$PP2, "",
      "Missense mechanism gene panel",
      scores$PP2_rationale %||% "",
      get_criterion_strength("PP2", criteria_meta),
      get_criterion_description("PP2", criteria_meta)
    ),
    build_criterion_evidence(
      "PP3", scores$PP3, scores$insilico_summary %||% "",
      sprintf(
        "REVEL>=%.2f; CADD>=%.0f; SpliceAI>=%.2f; AlphaMissense>=%.3f; >=%d tools",
        thresholds$revel_pp3 %||% 0.75,
        thresholds$cadd_pp3 %||% 20,
        thresholds$spliceai_pp3 %||% 0.5,
        thresholds$alphamissense_pp3 %||% 0.564,
        as.integer(thresholds$insilico_min_tools %||% 2),
        if (isTRUE(thresholds$insilico_revel_solo %||% TRUE)) "; REVEL solo allowed" else ""
      ),
      scores$PP3_rationale %||% "",
      get_criterion_strength("PP3", criteria_meta),
      get_criterion_description("PP3", criteria_meta)
    ),
    build_criterion_evidence(
      "PP4", scores$PP4, "",
      "Patient phenotype highly specific for gene/disease",
      scores$PP4_rationale %||% "",
      get_criterion_strength("PP4", criteria_meta),
      get_criterion_description("PP4", criteria_meta)
    ),
    build_criterion_evidence(
      "PP5", scores$PP5, "",
      "ClinVar pathogenic without conflict",
      scores$PP5_rationale %||% "",
      get_criterion_strength("PP5", criteria_meta),
      get_criterion_description("PP5", criteria_meta)
    ),
    build_criterion_evidence(
      "BP1", scores$BP1, "",
      "Truncating mechanism gene panel",
      scores$BP1_rationale %||% "",
      get_criterion_strength("BP1", criteria_meta),
      get_criterion_description("BP1", criteria_meta)
    ),
    build_criterion_evidence(
      "BP2", scores$BP2, "",
      "Observed in trans with pathogenic variant",
      scores$BP2_rationale %||% "",
      get_criterion_strength("BP2", criteria_meta),
      get_criterion_description("BP2", criteria_meta),
      review_required = include_all && !isTRUE(scores$BP2)
    ),
    build_criterion_evidence(
      "BP3", scores$BP3, "",
      "In-frame indel in non-functional region",
      scores$BP3_rationale %||% "",
      get_criterion_strength("BP3", criteria_meta),
      get_criterion_description("BP3", criteria_meta)
    ),
    build_criterion_evidence(
      "BP4", scores$BP4, scores$insilico_summary %||% "",
      sprintf(
        "REVEL<=%.2f; CADD<=%.0f; SpliceAI<=%.2f; AlphaMissense<=%.3f; >=%d tools",
        thresholds$revel_bp4 %||% 0.15,
        thresholds$cadd_bp4 %||% 10,
        thresholds$spliceai_bp4 %||% 0.1,
        thresholds$alphamissense_bp4 %||% 0.34,
        as.integer(thresholds$insilico_min_tools %||% 2),
        if (isTRUE(thresholds$insilico_revel_solo %||% TRUE)) "; REVEL solo allowed" else ""
      ),
      scores$BP4_rationale %||% "",
      get_criterion_strength("BP4", criteria_meta),
      get_criterion_description("BP4", criteria_meta)
    ),
    build_criterion_evidence(
      "BP6", scores$BP6, "",
      "ClinVar benign without conflict",
      scores$BP6_rationale %||% "",
      get_criterion_strength("BP6", criteria_meta),
      get_criterion_description("BP6", criteria_meta)
    ),
    build_criterion_evidence(
      "BP5", scores$BP5, "",
      "Reputable source reports variant as benign (non-ClinVar)",
      scores$BP5_rationale %||% "",
      get_criterion_strength("BP5", criteria_meta),
      get_criterion_description("BP5", criteria_meta),
      review_required = include_all && !isTRUE(scores$BP5)
    ),
    build_criterion_evidence(
      "BP7", scores$BP7, "",
      "Synonymous with no splice impact",
      scores$BP7_rationale %||% "",
      get_criterion_strength("BP7", criteria_meta),
      get_criterion_description("BP7", criteria_meta)
    )
  )

  tbl <- do.call(rbind, rows)
  if (identical(evidence_scope, "full")) {
    # keep all 28 rows
  } else if (identical(evidence_scope, "automated")) {
    tbl <- tbl[tbl$criterion %in% AUTOMATED_ACMG_CRITERIA, , drop = FALSE]
  } else {
    tbl <- tbl[tbl$criterion %in% c(AUTOMATED_ACMG_CRITERIA, CONTEXT_ASSISTED_CRITERIA), , drop = FALSE]
  }
  order_evidence_table(tbl)
}

evidence_table_to_json <- function(evidence_df, evidence_scope = c("triggered", "automated", "full")) {
  evidence_scope <- match.arg(evidence_scope)
  if (nrow(evidence_df) == 0) return("[]")
  export_df <- if (identical(evidence_scope, "full")) {
    evidence_df
  } else if (identical(evidence_scope, "automated")) {
    evidence_df[evidence_df$criterion %in% AUTOMATED_ACMG_CRITERIA, , drop = FALSE]
  } else {
    evidence_df[evidence_df$triggered, , drop = FALSE]
  }
  if (nrow(export_df) == 0) return("[]")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(jsonlite::toJSON(export_df, dataframe = "rows", auto_unbox = TRUE))
  }
  rows <- vapply(seq_len(nrow(export_df)), function(i) {
    sprintf(
      '{"criterion":"%s","status":"%s","observed_value":"%s","threshold":"%s","reason":"%s","strength":"%s"}',
      export_df$criterion[i], export_df$status[i],
      gsub('"', "'", export_df$observed_value[i]),
      gsub('"', "'", export_df$threshold[i]),
      gsub('"', "'", export_df$reason[i]),
      export_df$strength[i]
    )
  }, FUN.VALUE = character(1))
  paste0("[", paste(rows, collapse = ","), "]")
}

compute_confidence_score <- function(scores, thresholds = NULL) {
  if (is.null(thresholds)) {
    thresholds <- list(
      confidence_base = 50, confidence_per_criterion = 8, confidence_conflict_penalty = 15
    )
  }
  base <- thresholds$confidence_base %||% 50
  per_crit <- thresholds$confidence_per_criterion %||% 8
  conflict_pen <- thresholds$confidence_conflict_penalty %||% 15

  path_codes <- c(
    "PVS1", "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5"
  )
  benign_codes <- c("BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7")

  path_n <- sum(vapply(path_codes, function(c) isTRUE(scores[[c]]), logical(1L)))
  benign_n <- sum(vapply(benign_codes, function(c) isTRUE(scores[[c]]), logical(1L)))
  total_met <- path_n + benign_n

  score <- base + total_met * per_crit
  if (path_n > 0 && benign_n > 0) score <- score - conflict_pen
  if (isTRUE(scores$BA1)) score <- max(score - 20, 0)

  missing_penalty <- 0L
  if (is.na(scalar_num(scores$max_population_af %||% NA_real_))) missing_penalty <- missing_penalty + 5L
  if (!nzchar(scalar_chr(scores$insilico_summary %||% "", default = ""))) missing_penalty <- missing_penalty + 5L
  score <- score - missing_penalty

  score <- max(0, min(100, round(score)))
  list(
    confidence_score = score,
    confidence_label = if (score >= 80) "High" else if (score >= 60) "Moderate" else "Low",
    pathogenic_criteria_count = path_n,
    benign_criteria_count = benign_n,
    has_conflicting_evidence = path_n > 0 && benign_n > 0
  )
}

format_evidence_summary_text <- function(evidence_df) {
  triggered <- evidence_df[evidence_df$triggered, , drop = FALSE]
  if (nrow(triggered) == 0) return("No automated ACMG criteria triggered.")
  paste(
    vapply(seq_len(nrow(triggered)), function(i) {
      sprintf(
        "%s (%s): %s | Observed: %s | Threshold: %s",
        triggered$criterion[i],
        triggered$strength[i],
        triggered$reason[i],
        triggered$observed_value[i],
        triggered$threshold[i]
      )
    }, FUN.VALUE = character(1)),
    collapse = "\n"
  )
}

build_prediction_scores_row <- function(row) {
  data.frame(
    revel = row$revel_score %||% row$REVEL %||% NA_real_,
    cadd = row$cadd %||% NA_real_,
    spliceai_max = row$spliceai_max %||% NA_real_,
    alphamissense = row$alphamissense_score %||% NA_real_,
    polyphen = row$polyphen %||% NA_character_,
    polyphen_score = row$polyphen_score %||% NA_real_,
    sift = row$sift %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

prediction_scores_to_text <- function(row) {
  ps <- build_prediction_scores_row(row)
  parts <- c()
  if (!is.na(ps$revel)) parts <- c(parts, sprintf("REVEL=%.3f", ps$revel))
  if (!is.na(ps$cadd)) parts <- c(parts, sprintf("CADD=%.1f", ps$cadd))
  if (!is.na(ps$spliceai_max)) parts <- c(parts, sprintf("SpliceAI=%.3f", ps$spliceai_max))
  if (!is.na(ps$alphamissense)) parts <- c(parts, sprintf("AlphaMissense=%.3f", ps$alphamissense))
  if (!is.na(ps$polyphen_score)) {
    parts <- c(parts, sprintf("PolyPhen=%.3f", ps$polyphen_score))
  } else if (!is.na(ps$polyphen) && nzchar(ps$polyphen)) {
    parts <- c(parts, sprintf("PolyPhen=%s", ps$polyphen))
  }
  if (!is.na(ps$sift) && nzchar(ps$sift)) parts <- c(parts, sprintf("SIFT=%s", ps$sift))
  if (length(parts) == 0) "No prediction scores available" else paste(parts, collapse = "; ")
}
