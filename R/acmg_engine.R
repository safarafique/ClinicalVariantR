#' ClinicalVariantR engine: structured evidence and streaming integration.

acmg_pro_to_report <- function(
    scored_df,
    mode = c("full", "rapid"),
    session_id = NA_character_,
    run_metadata = NULL) {

  mode <- match.arg(mode)
  run_ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  pipeline_mode <- if (mode == "full") {
    "Group A (Prediction + curation)"
  } else {
    "Group B/C (Automated prediction)"
  }

  meta_row <- if (!is.null(run_metadata)) metadata_to_single_row(run_metadata) else NULL

  data.frame(
    variant_id = scored_df$variant_id,
    chrom = scored_df$chrom,
    pos = scored_df$pos,
    ref = scored_df$ref,
    alt = scored_df$alt,
    gene = scored_df$gene,
    consequence = scored_df$consequence,
    annotation_source = scored_df$annotation_source %||% NA_character_,
    genome_build_hint = scored_df$genome_build_hint %||% NA_character_,
    gnomad_af = scored_df$max_population_af %||% scored_df$gnomad_af %||% scored_df$population_af,
    revel_score = scored_df$revel_score %||% scored_df$REVEL,
    cadd_score = scored_df$cadd %||% NA_real_,
    spliceai_max = scored_df$spliceai_max %||% NA_real_,
    alphamissense_score = scored_df$alphamissense_score %||% NA_real_,
    prediction_scores = scored_df$prediction_scores %||% "",
    clinvar_classification = scored_df$clinvar_classification %||% scored_df$ClinVar,
    criteria_met = scored_df$criteria_met,
    criteria_rationale = scored_df$criteria_rationale,
    evidence_summary = scored_df$evidence_summary %||% "",
    evidence_json = scored_df$evidence_json %||% "[]",
    classification = scored_df$classification,
    confidence_score = scored_df$confidence_score %||% NA_integer_,
    confidence_label = scored_df$confidence_label %||% NA_character_,
    evidence_strength = scored_df$evidence_strength %||% NA_character_,
    pathogenic_evidence_count = scored_df$pathogenic_evidence_count %||% NA_integer_,
    benign_evidence_count = scored_df$benign_evidence_count %||% NA_integer_,
    prediction_limitations = scored_df$prediction_limitations %||% "",
    disease_profile = scored_df$disease_profile %||% DEFAULT_PROFILE_ID,
    pipeline_mode = pipeline_mode,
    classified_at = run_ts,
    analyst_session = session_id,
    engine = ACMG_PRO_ENGINE,
    app_version = if (!is.null(meta_row)) meta_row$app_version else APP_VERSION,
    acmg_guideline_version = if (!is.null(meta_row)) meta_row$acmg_guideline_version else ACMG_GUIDELINE_VERSION,
    input_vcf_checksum = if (!is.null(meta_row)) meta_row$input_vcf_checksum else NA_character_,
    stringsAsFactors = FALSE
  )
}

run_acmg_pro_chunk <- function(
    variants_df,
    mode = c("full", "rapid"),
    manual_inputs = list(),
    manual_by_variant = list(),
    clinical_context = NULL,
    pedigree_context = NULL,
    refs = NULL,
    session_id = NA_character_,
    profile_id = DEFAULT_PROFILE_ID,
    run_metadata = NULL,
    write_audit = TRUE) {

  mode <- match.arg(mode)
  if (nrow(variants_df) == 0) return(empty_report())

  scored <- score_variants_table(
    variants_df,
    manual_inputs = if (mode == "full") manual_inputs else list(),
    manual_by_variant = if (mode == "full") manual_by_variant else list(),
    clinical_context = clinical_context,
    pedigree_context = pedigree_context,
    refs = refs,
    profile_id = profile_id,
    evidence_scope = if (mode == "full") "full" else "automated"
  )
  report <- acmg_pro_to_report(scored, mode = mode, session_id = session_id, run_metadata = run_metadata)

  if (isTRUE(write_audit)) {
    audit_batch <- build_audit_entries_from_report(report, session_id = session_id)
    append_audit_log(audit_batch)
  }

  report
}

acmg_pro_criteria_summary <- function(scored_df) {
  data.frame(
    criterion = FULL_ACMG_CRITERIA,
    count = vapply(FULL_ACMG_CRITERIA, function(code) {
      if (code %in% names(scored_df)) sum(scored_df[[code]], na.rm = TRUE) else 0L
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}

parse_evidence_json <- function(json_text) {
  if (is.na(json_text) || !nzchar(json_text) || json_text == "[]") {
    return(data.frame(
      criterion = character(), triggered = logical(), status = character(),
      observed_value = character(), threshold = character(),
      reason = character(), strength = character(), description = character(),
      stringsAsFactors = FALSE
    ))
  }
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    df <- tryCatch(
      jsonlite::fromJSON(json_text, simplifyDataFrame = TRUE),
      error = function(e) NULL
    )
    if (!is.null(df) && nrow(df) > 0) return(order_evidence_table(df))
  }
  data.frame(
    criterion = character(), triggered = logical(), status = character(),
    observed_value = character(), threshold = character(),
    reason = character(), strength = character(), description = character(),
    stringsAsFactors = FALSE
  )
}
