#' Parse one variant from VCF and re-run ClinicalVariantR scoring (prediction curation).
#' @noRd

parse_single_variant_from_vcf <- function(vcf_path, chrom, pos, ref, alt, pass_only = FALSE) {
  if (!file.exists(vcf_path)) stop("VCF not found: ", vcf_path)

  con <- if (grepl("\\.gz$", vcf_path, ignore.case = TRUE)) gzfile(vcf_path, "rt") else file(vcf_path, "rt")
  on.exit(close(con), add = TRUE)

  header_cols <- NULL
  buffer_n <- if (exists("VCF_LINE_BUFFER", inherits = TRUE)) VCF_LINE_BUFFER else 50000L
  repeat {
    lines <- readLines(con, n = buffer_n, warn = FALSE)
    if (length(lines) == 0L) break
    for (line in lines) {
      if (grepl("^#CHROM\t", line)) {
        header_cols <- strsplit(sub("^#", "", line), "\t")[[1]]
        next
      }
      if (grepl("^#", line)) next
      parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(parts) < 8L) next
      if (isTRUE(pass_only) && !(parts[[7L]] %in% c("PASS", "."))) next

      parsed <- parse_vcf_line(line, header_cols)
      if (is.null(parsed)) next
      hit <- match_variant_rows(parsed, chrom, pos, ref, alt)
      if (nrow(hit) >= 1L) return(hit[1L, , drop = FALSE])
    }
  }
  NULL
}

rescore_variant_with_manual <- function(
    variant_row,
    manual_inputs = list(),
    clinical_context = NULL,
    pedigree_context = NULL,
    profile_id = DEFAULT_PROFILE_ID,
    refs = NULL,
    evidence_scope = "full") {

  if (is.null(variant_row) || nrow(variant_row) == 0L) {
    stop("Variant row missing for re-score.")
  }
  scored <- score_variants_table(
    variant_row,
    manual_inputs = normalize_manual_inputs(manual_inputs),
    clinical_context = clinical_context,
    pedigree_context = pedigree_context,
    profile_id = profile_id,
    refs = refs,
    evidence_scope = evidence_scope
  )
  if (nrow(scored) < 1L) stop("Re-score produced no output.")
  scored[1, , drop = FALSE]
}

rescore_variant_to_report_row <- function(
    scored_row,
    mode = c("full", "rapid"),
    session_id = NA_character_,
    run_metadata = NULL) {

  mode <- match.arg(mode)
  report <- acmg_pro_to_report(scored_row, mode = mode, session_id = session_id, run_metadata = run_metadata)
  if (nrow(report) != 1L) report[1, , drop = FALSE] else report
}

patch_report_row <- function(report_df, new_row) {
  if (is.null(report_df) || nrow(report_df) == 0L) return(new_row)
  key <- new_row$variant_id[1]
  if (!"variant_id" %in% names(report_df)) {
    report_df$variant_id <- variant_key_chr_pos_ref_alt(
      report_df$chrom, report_df$pos, report_df$ref, report_df$alt
    )
  }
  idx <- which(report_df$variant_id == key)
  for (col in names(new_row)) {
    if (!col %in% names(report_df)) report_df[[col]] <- NA
  }
  if (length(idx) == 0L) {
    return(rbind(report_df, new_row))
  }
  report_df[idx[1L], names(new_row)] <- new_row[1, names(new_row), drop = FALSE]
  report_df
}

write_report_csv <- function(report_df, csv_path) {
  utils::write.table(report_df, file = csv_path, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
  invisible(csv_path)
}

refresh_report_views <- function(full_df, selected_category) {
  if (is.null(full_df) || nrow(full_df) == 0L) {
    return(list(full = full_df, filtered = full_df))
  }
  if (is.null(selected_category) || !nzchar(selected_category)) {
    return(list(full = full_df, filtered = full_df[0, , drop = FALSE]))
  }
  sub <- full_df[full_df$classification == selected_category, , drop = FALSE]
  list(full = full_df, filtered = sub)
}
