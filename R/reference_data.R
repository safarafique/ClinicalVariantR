#' Load reference annotation placeholders for gnomAD, ClinVar, and REVEL.
#'
#' Production deployments should replace placeholder TSV files with curated
#' reference snapshots (gnomAD v4.1 allele frequencies, ClinVar submissions,
#' REVEL transcript scores).
load_reference_data <- function(paths = REFERENCE_PATHS) {
  refs <- list(
    gnomad = data.table::fread(paths$gnomad_v41, showProgress = FALSE),
    clinvar = data.table::fread(paths$clinvar, showProgress = FALSE),
    revel = data.table::fread(paths$revel, showProgress = FALSE)
  )

  # Normalize join keys
  refs$gnomad[, variant_key := paste(chrom, pos, ref, alt, sep = ":")]
  refs$clinvar[, variant_key := paste(chrom, pos, ref, alt, sep = ":")]
  refs$revel[, variant_key := paste(chrom, pos, ref, alt, sep = ":")]

  refs
}

#' Annotate variant table with reference placeholders via data.table joins.
annotate_variants <- function(variants_dt, refs) {
  dt <- data.table::as.data.table(variants_dt)
  if (!"variant_key" %in% names(dt)) {
    dt[, variant_key := paste(chrom, pos, ref, alt, sep = ":")]
  }

  dt <- merge(
    dt,
    refs$gnomad[, .(variant_key, gnomad_af = AF, gnomad_popmax = popmax_AF)],
    by = "variant_key",
    all.x = TRUE
  )
  dt <- merge(
    dt,
    refs$clinvar[, .(variant_key, clinvar_classification = clinical_significance,
                     clinvar_review_status = review_status)],
    by = "variant_key",
    all.x = TRUE
  )
  dt <- merge(
    dt,
    refs$revel[, .(variant_key, revel_score = REVEL)],
    by = "variant_key",
    all.x = TRUE
  )

  if ("AF" %in% names(dt)) {
    dt[, gnomad_af := data.table::fifelse(is.na(gnomad_af), AF, gnomad_af)]
  }
  if ("REVEL" %in% names(dt)) {
    dt[, revel_score := data.table::fifelse(is.na(revel_score), REVEL, revel_score)]
  }
  if ("ClinVar" %in% names(dt)) {
    dt[, clinvar_classification := data.table::fifelse(
      is.na(clinvar_classification), ClinVar, clinvar_classification
    )]
  }

  as.data.frame(dt)
}
