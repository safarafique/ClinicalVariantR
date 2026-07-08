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

coalesce_dt_col <- function(dt, target, sources) {
  existing <- intersect(sources, names(dt))
  if (length(existing) == 0L) {
    if (!target %in% names(dt)) dt[[target]] <- NA
    return(dt)
  }
  if (!target %in% names(dt)) {
    dt[[target]] <- dt[[existing[1]]]
    existing <- existing[-1]
  }
  for (src in existing) {
    dt[[target]] <- data.table::fifelse(is.na(dt[[target]]), dt[[src]], dt[[target]])
  }
  dt
}

#' Annotate variant table with reference placeholders via data.table joins.
annotate_variants <- function(variants_dt, refs) {
  dt <- data.table::as.data.table(variants_dt)
  if (!"variant_key" %in% names(dt)) {
    dt[, variant_key := paste(chrom, pos, ref, alt, sep = ":")]
  }

  dt <- merge(
    dt,
    refs$gnomad[, .(variant_key, ref_gnomad_af = AF, ref_gnomad_popmax = popmax_AF)],
    by = "variant_key",
    all.x = TRUE
  )
  dt <- merge(
    dt,
    refs$clinvar[, .(variant_key, ref_clinvar_classification = clinical_significance,
                     ref_clinvar_review_status = review_status)],
    by = "variant_key",
    all.x = TRUE
  )
  dt <- merge(
    dt,
    refs$revel[, .(variant_key, ref_revel_score = REVEL)],
    by = "variant_key",
    all.x = TRUE
  )

  dt <- coalesce_dt_col(dt, "gnomad_af", c("gnomad_af", "ref_gnomad_af", "population_af", "AF"))
  dt <- coalesce_dt_col(dt, "popmax_af", c("popmax_af", "ref_gnomad_popmax"))
  dt <- coalesce_dt_col(dt, "revel_score", c("revel_score", "ref_revel_score", "REVEL"))
  dt <- coalesce_dt_col(
    dt, "clinvar_classification",
    c("clinvar_classification", "ref_clinvar_classification", "ClinVar")
  )

  drop_cols <- intersect(
    c("ref_gnomad_af", "ref_revel_score",
      "ref_clinvar_classification", "ref_clinvar_review_status"),
    names(dt)
  )
  if (length(drop_cols) > 0) {
    dt[, (drop_cols) := NULL]
  }

  as.data.frame(dt)
}

dedupe_variants_by_key <- function(variants_df) {
  if (is.null(variants_df) || nrow(variants_df) == 0L) return(variants_df)
  dt <- data.table::as.data.table(variants_df)
  if (!"variant_key" %in% names(dt)) {
    dt[, variant_key := paste(chrom, pos, ref, alt, sep = ":")]
  }
  if (anyDuplicated(dt$variant_key) > 0L) {
    dt <- dt[!duplicated(variant_key)]
  }
  as.data.frame(dt)
}
