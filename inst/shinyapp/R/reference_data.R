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

  # Normalize join keys (chr-prefixed for cross-source matching)
  for (ref_name in names(refs)) {
    if ("chrom" %in% names(refs[[ref_name]])) {
      data.table::set(refs[[ref_name]], j = "chrom", value = normalize_chrom(refs[[ref_name]]$chrom))
    }
    data.table::set(
      refs[[ref_name]],
      j = "variant_key",
      value = variant_key_chr_pos_ref_alt(
        refs[[ref_name]]$chrom,
        refs[[ref_name]]$pos,
        refs[[ref_name]]$ref,
        refs[[ref_name]]$alt
      )
    )
  }

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
  if ("chrom" %in% names(dt)) {
    data.table::set(dt, j = "chrom", value = normalize_chrom(dt$chrom))
  }
  data.table::set(dt, j = "variant_key", value = variant_key_chr_pos_ref_alt(dt$chrom, dt$pos, dt$ref, dt$alt))

  dt <- merge(
    dt,
    data.frame(
      variant_key = refs$gnomad$variant_key,
      ref_gnomad_af = refs$gnomad$AF,
      ref_gnomad_popmax = refs$gnomad$popmax_AF,
      stringsAsFactors = FALSE
    ),
    by = "variant_key",
    all.x = TRUE
  )
  dt <- merge(
    dt,
    data.frame(
      variant_key = refs$clinvar$variant_key,
      ref_clinvar_classification = refs$clinvar$clinical_significance,
      ref_clinvar_review_status = refs$clinvar$review_status,
      stringsAsFactors = FALSE
    ),
    by = "variant_key",
    all.x = TRUE
  )
  dt <- merge(
    dt,
    data.frame(
      variant_key = refs$revel$variant_key,
      ref_revel_score = refs$revel$REVEL,
      stringsAsFactors = FALSE
    ),
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
    dt[, drop_cols] <- NULL
  }

  as.data.frame(dt)
}

dedupe_variants_by_key <- function(variants_df) {
  if (is.null(variants_df) || nrow(variants_df) == 0L) return(variants_df)
  dt <- data.table::as.data.table(variants_df)
  if ("chrom" %in% names(dt)) {
    data.table::set(dt, j = "chrom", value = normalize_chrom(dt$chrom))
  }
  data.table::set(dt, j = "variant_key", value = variant_key_chr_pos_ref_alt(dt$chrom, dt$pos, dt$ref, dt$alt))
  if (anyDuplicated(dt$variant_key) > 0L) {
    dt <- dt[!duplicated(dt$variant_key)]
  }
  as.data.frame(dt)
}
