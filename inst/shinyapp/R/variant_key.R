#' Canonical chromosome prefix for variant identification joins.
#'
#' Reference databases and VCF exports may use `17` or `chr17`; this normalizes
#' to a single `chr`-prefixed form (with `chrM` for mitochondrial contigs).
#' @noRd
normalize_chrom <- function(chrom) {
  chrom <- trimws(as.character(chrom))
  if (length(chrom) == 0L) return(chrom)

  out <- chrom
  invalid <- is.na(chrom) | !nzchar(chrom) | chrom == "."
  lower <- tolower(chrom)

  mitochondrial <- !invalid & lower %in% c("mt", "chrm", "m")
  out[mitochondrial] <- "chrM"

  prefixed <- !invalid & !mitochondrial & grepl("^chr", lower, perl = TRUE)
  prefixed_idx <- which(prefixed)
  suffix <- sub("^chr", "", chrom[prefixed_idx], ignore.case = TRUE)
  out[prefixed_idx] <- paste0("chr", suffix)
  out[prefixed_idx[toupper(suffix) == "M"]] <- "chrM"

  unprefixed <- !invalid & !mitochondrial & !prefixed
  out[unprefixed] <- paste0("chr", chrom[unprefixed])
  out
}

#' First ALT allele from a VCF ALT field (comma-separated multi-allelic).
#' @noRd
first_alt_allele <- function(alt) {
  alt <- as.character(alt)
  if (length(alt) == 0L || is.na(alt) || !nzchar(alt)) return(".")
  parts <- strsplit(alt, ",", fixed = TRUE)[[1L]]
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) "." else parts[[1L]]
}

#' Split a VCF ALT field into individual allele strings.
#' @noRd
split_vcf_alt_alleles <- function(alt_raw) {
  if (is.na(alt_raw) || !nzchar(alt_raw) || alt_raw == ".") return(".")
  alts <- strsplit(as.character(alt_raw), ",", fixed = TRUE)[[1L]]
  alts <- alts[nzchar(alts)]
  if (length(alts) == 0L) "." else alts
}

#' Build normalized genomic variant key: chrom:pos:ref:alt.
#' @noRd
variant_key_chr_pos_ref_alt <- function(chrom, pos, ref, alt) {
  chrom <- normalize_chrom(chrom)
  pos <- as.character(pos)
  ref <- as.character(ref)
  alt1 <- vapply(
    strsplit(as.character(alt), ",", fixed = TRUE),
    function(a) if (length(a) == 0L || !nzchar(a[[1L]])) "." else a[[1L]],
    character(1L)
  )
  paste(chrom, pos, ref, alt1, sep = ":")
}

#' Alias for benchmark and reference loaders.
#' @noRd
variant_key_from_parts <- function(chrom, pos, ref, alt) {
  variant_key_chr_pos_ref_alt(chrom, pos, ref, alt)
}

#' Alias used by annotation and join pipelines.
#' @noRd
normalize_variant_key <- function(chrom, pos, ref, alt) {
  variant_key_chr_pos_ref_alt(chrom, pos, ref, alt)
}

#' Compare two variant coordinates using normalized keys.
#' @noRd
variant_coords_equal <- function(chrom1, pos1, ref1, alt1, chrom2, pos2, ref2, alt2) {
  identical(
    variant_key_chr_pos_ref_alt(chrom1, pos1, ref1, alt1),
    variant_key_chr_pos_ref_alt(chrom2, pos2, ref2, alt2)
  )
}

#' Select rows from a parsed variant table matching normalized coordinates.
#' @noRd
match_variant_rows <- function(df, chrom, pos, ref, alt) {
  if (is.null(df) || nrow(df) == 0L) return(df[0, , drop = FALSE])
  target <- variant_key_chr_pos_ref_alt(chrom, pos, ref, alt)
  keys <- variant_key_chr_pos_ref_alt(df$chrom, df$pos, df$ref, df$alt)
  df[keys == target, , drop = FALSE]
}
