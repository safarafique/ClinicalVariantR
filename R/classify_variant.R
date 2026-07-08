#' Canonical list of all 28 ACMG/AMP criteria (Richards et al., 2015).
FULL_ACMG_CRITERIA <- c(
  "PVS1", "PS1", "PS2", "PS3", "PS4",
  "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
  "PP1", "PP2", "PP3", "PP4", "PP5",
  "BA1", "BS1", "BS2", "BS3", "BS4",
  "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"
)

#' Criteria automated from VCF fields in Group B / Group C (rapid mode).
AUTOMATED_ACMG_CRITERIA <- c(
  "PVS1", "PS1", "PS4",
  "PM1", "PM2", "PM4", "PM5",
  "PP2", "PP3", "PP5",
  "BA1", "BS1", "BS2",
  "BP1", "BP3", "BP4", "BP6", "BP7"
)

CONTEXT_ASSISTED_CRITERIA <- c("PS2", "PP1", "PP4")
MANUAL_ONLY_CRITERIA <- c("PS3", "PM3", "PM6", "BS3", "BS4", "BP2", "BP5")

#' ACMG/AMP evidence combining rules (Richards et al., 2015; updated framework).
#'
#' @param evidence List with counts: PVS, PS, PM, PP, BA, BS, BP
#' @return Character classification tier.
safe_int <- function(x, default = 0L) {
  if (is.null(x) || length(x) == 0L) return(default)
  val <- suppressWarnings(as.integer(x)[1L])
  if (length(val) == 0L || is.na(val)) default else val
}

combine_acmg_evidence <- function(evidence) {
  pvs <- safe_int(evidence$PVS)
  ps  <- safe_int(evidence$PS)
  pm  <- safe_int(evidence$PM)
  pp  <- safe_int(evidence$PP)
  ba  <- isTRUE(evidence$BA)
  bs  <- safe_int(evidence$BS)
  bp  <- safe_int(evidence$BP)

  path_count <- pvs + ps + pm + pp
  benign_count <- bs + bp

  if (isTRUE(ba)) {
    return("Benign")
  }

  if (path_count > 0L && benign_count > 0L) {
    return("VUS")
  }

  if (bs >= 2L) {
    return("Benign")
  }
  if ((bs == 1L && bp >= 1L) || bp >= 2L) {
    return("Likely Benign")
  }

  pathogenic <- (
    (pvs >= 1L && (ps >= 1L || pm >= 2L || (pm >= 1L && pp >= 1L) || pp >= 2L)) ||
    ps >= 2L ||
    (ps >= 1L && pm >= 3L) ||
    (ps >= 1L && pm >= 2L && pp >= 2L) ||
    (ps >= 1L && pm >= 1L && pp >= 4L)
  )
  if (pathogenic) return("Pathogenic")

  likely_pathogenic <- (
    (pvs >= 1L && pm >= 1L) ||
    (ps >= 1L && pm >= 1L && pm <= 2L) ||
    (ps >= 1L && pp >= 2L) ||
    pm >= 3L ||
    (pm >= 2L && pp >= 2L) ||
    (pm >= 1L && pp >= 4L)
  )
  if (likely_pathogenic) return("Likely Pathogenic")

  "VUS"
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

#' Coerce annotation fields to a single character value (avoids if()/&& length errors).
scalar_chr <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0L) return(default)
  if (is.list(x) && !is.data.frame(x)) x <- x[[1L]]
  val <- as.character(x)[1L]
  if (length(val) == 0L || is.na(val) || !nzchar(val)) return(default)
  val
}

scalar_num <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0L) return(default)
  if (is.list(x) && !is.data.frame(x)) x <- x[[1L]]
  val <- suppressWarnings(as.numeric(x)[1L])
  if (length(val) == 0L || is.na(val)) return(default)
  val
}

scalar_lgl <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L) return(default)
  if (is.list(x) && !is.data.frame(x)) x <- x[[1L]]
  isTRUE(x[1L])
}

#' Ensure variant scoring receives exactly one row with scalar columns.
normalize_variant_row_input <- function(row) {
  row <- as.data.frame(row, stringsAsFactors = FALSE)
  if (nrow(row) < 1L) stop("Empty variant row.")
  if (nrow(row) > 1L) row <- row[1L, , drop = FALSE]
  for (nm in names(row)) {
    val <- row[[nm]]
    if (is.list(val) && !is.data.frame(val)) {
      row[[nm]] <- if (length(val) > 0L) val[[1L]] else NA
    } else if (length(val) > 1L) {
      row[[nm]] <- val[[1L]]
    }
  }
  row
}

#' Summarize met criteria into strength counts.
summarize_evidence <- function(criteria_df) {
  met <- criteria_df[!is.na(criteria_df$met) & criteria_df$met, , drop = FALSE]
  strengths <- met$strength

  list(
    PVS = sum(strengths == "PVS", na.rm = TRUE),
    PS  = sum(strengths == "PS", na.rm = TRUE),
    PM  = sum(strengths == "PM", na.rm = TRUE),
    PP  = sum(strengths == "PP", na.rm = TRUE),
    BA  = any(met$criterion == "BA1", na.rm = TRUE),
    BS  = sum(strengths == "BS", na.rm = TRUE),
    BP  = sum(strengths == "BP", na.rm = TRUE),
    criteria_met = paste(met$criterion, collapse = ";"),
    criteria_strength = paste(strengths, collapse = ";")
  )
}
