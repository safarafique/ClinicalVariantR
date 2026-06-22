#' ACMG/AMP evidence combining rules (Richards et al., 2015; updated framework).
#'
#' @param evidence List with counts: PVS, PS, PM, PP, BA, BS, BP
#' @return Character classification tier.
combine_acmg_evidence <- function(evidence) {
  pvs <- evidence$PVS %||% 0L
  ps  <- evidence$PS %||% 0L
  pm  <- evidence$PM %||% 0L
  pp  <- evidence$PP %||% 0L
  ba  <- evidence$BA %||% FALSE
  bs  <- evidence$BS %||% 0L
  bp  <- evidence$BP %||% 0L

  if (isTRUE(ba) || bs >= 2L) {
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

#' Summarize met criteria into strength counts.
summarize_evidence <- function(criteria_df) {
  met <- criteria_df[criteria_df$met, , drop = FALSE]
  strengths <- met$strength

  list(
    PVS = sum(strengths == "PVS"),
    PS  = sum(strengths == "PS"),
    PM  = sum(strengths == "PM"),
    PP  = sum(strengths == "PP"),
    BA  = any(met$criterion == "BA1" & met$met),
    BS  = sum(strengths == "BS"),
    BP  = sum(strengths == "BP"),
    criteria_met = paste(met$criterion, collapse = ";"),
    criteria_strength = paste(strengths, collapse = ";")
  )
}
