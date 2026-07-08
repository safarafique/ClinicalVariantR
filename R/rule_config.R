#' Load configurable ACMG thresholds, criteria metadata, and disease profiles.

CONFIG_PATHS <- list(
  thresholds = file.path("config", "acmg_thresholds.csv"),
  criteria   = file.path("config", "acmg_criteria.csv"),
  profiles   = file.path("config", "disease_profiles.csv")
)

DEFAULT_PROFILE_ID <- "general_germline"

load_acmg_thresholds_table <- function(path = CONFIG_PATHS$thresholds) {
  if (!file.exists(path)) return(data.frame())
  read.csv(path, stringsAsFactors = FALSE)
}

load_acmg_criteria_table <- function(path = CONFIG_PATHS$criteria) {
  if (!file.exists(path)) return(data.frame())
  read.csv(path, stringsAsFactors = FALSE)
}

load_disease_profiles_table <- function(path = CONFIG_PATHS$profiles) {
  if (!file.exists(path)) return(data.frame())
  read.csv(path, stringsAsFactors = FALSE)
}

thresholds_table_to_list <- function(tbl) {
  fallback <- list(
    ba1_af = 0.05, bs1_af = 0.01, pm2_af = 0.01, pm2_strict_af = 0.001,
    revel_pp3 = 0.75, revel_bp4 = 0.15, cadd_pp3 = 20, cadd_bp4 = 10,
    spliceai_pp3 = 0.5, spliceai_bp4 = 0.1,
    alphamissense_pp3 = 0.564, alphamissense_bp4 = 0.34,
    polyphen_damaging = 0.85, polyphen_benign = 0.15,
    insilico_min_tools = 2, bs2_af = 0.001,
    confidence_base = 50, confidence_per_criterion = 8, confidence_conflict_penalty = 15
  )
  if (nrow(tbl) == 0L) return(fallback)
  keys <- tbl$threshold_key
  vals <- suppressWarnings(as.numeric(tbl$value))
  out <- as.list(setNames(vals, keys))
  out[!is.na(vals)]
}

load_rule_config <- function(profile_id = DEFAULT_PROFILE_ID) {
  thresh_tbl <- load_acmg_thresholds_table()
  criteria_tbl <- load_acmg_criteria_table()
  profiles_tbl <- load_disease_profiles_table()

  thresholds <- thresholds_table_to_list(thresh_tbl)

  if (nrow(profiles_tbl) > 0 && profile_id %in% profiles_tbl$profile_id) {
    prof <- profiles_tbl[profiles_tbl$profile_id == profile_id, , drop = FALSE][1, , drop = FALSE]
    for (col in c("ba1_af", "bs1_af", "pm2_af", "pm2_strict_af")) {
      if (col %in% names(prof) && !is.na(prof[[col]])) {
        num <- suppressWarnings(as.numeric(prof[[col]]))
        if (!is.na(num)) thresholds[[col]] <- num
      }
    }
  }

  criteria_meta <- if (nrow(criteria_tbl) > 0) {
    stats::setNames(
      lapply(seq_len(nrow(criteria_tbl)), function(i) as.list(criteria_tbl[i, , drop = FALSE])),
      criteria_tbl$criterion
    )
  } else {
    list()
  }

  list(
    profile_id = profile_id,
    profile_name = if (nrow(profiles_tbl) > 0 && profile_id %in% profiles_tbl$profile_id) {
      profiles_tbl$profile_name[profiles_tbl$profile_id == profile_id][1]
    } else {
      profile_id
    },
    thresholds = thresholds,
    criteria_meta = criteria_meta,
    thresholds_table = thresh_tbl,
    criteria_table = criteria_tbl,
    profiles_table = profiles_tbl
  )
}

get_criterion_strength <- function(criterion, criteria_meta = list()) {
  meta <- criteria_meta[[criterion]]
  if (is.null(meta) || is.null(meta$strength)) return("Unknown")
  meta$strength
}

get_criterion_description <- function(criterion, criteria_meta = list()) {
  meta <- criteria_meta[[criterion]]
  if (is.null(meta) || is.null(meta$description)) return("")
  meta$description
}
