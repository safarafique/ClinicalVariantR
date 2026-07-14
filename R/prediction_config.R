#' ACMGamp clinical prediction mode — stricter rules and report metadata.

PREDICTION_MODE <- TRUE

PREDICTION_SETTINGS <- list(
  require_lof_panel_for_pvs1 = TRUE,
  require_clinvar_no_conflict_for_pp5_bp6 = TRUE,
  min_clinvar_stars_pp5_bp6 = 0L,
  use_popmax_for_pm2_ba1 = TRUE,
  insilico_min_tools = 1L,
  insilico_revel_solo_pathogenic = 0.75,
  insilico_revel_solo_benign = 0.15,
  apply_single_tool_pp3 = TRUE,
  elevate_moderate_pathogenic_vus = TRUE,
  elevate_pvs1_with_supporting = TRUE,
  elevate_pvs1_lof_alone = TRUE,
  report_disclaimer = paste(
    "ACMGamp prediction output is decision support.",
    "Classifications require expert review before clinical use.",
    "Prediction sensitivity rules may upgrade strong VUS to Likely Pathogenic.",
    "See prediction_limitations column per variant."
  )
)

is_prediction_mode <- function() {
  isTRUE(PREDICTION_MODE)
}

prediction_threshold_overrides <- function() {
  if (!is_prediction_mode()) return(list())
  list(
    insilico_min_tools = PREDICTION_SETTINGS$insilico_min_tools,
    insilico_revel_solo = TRUE,
    revel_pp3 = PREDICTION_SETTINGS$insilico_revel_solo_pathogenic,
    revel_bp4 = PREDICTION_SETTINGS$insilico_revel_solo_benign
  )
}

reference_file_status <- function(path) {
  if (!file.exists(path)) {
    return(list(path = path, rows = 0L, placeholder = TRUE, label = "missing"))
  }
  lines <- readLines(path, warn = FALSE)
  data_lines <- lines[!grepl("^#", lines) & nzchar(trimws(lines))]
  n <- length(data_lines)
  if (n > 0L) n <- n - 1L
  placeholder <- n < 100L
  list(
    path = path,
    rows = max(0L, n),
    placeholder = placeholder,
    label = if (placeholder) "placeholder/small" else "loaded"
  )
}

check_reference_readiness <- function(paths = REFERENCE_PATHS) {
  status <- lapply(paths, reference_file_status)
  names(status) <- names(paths)
  ready <- all(!vapply(status, function(x) x$placeholder, logical(1L)))
  list(ready = ready, status = status)
}

build_prediction_limitations <- function(scores, row = NULL) {
  notes <- character()
  ref_chk <- check_reference_readiness()
  if (!ref_chk$ready) {
    notes <- c(notes, "Reference DBs are placeholders or incomplete; AF/ClinVar/REVEL may rely on VCF annotation only.")
  }
  if (isTRUE(scores$PVS1) && grepl("not in LoF mechanism panel", scores$PVS1_rationale %||% "", fixed = TRUE)) {
    notes <- c(notes, "PVS1 not applied in prediction mode without LoF gene panel support.")
  }
  if (!isTRUE(scores$PP5) && grepl("conflict", scores$PP5_rationale %||% "", ignore.case = TRUE)) {
    notes <- c(notes, "ClinVar conflict; PP5/BP6 withheld.")
  }
  if (isFALSE(scores$af_known %||% NA) ||
      is.na(scalar_num(scores$max_population_af %||% NA_real_))) {
    notes <- c(notes, "Population AF missing; PM2/BA1/BS rules may be incomplete.")
  }
  if (!nzchar(scalar_chr(scores$insilico_summary %||% "", default = ""))) {
    notes <- c(notes, "In silico scores sparse; PP3/BP4 may not trigger.")
  }
  if (is.null(row)) {
    return(if (length(notes) == 0L) "None identified." else paste(notes, collapse = " "))
  }
  if (length(notes) == 0L) "None identified." else paste(notes, collapse = " ")
}

#' Re-check VUS after automated scoring; upgrade when pathogenic evidence is strong and unconflicted.
apply_prediction_classification_refinement <- function(scores, evidence) {
  if (!is_prediction_mode()) return(scores)
  cls <- scalar_chr(scores$classification %||% "", default = "VUS")
  if (!identical(cls, "VUS")) return(scores)

  benign_n <- evidence$BS + evidence$BP + if (isTRUE(evidence$BA)) 1L else 0L
  path_n <- evidence$PVS + evidence$PS + evidence$PM + evidence$PP
  if (benign_n > 0L || isTRUE(scores$BA1)) return(scores)

  upgraded <- FALSE
  note <- ""

  if (evidence$PVS >= 1L && evidence$PM >= 1L) {
    upgraded <- TRUE
    cls <- "Likely Pathogenic"
    note <- "ACMG LP rule: PVS1 plus moderate pathogenic evidence."
  } else if (evidence$PS >= 1L && evidence$PM >= 1L) {
    upgraded <- TRUE
    cls <- "Likely Pathogenic"
    note <- "ACMG LP rule: strong plus moderate pathogenic evidence."
  } else if (evidence$PM >= 3L) {
    upgraded <- TRUE
    cls <- "Likely Pathogenic"
    note <- "ACMG LP rule: three or more moderate pathogenic criteria."
  } else if (evidence$PM >= 2L && evidence$PP >= 2L) {
    upgraded <- TRUE
    cls <- "Likely Pathogenic"
    note <- "ACMG LP rule: two moderate plus two supporting pathogenic criteria."
  } else if (isTRUE(PREDICTION_SETTINGS$elevate_pvs1_with_supporting) &&
             evidence$PVS >= 1L && evidence$PP >= 1L) {
    upgraded <- TRUE
    cls <- "Likely Pathogenic"
    note <- "Prediction sensitivity: PVS1 plus supporting pathogenic evidence."
  } else if (isTRUE(PREDICTION_SETTINGS$elevate_pvs1_lof_alone) && evidence$PVS >= 1L) {
    upgraded <- TRUE
    cls <- "Likely Pathogenic"
    note <- "Prediction sensitivity: PVS1 LoF in curated panel gene."
  } else if (isTRUE(PREDICTION_SETTINGS$elevate_moderate_pathogenic_vus) && path_n >= 2L) {
    upgraded <- TRUE
    cls <- "Likely Pathogenic"
    note <- "Prediction sensitivity: two or more pathogenic criteria without benign conflict."
  }

  if (!upgraded) return(scores)

  scores$classification <- cls
  lim <- scalar_chr(scores$prediction_limitations %||% "", default = "")
  scores$prediction_limitations <- paste(trimws(c(lim, note)), collapse = " ")
  scores
}

compute_evidence_strength <- function(scores) {
  evidence <- criteria_to_evidence(scores)
  path_n <- evidence$PVS + evidence$PS + evidence$PM + evidence$PP
  benign_n <- evidence$BS + evidence$BP + if (isTRUE(evidence$BA)) 1L else 0L
  if (path_n >= 3L || (evidence$PVS >= 1L && path_n >= 2L)) {
    return(list(strength = "Strong pathogenic", path_count = path_n, benign_count = benign_n))
  }
  if (path_n >= 1L && benign_n >= 1L) {
    return(list(strength = "Conflicting", path_count = path_n, benign_count = benign_n))
  }
  if (path_n >= 1L) {
    return(list(strength = "Moderate pathogenic", path_count = path_n, benign_count = benign_n))
  }
  if (isTRUE(evidence$BA) || benign_n >= 2L) {
    return(list(strength = "Strong benign", path_count = path_n, benign_count = benign_n))
  }
  if (benign_n >= 1L) {
    return(list(strength = "Moderate benign", path_count = path_n, benign_count = benign_n))
  }
  list(strength = "Insufficient", path_count = path_n, benign_count = benign_n)
}

MANUAL_EVIDENCE_KEYS <- c(
  "PS3_functional", "PP4_phenotype", "PS4_case_control", "PS2_de_novo",
  "PM6_de_novo", "PP1_segregation", "PP2_missense_mechanism"
)

empty_manual_inputs <- function() {
  stats::setNames(rep(list(FALSE), length(MANUAL_EVIDENCE_KEYS)), MANUAL_EVIDENCE_KEYS)
}

normalize_manual_inputs <- function(x) {
  base <- empty_manual_inputs()
  if (is.null(x) || length(x) == 0L) return(base)
  for (k in MANUAL_EVIDENCE_KEYS) {
    if (!is.null(x[[k]])) base[[k]] <- isTRUE(x[[k]])
  }
  base
}

manual_inputs_for_variant <- function(manual_by_variant, variant_id) {
  if (is.null(manual_by_variant) || length(manual_by_variant) == 0L) {
    return(empty_manual_inputs())
  }
  vid <- scalar_chr(variant_id, default = "")
  if (!nzchar(vid) || is.null(manual_by_variant[[vid]])) {
    return(empty_manual_inputs())
  }
  normalize_manual_inputs(manual_by_variant[[vid]])
}
