#' Reproducibility metadata for ClinicalVariantR runs.

file_checksum_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  tryCatch(
    as.character(tools::md5sum(path)),
    error = function(e) NA_character_
  )
}

capture_package_versions <- function() {
  pkgs <- c("shiny", "bslib", "DT", "data.table", "readr")
  installed <- rownames(installed.packages())
  vers <- vapply(pkgs, function(p) {
    if (p %in% installed) as.character(packageVersion(p)) else "not_installed"
  }, FUN.VALUE = character(1))
  paste(names(vers), vers, sep = "=", collapse = "; ")
}

build_run_metadata <- function(
    vcf_path = NA_character_,
    profile_id = DEFAULT_PROFILE_ID,
    rule_config = NULL,
    session_id = NA_character_,
    mode = "rapid") {

  if (is.null(rule_config)) rule_config <- load_rule_config(profile_id)

  list(
    app_version = APP_VERSION,
    engine = ACMG_PRO_ENGINE,
    acmg_guideline_version = ACMG_GUIDELINE_VERSION,
    disease_profile = rule_config$profile_id,
    disease_profile_name = rule_config$profile_name,
    thresholds_snapshot = paste(
      names(rule_config$thresholds),
      unlist(rule_config$thresholds),
      sep = "=",
      collapse = "; "
    ),
    input_vcf = basename(vcf_path),
    input_vcf_checksum = file_checksum_md5(vcf_path),
    r_packages = capture_package_versions(),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    session_id = session_id,
    pipeline_mode = mode,
    run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
    reference_paths = paste(names(REFERENCE_PATHS), unlist(REFERENCE_PATHS), sep = "=", collapse = "; "),
    prediction_mode = is_prediction_mode(),
    reference_ready = check_reference_readiness()$ready
  )
}

write_run_metadata_json <- function(metadata, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    txt <- jsonlite::toJSON(metadata, pretty = TRUE, auto_unbox = TRUE)
  } else {
    lines <- vapply(names(metadata), function(nm) {
      sprintf('  "%s": "%s"', nm, gsub('"', '\\"', as.character(metadata[[nm]])))
    }, FUN.VALUE = character(1))
    txt <- paste0("{\n", paste(lines, collapse = ",\n"), "\n}")
  }
  writeLines(txt, output_path, useBytes = TRUE)
  invisible(output_path)
}

metadata_to_single_row <- function(metadata) {
  data.frame(
    app_version = metadata$app_version,
    engine = metadata$engine,
    acmg_guideline_version = metadata$acmg_guideline_version,
    disease_profile = metadata$disease_profile,
    thresholds_snapshot = metadata$thresholds_snapshot,
    input_vcf = metadata$input_vcf,
    input_vcf_checksum = metadata$input_vcf_checksum,
    r_version = metadata$r_version,
    run_timestamp = metadata$run_timestamp,
    stringsAsFactors = FALSE
  )
}
