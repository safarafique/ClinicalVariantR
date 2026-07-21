#' App requirements for VCF and companion files.
vcf_app_requirements <- function(mode = c("full", "rapid")) {
  mode <- match.arg(mode)
  list(
    required_columns = c("CHROM", "POS", "REF", "ALT"),
    standard_columns = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"),
  info_fields = data.frame(
      field = c("CSQ", "ANN", "GENE", "AF", "REVEL"),
      tier = c("required", "required", "recommended", "recommended", "recommended"),
      used_for = c(
        "VEP consequence (PVS1, PM4, BP1, BP7, PP3)",
        "SnpEff consequence (PVS1, PM4, BP1, BP7, PP3)",
        "Gene symbol annotation",
        "Population allele frequency (BA1, BS1, PM2); gnomAD join if missing",
        "Computational pathogenicity (PP3, BP4); reference join if missing"
      ),
      stringsAsFactors = FALSE
    ),
    mode = mode,
    notes = if (mode == "full") {
      c(
        "Group A also requires Clinical Logs CSV and Pedigree CSV before analysis.",
        "Manual ACMG criteria (PS2-PS4, PP1-PP2, PP4, PM6) need curator input regardless of VCF content."
      )
    } else {
      c(
        "Group B uses automated rules only; clinical phenotype is not evaluated.",
        "PS4 is applied automatically when a variant matches the case-control reference DB.",
        "ClinVar and gnomAD are joined from reference files when absent in VCF."
      )
    }
  )
}

CLINICAL_REQUIRED_COLUMNS <- c("sample_id", "phenotype", "cml_phase", "tki_response")
PEDIGREE_REQUIRED_COLUMNS <- c("sample_id", "relation", "affected_status")

open_vcf_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "r") else file(path, "r")
}

#' Read VCF meta-header: column names, declared INFO/FORMAT fields, sample count.
#' When max_scan_variants is set, stops after that many variant rows (fast upload path).
read_vcf_header_meta <- function(vcf_path, max_data_rows = 100L, max_scan_variants = NULL,
                                 max_header_lines = 100000L) {
  if (!file.exists(vcf_path)) stop("VCF file not found.")

  if (grepl("\\.bcf$", vcf_path, ignore.case = TRUE)) {
    stop("BCF format detected. Please convert to VCF (.vcf or .vcf.gz) before upload.")
  }

  con <- open_vcf_connection(vcf_path)
  on.exit(close(con), add = TRUE)

  info_meta <- character()
  format_meta <- character()
  header_cols <- NULL
  sample_cols <- character()
  data_rows <- character()
  variant_count <- 0L
  count_truncated <- FALSE
  line_count <- 0L

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0) break
    line_count <- line_count + 1L

    if (is.null(header_cols) && line_count > max_header_lines) {
      stop(sprintf(
        "Invalid VCF: no #CHROM header found within the first %s lines.",
        format(max_header_lines, big.mark = ",")
      ))
    }

    if (grepl("^##INFO=<ID=([^,]+)", line, perl = TRUE)) {
      id <- sub("^##INFO=<ID=([^,]+).*$", "\\1", line)
      info_meta <- c(info_meta, id)
    } else if (grepl("^##FORMAT=<ID=([^,]+)", line, perl = TRUE)) {
      id <- sub("^##FORMAT=<ID=([^,]+).*$", "\\1", line)
      format_meta <- c(format_meta, id)
    } else if (grepl("^#CHROM\t", line)) {
      header_cols <- strsplit(sub("^#", "", line), "\t")[[1]]
      if (length(header_cols) > 8) {
        sample_cols <- header_cols[9:length(header_cols)]
      }
    } else if (!grepl("^#", line)) {
      variant_count <- variant_count + 1L
      if (length(data_rows) < max_data_rows) {
        data_rows[length(data_rows) + 1L] <- line
      }
      if (!is.null(max_scan_variants) && is.finite(max_scan_variants) &&
          variant_count >= max_scan_variants &&
          length(data_rows) >= min(max_data_rows, variant_count)) {
        count_truncated <- TRUE
        break
      }
    }
  }

  list(
    columns = header_cols,
    info_declared = unique(info_meta),
    format_declared = unique(format_meta),
    sample_columns = sample_cols,
    is_multi_sample = length(sample_cols) > 1,
    variant_count = variant_count,
    count_truncated = count_truncated,
    sample_rows = data_rows
  )
}

format_variant_count_label <- function(count, truncated = FALSE) {
  count <- as.integer(count)
  if (isTRUE(truncated)) {
    return(sprintf("At least %s variant row(s) (fast scan)", format(count, big.mark = ",")))
  }
  sprintf("%s variant row(s)", format(count, big.mark = ","))
}

#' Scan sample variant rows for populated INFO keys.
scan_info_fields_in_rows <- function(sample_rows, header_cols, fields) {
  if (length(sample_rows) == 0 || is.null(header_cols)) {
    return(setNames(rep(0L, length(fields)), fields))
  }

  info_idx <- match("INFO", header_cols)
  if (is.na(info_idx)) {
    return(setNames(rep(0L, length(fields)), fields))
  }

  mat <- strsplit(sample_rows, "\t")
  info_vals <- vapply(mat, function(x) x[info_idx], FUN.VALUE = character(1))
  counts <- integer(length(fields))
  names(counts) <- fields

  for (fld in fields) {
    counts[fld] <- sum(vapply(info_vals, function(s) {
      if (is.na(s) || s == ".") return(FALSE)
      pattern <- paste0("(^|;)", fld, "=")
      grepl(pattern, s, perl = TRUE)
    }, FUN.VALUE = logical(1)))
  }
  counts
}

#' Extract VEP CSQ subfield names from ##INFO CSQ header line.
extract_csq_format_fields <- function(vcf_path) {
  if (!file.exists(vcf_path)) return(character())
  con <- open_vcf_connection(vcf_path)
  on.exit(close(con), add = TRUE)
  fields <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (grepl("^##INFO=<ID=CSQ", line)) {
      m <- regexpr("Format: ([^\"]+)", line, perl = TRUE)
      if (m[1L] > 0L) {
        fmt <- sub('.*Format: ([^"]+).*', "\\1", line)
        fields <- strsplit(fmt, "|", fixed = TRUE)[[1L]]
      }
      break
    }
    if (grepl("^#CHROM\t", line)) break
  }
  fields
}

detect_csq_insilico_in_rows <- function(sample_rows, header_cols) {
  if (length(sample_rows) == 0) return(character())
  info_idx <- match("INFO", header_cols)
  if (is.na(info_idx)) return(character())
  found <- character()
  mat <- strsplit(sample_rows, "\t")
  info_vals <- vapply(mat, function(x) x[info_idx], FUN.VALUE = character(1))
  csq_vals <- vapply(info_vals, function(s) {
    if (is.na(s) || !grepl("CSQ=", s, fixed = TRUE)) return("")
    strsplit(sub(".*CSQ=", "", s), ";")[[1L]][1L]
  }, FUN.VALUE = character(1))
  csq_vals <- csq_vals[nzchar(csq_vals)]
  if (length(csq_vals) == 0) return(character())
  blob <- paste(csq_vals, collapse = " ")
  if (grepl("deleterious|tolerated", blob, ignore.case = TRUE)) found <- c(found, "SIFT (CSQ)")
  if (grepl("damaging|benign", blob, ignore.case = TRUE)) found <- c(found, "PolyPhen (CSQ)")
  if (grepl("REVEL", blob)) found <- c(found, "REVEL (CSQ)")
  if (grepl("CADD", blob)) found <- c(found, "CADD (CSQ)")
  unique(found)
}

#' Detect SnpEff-style ANN in INFO (Consequence subfield).
detect_ann_consequence <- function(sample_rows, header_cols) {
  if (length(sample_rows) == 0) return(FALSE)
  info_idx <- match("INFO", header_cols)
  if (is.na(info_idx)) return(FALSE)

  mat <- strsplit(sample_rows, "\t")
  info_vals <- vapply(mat, function(x) x[info_idx], FUN.VALUE = character(1))

  any(vapply(info_vals, function(s) {
    if (is.na(s) || s == ".") return(FALSE)
    if (grepl("ANN=", s, fixed = TRUE)) {
      ann_val <- sub(".*ANN=", "", s)
      ann_val <- strsplit(ann_val, ";")[[1]][1]
      return(grepl("missense|synonymous|stop_gained|frameshift|splice", ann_val, ignore.case = TRUE))
    }
    FALSE
  }, FUN.VALUE = logical(1)))
}

#' Detect VEP-style CSQ in INFO (Consequence subfield).
detect_csq_consequence <- function(sample_rows, header_cols) {
  if (length(sample_rows) == 0) return(FALSE)
  info_idx <- match("INFO", header_cols)
  if (is.na(info_idx)) return(FALSE)

  mat <- strsplit(sample_rows, "\t")
  info_vals <- vapply(mat, function(x) x[info_idx], FUN.VALUE = character(1))

  any(vapply(info_vals, function(s) {
    if (is.na(s) || s == ".") return(FALSE)
    if (grepl("CSQ=", s, fixed = TRUE)) {
      csq_val <- sub(".*CSQ=", "", s)
      csq_val <- strsplit(csq_val, ";")[[1]][1]
      return(grepl("missense|synonymous|stop_gained|frameshift|splice", csq_val, ignore.case = TRUE))
    }
    grepl("consequence=", s, ignore.case = TRUE)
  }, FUN.VALUE = logical(1)))
}

#' Validate uploaded VCF against app requirements.
new_vcf_check_recorder <- function() {
  checks_env <- new.env(parent = emptyenv())
  checks_env$items <- list()

  add <- function(category, requirement, status, detail) {
    items <- checks_env$items
    items[[length(items) + 1L]] <- list(
      category = category, requirement = requirement, status = status, detail = detail
    )
    checks_env$items <- items
    invisible(NULL)
  }

  to_df <- function() {
    checks <- checks_env$items
    if (length(checks) == 0) {
      return(data.frame(
        category = character(), requirement = character(),
        status = character(), detail = character(), stringsAsFactors = FALSE
      ))
    }
    df <- do.call(rbind, lapply(checks, as.data.frame, stringsAsFactors = FALSE))
    rownames(df) <- NULL
    df
  }

  list(add = add, to_df = to_df, items = function() checks_env$items)
}

unreadable_vcf_validation_result <- function(mode, checks_df) {
  list(
    valid = FALSE,
    can_analyze = FALSE,
    mode = mode,
    checks = checks_df,
    columns = character(),
    info_declared = character(),
    variant_count = 0L,
    summary = "VCF validation failed: file could not be read."
  )
}

detect_vcf_validation_annotation <- function(meta, vcf_path) {
  has_csq_consequence <- detect_csq_consequence(meta$sample_rows, meta$columns)
  has_ann_consequence <- detect_ann_consequence(meta$sample_rows, meta$columns)
  has_annovar <- any(c("Func.refGene", "Gene.refGene", "ExonicFunc.refGene") %in% meta$info_declared) ||
    any(grepl("Func\\.refGene=", meta$sample_rows, ignore.case = TRUE))

  csq_format_fields <- extract_csq_format_fields(vcf_path)
  csq_insilico <- detect_csq_insilico_in_rows(meta$sample_rows, meta$columns)
  csq_plugins <- character()
  if (length(csq_format_fields) > 0L) {
    if ("REVEL" %in% csq_format_fields) csq_plugins <- c(csq_plugins, "REVEL (CSQ)")
    if ("CADD_PHRED" %in% csq_format_fields || "CADD" %in% csq_format_fields) {
      csq_plugins <- c(csq_plugins, "CADD (CSQ)")
    }
    if (any(grepl("SpliceAI", csq_format_fields))) csq_plugins <- c(csq_plugins, "SpliceAI (CSQ)")
    if (any(grepl("am_pathogenicity|AlphaMissense", csq_format_fields))) {
      csq_plugins <- c(csq_plugins, "AlphaMissense (CSQ)")
    }
    if ("SIFT" %in% csq_format_fields) csq_plugins <- c(csq_plugins, "SIFT (CSQ)")
    if ("PolyPhen" %in% csq_format_fields) csq_plugins <- c(csq_plugins, "PolyPhen (CSQ)")
  }

  ann_info <- detect_vcf_annotation(meta$info_declared, character())
  if (has_csq_consequence) ann_info <- list(source = "CSQ", build_hint = "GRCh38", label = "VEP CSQ")
  if (has_ann_consequence && !has_csq_consequence) {
    ann_info <- list(source = "ANN", build_hint = "GRCh37/hg19", label = "SnpEff ANN")
  }
  if (has_annovar && !has_csq_consequence && !has_ann_consequence) {
    ann_info <- list(source = "ANNOVAR", build_hint = "unknown", label = "ANNOVAR")
  }

  list(
    has_csq_consequence = has_csq_consequence,
    has_ann_consequence = has_ann_consequence,
    has_annovar = has_annovar,
    has_consequence_annotation = has_csq_consequence || has_ann_consequence || has_annovar,
    csq_plugins = unique(c(csq_plugins, csq_insilico)),
    ann_info = ann_info
  )
}

validate_vcf <- function(vcf_path, mode = c("full", "rapid"), sample_rows = 100L) {
  mode <- match.arg(mode)
  reqs <- vcf_app_requirements(mode)
  checks <- new_vcf_check_recorder()
  add_check <- checks$add

  meta <- tryCatch(
    read_vcf_header_meta(
      vcf_path,
      max_data_rows = sample_rows,
      max_scan_variants = max(sample_rows, 1L) + 1L
    ),
    error = function(e) {
      add_check("File", "Readable VCF", "FAIL", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(meta)) {
    return(unreadable_vcf_validation_result(mode, checks$to_df()))
  }

  count_label <- format_variant_count_label(meta$variant_count, meta$count_truncated)
  add_check("File", "VCF format", "PASS", count_label)

  if (meta$variant_count == 0) {
    add_check("Content", "Variant records", "FAIL", "VCF contains no variant data rows.")
  } else {
    add_check("Content", "Variant records", "PASS", count_label)
  }

  missing_cols <- setdiff(reqs$required_columns, meta$columns)
  if (length(missing_cols) > 0) {
    add_check("Columns", "Required VCF columns", "FAIL", paste("Missing:", paste(missing_cols, collapse = ", ")))
  } else {
    add_check("Columns", "Required VCF columns", "PASS", paste(reqs$required_columns, collapse = ", "))
  }

  if (!is.null(meta$columns)) {
    add_check(
      "Columns", "VCF header columns", "PASS",
      paste(meta$columns, collapse = " | ")
    )
    extra_std <- setdiff(reqs$standard_columns, meta$columns)
    if (length(extra_std) > 0) {
      add_check(
        "Columns", "Standard columns", "WARN",
        paste("Missing standard columns:", paste(extra_std, collapse = ", "))
      )
    }
  }

  if (is.null(meta$columns) || !"INFO" %in% meta$columns) {
    add_check("Columns", "INFO column", "FAIL", "INFO column required for annotation fields.")
  } else {
    add_check("Columns", "INFO column", "PASS", "INFO column present.")
  }

  info_fields <- reqs$info_fields$field
  populated <- scan_info_fields_in_rows(meta$sample_rows, meta$columns, info_fields)
  annotation <- detect_vcf_validation_annotation(meta, vcf_path)
  has_csq_consequence <- annotation$has_csq_consequence
  has_ann_consequence <- annotation$has_ann_consequence
  has_consequence_annotation <- annotation$has_consequence_annotation
  csq_plugins <- annotation$csq_plugins
  ann_info <- annotation$ann_info

  for (i in seq_len(nrow(reqs$info_fields))) {
    fld <- reqs$info_fields$field[i]
    tier <- reqs$info_fields$tier[i]
    used <- reqs$info_fields$used_for[i]
    declared <- fld %in% meta$info_declared
    pop_n <- populated[[fld]]
    sample_n <- length(meta$sample_rows)

    if (fld == "CSQ") {
      if (has_csq_consequence) {
        add_check("INFO", "VEP CSQ / consequence", "PASS", sprintf("Consequence detected in sample rows. %s", used))
      } else if (has_ann_consequence) {
        add_check("INFO", "VEP CSQ / consequence", "PASS",
          sprintf("SnpEff ANN used instead (%s). %s", ann_info$label, used))
      } else if (declared) {
        add_check("INFO", "VEP CSQ / consequence", "WARN", "CSQ declared in header but not populated in sampled rows.")
      } else {
        add_check("INFO", "VEP CSQ / consequence", "WARN", paste("No CSQ annotation.", used))
      }
      next
    }

    if (fld == "ANN") {
      if (has_ann_consequence) {
        add_check("INFO", "SnpEff ANN / consequence", "PASS", sprintf("Consequence detected in sample rows. %s", used))
      } else if (has_csq_consequence) {
        add_check("INFO", "SnpEff ANN / consequence", "PASS", "VEP CSQ present; ANN not required.")
      } else if (declared) {
        add_check("INFO", "SnpEff ANN / consequence", "WARN", "ANN declared in header but not populated in sampled rows.")
      } else {
        add_check("INFO", "SnpEff ANN / consequence", "WARN", paste("No ANN annotation.", used))
      }
      next
    }

    if (fld == "GENE") {
      if (pop_n > 0) {
        add_check("INFO", "GENE", "PASS", sprintf("Present in %d/%d sampled rows. %s", pop_n, sample_n, used))
      } else if (has_ann_consequence || has_csq_consequence) {
        add_check("INFO", "GENE", "PASS", sprintf("Gene symbol read from %s annotation field.", ann_info$label))
      } else {
        add_check("INFO", "GENE", "WARN", paste("No gene annotation found.", used))
      }
      next
    }

    if (fld == "AF") {
      af_equiv <- c("AF", "1000Gp3_AF", "ESP6500_MAF", "SAMPLES_AF", "gnomAD_AF", "gnomAD_genome_ALL")
      present_af <- intersect(af_equiv, meta$info_declared)
      if (length(present_af) > 0) {
        add_check("INFO", "AF / population frequency", "PASS",
          sprintf("Population AF available via: %s. %s", paste(present_af, collapse = ", "), used))
      } else {
        add_check("INFO", "AF / population frequency", "WARN",
          paste("No population AF field found. Reference database join will be attempted.", used))
      }
      next
    }

    if (fld == "REVEL") {
      pred_equiv <- c("REVEL", "CADD_PHRED", "CADD", "SpliceAI", "MAX_SPLICEAI",
                      "am_pathogenicity", "AlphaMissense", "ESP6500_PH")
      present_pred <- intersect(pred_equiv, meta$info_declared)
      present_pred <- unique(c(present_pred, csq_plugins))
      if (length(present_pred) > 0) {
        add_check("INFO", "REVEL / prediction scores", "PASS",
          sprintf("Computational evidence available via: %s. %s", paste(present_pred, collapse = ", "), used))
      } else {
        add_check("INFO", "REVEL / prediction scores", "WARN",
          paste(
            "No REVEL/CADD/SpliceAI/AlphaMissense/SIFT/PolyPhen found.",
            "Re-annotate with VEP plugins for stronger PP3/BP4, or analysis will use population/consequence rules only."
          ))
      }
      next
    }

    ref_join_note <- if (fld %in% c("AF", "REVEL")) " Reference database join will be attempted." else ""

    if (pop_n > 0) {
      add_check("INFO", fld, "PASS", sprintf("Present in %d/%d sampled rows. %s", pop_n, sample_n, used))
    } else if (declared) {
      status <- if (tier == "required") "WARN" else "WARN"
      add_check("INFO", fld, status, sprintf("Declared in header but empty in sampled rows.%s", ref_join_note))
    } else {
      status <- if (tier == "required") "FAIL" else "WARN"
      add_check("INFO", fld, status, sprintf("Not found in VCF header.%s", ref_join_note))
    }
  }

  if (!has_consequence_annotation) {
    add_check(
      "INFO", "Consequence annotation", "FAIL",
      "No VEP CSQ or SnpEff ANN consequence found. Upload a VEP- or SnpEff-annotated VCF."
    )
  } else {
    add_check(
      "INFO", "ClinicalVariantR engine", "PASS",
      sprintf("%s detected - InterVar-style automated criteria enabled.", ann_info$label)
    )
  }

  if (meta$is_multi_sample) {
    add_check(
      "Samples", "Multi-sample VCF", "WARN",
      sprintf("%d sample columns detected. App uses variant-level INFO; trio phasing (PM3, PS2) is not automated.", length(meta$sample_columns))
    )
  } else if (length(meta$sample_columns) == 1) {
    add_check("Samples", "Single-sample VCF", "PASS", paste("Sample:", meta$sample_columns[1]))
  } else {
    add_check("Samples", "Genotype columns", "WARN", "No sample columns; variant-only VCF.")
  }

  check_items <- checks$items()
  fail_n <- sum(vapply(check_items, function(x) x$status == "FAIL", logical(1)))
  warn_n <- sum(vapply(check_items, function(x) x$status == "WARN", logical(1)))
  checks_df <- checks$to_df()

  can_analyze <- fail_n == 0 && meta$variant_count > 0 && length(missing_cols) == 0

  missing_vcf_columns <- missing_cols
  missing_info_fields <- character()
  if (!has_consequence_annotation) {
    missing_info_fields <- c(missing_info_fields, "CSQ or ANN (consequence)")
  }
  for (fld in info_fields) {
    if (fld %in% c("CSQ", "ANN")) next
    tier <- reqs$info_fields$tier[reqs$info_fields$field == fld]
    pop_n <- populated[[fld]]
    declared <- fld %in% meta$info_declared
    if (pop_n == 0 && !declared && tier == "required") {
      missing_info_fields <- c(missing_info_fields, fld)
    }
  }
  if (is.null(meta$columns) || !"INFO" %in% meta$columns) {
    missing_vcf_columns <- unique(c(missing_vcf_columns, "INFO"))
  }

  missing_items <- c(
    if (length(missing_vcf_columns) > 0) paste0("VCF column: ", missing_vcf_columns),
    if (length(missing_info_fields) > 0) paste0("INFO field: ", missing_info_fields)
  )
  if (meta$variant_count == 0) {
    missing_items <- c(missing_items, "Variant data rows (VCF is empty)")
  }

  summary <- if (!can_analyze) {
    if (length(missing_items) > 0) {
      paste0("Not ready. Missing: ", paste(missing_items, collapse = "; "))
    } else {
      sprintf("Not ready for analysis: %d failure(s), %d warning(s).", fail_n, warn_n)
    }
  } else if (warn_n > 0) {
    sprintf("Ready with caveats: %d warning(s). You may click Run Analysis.", warn_n)
  } else {
    "Ready for analysis. Click Run Analysis."
  }

  list(
    valid = fail_n == 0,
    can_analyze = can_analyze,
    readiness = if (can_analyze) "READY" else "NOT_READY",
    mode = mode,
    checks = checks_df,
    columns = meta$columns %||% character(),
    missing_vcf_columns = missing_vcf_columns,
    missing_info_fields = missing_info_fields,
    missing_items = missing_items,
    info_declared = meta$info_declared,
    format_declared = meta$format_declared,
    sample_columns = meta$sample_columns,
    variant_count = meta$variant_count,
    annotation = ann_info,
    summary = summary
  )
}

validate_clinical_csv <- function(path) {
  checks <- data.frame(
    category = "Clinical CSV", requirement = character(),
    status = character(), detail = character(), stringsAsFactors = FALSE
  )
  if (!file.exists(path)) {
    return(list(
      valid = FALSE, can_analyze = FALSE, readiness = "NOT_READY",
      missing_columns = CLINICAL_REQUIRED_COLUMNS,
      missing_items = "Clinical Logs CSV file not uploaded",
      checks = checks, summary = "Clinical Logs CSV not uploaded."
    ))
  }

  cols <- names(readr::read_csv(path, n_max = 0, show_col_types = FALSE))
  missing <- setdiff(CLINICAL_REQUIRED_COLUMNS, cols)

  if (length(missing) == 0) {
    checks <- rbind(checks, data.frame(
      category = "Clinical CSV", requirement = "Required columns",
      status = "PASS", detail = paste(CLINICAL_REQUIRED_COLUMNS, collapse = ", ")
    ))
    list(
      valid = TRUE, can_analyze = TRUE, readiness = "READY",
      missing_columns = character(), missing_items = character(),
      checks = checks, columns = cols,
      summary = "Clinical logs CSV meets Group A requirements."
    )
  } else {
    checks <- rbind(checks, data.frame(
      category = "Clinical CSV", requirement = "Required columns",
      status = "FAIL", detail = paste("Missing:", paste(missing, collapse = ", "))
    ))
    list(
      valid = FALSE, can_analyze = FALSE, readiness = "NOT_READY",
      missing_columns = missing, missing_items = paste0("Clinical CSV column: ", missing),
      checks = checks, columns = cols,
      summary = paste("Clinical CSV missing columns:", paste(missing, collapse = ", "))
    )
  }
}

validate_pedigree_csv <- function(path) {
  checks <- data.frame(
    category = "Pedigree CSV", requirement = character(),
    status = character(), detail = character(), stringsAsFactors = FALSE
  )
  if (!file.exists(path)) {
    return(list(
      valid = FALSE, can_analyze = FALSE, readiness = "NOT_READY",
      missing_columns = PEDIGREE_REQUIRED_COLUMNS,
      missing_items = "Pedigree CSV file not uploaded",
      checks = checks, summary = "Pedigree CSV not uploaded."
    ))
  }

  cols <- names(readr::read_csv(path, n_max = 0, show_col_types = FALSE))
  missing <- setdiff(PEDIGREE_REQUIRED_COLUMNS, cols)

  if (length(missing) == 0) {
    checks <- rbind(checks, data.frame(
      category = "Pedigree CSV", requirement = "Required columns",
      status = "PASS", detail = paste(PEDIGREE_REQUIRED_COLUMNS, collapse = ", ")
    ))
    list(
      valid = TRUE, can_analyze = TRUE, readiness = "READY",
      missing_columns = character(), missing_items = character(),
      checks = checks, columns = cols,
      summary = "Pedigree CSV meets Group A requirements."
    )
  } else {
    checks <- rbind(checks, data.frame(
      category = "Pedigree CSV", requirement = "Required columns",
      status = "FAIL", detail = paste("Missing:", paste(missing, collapse = ", "))
    ))
    list(
      valid = FALSE, can_analyze = FALSE, readiness = "NOT_READY",
      missing_columns = missing, missing_items = paste0("Pedigree CSV column: ", missing),
      checks = checks, columns = cols,
      summary = paste("Pedigree CSV missing columns:", paste(missing, collapse = ", "))
    )
  }
}

#' Combined Group A readiness across VCF + clinical + pedigree.
validate_group_a_inputs <- function(vcf_path, clinical_path, pedigree_path, vcf_val = NULL) {
  if (is.null(vcf_val)) {
    vcf_val <- validate_vcf(vcf_path, mode = "full")
  }
  clin_val <- validate_clinical_csv(clinical_path)
  ped_val <- validate_pedigree_csv(pedigree_path)

  checks <- rbind(vcf_val$checks, clin_val$checks, ped_val$checks)
  fail_n <- sum(checks$status == "FAIL")
  warn_n <- sum(checks$status == "WARN")

  can_analyze <- vcf_val$can_analyze && clin_val$can_analyze && ped_val$can_analyze

  missing_items <- c(
    vcf_val$missing_items %||% character(),
    clin_val$missing_items %||% character(),
    ped_val$missing_items %||% character()
  )
  missing_items <- unique(unlist(missing_items))

  summary <- if (!can_analyze) {
    if (length(missing_items) > 0) {
      paste0("Not ready. Missing: ", paste(missing_items, collapse = "; "))
    } else {
      sprintf("Group A not ready: %d failure(s), %d warning(s).", fail_n, warn_n)
    }
  } else if (warn_n > 0) {
    sprintf("Ready with %d warning(s). Click Run Analysis.", warn_n)
  } else {
    "All inputs validated. Ready for analysis - click Run Analysis."
  }

  list(
    valid = fail_n == 0,
    can_analyze = can_analyze,
    readiness = if (can_analyze) "READY" else "NOT_READY",
    vcf = vcf_val,
    clinical = clin_val,
    pedigree = ped_val,
    checks = checks,
    missing_items = missing_items,
    summary = summary
  )
}
