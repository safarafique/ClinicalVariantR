#' Benchmark ClinicalVariantR against reference ACMG classifications.

normalize_acmg_class <- function(x) {
  vapply(x, function(val) {
    if (is.na(val) || !nzchar(val)) return(NA_character_)
    val <- gsub("_", " ", tolower(trimws(val)))
    val <- gsub("benign/likely benign", "likely benign", val, fixed = TRUE)
    val <- gsub("pathogenic/likely pathogenic", "likely pathogenic", val, fixed = TRUE)
    if (grepl("vus|uncertain|conflict", val)) return("VUS")
    tools::toTitleCase(val)
  }, FUN.VALUE = character(1))
}

collapse_to_tier <- function(x) {
  x <- normalize_acmg_class(x)
  out <- rep("other", length(x))
  out[x %in% c("Pathogenic", "Likely Pathogenic")] <- "pathogenic"
  out[x %in% c("Benign", "Likely Benign")] <- "benign"
  out[x %in% c("VUS")] <- "vus"
  out
}

parse_testing_vcf <- function(vcf_path, pass_only = TRUE) {
  con <- if (grepl("\\.gz$", vcf_path, ignore.case = TRUE)) gzfile(vcf_path, "rt") else file(vcf_path, "rt")
  on.exit(close(con), add = TRUE)

  rows <- list()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (grepl("^#", line)) next
    parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 8L) next
    if (isTRUE(pass_only)) {
      filt <- parts[[7L]]
      if (!(filt %in% c("PASS", "."))) next
    }
    alts <- split_vcf_alt_alleles(parts[[5L]])
    for (alt in alts) {
      row <- parse_variant_from_vcf_fields(
        chrom = parts[[1L]], pos = parts[[2L]], ref = parts[[4L]],
        alt = alt,
        qual = scalar_num(parts[[6L]]),
        filter = parts[[7L]],
        info = parts[[8L]]
      )
      rows[[length(rows) + 1L]] <- row
    }
  }
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

load_reference_acmg_tsv <- function(tsv_path) {
  df <- read.delim(tsv_path, stringsAsFactors = FALSE, check.names = FALSE)
  df$variant_key <- variant_key_from_parts(df$chr, df$pos, df$ref, df$alt)
  df$ref_class <- normalize_acmg_class(df$acmg_classification_base)
  if ("pathogenicity_classification_combined_base" %in% names(df)) {
    alt_class <- normalize_acmg_class(df$pathogenicity_classification_combined_base)
    df$ref_class <- ifelse(is.na(df$ref_class) | !nzchar(df$ref_class), alt_class, df$ref_class)
  }
  df$ref_tier <- collapse_to_tier(df$ref_class)
  df$ref_criteria <- df$acmg_criteria_base %||% ""
  df
}

benchmark_one_sample <- function(vcf_path, tsv_path, profile_id = DEFAULT_PROFILE_ID) {
  sample_id <- tools::file_path_sans_ext(basename(vcf_path))
  t_start <- Sys.time()

  variants <- parse_testing_vcf(vcf_path, pass_only = TRUE)
  parse_secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  if (nrow(variants) == 0L) {
    return(list(sample_id = sample_id, error = "No PASS variants parsed", metrics = NULL))
  }

  t_score <- Sys.time()
  scored <- score_variants_table(variants, profile_id = profile_id)
  score_secs <- as.numeric(difftime(Sys.time(), t_score, units = "secs"))

  scored$variant_key <- variant_key_from_parts(scored$chrom, scored$pos, scored$ref, scored$alt)
  scored$pred_class <- normalize_acmg_class(scored$classification)
  scored$pred_tier <- collapse_to_tier(scored$pred_class)

  ref <- load_reference_acmg_tsv(tsv_path)

  # One reference row per variant key (prefer MTHFR missense if duplicates)
  ref_u <- ref[!duplicated(ref$variant_key), c("variant_key", "ref_class", "ref_tier", "ref_criteria", "gene_symbol_base"), drop = FALSE]
  names(ref_u)[names(ref_u) == "gene_symbol_base"] <- "ref_gene"

  cmp <- merge(
    scored[, c("variant_key", "gene", "consequence", "pred_class", "pred_tier",
               "classification", "criteria_met", "confidence_score")],
    ref_u,
    by = "variant_key",
    all = FALSE
  )

  total_secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  if (nrow(cmp) == 0L) {
    return(list(
      sample_id = sample_id,
      error = "No overlapping variants between VCF and reference TSV",
      n_vcf = nrow(variants),
      n_ref = nrow(ref_u),
      metrics = NULL,
      timing = list(parse_secs = parse_secs, score_secs = score_secs, total_secs = total_secs)
    ))
  }

  cmp$exact_match <- cmp$pred_class == cmp$ref_class
  cmp$tier_match <- cmp$pred_tier == cmp$ref_tier

  metrics <- list(
    sample_id = sample_id,
    n_vcf_variants = nrow(variants),
    n_reference_variants = nrow(ref_u),
    n_compared = nrow(cmp),
    exact_accuracy = mean(cmp$exact_match, na.rm = TRUE),
    tier_accuracy = mean(cmp$tier_match, na.rm = TRUE),
    parse_seconds = parse_secs,
    score_seconds = score_secs,
    total_seconds = total_secs,
    variants_per_second = nrow(variants) / max(total_secs, 0.001)
  )

  # Per-tier counts
  tab <- table(Predicted = cmp$pred_tier, Reference = cmp$ref_tier)
  metrics$confusion_tier <- tab

  # Pathogenic vs benign sensitivity/specificity (binary collapse)
  path_pred <- cmp$pred_tier == "pathogenic"
  path_ref <- cmp$ref_tier == "pathogenic"
  ben_ref <- cmp$ref_tier == "benign"
  if (any(path_ref)) {
    metrics$sensitivity_pathogenic <- mean(path_pred[path_ref], na.rm = TRUE)
  }
  if (any(ben_ref)) {
    metrics$specificity_benign <- mean(!path_pred[ben_ref], na.rm = TRUE)
  }

  list(
    sample_id = sample_id,
    error = NULL,
    metrics = metrics,
    comparison = cmp,
    timing = list(parse_secs = parse_secs, score_secs = score_secs, total_secs = total_secs)
  )
}

find_reference_tsv <- function(vcf_path, folder_path) {
  sample <- tools::file_path_sans_ext(basename(vcf_path))
  candidates <- c(
    file.path(folder_path, paste0(sample, ".acmg.tsv")),
    file.path(folder_path, paste0(sample, ".acmg.tsv.gz"))
  )
  id_match <- regmatches(sample, regexpr("^[0-9]+", sample))
  if (length(id_match) > 0 && nzchar(id_match[[1]])) {
    candidates <- c(
      candidates,
      file.path(folder_path, paste0(id_match[[1]], ".acmg.tsv")),
      file.path(folder_path, paste0(id_match[[1]], ".acmg.tsv.gz"))
    )
  }
  hits <- candidates[file.exists(candidates)]
  if (length(hits) > 0) hits[[1]] else NA_character_
}

benchmark_testing_folder <- function(
    folder_path,
    profile_id = DEFAULT_PROFILE_ID,
    output_dir = NULL) {

  vcf_files <- list.files(folder_path, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  vcf_files <- vcf_files[!grepl("\\.clean\\.vcf$", vcf_files, ignore.case = TRUE)]

  results <- lapply(vcf_files, function(vcf) {
    sample <- tools::file_path_sans_ext(basename(vcf))
    tsv <- find_reference_tsv(vcf, folder_path)
    if (is.na(tsv)) {
      return(list(sample_id = sample, error = paste("Missing reference TSV for", sample), metrics = NULL))
    }
    benchmark_one_sample(vcf, tsv, profile_id = profile_id)
  })

  summary_rows <- lapply(results, function(r) {
    if (is.null(r$metrics)) {
      return(data.frame(
        sample_id = r$sample_id,
        error = r$error %||% "unknown",
        n_compared = 0L,
        exact_accuracy = NA_real_,
        tier_accuracy = NA_real_,
        total_seconds = r$timing$total_secs %||% NA_real_,
        variants_per_second = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    m <- r$metrics
    data.frame(
      sample_id = m$sample_id,
      error = NA_character_,
      n_vcf_variants = m$n_vcf_variants,
      n_compared = m$n_compared,
      exact_accuracy = round(m$exact_accuracy, 4),
      tier_accuracy = round(m$tier_accuracy, 4),
      sensitivity_pathogenic = round(m$sensitivity_pathogenic %||% NA_real_, 4),
      specificity_benign = round(m$specificity_benign %||% NA_real_, 4),
      parse_seconds = round(m$parse_seconds, 2),
      score_seconds = round(m$score_seconds, 2),
      total_seconds = round(m$total_seconds, 2),
      variants_per_second = round(m$variants_per_second, 1),
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(summary_df, file.path(output_dir, "benchmark_summary.csv"), row.names = FALSE)
    for (r in results) {
      if (!is.null(r$comparison) && nrow(r$comparison) > 0) {
        write.csv(
          r$comparison,
          file.path(output_dir, paste0(r$sample_id, ".benchmark_comparison.csv")),
          row.names = FALSE
        )
      }
    }
  }

  list(summary = summary_df, results = results)
}
