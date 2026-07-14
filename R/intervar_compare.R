#' Compare ACMGamp classifications with InterVar or InterVar-style reference files.

ACMG_CLASS_LABELS <- c(
  "Pathogenic", "Likely Pathogenic", "VUS", "Likely Benign", "Benign"
)

normalize_acmg_class_label <- function(x) {
  vapply(x, function(val) {
    if (is.na(val) || !nzchar(val)) return(NA_character_)
    val <- gsub("_", " ", tolower(trimws(val)))
    val <- gsub("benign/likely benign", "likely benign", val, fixed = TRUE)
    val <- gsub("pathogenic/likely pathogenic", "likely pathogenic", val, fixed = TRUE)
    if (grepl("vus|uncertain|conflict", val)) return("VUS")
    tools::toTitleCase(val)
  }, FUN.VALUE = character(1))
}

collapse_acmg_tier <- function(x) {
  x <- normalize_acmg_class_label(x)
  out <- rep("other", length(x))
  out[x %in% c("Pathogenic", "Likely Pathogenic")] <- "pathogenic"
  out[x %in% c("Benign", "Likely Benign")] <- "benign"
  out[x %in% c("Vus")] <- "vus"
  out
}

pick_column <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) return(NA_character_)
  hit[[1]]
}

load_acmgamp_report_csv <- function(path) {
  if (!file.exists(path)) stop("ACMGamp report not found: ", path)
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  chrom_col <- pick_column(df, c("chrom", "chr", "CHROM", "Chr"))
  pos_col <- pick_column(df, c("pos", "POS", "Start", "start"))
  ref_col <- pick_column(df, c("ref", "REF", "Ref"))
  alt_col <- pick_column(df, c("alt", "ALT", "Alt"))
  class_col <- pick_column(df, c("classification", "Classification"))
  criteria_col <- pick_column(df, c("criteria_met", "criteria", "Criteria"))
  if (any(is.na(c(chrom_col, pos_col, ref_col, alt_col, class_col)))) {
    stop("ACMGamp CSV missing required columns (chrom/pos/ref/alt/classification).")
  }
  out <- data.frame(
    variant_key = variant_key_chr_pos_ref_alt(df[[chrom_col]], df[[pos_col]], df[[ref_col]], df[[alt_col]]),
    gene = if ("gene" %in% names(df)) df$gene else NA_character_,
    classification = normalize_acmg_class_label(df[[class_col]]),
    criteria = if (!is.na(criteria_col)) df[[criteria_col]] else "",
    source = "ACMGamp",
    stringsAsFactors = FALSE
  )
  out$tier <- collapse_acmg_tier(out$classification)
  out
}

load_intervar_reference_tsv <- function(path) {
  if (!file.exists(path)) stop("Reference TSV not found: ", path)
  df <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  chrom_col <- pick_column(df, c("chr", "chrom", "CHROM", "#Chr"))
  pos_col <- pick_column(df, c("pos", "POS", "Start", "start"))
  ref_col <- pick_column(df, c("ref", "REF", "Ref"))
  alt_col <- pick_column(df, c("alt", "ALT", "Alt"))
  class_col <- pick_column(df, c(
    "acmg_classification_base", "pathogenicity_classification_combined_base",
    "classification", "InterVar"
  ))
  criteria_col <- pick_column(df, c("acmg_criteria_base", "criteria", "acmg_criteria"))
  gene_col <- pick_column(df, c("gene_symbol_base", "gene", "Ref.Gene"))
  if (any(is.na(c(chrom_col, pos_col, ref_col, alt_col)))) {
    stop("Reference TSV missing chr/pos/ref/alt columns.")
  }
  out <- data.frame(
    variant_key = variant_key_chr_pos_ref_alt(df[[chrom_col]], df[[pos_col]], df[[ref_col]], df[[alt_col]]),
    gene = if (!is.na(gene_col)) df[[gene_col]] else NA_character_,
    classification = normalize_acmg_class_label(if (!is.na(class_col)) df[[class_col]] else NA_character_),
    criteria = if (!is.na(criteria_col)) df[[criteria_col]] else "",
    source = "Reference",
    stringsAsFactors = FALSE
  )
  out$tier <- collapse_acmg_tier(out$classification)
  out[!duplicated(out$variant_key), , drop = FALSE]
}

parse_intervar_evidence_field <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(list(classification = NA_character_, criteria = ""))
  }
  class_match <- regmatches(x, regexpr("InterVar:\\s*[^PVS1]+", x, perl = TRUE))
  cls <- if (length(class_match) == 0L) {
    NA_character_
  } else {
    gsub("^InterVar:\\s*", "", class_match[[1]])
  }
  cls <- gsub("\\s+PVS1.*$", "", cls)
  criteria <- sub("^.*?PVS1", "PVS1", x)
  if (!grepl("^PVS1", criteria)) criteria <- ""
  list(classification = cls, criteria = criteria)
}

load_intervar_output <- function(path) {
  if (!file.exists(path)) stop("InterVar output not found: ", path)
  first <- readLines(path, n = 1L, warn = FALSE)
  if (length(first) == 0L) stop("InterVar file is empty: ", path)
  sep <- if (grepl("\t", first)) "\t" else ","
  df <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, sep = sep)
  if (nrow(df) == 0L) stop("No rows in InterVar file: ", path)

  chrom_col <- pick_column(df, c("#Chr", "Chr", "chr", "CHROM", "chrom"))
  pos_col <- pick_column(df, c("Start", "start", "pos", "POS"))
  ref_col <- pick_column(df, c("Ref", "ref", "REF"))
  alt_col <- pick_column(df, c("Alt", "alt", "ALT"))
  gene_col <- pick_column(df, c("Ref.Gene", "Gene.refGene", "gene", "Gene"))
  evidence_col <- pick_column(df, c(
    "InterVar..InterVar.and.Evidence",
    "InterVar: InterVar and Evidence",
    "InterVar_and_Evidence",
    "InterVar"
  ))
  if (any(is.na(c(chrom_col, pos_col, ref_col, alt_col)))) {
    stop("InterVar file missing Chr/Start/Ref/Alt columns. Got: ", paste(names(df), collapse = ", "))
  }

  parsed <- if (!is.na(evidence_col)) {
    lapply(df[[evidence_col]], parse_intervar_evidence_field)
  } else {
    rep(list(list(classification = NA_character_, criteria = "")), nrow(df))
  }

  out <- data.frame(
    variant_key = variant_key_chr_pos_ref_alt(df[[chrom_col]], df[[pos_col]], df[[ref_col]], df[[alt_col]]),
    gene = if (!is.na(gene_col)) df[[gene_col]] else NA_character_,
    classification = normalize_acmg_class_label(vapply(parsed, `[[`, "", "classification")),
    criteria = vapply(parsed, `[[`, "", "criteria"),
    source = "InterVar",
    stringsAsFactors = FALSE
  )
  out$tier <- collapse_acmg_tier(out$classification)
  out[!duplicated(out$variant_key), , drop = FALSE]
}

classification_count_table <- function(df, label = "classification") {
  counts <- setNames(rep(0L, length(ACMG_CLASS_LABELS)), ACMG_CLASS_LABELS)
  if (!is.null(df) && nrow(df) > 0L) {
    tab <- table(df[[label]])
    for (nm in names(tab)) {
      norm <- normalize_acmg_class_label(nm)
      if (norm %in% names(counts)) counts[[norm]] <- as.integer(tab[[nm]])
    }
  }
  data.frame(
    classification = ACMG_CLASS_LABELS,
    variant_count = as.integer(unlist(counts[ACMG_CLASS_LABELS], use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}

compare_two_classification_tables <- function(left_df, right_df, left_name = "ACMGamp", right_name = "Reference") {
  cmp <- merge(
    left_df[, c("variant_key", "gene", "classification", "tier", "criteria")],
    right_df[, c("variant_key", "classification", "tier", "criteria")],
    by = "variant_key",
    suffixes = c("_left", "_right"),
    all = FALSE
  )
  if (nrow(cmp) == 0L) {
    return(list(
      comparison = cmp,
      metrics = list(
        n_overlap = 0L,
        n_left = nrow(left_df),
        n_right = nrow(right_df),
        exact_accuracy = NA_real_,
        tier_accuracy = NA_real_
      ),
      summary_by_class = data.frame()
    ))
  }

  names(cmp)[names(cmp) == "classification_left"] <- paste0("class_", left_name)
  names(cmp)[names(cmp) == "classification_right"] <- paste0("class_", right_name)
  names(cmp)[names(cmp) == "tier_left"] <- paste0("tier_", left_name)
  names(cmp)[names(cmp) == "tier_right"] <- paste0("tier_", right_name)
  names(cmp)[names(cmp) == "criteria_left"] <- paste0("criteria_", left_name)
  names(cmp)[names(cmp) == "criteria_right"] <- paste0("criteria_", right_name)

  class_left <- cmp[[paste0("class_", left_name)]]
  class_right <- cmp[[paste0("class_", right_name)]]
  tier_left <- cmp[[paste0("tier_", left_name)]]
  tier_right <- cmp[[paste0("tier_", right_name)]]

  cmp$exact_match <- class_left == class_right
  cmp$tier_match <- tier_left == tier_right

  metrics <- list(
    n_overlap = nrow(cmp),
    n_left = nrow(left_df),
    n_right = nrow(right_df),
    exact_accuracy = round(100 * mean(cmp$exact_match), 2),
    tier_accuracy = round(100 * mean(cmp$tier_match), 2),
    left_name = left_name,
    right_name = right_name
  )

  summary_by_class <- do.call(rbind, lapply(ACMG_CLASS_LABELS, function(cl) {
    idx <- class_right == cl
    data.frame(
      classification = cl,
      n_reference = sum(idx),
      n_exact_match = sum(idx & cmp$exact_match),
      n_tier_match = sum(idx & cmp$tier_match),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary_by_class) <- NULL

  list(comparison = cmp, metrics = metrics, summary_by_class = summary_by_class)
}

write_comparison_outputs <- function(result, output_prefix) {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  cmp_path <- paste0(output_prefix, ".comparison.csv")
  summary_path <- paste0(output_prefix, ".summary_by_class.csv")
  metrics_path <- paste0(output_prefix, ".metrics.txt")

  if (nrow(result$comparison) > 0L) {
    utils::write.csv(result$comparison, cmp_path, row.names = FALSE)
  }
  if (nrow(result$summary_by_class) > 0L) {
    utils::write.csv(result$summary_by_class, summary_path, row.names = FALSE)
  }

  m <- result$metrics
  cat(
    paste0(
      "ACMGamp vs ", m$right_name, " comparison\n",
      "Overlap variants: ", m$n_overlap, "\n",
      "ACMGamp variants: ", m$n_left, "\n",
      "Reference variants: ", m$n_right, "\n",
      "Exact 5-class accuracy (%): ", m$exact_accuracy, "\n",
      "Tier accuracy (pathogenic/benign/vus) (%): ", m$tier_accuracy, "\n"
    ),
    file = metrics_path
  )

  list(
    comparison_csv = if (file.exists(cmp_path)) cmp_path else NA_character_,
    summary_csv = if (file.exists(summary_path)) summary_path else NA_character_,
    metrics_txt = metrics_path
  )
}
