#!/usr/bin/env Rscript
# Parse HMC-1 SnpEff ANN VCF (hg19), filter rare coding variants, run ACMG.

args <- commandArgs(trailingOnly = TRUE)
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", initial_options[grep("^--file=", initial_options)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
} else {
  normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
}

acgm_root <- file.path(project_root, "results")
if (dir.exists(file.path(project_root, "R"))) {
  acgm_root <- project_root
  data_root <- normalizePath(file.path(project_root, ".."), winslash = "/", mustWork = FALSE)
} else {
  data_root <- project_root
  acgm_root <- file.path(project_root, "ClinicalVariantR")
}

vcf_path <- if (length(args) >= 2L) {
  normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
} else {
  file.path(data_root, "results", "HMC-1.coding.vcf.gz")
}
out_dir <- if (length(args) >= 3L) {
  normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
} else {
  file.path(data_root, "results")
}

if (!file.exists(vcf_path)) {
  stop("VCF not found: ", vcf_path)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

setwd(acgm_root)
source("global.R")

AUDIT_LOG_PATH <- file.path(acgm_root, "logs", "analysis_log.csv")

AF_CUTOFF <- 0.01

parse_info_tag <- function(info, tag) {
  m <- regexpr(paste0("(?:^|;)", tag, "=([^;]+)"), info, perl = TRUE)
  if (m[1L] == -1L) return(NA_character_)
  sub(paste0("^.*(?:^|;)", tag, "="), "", regmatches(info, m)[[1L]])
}

parse_info_num <- function(info, tag) {
  raw <- parse_info_tag(info, tag)
  if (is.na(raw) || !nzchar(raw)) return(NA_real_)
  nums <- vapply(unlist(strsplit(raw, "[,|]")), scalar_num, numeric(1L))
  nums <- nums[!is.na(nums)]
  if (length(nums) == 0L) NA_real_ else max(nums)
}

decode_clinvar_sig <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  codes <- c(
    "0" = "Uncertain significance",
    "1" = "not provided",
    "2" = "Benign",
    "3" = "Likely benign",
    "4" = "Likely pathogenic",
    "5" = "Pathogenic"
  )
  parts <- unlist(strsplit(as.character(x), "[,|]"))
  mapped <- unname(codes[parts])
  mapped <- mapped[!is.na(mapped)]
  if (length(mapped) == 0L) as.character(x) else paste(unique(mapped), collapse = ";")
}

parse_ann_primary <- function(ann_string, alt = NA_character_) {
  if (is.na(ann_string) || !nzchar(ann_string)) {
    return(list(
      gene = NA_character_, consequence = NA_character_, impact = NA_character_,
      hgvs_c = NA_character_, hgvs_p = NA_character_, polyphen = NA_character_
    ))
  }

  entries <- strsplit(ann_string, ",", fixed = TRUE)[[1L]]
  parsed <- lapply(entries, function(entry) {
  parts <- strsplit(entry, "|", fixed = TRUE)[[1L]]
  if (length(parts) < 4L) return(NULL)
  list(
    allele = parts[[1L]],
    consequence = parts[[2L]],
    impact = parts[[3L]],
    gene = parts[[4L]],
    hgvs_c = if (length(parts) >= 10L) parts[[10L]] else NA_character_,
    hgvs_p = if (length(parts) >= 11L) parts[[11L]] else NA_character_
  )
  })
  parsed <- Filter(Negate(is.null), parsed)
  if (length(parsed) == 0L) {
    return(list(
      gene = NA_character_, consequence = NA_character_, impact = NA_character_,
      hgvs_c = NA_character_, hgvs_p = NA_character_, polyphen = NA_character_
    ))
  }

  cons <- vapply(parsed, `[[`, "", "consequence")
  rank <- ifelse(grepl("stop_gained|frameshift|splice_donor|splice_acceptor|stop_lost|start_lost", cons), 1L,
    ifelse(grepl("missense_variant", cons), 2L,
      ifelse(grepl("inframe_", cons), 3L, 4L)))
  best <- parsed[[which.min(rank)]]

  list(
    gene = best$gene,
    consequence = best$consequence,
    impact = best$impact,
    hgvs_c = best$hgvs_c,
    hgvs_p = best$hgvs_p,
    polyphen = NA_character_
  )
}

parse_vcf_ann <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)

  header_cols <- NULL
  rows <- list()

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (grepl("^#CHROM\t", line)) {
      header_cols <- strsplit(sub("^#", "", line), "\t", fixed = TRUE)[[1L]]
      next
    }
    if (grepl("^#", line)) next

    parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 8L) next

    chrom <- parts[[1L]]
    pos <- as.integer(parts[[2L]])
    ref <- parts[[4L]]
    alt <- parts[[5L]]
    info <- parts[[8L]]
    sample_col <- if (length(parts) >= 10L) parts[[10L]] else NA_character_

    ann <- parse_info_tag(info, "ANN")
    ann_fields <- parse_ann_primary(ann, alt)
    af_1kg <- parse_info_num(info, "1000Gp3_AF")
    af_esp <- parse_info_num(info, "ESP6500_MAF")
    af_use <- if (!is.na(af_1kg)) af_1kg else if (!is.na(af_esp)) af_esp / 100 else NA_real_
    ph <- parse_info_tag(info, "ESP6500_PH")
    polyphen_score <- NA_real_
    if (!is.na(ph) && grepl(":", ph)) {
      polyphen_score <- scalar_num(sub(".*:", "", ph))
    }
    clinvar <- decode_clinvar_sig(parse_info_tag(info, "CLINVAR_CLNSIG"))

    gt <- NA_character_
    if (!is.na(sample_col) && grepl(":", sample_col)) {
      gt <- strsplit(sample_col, ":", fixed = TRUE)[[1L]][[1L]]
    }

    rows[[length(rows) + 1L]] <- data.frame(
      variant_id = paste(chrom, pos, ref, alt, sep = ":"),
      chrom = chrom,
      pos = pos,
      ref = ref,
      alt = alt,
      gene = ann_fields$gene,
      consequence = ann_fields$consequence,
      impact = ann_fields$impact,
      hgvs_c = ann_fields$hgvs_c,
      hgvs_p = ann_fields$hgvs_p,
      genotype = gt,
      af_1000g = af_1kg,
      af_esp6500 = af_esp,
      population_af = af_use,
      polyphen = ph,
      polyphen_score = polyphen_score,
      clinvar_classification = clinvar,
      AF = af_use,
      REVEL = NA_real_,
      ClinVar = clinvar,
      gnomad_af = af_use,
      revel_score = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) {
    stop("No variants parsed from: ", path)
  }

  do.call(rbind, rows)
}

message("Parsing: ", vcf_path)
variants <- parse_vcf_ann(vcf_path)
message("Parsed variants: ", nrow(variants))

parsed_path <- file.path(out_dir, "HMC-1.parsed_variants.csv")
write.csv(variants, parsed_path, row.names = FALSE)
message("Wrote: ", parsed_path)

message("Scoring detailed ACMG criteria with transparent evidence...")
lof_panel <- file.path(acgm_root, "data", "gene_panels", "lof_disease_mechanism_genes.csv")
detailed <- score_variants_table(variants, lof_panel_path = lof_panel, profile_id = "hematologic_predisposition")
evidence_report <- acmg_pro_to_report(detailed, mode = "rapid", session_id = "HMC1-CLI")

detailed_path <- file.path(out_dir, "HMC-1.acmg_detailed.csv")
write.csv(detailed, detailed_path, row.names = FALSE)
message("Wrote: ", detailed_path)

evidence_report_path <- file.path(out_dir, "HMC-1.evidence_report.csv")
write.csv(evidence_report, evidence_report_path, row.names = FALSE)
message("Wrote: ", evidence_report_path)

metadata <- build_run_metadata(vcf_path = vcf_path, profile_id = "hematologic_predisposition", mode = "rapid", session_id = "HMC1-CLI")
write_run_metadata_json(metadata, file.path(out_dir, "HMC-1.metadata.json"))
message("Wrote: ", file.path(out_dir, "HMC-1.metadata.json"))

rare_detailed <- detailed[is.na(detailed$max_population_af) | detailed$max_population_af <= AF_CUTOFF, , drop = FALSE]
message("Rare variants (AF <= ", AF_CUTOFF, " or missing): ", nrow(rare_detailed))

rare_path <- file.path(out_dir, "HMC-1.rare_coding_variants.csv")
write.csv(rare_detailed, rare_path, row.names = FALSE)
message("Wrote: ", rare_path)

acmg_path <- file.path(out_dir, "HMC-1.acmg_results.csv")
write.csv(rare_detailed, acmg_path, row.names = FALSE)
message("Wrote: ", acmg_path)

priority <- rare_detailed[order(
  match(rare_detailed$classification, c("Pathogenic", "Likely Pathogenic", "VUS", "Likely Benign", "Benign")),
  rare_detailed$max_population_af
), , drop = FALSE]

top <- priority[priority$classification %in% c("Pathogenic", "Likely Pathogenic", "VUS"), , drop = FALSE]
if (nrow(top) > 200L) top <- top[seq_len(200L), , drop = FALSE]

top_path <- file.path(out_dir, "HMC-1.top_candidates.csv")
write.csv(top, top_path, row.names = FALSE)
message("Wrote: ", top_path)

criteria_summary <- acmg_pro_criteria_summary(rare_detailed)
summary_path <- file.path(out_dir, "HMC-1.acmg_criteria_summary.csv")
write.csv(criteria_summary, summary_path, row.names = FALSE)
message("Wrote: ", summary_path)

cat("\nSummary\n")
cat("-------\n")
cat("Input variants:     ", nrow(variants), "\n", sep = "")
cat("Rare variants:      ", nrow(rare_detailed), "\n", sep = "")
cat("ACMG classified:    ", nrow(rare_detailed), "\n", sep = "")
print(table(rare_detailed$classification, useNA = "ifany"))
cat("\nCriteria triggered (rare set):\n")
print(criteria_summary)
cat("\nDone.\n")
