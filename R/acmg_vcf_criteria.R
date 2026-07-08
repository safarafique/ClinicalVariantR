#' Automated ACMG criteria from VCF-derived variant tables.
#'
#' Supports population AF, SnpEff ANN / VEP CSQ consequence, in silico scores,
#' ClinVar fields, and Group A clinical/pedigree context (28-criteria reporting).

DEFAULT_ACMG_THRESHOLDS <- list(
  ba1_af = 0.05,
  bs1_af = 0.01,
  pm2_af = 0.01,
  pm2_strict_af = 0.001,
  bs2_af = 0.001,
  revel_pp3 = 0.75,
  revel_bp4 = 0.15,
  cadd_pp3 = 20,
  cadd_bp4 = 10,
  spliceai_pp3 = 0.5,
  spliceai_bp4 = 0.1,
  alphamissense_pp3 = 0.564,
  alphamissense_bp4 = 0.34,
  polyphen_damaging = 0.85,
  polyphen_benign = 0.15,
  insilico_min_tools = 2,
  insilico_revel_solo = TRUE,
  confidence_base = 50,
  confidence_per_criterion = 8,
  confidence_conflict_penalty = 15
)

threshold_num <- function(thresholds, key, default) {
  val <- thresholds[[key]]
  if (is.null(val) || length(val) == 0L) return(default)
  num <- suppressWarnings(as.numeric(val)[1L])
  if (length(num) == 0L || is.na(num)) default else num
}

resolve_thresholds <- function(thresholds = NULL, profile_id = DEFAULT_PROFILE_ID) {
  if (!is.null(thresholds)) {
    out <- modifyList(DEFAULT_ACMG_THRESHOLDS, thresholds)
  } else {
    cfg <- load_rule_config(profile_id)
    out <- modifyList(DEFAULT_ACMG_THRESHOLDS, cfg$thresholds)
  }
  for (key in names(DEFAULT_ACMG_THRESHOLDS)) {
    out[[key]] <- threshold_num(out, key, DEFAULT_ACMG_THRESHOLDS[[key]])
  }
  if (is_prediction_mode()) {
    out <- modifyList(out, prediction_threshold_overrides())
  }
  out
}

CRITERION_FLAG_FIELDS <- c(
  "BA1", "BS1", "BS2", "BS3", "BS4",
  "PM2", "PVS1", "PS1", "PS2", "PS3", "PS4",
  "PM1", "PM3", "PM4", "PM5", "PM6",
  "PP1", "PP2", "PP3", "PP4", "PP5",
  "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
  "gwas_supplementary"
)

normalize_criterion_flags <- function(scores) {
  for (fld in CRITERION_FLAG_FIELDS) {
    if (fld %in% names(scores)) scores[[fld]] <- isTRUE(scores[[fld]])
  }
  scores
}

LOF_CONSEQUENCES <- c(
  "stop_gained", "nonsense", "frameshift", "frameshift_variant",
  "splice_donor", "splice_donor_variant", "splice_acceptor", "splice_acceptor_variant",
  "start_lost", "stop_lost", "stopgain"
)

normalize_consequence <- function(consequence) {
  if (is.null(consequence) || length(consequence) == 0L) return("")
  x <- tolower(paste(as.character(consequence), collapse = "&"))
  x <- gsub("nonframeshift", "inframe", x, fixed = TRUE)
  x <- gsub("stoploss", "stop_lost", x, fixed = TRUE)
  x
}

load_gene_panel <- function(path) {
  if (!file.exists(path)) return(character())
  unique(trimws(read.csv(path, stringsAsFactors = FALSE)$gene))
}

load_lof_gene_panel <- function(path = file.path("data", "gene_panels", "lof_disease_mechanism_genes.csv")) {
  load_gene_panel(path)
}

PM4_CONSEQUENCES <- c(
  "inframe_insertion", "inframe_deletion", "stop_lost"
)

BP3_CONSEQUENCES <- c(
  "inframe_insertion", "inframe_deletion"
)


best_population_af <- function(af_1000g = NA_real_, af_esp6500 = NA_real_,
                               gnomad_af = NA_real_, population_af = NA_real_) {
  af_1000g <- scalar_num(af_1000g)
  af_esp6500 <- scalar_num(af_esp6500)
  gnomad_af <- scalar_num(gnomad_af)
  population_af <- scalar_num(population_af)
  vals <- c(
    af_1000g,
    if (!is.na(af_esp6500)) af_esp6500 / 100 else NA_real_,
    gnomad_af,
    population_af
  )
  vals <- suppressWarnings(as.numeric(vals))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) return(NA_real_)
  max(vals)
}

VEP_POPMAX_AF_FIELDS <- c(
  "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
  "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_ASJ_AF", "gnomADe_EAS_AF",
  "gnomADe_FIN_AF", "gnomADe_MID_AF", "gnomADe_NFE_AF", "gnomADe_REMAINING_AF", "gnomADe_SAS_AF",
  "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF", "gnomADg_ASJ_AF",
  "gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF", "gnomADg_NFE_AF",
  "gnomADg_REMAINING_AF", "gnomADg_SAS_AF"
)

best_popmax_af <- function(popmax_af = NA_real_, row = NULL) {
  vals <- suppressWarnings(as.numeric(popmax_af))
  vals <- vals[!is.na(vals)]
  if (!is.null(row)) {
    for (fld in VEP_POPMAX_AF_FIELDS) {
      if (fld %in% names(row)) {
        v <- scalar_num(row[[fld]])
        if (!is.na(v)) vals <- c(vals, v)
      }
    }
  }
  if (length(vals) == 0L) return(NA_real_)
  max(vals)
}

use_popmax_for_frequency_rules <- function() {
  is_prediction_mode() && isTRUE(PREDICTION_SETTINGS$use_popmax_for_pm2_ba1)
}

resolve_scoring_population_af <- function(row) {
  global_af <- best_population_af(
    af_1000g = scalar_num(row$af_1000g %||% NA_real_),
    af_esp6500 = scalar_num(row$af_esp6500 %||% NA_real_),
    gnomad_af = scalar_num(row$gnomad_af %||% NA_real_),
    population_af = scalar_num(row$population_af %||% NA_real_)
  )
  popmax_af <- best_popmax_af(
    popmax_af = scalar_num(row$popmax_af %||% row$ref_gnomad_popmax %||% NA_real_),
    row = row
  )
  use_popmax <- use_popmax_for_frequency_rules()
  effective <- if (use_popmax && !is.na(popmax_af)) popmax_af else global_af
  list(
    effective = effective,
    global = global_af,
    popmax = popmax_af,
    af_known = !is.na(effective)
  )
}

is_lof_consequence <- function(consequence) {
  cons <- normalize_consequence(consequence)
  if (!nzchar(cons)) return(FALSE)
  any(vapply(LOF_CONSEQUENCES, function(x) grepl(x, cons, fixed = TRUE), logical(1L)))
}

is_missense_consequence <- function(consequence) {
  any(grepl("missense", normalize_consequence(consequence), fixed = TRUE))
}

is_synonymous_consequence <- function(consequence) {
  cons <- normalize_consequence(consequence)
  any(grepl("synonymous", cons, fixed = TRUE)) || any(grepl("coding-synon", cons, fixed = TRUE))
}

is_pm4_consequence <- function(consequence) {
  cons <- normalize_consequence(consequence)
  if (!nzchar(cons)) return(FALSE)
  any(vapply(PM4_CONSEQUENCES, function(x) grepl(x, cons, fixed = TRUE), logical(1L)))
}

is_bp3_consequence <- function(consequence) {
  cons <- normalize_consequence(consequence)
  if (!nzchar(cons)) return(FALSE)
  any(vapply(BP3_CONSEQUENCES, function(x) grepl(x, cons, fixed = TRUE), logical(1L)))
}

parse_polyphen_call <- function(polyphen_text, polyphen_score = NA_real_) {
  polyphen_text <- scalar_chr(polyphen_text, default = "")
  polyphen_score <- scalar_num(polyphen_score)
  text <- tolower(paste0(polyphen_text, collapse = " "))
  score <- polyphen_score
  if (!is.na(score) && score >= 0) {
    if (score >= 0.85) return("damaging")
    if (score <= 0.15) return("benign")
  }
  if (grepl("probably_damaging|possibly_damaging|damaging", text)) return("damaging")
  if (grepl("benign", text)) return("benign")
  "unknown"
}

parse_sift_call <- function(sift_text) {
  sift_text <- scalar_chr(sift_text, default = "")
  text <- tolower(paste0(sift_text, collapse = " "))
  if (grepl("deleterious", text)) return("deleterious")
  if (grepl("tolerated", text)) return("tolerated")
  "unknown"
}

clinvar_has_conflict <- function(clinvar_text) {
  clinvar_text <- scalar_chr(clinvar_text, default = "")
  if (!nzchar(clinvar_text)) return(FALSE)
  if (grepl("conflict|conflicting", clinvar_text, ignore.case = TRUE)) return(TRUE)
  has_path <- isTRUE(grepl("pathogenic", clinvar_text, ignore.case = TRUE))
  has_benign <- isTRUE(grepl("benign", clinvar_text, ignore.case = TRUE))
  isTRUE(has_path && has_benign)
}

clinvar_is_pathogenic <- function(clinvar_text) {
  clinvar_text <- scalar_chr(clinvar_text, default = "")
  if (!nzchar(clinvar_text)) return(FALSE)
  if (clinvar_has_conflict(clinvar_text)) return(FALSE)
  grepl("pathogenic", clinvar_text, ignore.case = TRUE) &&
    !grepl("benign|conflict|uncertain", clinvar_text, ignore.case = TRUE)
}

clinvar_is_benign <- function(clinvar_text) {
  clinvar_text <- scalar_chr(clinvar_text, default = "")
  if (!nzchar(clinvar_text)) return(FALSE)
  if (clinvar_has_conflict(clinvar_text)) return(FALSE)
  grepl("likely benign|likely_benign|\\bbenign\\b", clinvar_text, ignore.case = TRUE) &&
    !grepl("pathogenic|conflict", clinvar_text, ignore.case = TRUE)
}

score_population_criteria <- function(af, thresholds = DEFAULT_ACMG_THRESHOLDS, af_known = TRUE) {
  af <- scalar_num(af, default = NA_real_)
  ba1_af <- threshold_num(thresholds, "ba1_af", 0.05)
  bs1_af <- threshold_num(thresholds, "bs1_af", 0.01)
  bs2_af <- threshold_num(thresholds, "bs2_af", 0.001)
  pm2_af <- threshold_num(thresholds, "pm2_af", 0.01)
  pm2_cutoff <- threshold_num(thresholds, "pm2_strict_af", bs2_af)
  out <- list(
    BA1 = FALSE, BS1 = FALSE, BS2 = FALSE, PM2 = FALSE,
    BA1_rationale = "", BS1_rationale = "", BS2_rationale = "", PM2_rationale = ""
  )

  if (!isTRUE(af_known) || is.na(af)) {
    out$PM2_rationale <- paste(
      "Population allele frequency unavailable;",
      if (is_prediction_mode()) {
        "PM2 withheld in prediction mode (verify gnomAD/Popmax before clinical use)."
      } else {
        "PM2/BA1/BS rules skipped until AF is annotated."
      }
    )
    return(out)
  }

  if (!is.na(af) && isTRUE(af > ba1_af)) {
    out$BA1 <- TRUE
    out$BA1_rationale <- sprintf(
      "Population MAX_AF = %.8f. Threshold: >%.4f (%.0f%%). Evidence: BA1 (Stand-alone).",
      af, ba1_af, ba1_af * 100
    )
    return(out)
  }

  if (!is.na(af) && isTRUE(af > bs1_af)) {
    out$BS1 <- TRUE
    out$BS1_rationale <- sprintf(
      "Population MAX_AF = %.8f. Threshold: >%.4f (%.0f%%). Evidence: BS1 (Strong).",
      af, bs1_af, bs1_af * 100
    )
  }

  if (!is.na(af) && isTRUE(af >= bs2_af && af < bs1_af)) {
    out$BS2 <- TRUE
    out$BS2_rationale <- sprintf(
      "Population MAX_AF = %.8f is below disorder threshold but observed in population databases. Evidence: BS2 (Strong).",
      af
    )
  }

  if (!isTRUE(out$BS2) && (af == 0 || isTRUE(af < pm2_cutoff))) {
    out$PM2 <- TRUE
    if (af == 0) {
      out$PM2_rationale <- sprintf(
        "Population MAX_AF = 0. Threshold: <%.4f. Evidence: PM2 (Moderate).",
        pm2_af
      )
    } else {
      popmax_note <- if (use_popmax_for_frequency_rules()) " (popmax-aware)" else ""
      out$PM2_rationale <- sprintf(
        "Population MAX_AF%s = %.8f. Threshold: <%.4f. Evidence: PM2 (Moderate).",
        popmax_note, af, pm2_af
      )
    }
  }

  out
}

score_gene_mechanism_criteria <- function(consequence, gene,
                                          pp2_genes = character(),
                                          bp1_genes = character()) {
  out <- list(PP2 = FALSE, BP1 = FALSE, PP2_rationale = "", BP1_rationale = "")
  gene <- scalar_chr(gene, default = NA_character_)
  if (is_missense_consequence(consequence) && !is.na(gene) && gene %in% pp2_genes) {
    out$PP2 <- TRUE
    out$PP2_rationale <- sprintf(
      "Missense in %s; gene panel indicates missense is a common disease mechanism.", gene
    )
  }
  if (is_missense_consequence(consequence) && !is.na(gene) && gene %in% bp1_genes) {
    out$BP1 <- TRUE
    out$BP1_rationale <- sprintf(
      "Missense in %s; truncating variants are the primary disease mechanism.", gene
    )
  }
  out
}

score_bp7_criteria <- function(consequence, sift = NA_character_, polyphen = NA_character_,
                               polyphen_score = NA_real_) {
  consequence <- scalar_chr(consequence, default = "")
  out <- list(BP7 = FALSE, BP7_rationale = "")
  if (!is_synonymous_consequence(consequence)) return(out)
  sift_call <- parse_sift_call(sift)
  ph_call <- parse_polyphen_call(polyphen, polyphen_score)
  if (sift_call %in% c("tolerated", "unknown") && ph_call %in% c("benign", "unknown")) {
    out$BP7 <- TRUE
    out$BP7_rationale <- "Synonymous variant with no strong in silico splice/damage signal (BP7)."
  }
  out
}

score_consequence_criteria <- function(consequence, gene, lof_genes = character(),
                                       is_protein_coding = FALSE,
                                       is_canonical_transcript = FALSE,
                                       impact = NA_character_,
                                       ref = NA_character_, alt = NA_character_,
                                       exon = NA_character_, biotype = NA_character_) {
  consequence <- scalar_chr(consequence, default = "")
  gene <- scalar_chr(gene, default = NA_character_)
  ref <- scalar_chr(ref, default = "")
  alt <- scalar_chr(alt, default = "")
  is_protein_coding <- scalar_lgl(is_protein_coding)
  is_canonical_transcript <- scalar_lgl(is_canonical_transcript)
  out <- list(
    PVS1 = FALSE, PM4 = FALSE, BP3 = FALSE,
    PVS1_rationale = "", PM4_rationale = "", BP3_rationale = ""
  )

  pvs1 <- score_autopvs1(
    consequence = consequence,
    gene = gene,
    lof_genes = lof_genes,
    is_protein_coding = is_protein_coding,
    is_canonical_transcript = is_canonical_transcript,
    impact = impact,
    exon = exon,
    biotype = biotype
  )
  out$PVS1 <- isTRUE(pvs1$PVS1)
  out$PVS1_rationale <- pvs1$PVS1_rationale %||% ""

  if (is_pm4_consequence(consequence)) {
    out$PM4 <- TRUE
    out$PM4_rationale <- sprintf("Protein-length altering consequence: %s.", consequence)
  }

  if (is_bp3_consequence(consequence)) {
    if (is_homopolymer_indel(ref, alt)) {
      out$BP3 <- TRUE
      out$BP3_rationale <- sprintf(
        "In-frame indel (%s) in homopolymer/repetitive sequence (BP3).", consequence
      )
    } else {
      out$BP3_rationale <- sprintf(
        "In-frame indel (%s); not in annotated repetitive region (BP3 not applied).", consequence
      )
    }
  }

  out
}

score_insilico_criteria <- function(
    revel = NA_real_,
    polyphen = NA_character_,
    polyphen_score = NA_real_,
    sift = NA_character_,
    cadd = NA_real_,
    spliceai_max = NA_real_,
    alphamissense = NA_real_,
    thresholds = DEFAULT_ACMG_THRESHOLDS) {

  revel <- scalar_num(revel)
  polyphen <- scalar_chr(polyphen, default = NA_character_)
  polyphen_score <- scalar_num(polyphen_score)
  sift <- scalar_chr(sift, default = NA_character_)
  cadd <- scalar_num(cadd)
  spliceai_max <- scalar_num(spliceai_max)
  alphamissense <- scalar_num(alphamissense)

  out <- list(
    PP3 = FALSE, BP4 = FALSE,
    PP3_rationale = "", BP4_rationale = "",
    insilico_summary = ""
  )

  ph_call <- parse_polyphen_call(polyphen, polyphen_score)
  sift_call <- parse_sift_call(sift)
  tool_notes <- character()
  damaging_hits <- 0L
  benign_hits <- 0L
  min_tools <- as.integer(thresholds$insilico_min_tools %||% 2)
  revel_solo <- isTRUE(thresholds$insilico_revel_solo %||% TRUE)

  add_tool <- function(name, value_text, damaging = FALSE, benign = FALSE) {
    tool_notes <<- c(tool_notes, sprintf("%s=%s", name, value_text))
    if (isTRUE(damaging)) damaging_hits <<- damaging_hits + 1L
    if (isTRUE(benign)) benign_hits <<- benign_hits + 1L
  }

  if (!is.na(revel)) {
    revel_pp3 <- threshold_num(thresholds, "revel_pp3", 0.75)
    revel_bp4 <- threshold_num(thresholds, "revel_bp4", 0.15)
    add_tool("REVEL", sprintf("%.3f", revel),
             damaging = isTRUE(revel >= revel_pp3),
             benign = isTRUE(revel <= revel_bp4))
    if (!isTRUE(out$PP3) && isTRUE(revel >= revel_pp3)) {
      out$PP3 <- TRUE
      out$PP3_rationale <- sprintf(
        "REVEL score = %.3f. Threshold: >=%.2f. Evidence: PP3 (Supporting).",
        revel, revel_pp3
      )
    } else if (!isTRUE(out$BP4) && isTRUE(revel <= revel_bp4)) {
      out$BP4 <- TRUE
      out$BP4_rationale <- sprintf(
        "REVEL score = %.3f. Threshold: <=%.2f. Evidence: BP4 (Supporting).",
        revel, revel_bp4
      )
    }
  }

  if (!is.na(cadd)) {
    cadd_pp3 <- threshold_num(thresholds, "cadd_pp3", 20)
    cadd_bp4 <- threshold_num(thresholds, "cadd_bp4", 10)
    if (!isTRUE(out$PP3) && isTRUE(cadd >= cadd_pp3)) {
      out$PP3 <- TRUE
      out$PP3_rationale <- sprintf(
        "CADD Phred = %.1f. Threshold: >=%.0f. Evidence: PP3 (Supporting).",
        cadd, cadd_pp3
      )
    } else if (!isTRUE(out$BP4) && isTRUE(cadd <= cadd_bp4)) {
      out$BP4 <- TRUE
      out$BP4_rationale <- sprintf(
        "CADD Phred = %.1f. Threshold: <=%.0f. Evidence: BP4 (Supporting).",
        cadd, cadd_bp4
      )
    }
    add_tool("CADD", sprintf("%.1f", cadd),
             damaging = isTRUE(cadd >= cadd_pp3),
             benign = isTRUE(cadd <= cadd_bp4))
  }

  if (!is.na(spliceai_max)) {
    spliceai_pp3 <- threshold_num(thresholds, "spliceai_pp3", 0.5)
    spliceai_bp4 <- threshold_num(thresholds, "spliceai_bp4", 0.1)
    if (!isTRUE(out$PP3) && isTRUE(spliceai_max >= spliceai_pp3)) {
      out$PP3 <- TRUE
      out$PP3_rationale <- sprintf(
        "SpliceAI max delta = %.3f. Threshold: >=%.2f. Evidence: PP3 (Supporting).",
        spliceai_max, spliceai_pp3
      )
    } else if (!isTRUE(out$BP4) && isTRUE(spliceai_max <= spliceai_bp4)) {
      out$BP4 <- TRUE
      out$BP4_rationale <- sprintf(
        "SpliceAI max delta = %.3f. Threshold: <=%.2f. Evidence: BP4 (Supporting).",
        spliceai_max, spliceai_bp4
      )
    }
    add_tool("SpliceAI", sprintf("%.3f", spliceai_max),
             damaging = isTRUE(spliceai_max >= spliceai_pp3),
             benign = isTRUE(spliceai_max <= spliceai_bp4))
  }

  if (!is.na(alphamissense)) {
    alphamissense_pp3 <- threshold_num(thresholds, "alphamissense_pp3", 0.564)
    alphamissense_bp4 <- threshold_num(thresholds, "alphamissense_bp4", 0.34)
    if (!isTRUE(out$PP3) && isTRUE(alphamissense >= alphamissense_pp3)) {
      out$PP3 <- TRUE
      out$PP3_rationale <- sprintf(
        "AlphaMissense score = %.3f. Threshold: >=%.3f. Evidence: PP3 (Supporting).",
        alphamissense, alphamissense_pp3
      )
    } else if (!isTRUE(out$BP4) && isTRUE(alphamissense <= alphamissense_bp4)) {
      out$BP4 <- TRUE
      out$BP4_rationale <- sprintf(
        "AlphaMissense score = %.3f. Threshold: <=%.3f. Evidence: BP4 (Supporting).",
        alphamissense, alphamissense_bp4
      )
    }
    add_tool("AlphaMissense", sprintf("%.3f", alphamissense),
             damaging = isTRUE(alphamissense >= alphamissense_pp3),
             benign = isTRUE(alphamissense <= alphamissense_bp4))
  }

  if (ph_call == "damaging") {
    damaging_hits <- damaging_hits + 1L
    add_tool("PolyPhen", "damaging", damaging = TRUE)
  } else if (ph_call == "benign") {
    benign_hits <- benign_hits + 1L
    add_tool("PolyPhen", "benign", benign = TRUE)
  }

  if (sift_call == "deleterious") {
    damaging_hits <- damaging_hits + 1L
    add_tool("SIFT", "deleterious", damaging = TRUE)
  } else if (sift_call == "tolerated") {
    benign_hits <- benign_hits + 1L
    add_tool("SIFT", "tolerated", benign = TRUE)
  }

  effective_min <- if (isTRUE(revel_solo)) 1L else min_tools

  if (!isTRUE(out$PP3) && damaging_hits >= effective_min) {
    out$PP3 <- TRUE
    out$PP3_rationale <- sprintf(
      "%d in silico tools support damaging effect (threshold: >=%d). Evidence: PP3 (Supporting). Tools: %s.",
      damaging_hits, effective_min, paste(tool_notes, collapse = "; ")
    )
  }
  if (!isTRUE(out$BP4) && benign_hits >= effective_min) {
    out$BP4 <- TRUE
    out$BP4_rationale <- sprintf(
      "%d in silico tools support benign effect (threshold: >=%d). Evidence: BP4 (Supporting). Tools: %s.",
      benign_hits, effective_min, paste(tool_notes, collapse = "; ")
    )
  }

  if (!isTRUE(out$PP3) && !isTRUE(out$BP4)) {
    if (damaging_hits == 1L) {
      out$PP3_rationale <- sprintf("Single-tool damaging signal (partial PP3 review): %s.", paste(tool_notes, collapse = "; "))
    } else if (benign_hits == 1L) {
      out$BP4_rationale <- sprintf("Single-tool benign signal (partial BP4 review): %s.", paste(tool_notes, collapse = "; "))
    }
  }

  out$insilico_summary <- if (length(tool_notes) > 0) paste(tool_notes, collapse = "; ") else ""
  out
}

score_clinvar_criteria <- function(clinvar_text) {
  clinvar_text <- scalar_chr(clinvar_text, default = "")
  out <- list(PP5 = FALSE, BP6 = FALSE, PP5_rationale = "", BP6_rationale = "")

  if (clinvar_has_conflict(clinvar_text)) {
    out$PP5_rationale <- "ClinVar entry has conflict/uncertain status; PP5/BP6 not applied."
    out$BP6_rationale <- out$PP5_rationale
    return(out)
  }

  if (clinvar_is_pathogenic(clinvar_text)) {
    out$PP5 <- TRUE
    out$PP5_rationale <- sprintf("ClinVar reports pathogenic/likely pathogenic: %s.", clinvar_text)
  }

  if (clinvar_is_benign(clinvar_text)) {
    out$BP6 <- TRUE
    out$BP6_rationale <- sprintf("ClinVar reports benign/likely benign: %s.", clinvar_text)
  }

  out
}

criteria_to_evidence <- function(scores) {
  flag_count <- function(code) if (isTRUE(scores[[code]])) 1L else 0L
  bs_count <- safe_int(scores$BS_count %||% 0L)
  if (isTRUE(scores$BS1)) bs_count <- bs_count + 1L
  if (isTRUE(scores$BS2)) bs_count <- bs_count + 1L
  ps_count <- safe_int(scores$PS_count %||% 0L)
  if (isTRUE(scores$PS1)) ps_count <- ps_count + 1L
  if (isTRUE(scores$PS2)) ps_count <- ps_count + 1L
  if (isTRUE(scores$PS3)) ps_count <- ps_count + 1L
  if (isTRUE(scores$PS4)) ps_count <- ps_count + 1L
  pp_count <- sum(c(
    flag_count("PP2"), flag_count("PP3"), flag_count("PP5"),
    flag_count("PP1"), flag_count("PP4"),
    safe_int(scores$PP_manual %||% 0L)
  ))
  list(
    PVS = flag_count("PVS1"),
    PS  = ps_count,
    PM  = sum(c(
      flag_count("PM1"), flag_count("PM2"), flag_count("PM3"),
      flag_count("PM4"), flag_count("PM5"), flag_count("PM6"),
      safe_int(scores$PM_manual %||% 0L)
    )),
    PP  = pp_count,
    BA  = isTRUE(scores$BA1),
    BS  = bs_count + flag_count("BS3") + flag_count("BS4"),
    BP  = sum(c(
      flag_count("BP1"), flag_count("BP2"), flag_count("BP3"), flag_count("BP4"),
      flag_count("BP5"), flag_count("BP6"), flag_count("BP7"),
      safe_int(scores$BP_manual %||% 0L)
    ))
  )
}

collect_met_criteria <- function(scores) {
  met <- FULL_ACMG_CRITERIA[vapply(FULL_ACMG_CRITERIA, function(code) isTRUE(scores[[code]]), logical(1L))]
  paste(met, collapse = ";")
}

collect_criteria_rationale <- function(scores) {
  fields <- paste0(FULL_ACMG_CRITERIA, "_rationale")
  vals <- unlist(scores[fields], use.names = FALSE)
  vals <- vals[nzchar(vals %||% "")]
  paste(vals, collapse = " | ")
}

score_variant_row <- function(row, lof_genes = character(),
                              pp2_genes = character(), bp1_genes = character(),
                              pm1_genes = character(),
                              clinvar_protein_db = NULL,
                              ps4_db = NULL,
                              pm1_domains = NULL,
                              clinical_context = NULL,
                              pedigree_context = NULL,
                              refs = NULL,
                              thresholds = NULL,
                              profile_id = DEFAULT_PROFILE_ID,
                              criteria_meta = list(),
                              evidence_scope = c("triggered", "automated", "full")) {
  evidence_scope <- match.arg(evidence_scope)
  thresholds <- resolve_thresholds(thresholds, profile_id)
  if (length(criteria_meta) == 0L) {
    criteria_meta <- load_rule_config(profile_id)$criteria_meta
  }
  row <- normalize_variant_row_input(row)
  if (is.null(clinvar_protein_db)) clinvar_protein_db <- load_clinvar_protein_db()
  if (is.null(ps4_db)) ps4_db <- load_ps4_case_control_db()
  gwas_db <- load_gwas_supplementary_db()
  if (is.null(pm1_domains)) pm1_domains <- load_pm1_critical_domains()

  af_info <- resolve_scoring_population_af(row)
  af <- af_info$effective
  pop <- score_population_criteria(af, thresholds, af_known = af_info$af_known)
  variant_gene <- scalar_chr(row$gene %||% "")
  variant_consequence <- scalar_chr(row$consequence %||% "")
  variant_ref <- scalar_chr(row$ref %||% "")
  variant_alt <- scalar_chr(row$alt %||% "")
  cons <- score_consequence_criteria(
    variant_consequence, variant_gene, lof_genes,
    is_protein_coding = scalar_lgl(row$is_protein_coding),
    is_canonical_transcript = scalar_lgl(row$is_canonical_transcript),
    impact = scalar_chr(row$impact %||% ""),
    ref = variant_ref, alt = variant_alt,
    exon = scalar_chr(row$exon %||% ""),
    biotype = scalar_chr(row$biotype %||% "")
  )
  csq_catalog <- parse_csq_pathogenic_catalog(
    scalar_chr(row$info_csq %||% ""), variant_ref, variant_alt
  )
  pathogenic_protein <- score_ps1_pm5_pm1_criteria(
    gene = variant_gene,
    consequence = variant_consequence,
    hgvs_p = scalar_chr(row$hgvs_p %||% ""),
    amino_acids = scalar_chr(row$amino_acids %||% ""),
    protein_position = scalar_chr(row$protein_position %||% ""),
    clinvar_db = clinvar_protein_db,
    csq_catalog = csq_catalog,
    pm1_genes = pm1_genes,
    pm1_domains = pm1_domains
  )
  ps4 <- score_ps4_criteria(row, ps4_db = ps4_db, gwas_db = gwas_db, thresholds = thresholds)
  gene_rules <- score_gene_mechanism_criteria(
    variant_consequence, variant_gene, pp2_genes, bp1_genes
  )
  bp7 <- score_bp7_criteria(
    variant_consequence,
    scalar_chr(row$sift %||% "", default = NA_character_),
    scalar_chr(row$polyphen %||% "", default = NA_character_),
    scalar_num(row$polyphen_score %||% NA_real_)
  )
  ins <- score_insilico_criteria(
    revel = scalar_num(row$revel_score %||% row$REVEL %||% NA_real_),
    polyphen = scalar_chr(row$polyphen %||% "", default = NA_character_),
    polyphen_score = scalar_num(row$polyphen_score %||% NA_real_),
    sift = scalar_chr(row$sift %||% "", default = NA_character_),
    cadd = scalar_num(row$cadd %||% NA_real_),
    spliceai_max = scalar_num(row$spliceai_max %||% NA_real_),
    alphamissense = scalar_num(row$alphamissense_score %||% NA_real_),
    thresholds = thresholds
  )
  if (is_synonymous_consequence(variant_consequence) &&
      isTRUE(bp7$BP7) && !isTRUE(ins$PP3) && !isTRUE(ins$BP4)) {
    ins$BP4 <- TRUE
    ins$BP4_rationale <- paste(
      "Synonymous variant; in silico predictors do not support a damaging effect (BP4)."
    )
  }
  clin <- score_clinvar_criteria(
    scalar_chr(row$clinvar_classification %||% row$ClinVar %||% "", default = "")
  )

  clinical_ped <- if (identical(evidence_scope, "full")) {
    score_clinical_pedigree_criteria(
      gene = variant_gene,
      clinical_context = clinical_context,
      pedigree_context = pedigree_context,
      sample_genotypes = row$sample_genotypes %||% "{}",
      clinvar_text = scalar_chr(row$clinvar_classification %||% row$ClinVar %||% "", default = "")
    )
  } else {
    empty_group_a_context_scores()
  }

  scores <- c(
    list(
      max_population_af = af,
      popmax_af = af_info$popmax,
      global_population_af = af_info$global,
      af_known = af_info$af_known,
      PS_count = 0L, BS_count = 0L,
      PM_manual = 0L, PP_manual = 0L, BP_manual = 0L,
      PS3 = FALSE, PM3 = FALSE, PM6 = FALSE, BS3 = FALSE, BS4 = FALSE, BP2 = FALSE, BP5 = FALSE,
      PS3_rationale = "", PM3_rationale = "", PM6_rationale = "",
      BS3_rationale = "", BS4_rationale = "", BP2_rationale = "", BP5_rationale = ""
    ),
    pop, cons, pathogenic_protein, ps4, gene_rules, bp7, ins, clin, clinical_ped
  )
  scores <- normalize_criterion_flags(scores)

  evidence <- criteria_to_evidence(scores)
  scores$classification <- combine_acmg_evidence(evidence)
  scores$criteria_met <- collect_met_criteria(scores)
  scores$criteria_rationale <- collect_criteria_rationale(scores)

  evidence_tbl <- build_variant_evidence_table(scores, thresholds, criteria_meta, evidence_scope = evidence_scope)
  conf <- compute_confidence_score(scores, thresholds)
  scores$evidence_json <- evidence_table_to_json(evidence_tbl, evidence_scope = evidence_scope)
  scores$evidence_summary <- format_evidence_summary_text(evidence_tbl)
  scores$prediction_scores <- prediction_scores_to_text(row)
  scores$confidence_score <- conf$confidence_score
  scores$confidence_label <- conf$confidence_label
  strength <- compute_evidence_strength(scores)
  scores$evidence_strength <- strength$strength
  scores$pathogenic_evidence_count <- strength$path_count
  scores$benign_evidence_count <- strength$benign_count
  scores$prediction_limitations <- build_prediction_limitations(scores, row)
  scores$disease_profile <- profile_id
  scores
}

apply_manual_evidence_overlay <- function(scores, manual_inputs = list(), thresholds = NULL,
                                          criteria_meta = list(), evidence_scope = c("triggered", "automated", "full")) {
  evidence_scope <- match.arg(evidence_scope)
  if (length(manual_inputs) == 0L) {
    if (evidence_scope %in% c("full", "automated")) {
      thresholds <- thresholds %||% resolve_thresholds(NULL, scores$disease_profile %||% DEFAULT_PROFILE_ID)
      if (length(criteria_meta) == 0L) {
        criteria_meta <- load_rule_config(scores$disease_profile %||% DEFAULT_PROFILE_ID)$criteria_meta
      }
      evidence_tbl <- build_variant_evidence_table(scores, thresholds, criteria_meta, evidence_scope = evidence_scope)
      scores$evidence_json <- evidence_table_to_json(evidence_tbl, evidence_scope = evidence_scope)
    }
    return(scores)
  }
  if (isTRUE(manual_inputs$PS2_de_novo) && !isTRUE(scores$PS2)) {
    scores$PS2 <- TRUE
    scores$PS2_rationale <- "Curator confirmed de novo occurrence (manual PS2)."
  }
  if (isTRUE(manual_inputs$PS3_functional)) {
    scores$PS3 <- TRUE
    scores$PS3_rationale <- "Curator confirmed functional study supports damaging effect (manual PS3)."
  }
  if (isTRUE(manual_inputs$PS4_case_control) && !isTRUE(scores$PS4)) {
    scores$PS4 <- TRUE
    scores$PS4_rationale <- "Curator confirmed case-control enrichment (manual PS4)."
  }
  if (isTRUE(manual_inputs$PP1_segregation) && !isTRUE(scores$PP1)) {
    scores$PP1 <- TRUE
    scores$PP1_rationale <- "Curator confirmed co-segregation (manual PP1)."
  }
  if (isTRUE(manual_inputs$PP4_phenotype) && !isTRUE(scores$PP4)) {
    scores$PP4 <- TRUE
    scores$PP4_rationale <- "Curator confirmed phenotype match (manual PP4)."
  }
  if (isTRUE(manual_inputs$PM6_de_novo)) {
    scores$PM6 <- TRUE
    scores$PM6_rationale <- "Curator confirmed assumed de novo without parental confirmation (manual PM6)."
  }
  if (isTRUE(manual_inputs$PP2_missense_mechanism) && !isTRUE(scores$PP2)) {
    scores$PP2 <- TRUE
    scores$PP2_rationale <- paste(scores$PP2_rationale, "Curator confirmed missense mechanism.")
  }
  scores <- normalize_criterion_flags(scores)
  evidence <- criteria_to_evidence(scores)
  scores$classification <- combine_acmg_evidence(evidence)
  scores$criteria_met <- collect_met_criteria(scores)
  scores$criteria_rationale <- collect_criteria_rationale(scores)
  thresholds <- thresholds %||% resolve_thresholds(NULL, scores$disease_profile %||% DEFAULT_PROFILE_ID)
  if (length(criteria_meta) == 0L) {
    criteria_meta <- load_rule_config(scores$disease_profile %||% DEFAULT_PROFILE_ID)$criteria_meta
  }
  evidence_tbl <- build_variant_evidence_table(scores, thresholds, criteria_meta, evidence_scope = evidence_scope)
  conf <- compute_confidence_score(scores, thresholds)
  scores$evidence_json <- evidence_table_to_json(evidence_tbl, evidence_scope = evidence_scope)
  scores$evidence_summary <- format_evidence_summary_text(evidence_tbl)
  scores$confidence_score <- conf$confidence_score
  scores$confidence_label <- conf$confidence_label
  strength <- compute_evidence_strength(scores)
  scores$evidence_strength <- strength$strength
  scores$pathogenic_evidence_count <- strength$path_count
  scores$benign_evidence_count <- strength$benign_count
  scores$prediction_limitations <- build_prediction_limitations(scores)
  scores
}

score_variants_table <- function(variants_df, lof_genes = NULL, pp2_genes = NULL, bp1_genes = NULL,
                                 thresholds = NULL,
                                 profile_id = DEFAULT_PROFILE_ID,
                                 manual_inputs = list(),
                                 manual_by_variant = list(),
                                 clinical_context = NULL,
                                 pedigree_context = NULL,
                                 evidence_scope = c("triggered", "automated", "full"),
                                 lof_panel_path = file.path("data", "gene_panels", "lof_disease_mechanism_genes.csv"),
                                 pp2_panel_path = file.path("data", "gene_panels", "pp2_missense_mechanism_genes.csv"),
                                 bp1_panel_path = file.path("data", "gene_panels", "bp1_truncating_mechanism_genes.csv"),
                                 pm1_panel_path = PM1_HOTSPOT_PANEL_PATH,
                                 clinvar_protein_db_path = CLINVAR_PROTEIN_DB_PATH,
                                 ps4_db_path = PS4_CASE_CONTROL_DB_PATH,
                                 refs = NULL) {
  if (is.null(lof_genes)) lof_genes <- load_lof_gene_panel(lof_panel_path)
  if (is.null(pp2_genes)) pp2_genes <- load_gene_panel(pp2_panel_path)
  if (is.null(bp1_genes)) bp1_genes <- load_gene_panel(bp1_panel_path)
  pm1_genes <- load_pm1_hotspot_genes(pm1_panel_path)
  clinvar_protein_db <- load_clinvar_protein_db(clinvar_protein_db_path)
  ps4_db <- load_ps4_case_control_db(ps4_db_path)
  pm1_domains <- load_pm1_critical_domains()

  if (is.null(refs)) {
    refs <- tryCatch(load_reference_data(), error = function(e) NULL)
  }
  if (!is.null(refs)) {
    variants_df <- dedupe_variants_by_key(annotate_variants(variants_df, refs))
  }
  thresholds <- resolve_thresholds(thresholds, profile_id)
  criteria_meta <- load_rule_config(profile_id)$criteria_meta

  evidence_scope <- match.arg(evidence_scope)
  if (identical(evidence_scope, "full") &&
      is.null(clinical_context) && is.null(pedigree_context) &&
      length(manual_inputs) == 0L) {
    evidence_scope <- "automated"
  }

  scored <- lapply(seq_len(nrow(variants_df)), function(i) {
    tryCatch({
      s <- score_variant_row(
        variants_df[i, , drop = FALSE],
        lof_genes = lof_genes, pp2_genes = pp2_genes, bp1_genes = bp1_genes,
        pm1_genes = pm1_genes, pm1_domains = pm1_domains,
        clinvar_protein_db = clinvar_protein_db, ps4_db = ps4_db,
        clinical_context = clinical_context, pedigree_context = pedigree_context,
        thresholds = thresholds, profile_id = profile_id, criteria_meta = criteria_meta,
        evidence_scope = evidence_scope
      )
      vid <- variants_df$variant_id[i] %||% paste(
        variants_df$chrom[i], variants_df$pos[i], variants_df$ref[i], variants_df$alt[i], sep = ":"
      )
      apply_manual_evidence_overlay(
        s,
        if (length(manual_by_variant) > 0L) {
          manual_inputs_for_variant(manual_by_variant, vid)
        } else {
          normalize_manual_inputs(manual_inputs)
        },
        thresholds = thresholds, criteria_meta = criteria_meta,
        evidence_scope = evidence_scope
      )
    }, error = function(e) {
      vid <- variants_df$variant_id[i] %||% paste(
        variants_df$chrom[i], variants_df$pos[i], variants_df$ref[i], variants_df$alt[i], sep = ":"
      )
      stop(sprintf("Scoring failed for variant %s (chunk row %d): %s", vid, i, conditionMessage(e)), call. = FALSE)
    })
  })

  criteria_cols <- c(
    "max_population_af", "popmax_af", "global_population_af", "af_known",
    "BA1", "BS1", "BS2", "BS3", "BS4", "PM2",
    "PVS1", "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
    "gwas_supplementary", "gwas_supplementary_note",
    "classification", "criteria_met", "criteria_rationale",
    "evidence_json", "evidence_summary", "prediction_scores",
    "confidence_score", "confidence_label", "evidence_strength",
    "pathogenic_evidence_count", "benign_evidence_count", "prediction_limitations",
    "disease_profile",
    "insilico_summary",
    paste0(FULL_ACMG_CRITERIA, "_rationale")
  )

  criteria_df <- do.call(rbind, lapply(scored, function(x) {
    padded <- lapply(criteria_cols, function(col) {
      val <- x[[col]]
      if (is.null(val) || length(val) == 0L) {
        return(if (grepl("_rationale$", col)) "" else FALSE)
      }
      if (grepl("_rationale$", col) || col %in% c(
        "criteria_met", "criteria_rationale", "evidence_json", "evidence_summary",
        "prediction_scores", "insilico_summary", "classification", "confidence_label",
        "evidence_strength", "prediction_limitations", "disease_profile", "gwas_supplementary_note"
      )) {
        scalar_chr(val, default = "")
      } else if (col %in% c("max_population_af", "popmax_af", "global_population_af", "confidence_score", "pathogenic_evidence_count", "benign_evidence_count")) {
        scalar_num(val)
      } else if (col == "af_known") {
        scalar_lgl(val)
      } else {
        scalar_lgl(val)
      }
    })
    names(padded) <- criteria_cols
    as.data.frame(padded, stringsAsFactors = FALSE)
  }))
  rownames(criteria_df) <- NULL

  cbind(variants_df, criteria_df)
}
