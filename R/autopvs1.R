#' ClinGen AutoPVS1-style loss-of-function evaluation (simplified implementation).
#'
#' Reference: Abou Tayoun et al., Genetics in Medicine 2018; ClinGen SVI PVS1 specification.

CRITICAL_EXON_DOMAINS_PATH <- file.path("data", "reference", "critical_exon_domains.tsv")
.critical_exon_domains <- NULL

load_critical_exon_domains <- function(path = CRITICAL_EXON_DOMAINS_PATH) {
  if (!is.null(.critical_exon_domains) && identical(attr(.critical_exon_domains, "path"), path)) {
    return(.critical_exon_domains)
  }
  if (!file.exists(path)) {
    empty <- data.frame(gene = character(), start_exon = integer(), end_exon = integer(), stringsAsFactors = FALSE)
    attr(empty, "path") <- path
    .critical_exon_domains <<- empty
    return(empty)
  }
  df <- utils::read.delim(path, stringsAsFactors = FALSE, comment.char = "#")
  df$gene <- toupper(trimws(df$gene))
  df$start_exon <- as.integer(df$start_exon)
  df$end_exon <- as.integer(df$end_exon)
  attr(df, "path") <- path
  .critical_exon_domains <<- df
  df
}

NMD_ESCAPE_CONSEQUENCES <- c("stop_lost", "start_lost")

is_nmd_escape_consequence <- function(consequence) {
  cons <- normalize_consequence(consequence)
  any(vapply(NMD_ESCAPE_CONSEQUENCES, function(x) grepl(x, cons, fixed = TRUE), logical(1L)))
}

parse_exon_number <- function(exon_field) {
  exon_field <- scalar_chr(exon_field, default = "")
  if (!nzchar(exon_field)) return(NA_integer_)
  num <- suppressWarnings(as.integer(sub("/.*$", "", exon_field)))
  if (is.na(num)) suppressWarnings(as.integer(sub("^.*?/", "", exon_field)))
  else num
}

parse_coding_exon_count <- function(exon_field) {
  exon_field <- scalar_chr(exon_field, default = "")
  if (!grepl("/", exon_field, fixed = TRUE)) return(NA_integer_)
  suppressWarnings(as.integer(sub("^.*/", "", exon_field)))
}

is_critical_exon <- function(gene, exon_num, domains_df) {
  gene <- toupper(scalar_chr(gene, default = ""))
  if (!nzchar(gene) || is.na(exon_num) || nrow(domains_df) == 0L) return(FALSE)
  if (!all(c("gene", "start_exon", "end_exon") %in% names(domains_df))) return(FALSE)
  rows <- domains_df[toupper(domains_df$gene) == gene, , drop = FALSE]
  if (nrow(rows) == 0L) return(FALSE)
  any(vapply(seq_len(nrow(rows)), function(i) {
    exon_num >= rows$start_exon[i] && exon_num <= rows$end_exon[i]
  }, logical(1L)))
}

is_high_impact_lof <- function(consequence, impact = NA_character_) {
  if (!is_lof_consequence(consequence)) return(FALSE)
  impact_hi <- grepl("^HIGH$", toupper(scalar_chr(impact, default = "")))
  if (impact_hi) return(TRUE)
  cons <- normalize_consequence(consequence)
  grepl("stop_gained|frameshift|splice_donor|splice_acceptor|nonsense", cons)
}

lof_transcript_eligible <- function(
    is_protein_coding,
    is_canonical_transcript,
    consequence,
    impact = NA_character_) {
  if (isTRUE(is_protein_coding) && isTRUE(is_canonical_transcript)) return(TRUE)
  isTRUE(is_protein_coding) && is_high_impact_lof(consequence, impact)
}

score_autopvs1 <- function(
    consequence,
    gene,
    lof_genes = character(),
    is_protein_coding = FALSE,
    is_canonical_transcript = FALSE,
    impact = NA_character_,
    exon = NA_character_,
    biotype = NA_character_,
    domains_df = NULL) {

  out <- list(
    PVS1 = FALSE,
    PVS1_strength = NA_character_,
    PVS1_rationale = ""
  )

  if (!is_lof_consequence(consequence)) {
    out$PVS1_rationale <- "Not a predicted loss-of-function variant."
    return(out)
  }

  if (is_nmd_escape_consequence(consequence)) {
    out$PVS1_rationale <- "NMD-escape consequence (stop_loss/start_loss); PVS1 not applied (AutoPVS1)."
    return(out)
  }

  gene <- scalar_chr(gene, default = "")
  in_lof_panel <- nzchar(gene) && gene %in% lof_genes
  transcript_ok <- lof_transcript_eligible(
    is_protein_coding, is_canonical_transcript, consequence, impact
  )

  if (!in_lof_panel && is_prediction_mode()) {
    out$PVS1_rationale <- sprintf(
      "LoF in %s but gene not in LoF mechanism panel; PVS1 withheld (prediction mode).",
      ifelse(nzchar(gene), gene, "unknown gene")
    )
    return(out)
  }

  exon_num <- parse_exon_number(exon)
  exon_total <- parse_coding_exon_count(exon)
  if (is.null(domains_df)) domains_df <- load_critical_exon_domains()

  last_exon <- !is.na(exon_num) && !is.na(exon_total) && exon_num >= exon_total
  if (isTRUE(last_exon)) {
    if (in_lof_panel && grepl("stop_gained", normalize_consequence(consequence), fixed = TRUE)) {
      out$PVS1 <- TRUE
      out$PVS1_strength <- "Very Strong"
      out$PVS1_rationale <- sprintf(
        "AutoPVS1: stop_gained in terminal exon (%s/%s) of %s (LoF panel; NMD review recommended).",
        exon_num, exon_total, gene
      )
      return(out)
    }
    out$PVS1_rationale <- sprintf(
      "LoF in last exon (%s/%s); likely NMD escape — PVS1 not applied (AutoPVS1).",
      exon_num, exon_total
    )
    return(out)
  }

  critical <- is_critical_exon(gene, exon_num, domains_df)

  if (in_lof_panel && transcript_ok) {
    out$PVS1 <- TRUE
    out$PVS1_strength <- "Very Strong"
    canon_note <- if (isTRUE(is_canonical_transcript)) {
      "canonical transcript"
    } else {
      "protein-coding HIGH-impact LoF transcript (canonical fallback)"
    }
    out$PVS1_rationale <- sprintf(
      "AutoPVS1: LoF (%s) in %s with established LoF mechanism, %s%s.",
      consequence, gene, canon_note,
      if (critical) " in critical exon/domain" else ""
    )
    return(out)
  }

  if (transcript_ok && !is_prediction_mode()) {
    out$PVS1 <- TRUE
    out$PVS1_strength <- "Very Strong"
    out$PVS1_rationale <- sprintf(
      "LoF on protein-coding transcript in %s (PVS1; verify LoF mechanism).", gene
    )
    return(out)
  }

  out$PVS1_rationale <- "LoF variant fails AutoPVS1 checks (gene panel / transcript / NMD review)."
  out
}
