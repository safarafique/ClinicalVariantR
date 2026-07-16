#' PS1, PM5, and PM1 criteria using ClinVar protein-level reference data.
#'
#' InterVar-style logic:
#'   PS1 - same amino acid change as established pathogenic variant
#'   PM5 - missense at residue with different established pathogenic missense
#'   PM1 - missense in curated critical domain / hotspot gene panel

CLINVAR_PROTEIN_DB_PATH <- file.path("data", "reference", "clinvar_pathogenic_protein.tsv")
PM1_HOTSPOT_PANEL_PATH <- file.path("data", "gene_panels", "pm1_hotspot_genes.csv")
PM1_CRITICAL_DOMAINS_PATH <- file.path("data", "reference", "pm1_critical_domains.tsv")

.clinvar_protein_cache <- new.env(parent = emptyenv())
.clinvar_protein_cache$db <- NULL
.pm1_domains_cache <- new.env(parent = emptyenv())
.pm1_domains_cache$db <- NULL

normalize_aa_letter <- function(aa) {
  if (is.na(aa) || !nzchar(aa)) return(NA_character_)
  aa <- toupper(trimws(aa))
  if (aa %in% c("*", "X", "TER", "STOP")) return("*")
  aa
}

normalize_hgvs_protein <- function(hgvs_p) {
  hgvs_p <- scalar_chr(hgvs_p, default = "")
  if (!nzchar(hgvs_p)) return(NA_character_)
  x <- trimws(as.character(hgvs_p))
  x <- sub("^p\\.", "", x, ignore.case = TRUE)
  x <- gsub("\\*", "Ter", x)
  x <- gsub("Ter", "*", x)
  x
}

parse_protein_change <- function(hgvs_p = NA_character_, amino_acids = NA_character_,
                                 protein_position = NA_character_) {
  hgvs_p <- scalar_chr(hgvs_p, default = "")
  amino_acids <- scalar_chr(amino_acids, default = "")
  protein_position <- scalar_chr(protein_position, default = "")
  out <- list(
    protein_position = NA_integer_,
    ref_aa = NA_character_,
    alt_aa = NA_character_,
    hgvs_p_norm = normalize_hgvs_protein(hgvs_p)
  )

  if (!is.na(out$hgvs_p_norm) && grepl("^([A-Za-z\\*]+)(\\d+)([A-Za-z\\*]+)$", out$hgvs_p_norm, perl = TRUE)) {
    m <- regmatches(out$hgvs_p_norm, regexec("^([A-Za-z\\*]+)(\\d+)([A-Za-z\\*]+)$", out$hgvs_p_norm, perl = TRUE))[[1L]]
    out$ref_aa <- normalize_aa_letter(m[[2L]])
    out$protein_position <- as.integer(m[[3L]])
    out$alt_aa <- normalize_aa_letter(m[[4L]])
    return(out)
  }

  if (!is.na(amino_acids) && grepl("/", amino_acids, fixed = TRUE)) {
    parts <- strsplit(amino_acids, "/", fixed = TRUE)[[1L]]
    if (length(parts) >= 2L) {
      out$ref_aa <- normalize_aa_letter(parts[[1L]])
      out$alt_aa <- normalize_aa_letter(parts[[2L]])
    }
  }

  if (!is.na(protein_position) && nzchar(protein_position)) {
    pos <- scalar_int(sub("/.*$", "", protein_position))
    if (!is.na(pos)) out$protein_position <- pos
  }

  out
}

is_pathogenic_significance <- function(x) {
  x <- scalar_chr(x, default = "")
  if (!nzchar(x)) return(FALSE)
  grepl("pathogenic", x, ignore.case = TRUE) &&
    !grepl("conflict|benign|uncertain|vus", x, ignore.case = TRUE)
}

load_clinvar_protein_db <- function(path = CLINVAR_PROTEIN_DB_PATH) {
  if (!is.null(.clinvar_protein_cache$db) && identical(attr(.clinvar_protein_cache$db, "path"), path)) {
    return(.clinvar_protein_cache$db)
  }
  if (!file.exists(path)) {
    db <- data.frame(
      gene = character(), protein_position = integer(),
      ref_aa = character(), alt_aa = character(),
      hgvs_p = character(), clinical_significance = character(),
      stringsAsFactors = FALSE
    )
    attr(db, "path") <- path
    .clinvar_protein_cache$db <- db
    return(db)
  }

  db <- utils::read.delim(path, stringsAsFactors = FALSE, comment.char = "#")
  required <- c("gene", "protein_position", "ref_aa", "alt_aa", "clinical_significance")
  if (!all(required %in% names(db))) {
    stop("ClinVar protein DB missing columns: ", paste(setdiff(required, names(db)), collapse = ", "))
  }
  db$gene <- toupper(trimws(db$gene))
  db$protein_position <- as.integer(db$protein_position)
  db$ref_aa <- vapply(db$ref_aa, normalize_aa_letter, character(1L))
  db$alt_aa <- vapply(db$alt_aa, normalize_aa_letter, character(1L))
  if (!"hgvs_p" %in% names(db)) db$hgvs_p <- NA_character_
  db <- db[is_pathogenic_significance(db$clinical_significance), , drop = FALSE]
  attr(db, "path") <- path
  .clinvar_protein_cache$db <- db
  db
}

load_pm1_hotspot_genes <- function(path = PM1_HOTSPOT_PANEL_PATH) {
  if (!file.exists(path)) return(character())
  genes <- unique(toupper(trimws(utils::read.csv(path, stringsAsFactors = FALSE)$gene)))
  genes[nzchar(genes)]
}

load_pm1_critical_domains <- function(path = PM1_CRITICAL_DOMAINS_PATH) {
  if (!is.null(.pm1_domains_cache$db) && identical(attr(.pm1_domains_cache$db, "path"), path)) {
    return(.pm1_domains_cache$db)
  }
  if (!file.exists(path)) {
    empty <- data.frame(gene = character(), aa_start = integer(), aa_end = integer(), stringsAsFactors = FALSE)
    attr(empty, "path") <- path
    .pm1_domains_cache$db <- empty
    return(empty)
  }
  df <- utils::read.delim(path, stringsAsFactors = FALSE, comment.char = "#")
  df$gene <- toupper(trimws(df$gene))
  df$aa_start <- as.integer(df$aa_start)
  df$aa_end <- as.integer(df$aa_end)
  attr(df, "path") <- path
  .pm1_domains_cache$db <- df
  df
}

pm1_in_critical_domain <- function(gene, protein_position, domains = NULL) {
  if (is.na(protein_position) || is.na(gene) || !nzchar(gene)) return(FALSE)
  if (is.null(domains)) domains <- load_pm1_critical_domains()
  if (nrow(domains) == 0L) return(FALSE)
  pos <- as.integer(protein_position)
  hits <- domains[toupper(domains$gene) == toupper(gene), , drop = FALSE]
  if (nrow(hits) == 0L) return(FALSE)
  any(isTRUE(!is.na(hits$aa_start) & !is.na(hits$aa_end) & pos >= hits$aa_start & pos <= hits$aa_end))
}

score_ps1_pm5_pm1_criteria <- function(
    gene = NA_character_,
    consequence = NA_character_,
    hgvs_p = NA_character_,
    amino_acids = NA_character_,
    protein_position = NA_character_,
    clinvar_db = NULL,
    csq_catalog = NULL,
    pm1_genes = character(),
    pm1_domains = NULL) {

  out <- list(
    PS1 = FALSE, PM1 = FALSE, PM5 = FALSE,
    PS1_rationale = "", PM1_rationale = "", PM5_rationale = ""
  )

  gene <- scalar_chr(gene, default = "")
  consequence <- scalar_chr(consequence, default = "")

  if (is.null(clinvar_db)) clinvar_db <- load_clinvar_protein_db()
  if (!is.null(csq_catalog) && nrow(csq_catalog) > 0L) {
    clinvar_db <- rbind(clinvar_db, csq_catalog)
    clinvar_db <- clinvar_db[!duplicated(
      paste(clinvar_db$gene, clinvar_db$protein_position, clinvar_db$alt_aa, sep = ":")
    ), , drop = FALSE]
  }

  change <- parse_protein_change(hgvs_p, amino_acids, protein_position)
  gene_key <- toupper(trimws(gene))

  if (is_missense_consequence(consequence) && nzchar(gene_key)) {
    in_hotspot <- isTRUE(gene_key %in% toupper(pm1_genes))
    in_domain <- isTRUE(pm1_in_critical_domain(gene_key, change$protein_position, pm1_domains))
    if (isTRUE(in_hotspot || in_domain)) {
      out$PM1 <- TRUE
      out$PM1_rationale <- if (isTRUE(in_domain)) {
        sprintf(
          "Missense at amino acid %s in critical domain of %s (PM1).",
          change$protein_position %||% "?", gene
        )
      } else {
        sprintf("Missense in %s; gene is in curated PM1 hotspot panel.", gene)
      }
    }
  }

  if (!is_missense_consequence(consequence) || !nzchar(gene_key) ||
      is.na(change$protein_position) || is.na(change$alt_aa)) {
    return(out)
  }

  hits <- clinvar_db[
    clinvar_db$gene == gene_key &
      clinvar_db$protein_position == change$protein_position,
    , drop = FALSE
  ]
  if (nrow(hits) == 0L) return(out)

  same_protein <- hits[
    !is.na(hits$alt_aa) & hits$alt_aa == change$alt_aa,
    , drop = FALSE
  ]
  if (nrow(same_protein) > 0L) {
    ref_hgvs <- same_protein$hgvs_p[[1L]]
    out$PS1 <- TRUE
    out$PS1_rationale <- sprintf(
      "Same amino acid change (%s) as established pathogenic variant in %s (PS1). Reference: %s.",
      change$hgvs_p_norm %||% paste0(change$ref_aa, change$protein_position, change$alt_aa),
      gene, ref_hgvs %||% "ClinVar protein DB"
    )
    return(out)
  }

  diff_protein <- hits[
    !is.na(hits$alt_aa) & hits$alt_aa != change$alt_aa,
    , drop = FALSE
  ]
  if (nrow(diff_protein) > 0L) {
    ref_hgvs <- diff_protein$hgvs_p[[1L]]
    out$PM5 <- TRUE
    out$PM5_rationale <- sprintf(
      "Missense at residue %s (%s) where a different pathogenic missense is known (PM5). Example: %s.",
      change$protein_position, gene, ref_hgvs %||% paste0(
        diff_protein$ref_aa[[1L]], change$protein_position, diff_protein$alt_aa[[1L]]
      )
    )
  }

  out
}
