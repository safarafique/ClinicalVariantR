#' PS4 - variant-specific case-control enrichment (ACMG Strong).
#'
#' Professional germline interpretation policy:
#'   - PS4 is awarded ONLY from curated case-control / cohort frequency data
#'     demonstrating enrichment in affected vs controls.
#'   - GWAS Catalog associations are NEVER auto-awarded as PS4.
#'     They are stored separately as supplementary evidence for curator review.

PS4_CASE_CONTROL_DB_PATH <- file.path("data", "reference", "ps4_case_control_enrichment.tsv")
GWAS_SUPPLEMENTARY_DB_PATH <- file.path("data", "reference", "gwas_supplementary_evidence.tsv")

.ps4_case_control_db <- NULL
.gwas_supplementary_db <- NULL

DEFAULT_PS4_THRESHOLDS <- list(
  ps4_min_case_af = 0.001,
  ps4_min_enrichment_ratio = 5.0,
  ps4_max_control_af = 0.01
)

normalize_rsid <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  x <- toupper(trimws(as.character(x)))
  sub("^RS", "rs", x)
}

split_rsid_field <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  parts <- unlist(strsplit(x, "[,;&|]", perl = TRUE))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  parts <- vapply(parts, normalize_rsid, character(1L), USE.NAMES = FALSE)
  parts <- parts[grepl("^rs[0-9]+$", parts, ignore.case = TRUE)]
  unique(parts)
}

normalize_ps4_source <- function(x) {
  if (is.na(x) || !nzchar(x)) return("case_control")
  gsub("[^a-z0-9]+", "_", tolower(trimws(x)))
}

load_ps4_case_control_db <- function(path = PS4_CASE_CONTROL_DB_PATH) {
  if (!is.null(.ps4_case_control_db) && identical(attr(.ps4_case_control_db, "path"), path)) {
    return(.ps4_case_control_db)
  }
  empty <- data.frame(
    variant_key = character(), chrom = character(), pos = integer(),
    ref = character(), alt = character(), gene = character(), rsid = character(),
    hgvs_p = character(), case_af = numeric(), control_af = numeric(),
    disease = character(), source = character(), ps4_confirmed = logical(),
    stringsAsFactors = FALSE
  )
  if (!file.exists(path)) {
    attr(empty, "path") <- path
    .ps4_case_control_db <<- empty
    return(empty)
  }

  db <- utils::read.delim(path, stringsAsFactors = FALSE, comment.char = "#")
  if (nrow(db) == 0L) {
    attr(empty, "path") <- path
    .ps4_case_control_db <<- empty
    return(empty)
  }

  if ("chrom" %in% names(db)) db$chrom <- normalize_chrom(db$chrom)
  if (!"variant_key" %in% names(db) && all(c("chrom", "pos", "ref", "alt") %in% names(db))) {
    db$variant_key <- variant_key_chr_pos_ref_alt(db$chrom, db$pos, db$ref, db$alt)
  } else if ("variant_key" %in% names(db) && all(c("chrom", "pos", "ref", "alt") %in% names(db))) {
    db$variant_key <- variant_key_chr_pos_ref_alt(db$chrom, db$pos, db$ref, db$alt)
  }
  if (!"variant_key" %in% names(db)) db$variant_key <- NA_character_
  for (col in c("rsid", "gene", "hgvs_p", "disease", "source")) {
    if (!col %in% names(db)) db[[col]] <- NA_character_
  }
  for (col in c("case_af", "control_af")) {
    if (!col %in% names(db)) db[[col]] <- NA_real_
  }
  if (!"ps4_confirmed" %in% names(db)) {
    db$ps4_confirmed <- FALSE
  } else {
    db$ps4_confirmed <- tolower(as.character(db$ps4_confirmed)) %in% c("true", "yes", "1")
  }

  db$source_norm <- vapply(db$source, normalize_ps4_source, character(1L))
  db <- db[!grepl("gwas", db$source_norm, fixed = TRUE), , drop = FALSE]

  db$rsid <- vapply(db$rsid, normalize_rsid, character(1L))
  db$gene <- toupper(trimws(db$gene))
  db$hgvs_p_norm <- vapply(db$hgvs_p, function(x) {
    if (is.na(x) || !nzchar(x)) return(NA_character_)
    gsub("\\s+", "", sub("^p\\.", "", trimws(x), ignore.case = TRUE))
  }, character(1L))

  attr(db, "path") <- path
  .ps4_case_control_db <<- db
  db
}

load_gwas_supplementary_db <- function(path = GWAS_SUPPLEMENTARY_DB_PATH) {
  if (!is.null(.gwas_supplementary_db) && identical(attr(.gwas_supplementary_db, "path"), path)) {
    return(.gwas_supplementary_db)
  }
  empty <- data.frame(
    variant_key = character(), chrom = character(), pos = integer(),
    rsid = character(), gene = character(), trait = character(),
    pvalue = numeric(), odds_ratio = numeric(), maf = numeric(),
    source = character(), stringsAsFactors = FALSE
  )
  if (!file.exists(path)) {
    attr(empty, "path") <- path
    .gwas_supplementary_db <<- empty
    return(empty)
  }

  db <- utils::read.delim(path, stringsAsFactors = FALSE, comment.char = "#")
  if (nrow(db) == 0L) {
    attr(empty, "path") <- path
    .gwas_supplementary_db <<- empty
    return(empty)
  }

  if ("chrom" %in% names(db)) db$chrom <- normalize_chrom(db$chrom)
  if (!"variant_key" %in% names(db) && all(c("chrom", "pos") %in% names(db))) {
    ref_col <- if ("ref" %in% names(db)) db$ref else "."
    alt_col <- if ("alt" %in% names(db)) db$alt else "."
    db$variant_key <- variant_key_chr_pos_ref_alt(db$chrom, db$pos, ref_col, alt_col)
  } else if ("variant_key" %in% names(db) && all(c("chrom", "pos") %in% names(db))) {
    ref_col <- if ("ref" %in% names(db)) db$ref else "."
    alt_col <- if ("alt" %in% names(db)) db$alt else "."
    db$variant_key <- variant_key_chr_pos_ref_alt(db$chrom, db$pos, ref_col, alt_col)
  }
  if (!"rsid" %in% names(db)) db$rsid <- NA_character_
  if (!"gene" %in% names(db)) db$gene <- NA_character_
  if (!"trait" %in% names(db)) db$trait <- NA_character_
  if (!"pvalue" %in% names(db)) db$pvalue <- NA_real_
  if (!"odds_ratio" %in% names(db)) db$odds_ratio <- NA_real_
  if (!"maf" %in% names(db)) db$maf <- NA_real_
  if (!"source" %in% names(db)) db$source <- "GWAS_Catalog"

  db$rsid <- vapply(db$rsid, normalize_rsid, character(1L))
  db$gene <- toupper(trimws(db$gene))

  attr(db, "path") <- path
  .gwas_supplementary_db <<- db
  db
}

ps4_enrichment_passes <- function(case_af, control_af, ps4_confirmed = FALSE,
                                  thresholds = DEFAULT_PS4_THRESHOLDS) {
  if (isTRUE(ps4_confirmed)) return(TRUE)
  if (is.na(case_af)) return(FALSE)
  control_af <- if (is.na(control_af)) 0 else control_af
  min_case <- thresholds$ps4_min_case_af %||% 0.001
  max_control <- thresholds$ps4_max_control_af %||% 0.01
  ratio_min <- thresholds$ps4_min_enrichment_ratio %||% 5.0
  if (case_af < min_case) return(FALSE)
  if (control_af <= 0) return(case_af >= min_case)
  (case_af / control_af) >= ratio_min && control_af <= max_control
}

find_reference_hits <- function(row, db) {
  if (nrow(db) == 0L) return(db[0, , drop = FALSE])

  vkey <- variant_key_chr_pos_ref_alt(row$chrom, row$pos, row$ref, row$alt)
  if ("variant_key" %in% names(db)) {
    hits <- db[db$variant_key == vkey, , drop = FALSE]
    if (nrow(hits) > 0L) return(hits)
  }

  rsids <- split_rsid_field(row$rsids %||% row$rsid %||% "")
  if (length(rsids) > 0L && "rsid" %in% names(db)) {
    hits <- db[db$rsid %in% rsids, , drop = FALSE]
    if (nrow(hits) > 0L) return(hits)
  }

  gene_key <- toupper(trimws(scalar_chr(row$gene %||% "")))
  hgvs_norm <- {
    hgvs_p <- scalar_chr(row$hgvs_p %||% "", default = "")
    if (nzchar(hgvs_p)) {
      gsub("\\s+", "", sub("^p\\.", "", trimws(hgvs_p), ignore.case = TRUE))
    } else {
      NA_character_
    }
  }
  if (nzchar(gene_key) && !is.na(hgvs_norm) && nzchar(hgvs_norm) &&
      all(c("gene", "hgvs_p_norm") %in% names(db))) {
    hits <- db[db$gene == gene_key & db$hgvs_p_norm == hgvs_norm, , drop = FALSE]
    if (nrow(hits) > 0L) return(hits)
  }

  db[0, , drop = FALSE]
}

format_gwas_supplementary_note <- function(hit) {
  trait <- hit$trait[[1L]] %||% "associated trait"
  pval <- hit$pvalue[[1L]]
  or_val <- hit$odds_ratio[[1L]]
  src <- hit$source[[1L]] %||% "GWAS Catalog"
  parts <- c(
    sprintf("GWAS supplementary: %s", trait),
    if (!is.na(pval)) sprintf("p=%.2e", pval),
    if (!is.na(or_val)) sprintf("OR=%.2f", or_val),
    sprintf("source=%s", src)
  )
  paste0(
    paste(parts, collapse = "; "),
    ". Not applied as automated PS4 - requires variant-specific case-control enrichment review."
  )
}

lookup_gwas_supplementary_evidence <- function(row, gwas_db = NULL) {
  out <- list(
    gwas_supplementary = FALSE,
    gwas_supplementary_note = ""
  )
  if (is.null(gwas_db)) gwas_db <- load_gwas_supplementary_db()
  if (nrow(gwas_db) == 0L) return(out)

  hits <- find_reference_hits(row, gwas_db)
  if (nrow(hits) == 0L) return(out)

  out$gwas_supplementary <- TRUE
  out$gwas_supplementary_note <- format_gwas_supplementary_note(hits[1L, , drop = FALSE])
  out
}

score_ps4_criteria <- function(row, ps4_db = NULL, gwas_db = NULL, thresholds = NULL) {
  row <- normalize_variant_row_input(row)
  out <- list(
    PS4 = FALSE,
    PS4_rationale = "",
    gwas_supplementary = FALSE,
    gwas_supplementary_note = ""
  )

  gwas_hit <- lookup_gwas_supplementary_evidence(row, gwas_db = gwas_db)
  out$gwas_supplementary <- gwas_hit$gwas_supplementary
  out$gwas_supplementary_note <- gwas_hit$gwas_supplementary_note

  if (is.null(ps4_db)) ps4_db <- load_ps4_case_control_db()
  if (nrow(ps4_db) == 0L) {
    if (isTRUE(out$gwas_supplementary)) {
      out$PS4_rationale <- out$gwas_supplementary_note
    }
    return(out)
  }

  th <- modifyList(DEFAULT_PS4_THRESHOLDS, thresholds %||% list())
  hits <- find_reference_hits(row, ps4_db)
  if (nrow(hits) == 0L) {
    if (isTRUE(out$gwas_supplementary)) {
      out$PS4_rationale <- out$gwas_supplementary_note
    }
    return(out)
  }

  for (i in seq_len(nrow(hits))) {
    hit <- hits[i, , drop = FALSE]
    if (ps4_enrichment_passes(
      case_af = hit$case_af[[1L]],
      control_af = hit$control_af[[1L]],
      ps4_confirmed = isTRUE(hit$ps4_confirmed[[1L]]),
      thresholds = th
    )) {
      out$PS4 <- TRUE
      disease <- hit$disease[[1L]] %||% "reported disorder"
      source <- hit$source[[1L]] %||% "case-control reference"
      case_af <- hit$case_af[[1L]]
      control_af <- hit$control_af[[1L]]
      if (isTRUE(hit$ps4_confirmed[[1L]])) {
        out$PS4_rationale <- sprintf(
          "Variant-specific case-control enrichment confirmed for %s (PS4). Source: %s.",
          disease, source
        )
      } else {
        out$PS4_rationale <- sprintf(
          "Case AF = %.6f vs control AF = %.6f for %s (ratio >= %.1f). Variant-specific case-control evidence: PS4 (Strong). Source: %s.",
          case_af, control_af %||% 0, disease,
          th$ps4_min_enrichment_ratio %||% 5.0, source
        )
      }
      if (isTRUE(out$gwas_supplementary) && nzchar(out$gwas_supplementary_note)) {
        out$PS4_rationale <- paste(out$PS4_rationale, out$gwas_supplementary_note, sep = " | ")
      }
      return(out)
    }
  }

  if (isTRUE(out$gwas_supplementary)) {
    out$PS4_rationale <- out$gwas_supplementary_note
  }
  out
}
