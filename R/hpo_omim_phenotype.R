#' HPO and OMIM integration for PP4 and phenotype-driven prediction.

HPO_GENE_MAP_PATH <- file.path("data", "reference", "hpo_gene_associations.tsv")
OMIM_GENE_MAP_PATH <- file.path("data", "reference", "omim_gene_map.tsv")

.hpo_gene_db <- NULL
.omim_gene_db <- NULL

normalize_hpo_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^HP:", "HP:", x)
  if (grepl("^HP[0-9]+$", x)) {
    num <- sub("^HP", "", x)
    return(paste0("HP:", sprintf("%07d", as.integer(num))))
  }
  if (grepl("^[0-9]+$", x)) return(paste0("HP:", sprintf("%07d", as.integer(x))))
  x
}

load_hpo_gene_map <- function(path = HPO_GENE_MAP_PATH) {
  if (!is.null(.hpo_gene_db) && identical(attr(.hpo_gene_db, "path"), path)) {
    return(.hpo_gene_db)
  }
  if (!file.exists(path)) return(data.frame())
  df <- utils::read.delim(path, stringsAsFactors = FALSE, comment.char = "#")
  if (nrow(df) == 0L) return(df)
  df$hpo_id <- vapply(df$hpo_id, normalize_hpo_id, character(1L))
  df$gene <- toupper(trimws(df$gene))
  attr(df, "path") <- path
  .hpo_gene_db <<- df
  df
}

load_omim_gene_map <- function(path = OMIM_GENE_MAP_PATH) {
  if (!is.null(.omim_gene_db) && identical(attr(.omim_gene_db, "path"), path)) {
    return(.omim_gene_db)
  }
  if (!file.exists(path)) return(data.frame())
  df <- utils::read.delim(path, stringsAsFactors = FALSE, comment.char = "#")
  if (nrow(df) == 0L) return(df)
  df$gene <- toupper(trimws(df$gene))
  attr(df, "path") <- path
  .omim_gene_db <<- df
  df
}

extract_hpo_terms_from_clinical <- function(clinical_context) {
  if (is.null(clinical_context) || nrow(clinical_context) == 0L) return(character())
  terms <- character()
  if ("hpo_terms" %in% names(clinical_context)) {
    raw <- paste(clinical_context$hpo_terms, collapse = ";")
    parts <- unlist(strsplit(raw, "[;,|]+"))
    parts <- trimws(parts)
    terms <- c(terms, parts[nzchar(parts)])
  }
  if ("hpo_ids" %in% names(clinical_context)) {
    raw <- paste(clinical_context$hpo_ids, collapse = ";")
    parts <- unlist(strsplit(raw, "[;,|]+"))
    parts <- trimws(parts)
    terms <- c(terms, parts[nzchar(parts)])
  }
  unique(normalize_hpo_id(terms[nzchar(terms)]))
}

extract_phenotype_text <- function(clinical_context) {
  if (is.null(clinical_context) || nrow(clinical_context) == 0L) return("")
  cols <- intersect(c("phenotype", "cml_phase", "tki_response", "disease", "clinical_summary"), names(clinical_context))
  if (length(cols) == 0L) return("")
  paste(apply(clinical_context[, cols, drop = FALSE], 1L, paste, collapse = " "), collapse = "; ")
}

score_hpo_omim_pp4 <- function(gene, clinical_context = NULL,
                               hpo_map = NULL, omim_map = NULL) {
  gene <- toupper(scalar_chr(gene, default = ""))
  out <- list(PP4 = FALSE, PP4_rationale = "")
  if (!nzchar(gene)) return(out)

  if (is.null(hpo_map)) hpo_map <- load_hpo_gene_map()
  if (is.null(omim_map)) omim_map <- load_omim_gene_map()

  hpo_terms <- extract_hpo_terms_from_clinical(clinical_context)
  if (length(hpo_terms) > 0L && nrow(hpo_map) > 0L) {
    hits <- hpo_map[hpo_map$hpo_id %in% hpo_terms & hpo_map$gene == gene, , drop = FALSE]
    if (nrow(hits) > 0L) {
      out$PP4 <- TRUE
      out$PP4_rationale <- sprintf(
        "HPO term(s) %s match curated gene association for %s (%s).",
        paste(unique(hits$hpo_id), collapse = ", "),
        gene,
        paste(unique(hits$hpo_label), collapse = "; ")
      )
      return(out)
    }
  }

  pheno <- extract_phenotype_text(clinical_context)
  if (nzchar(pheno) && nrow(omim_map) > 0L) {
    om <- omim_map[omim_map$gene == gene, , drop = FALSE]
    if (nrow(om) > 0L) {
      for (i in seq_len(nrow(om))) {
        disease <- om$disease_name[i]
        keys <- unlist(strsplit(tolower(disease), "[^a-z0-9]+"))
        keys <- keys[nchar(keys) >= 4L]
        if (any(vapply(keys, function(k) grepl(k, tolower(pheno), fixed = TRUE), logical(1L)))) {
          out$PP4 <- TRUE
          out$PP4_rationale <- sprintf(
            "Patient phenotype matches OMIM disease '%s' (OMIM:%s) for gene %s (PP4).",
            disease, om$omim_id[i], gene
          )
          return(out)
        }
      }
    }
  }

  out
}
