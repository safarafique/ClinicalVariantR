#' PP4, PP1, PS2, and Group A context criteria from clinical logs and pedigree.
#' @noRd

GENE_PHENOTYPE_MAP_PATH <- file.path("data", "gene_panels", "gene_phenotype_map.csv")

load_gene_phenotype_map <- function(path = GENE_PHENOTYPE_MAP_PATH) {
  if (!file.exists(path)) return(data.frame())
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  df$gene <- toupper(trimws(df$gene))
  df
}

phenotype_matches_gene <- function(gene, phenotype_text, map_df) {
  gene <- scalar_chr(gene, default = "")
  if (!nzchar(gene) || nrow(map_df) == 0L) return(FALSE)
  row <- map_df[toupper(map_df$gene) == toupper(gene), , drop = FALSE]
  if (nrow(row) == 0L) return(FALSE)
  pheno <- paste(phenotype_text, collapse = " ")
  if (!nzchar(pheno)) return(FALSE)
  keys <- unlist(strsplit(row$phenotype_keywords[[1L]], "|", fixed = TRUE))
  any(vapply(keys, function(k) grepl(k, pheno, ignore.case = TRUE), logical(1L)))
}

pedigree_summary <- function(pedigree_context) {
  if (is.null(pedigree_context) || nrow(pedigree_context) == 0L) {
    return(list(
      n_affected = 0L, has_parent = FALSE, parents_unaffected = FALSE,
      proband_affected = FALSE, n_members = 0L
    ))
  }
  rel <- tolower(as.character(pedigree_context$relation %||% ""))
  aff <- tolower(as.character(pedigree_context$affected_status %||% ""))
  n_affected <- sum(grepl("yes|affected|1|true", aff, ignore.case = TRUE))
  parent_mask <- grepl("mother|father|parent|maternal|paternal", rel, ignore.case = TRUE)
  has_parent <- any(parent_mask)
  parents_unaffected <- has_parent && all(
    !grepl("yes|affected|1|true", aff[parent_mask], ignore.case = TRUE)
  )
  proband_affected <- any(
    grepl("proband|patient|self|index", rel, ignore.case = TRUE) &
      grepl("yes|affected|1|true", aff, ignore.case = TRUE)
  )
  list(
    n_affected = n_affected,
    has_parent = has_parent,
    parents_unaffected = parents_unaffected,
    proband_affected = proband_affected,
    n_members = length(rel)
  )
}

score_clinical_pedigree_criteria <- function(
    gene = NA_character_,
    clinical_context = NULL,
    pedigree_context = NULL,
    gene_phenotype_map = NULL,
    sample_genotypes = NULL,
    clinvar_text = NA_character_) {

  gene <- scalar_chr(gene, default = "")
  out <- list(
    PP4 = FALSE, PP1 = FALSE, PS2 = FALSE,
    PM3 = FALSE, PM6 = FALSE, BS3 = FALSE, BS4 = FALSE, BP2 = FALSE, BP5 = FALSE,
    PP4_rationale = "", PP1_rationale = "", PS2_rationale = "",
    PM3_rationale = "", PM6_rationale = "", BS3_rationale = "",
    BS4_rationale = "", BP2_rationale = "", BP5_rationale = ""
  )

  if (is.null(gene_phenotype_map)) gene_phenotype_map <- load_gene_phenotype_map()

  if (!is.null(clinical_context) && nrow(clinical_context) > 0L) {
    hpo_pp4 <- score_hpo_omim_pp4(gene, clinical_context)
    if (isTRUE(hpo_pp4$PP4)) {
      out$PP4 <- TRUE
      out$PP4_rationale <- hpo_pp4$PP4_rationale
    }

    pheno_cols <- intersect(c("phenotype", "cml_phase", "tki_response"), names(clinical_context))
    pheno_text <- if (length(pheno_cols) > 0L) {
      paste(apply(clinical_context[, pheno_cols, drop = FALSE], 1L, paste, collapse = " "), collapse = "; ")
    } else {
      ""
    }
    if (!isTRUE(out$PP4) && phenotype_matches_gene(gene, pheno_text, gene_phenotype_map)) {
      out$PP4 <- TRUE
      out$PP4_rationale <- sprintf(
        "Patient phenotype (%s) matches curated gene-disease association for %s (PP4).",
        pheno_text, gene
      )
    } else if (!isTRUE(out$PP4) &&
               isTRUE(any(grepl("accelerated|blast|advanced", clinical_context$cml_phase %||% "", ignore.case = TRUE), na.rm = TRUE)) &&
               toupper(gene) %in% c("BCR", "ABL1", "RUNX1", "GATA2")) {
      out$PP4 <- TRUE
      out$PP4_rationale <- sprintf(
        "Advanced CML phase with variant in %s - phenotype consistent with hematologic disease context (PP4).",
        gene
      )
    }
  }

  ped <- if (!is.null(pedigree_context) && nrow(pedigree_context) > 0L) {
    pedigree_summary(pedigree_context)
  } else {
    NULL
  }

  if (!is.null(ped)) {
    if (ped$n_affected >= 2L) {
      out$PP1 <- TRUE
      out$PP1_rationale <- sprintf(
        "Pedigree shows %d affected member(s) - co-segregation with disease supported (PP1).", ped$n_affected
      )
    }
    if (ped$proband_affected && ped$parents_unaffected) {
      out$PS2 <- TRUE
      out$PS2_rationale <- "Proband affected with unaffected parents - de novo occurrence supported (PS2; confirm phasing)."
    }
    if (ped$proband_affected && !ped$has_parent) {
      out$PM6_rationale <- paste(
        "Parental genotypes not listed in pedigree upload;",
        "PM6 may apply if de novo is assumed without parental confirmation (curator review)."
      )
    }
    if (ped$n_members >= 2L && ped$n_affected == 1L && ped$proband_affected) {
      out$BS4_rationale <- paste(
        "Pedigree lists additional relatives but only proband affected;",
        "review BS4 if variant fails to segregate with disease in other affected relatives."
      )
    }
    if (ped$n_affected >= 2L) {
      out$PM3_rationale <- paste(
        "Pedigree suggests familial disease;",
        "PM3 requires phased trans/compound-heterozygous configuration with a pathogenic allele."
      )
      out$BP2_rationale <- paste(
        "Pedigree available;",
        "BP2 requires observation in trans with a pathogenic variant (phased genotypes not in upload)."
      )
    }
  }

  trio <- score_trio_genotype_criteria(
    sample_genotypes = sample_genotypes,
    pedigree_context = pedigree_context,
    gene = gene,
    clinvar_text = clinvar_text
  )
  if (isTRUE(trio$PS2) && !isTRUE(out$PS2)) {
    out$PS2 <- TRUE
    out$PS2_rationale <- trio$PS2_rationale
  } else if (isTRUE(trio$PS2) && isTRUE(out$PS2)) {
    out$PS2_rationale <- paste(out$PS2_rationale, trio$PS2_rationale, sep = " | ")
  }
  if (nzchar(trio$PM3_rationale %||% "") &&
      (grepl("review PM3|compound het", trio$PM3_rationale, ignore.case = TRUE))) {
    out$PM3_rationale <- trio$PM3_rationale
  }
  if (nzchar(trio$BP2_rationale %||% "") &&
      grepl("review BP2|trans", trio$BP2_rationale, ignore.case = TRUE)) {
    out$BP2_rationale <- trio$BP2_rationale
  }
  if (!nzchar(out$PM3_rationale)) out$PM3_rationale <- trio$PM3_rationale
  if (!nzchar(out$BP2_rationale)) out$BP2_rationale <- trio$BP2_rationale
  if (!nzchar(out$PM3_rationale)) {
    out$PM3_rationale <- paste(
      "PM3 requires variant in trans with a pathogenic allele;",
      "phased genotypes or second variant not inferred from VCF/clinical/pedigree uploads."
    )
  }
  if (!nzchar(out$BS3_rationale)) {
    out$BS3_rationale <- "BS3 requires well-established functional studies showing no damaging effect (curator/literature review)."
  }
  if (!nzchar(out$BS4_rationale)) {
    out$BS4_rationale <- "BS4 requires lack of segregation in affected family members; review pedigree segregation manually."
  }
  if (!nzchar(out$BP2_rationale)) {
    out$BP2_rationale <- "BP2 requires observation in trans with a pathogenic variant; trio phasing not available from VCF alone."
  }
  if (!nzchar(out$BP5_rationale)) {
    out$BP5_rationale <- "BP5 requires reputable independent benign report beyond ClinVar (BP6 covers ClinVar submissions)."
  }

  out
}

empty_group_a_context_scores <- function() {
  list(
    PP4 = FALSE, PP1 = FALSE, PS2 = FALSE,
    PM3 = FALSE, PM6 = FALSE, BS3 = FALSE, BS4 = FALSE, BP2 = FALSE, BP5 = FALSE,
    PP4_rationale = "", PP1_rationale = "", PS2_rationale = "",
    PM3_rationale = "", PM6_rationale = "", BS3_rationale = "",
    BS4_rationale = "", BP2_rationale = "", BP5_rationale = ""
  )
}
