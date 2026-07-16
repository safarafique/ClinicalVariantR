#' Trio VCF genotype parsing for PS2, PM3, and BP2 prediction.

parse_vcf_genotypes <- function(format_str, sample_fields, sample_names) {
  if (length(sample_fields) == 0L || length(sample_names) == 0L) {
    return(list())
  }
  format_keys <- strsplit(format_str, ":", fixed = TRUE)[[1L]]
  gt_idx <- match("GT", format_keys)
  if (is.na(gt_idx)) return(list())

  out <- list()
  for (i in seq_along(sample_names)) {
    if (i > length(sample_fields)) break
    parts <- strsplit(sample_fields[i], ":", fixed = TRUE)[[1L]]
    if (length(parts) >= gt_idx) {
      gt <- parts[[gt_idx]]
      if (nzchar(gt) && gt != ".") out[[sample_names[i]]] <- gt
    }
  }
  out
}

normalize_gt <- function(gt) {
  gt <- as.character(gt %||% "")
  gt <- gsub("\\|", "/", gt)
  gt
}

gt_has_alt <- function(gt) {
  gt <- normalize_gt(gt)
  if (!nzchar(gt) || gt %in% c(".", "./.")) return(FALSE)
  alleles <- strsplit(gt, "/", fixed = TRUE)[[1L]]
  alleles <- vapply(alleles, scalar_int, integer(1L))
  any(!is.na(alleles) & alleles > 0L)
}

gt_is_hom_ref <- function(gt) {
  gt <- normalize_gt(gt)
  gt %in% c("0/0", "./.")
}

gt_is_het <- function(gt) {
  gt <- normalize_gt(gt)
  grepl("^0/1$|^1/0$", gt)
}

gt_is_hom_alt <- function(gt) {
  gt <- normalize_gt(gt)
  gt %in% c("1/1", "2/2")
}

map_pedigree_to_vcf_samples <- function(pedigree_context, vcf_sample_names) {
  if (is.null(pedigree_context) || nrow(pedigree_context) == 0L) {
    return(list(proband = NA_character_, mother = NA_character_, father = NA_character_))
  }
  vcf_sample_names <- as.character(vcf_sample_names)
  pick <- function(relation_pattern) {
    rel <- tolower(pedigree_context$relation %||% "")
    idx <- grep(relation_pattern, rel, ignore.case = TRUE)
    if (length(idx) == 0L) return(NA_character_)
    sid <- pedigree_context$sample_id[idx[1L]]
    if ("vcf_sample" %in% names(pedigree_context) && nzchar(pedigree_context$vcf_sample[idx[1L]] %||% "")) {
      return(as.character(pedigree_context$vcf_sample[idx[1L]]))
    }
    if (sid %in% vcf_sample_names) return(as.character(sid))
    hit <- vcf_sample_names[grepl(sid, vcf_sample_names, ignore.case = TRUE)]
    if (length(hit) > 0L) return(hit[[1L]])
    NA_character_
  }
  list(
    proband = pick("proband|patient|self|index"),
    mother = pick("mother|maternal"),
    father = pick("father|paternal")
  )
}

get_gt <- function(genotypes, sample_name) {
  if (is.null(genotypes) || length(genotypes) == 0L || is.na(sample_name)) return(NA_character_)
  genotypes[[sample_name]] %||% NA_character_
}

parse_sample_genotypes_field <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x) || x == "{}") return(list())
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(tryCatch(jsonlite::fromJSON(x), error = function(e) list()))
  }
  list()
}

score_trio_genotype_criteria <- function(
    sample_genotypes = NULL,
    pedigree_context = NULL,
    gene = NA_character_,
    clinvar_text = NA_character_) {

  out <- list(
    PS2 = FALSE, PM3 = FALSE, BP2 = FALSE,
    PS2_rationale = "", PM3_rationale = "", BP2_rationale = ""
  )

  if (is.character(sample_genotypes) && length(sample_genotypes) == 1L) {
    sample_genotypes <- parse_sample_genotypes_field(sample_genotypes)
  }
  if (is.null(sample_genotypes) || length(sample_genotypes) == 0L) {
    out$PS2_rationale <- "No FORMAT/GT data in VCF; PS2 requires trio genotypes."
    out$PM3_rationale <- "PM3 requires phased compound-heterozygous or trans configuration (GT + second variant)."
    out$BP2_rationale <- "BP2 requires trans observation with pathogenic variant (GT + phasing)."
    return(out)
  }

  vcf_names <- names(sample_genotypes)
  mapped <- map_pedigree_to_vcf_samples(pedigree_context, vcf_names)

  gt_prob <- get_gt(sample_genotypes, mapped$proband)
  gt_mom <- get_gt(sample_genotypes, mapped$mother)
  gt_dad <- get_gt(sample_genotypes, mapped$father)

  if (!is.na(gt_prob) && gt_has_alt(gt_prob)) {
    parents_present <- !is.na(gt_mom) || !is.na(gt_dad)
    parents_ref <- TRUE
    if (!is.na(gt_mom)) parents_ref <- parents_ref && gt_is_hom_ref(gt_mom)
    if (!is.na(gt_dad)) parents_ref <- parents_ref && gt_is_hom_ref(gt_dad)
    if (parents_present && parents_ref) {
      out$PS2 <- TRUE
      out$PS2_rationale <- sprintf(
        "Trio GT: proband %s with unaffected parents (%s/%s) - de novo supported (PS2; confirm phasing).",
        gt_prob,
        if (is.na(gt_mom)) "NA" else gt_mom,
        if (is.na(gt_dad)) "NA" else gt_dad
      )
    }
  }

  if (!is.na(gt_prob) && gt_is_het(gt_prob) && !is.na(gt_mom) && !is.na(gt_dad)) {
    if (gt_is_het(gt_mom) && gt_is_hom_ref(gt_dad)) {
      out$PM3_rationale <- sprintf(
        "Trio GT: proband %s, mother %s, father %s - review PM3/compound het with second pathogenic variant in %s.",
        gt_prob, gt_mom, gt_dad, gene
      )
    } else if (gt_is_het(gt_dad) && gt_is_hom_ref(gt_mom)) {
      out$PM3_rationale <- sprintf(
        "Trio GT: proband %s, mother %s, father %s - review PM3/compound het with second pathogenic variant in %s.",
        gt_prob, gt_mom, gt_dad, gene
      )
    }
  }

  if (!is.na(gt_prob) && gt_is_het(gt_prob) && clinvar_is_pathogenic(clinvar_text)) {
    if ((!is.na(gt_mom) && gt_is_het(gt_mom)) || (!is.na(gt_dad) && gt_is_het(gt_dad))) {
      out$BP2_rationale <- paste(
        "Trio GT with pathogenic ClinVar annotation;",
        "review BP2 if variant observed in trans with another pathogenic allele."
      )
    }
  }

  if (!nzchar(out$PM3_rationale)) {
    out$PM3_rationale <- "PM3 not met; requires trans pathogenic variant or compound het phasing."
  }
  if (!nzchar(out$BP2_rationale)) {
    out$BP2_rationale <- "BP2 not met from trio GT; requires trans pathogenic configuration."
  }

  out
}

serialize_sample_genotypes <- function(genotypes) {
  if (length(genotypes) == 0L) return("{}")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(jsonlite::toJSON(genotypes, auto_unbox = TRUE))
  }
  "{}"
}
