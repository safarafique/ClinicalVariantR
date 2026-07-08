#' Unified VCF parsing for VEP CSQ, SnpEff ANN, and ANNOVAR-style INFO.

VEP_CSQ_COLUMN_NAMES <- c(
  "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE",
  "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
  "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS",
  "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "MANE", "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL",
  "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "GENE_PHENO",
  "SIFT", "PolyPhen", "DOMAINS", "miRNA", "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
  "gnomADe_AF", "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_ASJ_AF", "gnomADe_EAS_AF",
  "gnomADe_FIN_AF", "gnomADe_MID_AF", "gnomADe_NFE_AF", "gnomADe_REMAINING_AF", "gnomADe_SAS_AF",
  "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF", "gnomADg_ASJ_AF",
  "gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF", "gnomADg_NFE_AF", "gnomADg_REMAINING_AF",
  "gnomADg_SAS_AF", "MAX_AF", "MAX_AF_POPS", "FREQS", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED",
  "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS"
)

detect_vcf_annotation <- function(info_declared = character(), sample_info = character()) {
  has_csq <- "CSQ" %in% info_declared ||
    any(grepl("(^|;)CSQ=", sample_info, perl = TRUE))
  has_ann <- "ANN" %in% info_declared ||
    any(grepl("(^|;)ANN=", sample_info, perl = TRUE))
  has_annovar <- any(c("Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene") %in% info_declared) ||
    any(grepl("(^|;)(Func\\.refGene|Gene\\.refGene|ExonicFunc\\.refGene)=", sample_info, perl = TRUE))

  if (has_csq) {
    list(source = "CSQ", build_hint = "GRCh38", label = "VEP CSQ")
  } else if (has_ann) {
    list(source = "ANN", build_hint = "GRCh37", label = "SnpEff ANN")
  } else if (has_annovar) {
    list(source = "ANNOVAR", build_hint = "unknown", label = "ANNOVAR")
  } else {
    list(source = "unknown", build_hint = "unknown", label = "No CSQ/ANN/ANNOVAR")
  }
}

parse_info_tag <- function(info, tag) {
  m <- regexpr(paste0("(?:^|;)", tag, "=([^;]+)"), info, perl = TRUE)
  if (m[1L] == -1L) return(NA_character_)
  sub(paste0("^.*(?:^|;)", tag, "="), "", regmatches(info, m)[[1L]])
}

parse_info_num_max <- function(info, tag, scale_esp = FALSE) {
  raw <- parse_info_tag(info, tag)
  if (is.na(raw) || !nzchar(raw)) return(NA_real_)
  nums <- suppressWarnings(as.numeric(unlist(strsplit(raw, "[,|:&]"))))
  nums <- nums[!is.na(nums)]
  if (length(nums) == 0L) return(NA_real_)
  val <- max(nums)
  if (isTRUE(scale_esp) && tag == "ESP6500_MAF") val <- val / 100
  val
}

parse_spliceai_max <- function(raw) {
  if (is.na(raw) || !nzchar(raw)) return(NA_real_)
  nums <- suppressWarnings(as.numeric(unlist(strsplit(raw, "[|,&:]"))))
  nums <- nums[!is.na(nums)]
  if (length(nums) == 0L) return(NA_real_)
  max(nums)
}

decode_clinvar_sig <- function(x) {
  x <- scalar_chr(x, default = "")
  if (!nzchar(x) || x == ".") return(NA_character_)
  if (grepl("pathogenic|benign|uncertain", x, ignore.case = TRUE)) {
    return(gsub("_", " ", as.character(x), fixed = TRUE))
  }
  codes <- c(
    "0" = "Uncertain significance", "1" = "not provided",
    "2" = "Benign", "3" = "Likely benign",
    "4" = "Likely pathogenic", "5" = "Pathogenic"
  )
  parts <- unlist(strsplit(as.character(x), "[,|&]"))
  mapped <- unname(codes[parts])
  mapped <- mapped[!is.na(mapped)]
  if (length(mapped) == 0L) as.character(x) else paste(unique(mapped), collapse = ";")
}

parse_vep_prediction_field <- function(raw) {
  raw <- scalar_chr(raw, default = "")
  if (!nzchar(raw) || raw == ".") {
    return(list(text = NA_character_, score = NA_real_))
  }
  parts <- unlist(strsplit(as.character(raw), "|", fixed = TRUE))
  text <- parts[[1L]]
  score <- NA_real_
  m <- regexpr("\\(([0-9.]+)\\)", text, perl = TRUE)
  if (m[1L] != -1L) {
    score <- suppressWarnings(as.numeric(sub(".*\\(([^)]+)\\).*", "\\1", regmatches(text, m)[[1L]])))
  }
  list(text = text, score = score)
}

csq_get_part <- function(parts, name, column_names = VEP_CSQ_COLUMN_NAMES) {
  idx <- match(name, column_names)
  if (is.na(idx) || length(parts) < idx) return(NA_character_)
  val <- parts[[idx]]
  if (is.na(val) || !nzchar(val) || val == ".") return(NA_character_)
  val
}

csq_numeric_max <- function(parts, names, column_names = VEP_CSQ_COLUMN_NAMES) {
  vals <- vapply(names, function(nm) {
    raw <- csq_get_part(parts, nm, column_names)
    if (is.na(raw)) return(NA_real_)
    suppressWarnings(as.numeric(raw))
  }, numeric(1))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) NA_real_ else max(vals)
}

csq_consequence_priority <- function(consequence) {
  cons <- tolower(consequence %||% "")
  if (grepl("stop_gained|frameshift|splice_donor|splice_acceptor|stop_lost|start_lost", cons)) return(1L)
  if (grepl("missense_variant", cons, fixed = TRUE)) return(2L)
  if (grepl("synonymous_variant", cons, fixed = TRUE)) return(3L)
  if (grepl("splice_", cons, fixed = TRUE)) return(4L)
  if (grepl("5_prime_utr|3_prime_utr|utr", cons, fixed = TRUE)) return(5L)
  if (grepl("intron", cons, fixed = TRUE)) return(6L)
  if (grepl("regulatory|upstream|downstream", cons, fixed = TRUE)) return(8L)
  7L
}

#' Match VEP CSQ Allele field to VCF ref/alt (handles insertions and deletions).
csq_allele_matches_vcf <- function(csq_allele, ref = NA_character_, alt = NA_character_) {
  csq_allele <- scalar_chr(csq_allele, default = "")
  ref <- scalar_chr(ref, default = "")
  alt <- scalar_chr(alt, default = "")
  if (!nzchar(csq_allele) || csq_allele == ".") return(TRUE)
  if (!nzchar(alt)) return(TRUE)
  if (identical(csq_allele, alt)) return(TRUE)
  if (csq_allele == "-") return(TRUE)

  ref <- ref %||% ""
  if (nzchar(ref) && nchar(alt) > nchar(ref)) {
    inserted <- substr(alt, nchar(ref) + 1L, nchar(alt))
    if (identical(csq_allele, inserted)) return(TRUE)
  }
  if (nzchar(ref) && nchar(ref) > nchar(alt) && (identical(csq_allele, alt) || csq_allele == "-")) {
    return(TRUE)
  }
  FALSE
}

csq_entry_priority <- function(entry) {
  score <- 0L
  biotype <- tolower(entry$biotype %||% "")
  cons <- tolower(entry$consequence %||% "")

  if (grepl("pseudogene|processed_pseudogene|lncrna|lincrna|mirna|snorna|snrna|rrna|misc_rna", biotype)) {
    score <- score - 3000L
  }
  if (grepl("upstream_gene_variant|downstream_gene_variant|intergenic_variant|regulatory_region", cons)) {
    score <- score - 2500L
  }
  if (grepl("intron_variant", cons, fixed = TRUE) && !grepl("splice", cons)) {
    score <- score - 800L
  }

  if (isTRUE(entry$mane_select)) score <- score + 2500L
  if (identical(entry$canonical, "YES")) score <- score + 1200L
  if (identical(biotype, "protein_coding")) score <- score + 1500L

  if (grepl("stop_gained|frameshift|splice_donor|splice_acceptor|stop_lost|start_lost", cons)) {
    score <- score + 2000L
  } else if (grepl("missense_variant", cons, fixed = TRUE)) {
    score <- score + 900L
  } else if (grepl("inframe", cons)) {
    score <- score + 500L
  }

  if (nzchar(scalar_chr(entry$gene %||% "", default = ""))) score <- score + 50L
  score <- score - csq_consequence_priority(entry$consequence) * 5L
  if (nzchar(scalar_chr(entry$impact %||% "", default = ""))) {
    score <- score + switch(toupper(scalar_chr(entry$impact)), 
      "HIGH" = 600L, "MODERATE" = 350L, "LOW" = 120L, "MODIFIER" = 0L, 0L)
  }
  score
}

parse_one_csq_entry <- function(parts, column_names = VEP_CSQ_COLUMN_NAMES) {
  sift <- parse_vep_prediction_field(csq_get_part(parts, "SIFT", column_names))
  polyphen <- parse_vep_prediction_field(csq_get_part(parts, "PolyPhen", column_names))
  pop_af <- csq_numeric_max(parts, c(
    "MAX_AF", "gnomADg_AF", "gnomADe_AF", "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF"
  ), column_names)
  popmax_af <- csq_numeric_max(parts, c(
    "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
    "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_ASJ_AF", "gnomADe_EAS_AF",
    "gnomADe_FIN_AF", "gnomADe_MID_AF", "gnomADe_NFE_AF", "gnomADe_REMAINING_AF", "gnomADe_SAS_AF",
    "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF", "gnomADg_ASJ_AF",
    "gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF", "gnomADg_NFE_AF",
    "gnomADg_REMAINING_AF", "gnomADg_SAS_AF"
  ), column_names)

  list(
    allele = csq_get_part(parts, "Allele", column_names),
    consequence = csq_get_part(parts, "Consequence", column_names),
    impact = csq_get_part(parts, "IMPACT", column_names),
    gene = csq_get_part(parts, "SYMBOL", column_names),
    biotype = csq_get_part(parts, "BIOTYPE", column_names),
    canonical = csq_get_part(parts, "CANONICAL", column_names),
    mane_select = {
      ms <- csq_get_part(parts, "MANE_SELECT", column_names)
      !is.na(ms) && nzchar(ms)
    },
    hgvs_c = csq_get_part(parts, "HGVSc", column_names),
    hgvs_p = csq_get_part(parts, "HGVSp", column_names),
    protein_position = csq_get_part(parts, "Protein_position", column_names),
    amino_acids = csq_get_part(parts, "Amino_acids", column_names),
    exon = csq_get_part(parts, "EXON", column_names),
    sift = sift$text,
    sift_score = sift$score,
    polyphen = polyphen$text,
    polyphen_score = polyphen$score,
    gnomad_af = pop_af,
    popmax_af = popmax_af,
    clin_sig = decode_clinvar_sig(csq_get_part(parts, "CLIN_SIG", column_names)),
    rsids = csq_get_part(parts, "Existing_variation", column_names),
    cadd = NA_real_,
    revel = NA_real_,
    spliceai_max = NA_real_,
    alphamissense = NA_real_
  )
}

parse_csq_entry <- function(csq_string, alt = NA_character_, ref = NA_character_,
                            column_names = VEP_CSQ_COLUMN_NAMES) {
  empty <- list(
    gene = NA_character_, consequence = NA_character_, impact = NA_character_,
    hgvs_c = NA_character_, hgvs_p = NA_character_,
    protein_position = NA_character_, amino_acids = NA_character_, exon = NA_character_,
    biotype = NA_character_, is_canonical_transcript = FALSE,
    is_protein_coding = FALSE,
    sift = NA_character_, sift_score = NA_real_,
    polyphen = NA_character_, polyphen_score = NA_real_,
    gnomad_af = NA_real_, popmax_af = NA_real_, clin_sig = NA_character_, rsids = NA_character_,
    cadd = NA_real_, revel = NA_real_,
    spliceai_max = NA_real_, alphamissense = NA_real_
  )
  if (is.na(csq_string) || !nzchar(csq_string)) return(empty)

  entries <- strsplit(csq_string, ",", fixed = TRUE)[[1L]]
  parsed <- lapply(entries, function(entry) {
    parts <- strsplit(entry, "|", fixed = TRUE)[[1L]]
    if (length(parts) < 4L) return(NULL)
    row <- parse_one_csq_entry(parts, column_names)
    if (!csq_allele_matches_vcf(row$allele, ref, alt)) {
      return(NULL)
    }
    row
  })
  parsed <- Filter(Negate(is.null), parsed)
  if (length(parsed) == 0L) return(empty)

  pri <- vapply(parsed, csq_entry_priority, integer(1L))
  best <- parsed[[which.max(pri)]]
  best$is_canonical_transcript <- identical(best$canonical, "YES") || isTRUE(best$mane_select)
  best$is_protein_coding <- identical(tolower(best$biotype %||% ""), "protein_coding")
  best$rsids <- best$rsids %||% NA_character_
  best$mane_select <- NULL
  best$canonical <- NULL
  best$allele <- NULL
  best
}

parse_csq_all_genes <- function(csq_string, column_names = VEP_CSQ_COLUMN_NAMES) {
  if (is.na(csq_string) || !nzchar(csq_string)) return(character())
  genes <- character()
  entries <- strsplit(csq_string, ",", fixed = TRUE)[[1L]]
  for (entry in entries) {
    parts <- strsplit(entry, "|", fixed = TRUE)[[1L]]
    if (length(parts) < 4L) next
    row <- parse_one_csq_entry(parts, column_names)
    if (!is.na(row$gene) && nzchar(row$gene) && row$gene != ".") {
      genes <- c(genes, row$gene)
    }
  }
  unique(trimws(genes))
}

extract_all_annotation_genes <- function(info) {
  if (is.na(info) || !nzchar(info)) return(character())
  genes <- character()

  ann_val <- parse_info_tag(info, "ANN")
  if (!is.na(ann_val) && nzchar(ann_val)) {
    entries <- strsplit(ann_val, ",", fixed = TRUE)[[1L]]
    for (entry in entries) {
      parts <- strsplit(entry, "|", fixed = TRUE)[[1L]]
      if (length(parts) >= 4L && nzchar(parts[[4L]]) && parts[[4L]] != ".") {
        genes <- c(genes, parts[[4L]])
      }
    }
  }

  csq_val <- parse_info_tag(info, "CSQ")
  if (!is.na(csq_val) && nzchar(csq_val)) {
    genes <- c(genes, parse_csq_all_genes(csq_val))
  }

  gene <- parse_info_tag(info, "GENE")
  if (!is.na(gene) && nzchar(gene) && gene != ".") genes <- c(genes, gene)
  gene <- parse_info_tag(info, "Gene.refGene")
  if (!is.na(gene) && nzchar(gene) && gene != ".") genes <- c(genes, gene)
  gene <- parse_info_tag(info, "Gene.knownGene")
  if (!is.na(gene) && nzchar(gene) && gene != ".") genes <- c(genes, gene)

  unique(trimws(genes[nzchar(genes)]))
}

parse_ann_entry <- function(ann_string, alt = NA_character_) {
  if (is.na(ann_string) || !nzchar(ann_string)) {
    return(list(gene = NA_character_, consequence = NA_character_, impact = NA_character_,
                hgvs_c = NA_character_, hgvs_p = NA_character_))
  }
  entries <- strsplit(ann_string, ",", fixed = TRUE)[[1L]]
  parsed <- lapply(entries, function(entry) {
    parts <- strsplit(entry, "|", fixed = TRUE)[[1L]]
    if (length(parts) < 4L) return(NULL)
    list(
      allele = parts[[1L]], consequence = parts[[2L]], impact = parts[[3L]], gene = parts[[4L]],
      hgvs_c = if (length(parts) >= 10L) parts[[10L]] else NA_character_,
      hgvs_p = if (length(parts) >= 11L) parts[[11L]] else NA_character_
    )
  })
  parsed <- Filter(Negate(is.null), parsed)
  if (length(parsed) == 0L) {
    return(list(gene = NA_character_, consequence = NA_character_, impact = NA_character_,
                hgvs_c = NA_character_, hgvs_p = NA_character_))
  }
  cons <- vapply(parsed, `[[`, "", "consequence")
  rank <- ifelse(grepl("stop_gained|frameshift|splice_donor|splice_acceptor|stop_lost|start_lost", cons), 1L,
    ifelse(grepl("missense_variant", cons), 2L,
      ifelse(grepl("inframe_", cons), 3L, 4L)))
  parsed[[which.min(rank)]]
}

parse_annovar_fields <- function(info) {
  gene <- parse_info_tag(info, "Gene.refGene")
  if (is.na(gene)) gene <- parse_info_tag(info, "Gene.knownGene")
  func <- parse_info_tag(info, "Func.refGene")
  exonic <- parse_info_tag(info, "ExonicFunc.refGene")
  consequence <- if (!is.na(exonic) && nzchar(exonic) && exonic != ".") exonic else func
  clin <- decode_clinvar_sig(parse_info_tag(info, "CLNSIG"))
  if (is.na(clin)) clin <- decode_clinvar_sig(parse_info_tag(info, "clinvar_CLNSIG"))

  sift_pred <- parse_info_tag(info, "SIFT_pred")
  polyphen_pred <- parse_info_tag(info, "Polyphen2_HVAR_pred")
  if (is.na(polyphen_pred)) polyphen_pred <- parse_info_tag(info, "Polyphen2_HDIV_pred")

  list(
    gene = gene,
    consequence = consequence,
    impact = NA_character_,
    hgvs_c = NA_character_,
    hgvs_p = NA_character_,
    sift = sift_pred,
    polyphen = polyphen_pred,
    gnomad_af = parse_info_num_max(info, "gnomAD_genome_ALL"),
    clin_sig = clin
  )
}

parse_polyphen_from_esp <- function(ph) {
  if (is.na(ph) || !nzchar(ph)) return(list(text = NA_character_, score = NA_real_))
  first <- strsplit(ph, ",")[[1L]][1L]
  score <- NA_real_
  if (grepl(":", first, fixed = TRUE)) {
    score <- suppressWarnings(as.numeric(sub(".*:", "", first)))
  }
  list(text = first, score = score)
}

population_af_from_annotation <- function(
    af_1000g = NA_real_,
    af_esp = NA_real_,
    csq_af = NA_real_,
    gnomad_info = NA_real_) {
  af_1000g <- scalar_num(af_1000g)
  af_esp <- scalar_num(af_esp)
  csq_af <- scalar_num(csq_af)
  gnomad_info <- scalar_num(gnomad_info)
  vals <- c(
    af_1000g,
    if (!is.na(af_esp)) af_esp / 100 else NA_real_,
    csq_af,
    gnomad_info
  )
  vals <- suppressWarnings(as.numeric(vals))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) NA_real_ else max(vals)
}

parse_variant_from_vcf_fields <- function(chrom, pos, ref, alt, qual = NA_real_, filter = ".", info = ".") {
  chrom <- scalar_chr(chrom)
  pos <- suppressWarnings(as.integer(scalar_chr(pos)))
  ref <- scalar_chr(ref)
  alt <- scalar_chr(alt)
  qual <- scalar_num(qual)
  ann <- detect_vcf_annotation(sample_info = scalar_chr(info, default = ""))
  gene <- NA_character_
  consequence <- NA_character_
  impact <- NA_character_
  hgvs_c <- NA_character_
  hgvs_p <- NA_character_
  protein_position <- NA_character_
  amino_acids <- NA_character_
  exon <- NA_character_
  biotype <- NA_character_
  is_canonical_transcript <- FALSE
  is_protein_coding <- FALSE
  rsids <- NA_character_
  sift <- NA_character_
  sift_score <- NA_real_
  polyphen <- NA_character_
  polyphen_score <- NA_real_
  clinvar <- NA_character_
  csq_af <- NA_real_
  csq_popmax <- NA_real_
  cadd <- NA_real_
  revel <- NA_real_
  spliceai_max <- NA_real_
  alphamissense <- NA_real_

  if (ann$source == "CSQ") {
    csq <- parse_info_tag(info, "CSQ")
    fields <- parse_csq_entry(csq, alt, ref)
    gene <- fields$gene
    consequence <- fields$consequence
    impact <- fields$impact
    hgvs_c <- fields$hgvs_c
    hgvs_p <- fields$hgvs_p
    protein_position <- fields$protein_position
    amino_acids <- fields$amino_acids
    exon <- fields$exon %||% NA_character_
    biotype <- fields$biotype
    is_canonical_transcript <- scalar_lgl(fields$is_canonical_transcript)
    is_protein_coding <- scalar_lgl(fields$is_protein_coding)
    rsids <- fields$rsids
    sift <- fields$sift
    sift_score <- fields$sift_score
    polyphen <- fields$polyphen
    polyphen_score <- fields$polyphen_score
    csq_af <- fields$gnomad_af
    csq_popmax <- fields$popmax_af
    clinvar <- fields$clin_sig
    cadd <- fields$cadd
    revel <- fields$revel
    spliceai_max <- fields$spliceai_max
    alphamissense <- fields$alphamissense
  } else if (ann$source == "ANN") {
    ann_val <- parse_info_tag(info, "ANN")
    fields <- parse_ann_entry(ann_val, alt)
    gene <- fields$gene
    consequence <- fields$consequence
    impact <- fields$impact
    hgvs_c <- fields$hgvs_c
    hgvs_p <- fields$hgvs_p
    ph <- parse_polyphen_from_esp(parse_info_tag(info, "ESP6500_PH"))
    polyphen <- ph$text
    polyphen_score <- ph$score
    clinvar <- decode_clinvar_sig(parse_info_tag(info, "CLINVAR_CLNSIG"))
  } else if (ann$source == "ANNOVAR") {
    fields <- parse_annovar_fields(info)
    gene <- fields$gene
    consequence <- fields$consequence
    impact <- fields$impact
    hgvs_c <- fields$hgvs_c
    hgvs_p <- fields$hgvs_p
    sift <- fields$sift
    polyphen <- fields$polyphen
    csq_af <- fields$gnomad_af
    clinvar <- fields$clin_sig
  }

  af_1000g <- parse_info_num_max(info, "1000Gp3_AF")
  af_esp <- parse_info_num_max(info, "ESP6500_MAF")
  gnomad_info <- parse_info_num_max(info, "gnomAD_AF")
  if (is.na(scalar_num(gnomad_info))) gnomad_info <- parse_info_num_max(info, "AF_gnomAD")

  population_af <- population_af_from_annotation(
    af_1000g = af_1000g,
    af_esp = af_esp,
    csq_af = scalar_num(csq_af),
    gnomad_info = gnomad_info
  )

  revel <- scalar_num(revel)
  if (is.na(revel)) revel <- scalar_num(parse_info_tag(info, "REVEL"))
  cadd <- scalar_num(cadd)
  if (is.na(cadd)) cadd <- scalar_num(parse_info_tag(info, "CADD_PHRED"))
  if (is.na(cadd)) cadd <- scalar_num(parse_info_tag(info, "CADD"))
  spliceai_max <- scalar_num(spliceai_max)
  if (is.na(spliceai_max)) {
    spliceai_max <- parse_spliceai_max(parse_info_tag(info, "SpliceAI"))
    if (is.na(scalar_num(spliceai_max))) spliceai_max <- parse_spliceai_max(parse_info_tag(info, "MAX_SPLICEAI"))
  }
  spliceai_max <- scalar_num(spliceai_max)
  alphamissense <- scalar_num(alphamissense)
  if (is.na(alphamissense)) {
    alphamissense <- scalar_num(parse_info_tag(info, "am_pathogenicity"))
    if (is.na(alphamissense)) alphamissense <- scalar_num(parse_info_tag(info, "AlphaMissense"))
  }

  gene <- scalar_chr(gene, default = NA_character_)
  consequence <- scalar_chr(consequence, default = NA_character_)
  hgvs_c <- scalar_chr(hgvs_c, default = NA_character_)
  hgvs_p <- scalar_chr(hgvs_p, default = NA_character_)
  protein_position <- scalar_chr(protein_position, default = NA_character_)
  amino_acids <- scalar_chr(amino_acids, default = NA_character_)
  exon <- scalar_chr(exon, default = NA_character_)
  biotype <- scalar_chr(biotype, default = NA_character_)
  rsids <- scalar_chr(rsids, default = NA_character_)
  sift <- scalar_chr(sift, default = NA_character_)
  polyphen <- scalar_chr(polyphen, default = NA_character_)
  clinvar <- scalar_chr(clinvar, default = NA_character_)
  sift_score <- scalar_num(sift_score)
  polyphen_score <- scalar_num(polyphen_score)
  csq_af <- scalar_num(csq_af)
  csq_popmax <- scalar_num(csq_popmax)
  impact <- scalar_chr(impact, default = NA_character_)

  all_genes <- extract_all_annotation_genes(info)
  if (length(all_genes) == 0L && nzchar(gene)) all_genes <- gene

  data.frame(
    variant_id = paste(chrom, pos, ref, alt, sep = ":"),
    chrom = chrom, pos = as.integer(pos), ref = ref, alt = alt,
    gene = gene, all_genes = paste(all_genes, collapse = ";"),
    consequence = consequence, impact = impact,
    hgvs_c = hgvs_c, hgvs_p = hgvs_p,
    protein_position = protein_position, amino_acids = amino_acids, exon = exon,
    biotype = biotype,
    is_canonical_transcript = is_canonical_transcript,
    is_protein_coding = is_protein_coding,
    rsids = rsids %||% NA_character_,
    annotation_source = ann$source,
    genome_build_hint = ann$build_hint,
    af_1000g = af_1000g, af_esp6500 = af_esp,
    population_af = population_af,
    gnomad_af = population_af,
    popmax_af = csq_popmax,
    sift = sift, sift_score = sift_score,
    polyphen = polyphen, polyphen_score = polyphen_score,
    revel_score = revel, REVEL = revel,
    cadd = cadd, spliceai_max = spliceai_max, alphamissense_score = alphamissense,
    clinvar_classification = clinvar, ClinVar = clinvar,
    info_csq = if (ann$source == "CSQ") parse_info_tag(info, "CSQ") else NA_character_,
    qual = qual, filter = filter,
    stringsAsFactors = FALSE
  )
}

#' Build PS1/PM5 lookup rows from VEP CSQ CLIN_SIG pathogenic entries in the same VCF record.
parse_csq_pathogenic_catalog <- function(csq_string, ref = NA_character_, alt = NA_character_) {
  empty <- data.frame(
    gene = character(), protein_position = integer(),
    ref_aa = character(), alt_aa = character(),
    hgvs_p = character(), clinical_significance = character(),
    stringsAsFactors = FALSE
  )
  if (is.na(csq_string) || !nzchar(csq_string)) return(empty)

  rows <- list()
  entries <- strsplit(csq_string, ",", fixed = TRUE)[[1L]]
  for (entry in entries) {
    parts <- strsplit(entry, "|", fixed = TRUE)[[1L]]
    if (length(parts) < 4L) next
    row <- parse_one_csq_entry(parts)
    if (!csq_allele_matches_vcf(row$allele, ref, alt)) next
    clin <- row$clin_sig %||% ""
    if (!is_pathogenic_significance(clin)) next
    change <- parse_protein_change(row$hgvs_p, row$amino_acids, row$protein_position)
    if (is.na(change$protein_position) || is.na(change$alt_aa)) next
    rows[[length(rows) + 1L]] <- data.frame(
      gene = toupper(trimws(row$gene %||% "")),
      protein_position = change$protein_position,
      ref_aa = change$ref_aa,
      alt_aa = change$alt_aa,
      hgvs_p = row$hgvs_p %||% "",
      clinical_significance = clin,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) return(empty)
  do.call(rbind, rows)
}

is_homopolymer_indel <- function(ref, alt) {
  ref <- toupper(as.character(ref %||% ""))
  alt <- toupper(as.character(alt %||% ""))
  if (!nzchar(ref) || !nzchar(alt)) return(FALSE)
  if (nchar(ref) == nchar(alt)) return(FALSE)
  inserted <- if (nchar(alt) > nchar(ref)) substr(alt, nchar(ref) + 1L, nchar(alt)) else ""
  deleted <- if (nchar(ref) > nchar(alt)) substr(ref, nchar(alt) + 1L, nchar(ref)) else ""
  seq <- if (nzchar(inserted)) inserted else deleted
  if (!nzchar(seq) || nchar(seq) < 4L) return(FALSE)
  length(unique(strsplit(seq, "")[[1L]])) == 1L
}
