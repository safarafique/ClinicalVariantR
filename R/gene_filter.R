#' Common fusion / alias symbols mapped to HGNC gene symbols for panel search.
#' @noRd
GENE_SYMBOL_ALIASES <- c(
  "BCR-ABL" = "ABL1,BCR",
  "BCR::ABL1" = "ABL1,BCR",
  "BCRABL1" = "ABL1,BCR",
  "PML-RARA" = "PML,RARA",
  "PML::RARA" = "PML,RARA"
)

expand_gene_aliases <- function(genes) {
  if (length(genes) == 0L) return(character())
  genes <- trimws(as.character(genes))
  genes <- genes[nzchar(genes)]
  out <- character()
  for (g in genes) {
    key <- toupper(gsub("\\s+", "", g))
    alias_keys <- toupper(gsub("\\s+", "", names(GENE_SYMBOL_ALIASES)))
    idx <- match(key, alias_keys)
    if (!is.na(idx)) {
      out <- c(out, unlist(strsplit(GENE_SYMBOL_ALIASES[[idx]], ",", fixed = TRUE)))
    } else {
      out <- c(out, g)
    }
  }
  unique(trimws(out[nzchar(out)]))
}

#' Parse user-entered gene symbols (comma, semicolon, space, or newline separated).
#' @noRd
parse_gene_filter <- function(text) {
  if (is.null(text) || length(text) == 0L) return(character())
  text <- paste(as.character(text), collapse = "\n")
  if (!nzchar(trimws(text))) return(character())
  parts <- unlist(strsplit(text, "[,;\\s]+", perl = TRUE))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  expand_gene_aliases(unique(parts))
}

variant_row_genes <- function(variants_df, row_idx) {
  genes <- character()
  if (!is.null(variants_df$gene)) {
    g <- variants_df$gene[[row_idx]]
    if (!is.na(g) && nzchar(g)) genes <- c(genes, g)
  }
  if (!is.null(variants_df$all_genes)) {
    extra <- unlist(strsplit(as.character(variants_df$all_genes[[row_idx]] %||% ""), ";", fixed = TRUE))
    extra <- trimws(extra)
    genes <- c(genes, extra[nzchar(extra)])
  }
  unique(genes)
}

variant_row_matches_genes <- function(variants_df, row_idx, genes_norm) {
  if (length(genes_norm) == 0L) return(TRUE)
  row_genes <- toupper(variant_row_genes(variants_df, row_idx))
  any(row_genes %in% genes_norm)
}

#' After filtering, set primary gene column to the panel symbol that matched.
#' @noRd
prioritize_panel_gene_annotation <- function(variants_df, genes) {
  if (is.null(variants_df) || nrow(variants_df) == 0L) return(variants_df)
  genes <- parse_gene_filter(genes)
  if (length(genes) == 0L) return(variants_df)
  genes_norm <- toupper(genes)
  gene_lookup <- setNames(genes, genes_norm)
  for (i in seq_len(nrow(variants_df))) {
    row_genes <- toupper(variant_row_genes(variants_df, i))
    hit <- intersect(row_genes, genes_norm)
    if (length(hit) > 0L) variants_df$gene[[i]] <- gene_lookup[[hit[[1L]]]]
  }
  variants_df
}

#' Keep variant rows whose SYMBOL/gene matches the requested gene list (case-insensitive).
#' Matches any gene in VEP CSQ / SnpEff ANN / ANNOVAR annotation, not only the primary transcript.
#' @noRd
filter_variants_by_genes <- function(variants_df, genes) {
  if (is.null(variants_df) || nrow(variants_df) == 0L) return(variants_df)
  genes <- parse_gene_filter(genes)
  if (length(genes) == 0L) return(variants_df)
  genes_norm <- toupper(genes)
  keep <- vapply(seq_len(nrow(variants_df)), function(i) {
    variant_row_matches_genes(variants_df, i, genes_norm)
  }, FUN.VALUE = logical(1))
  out <- variants_df[keep, , drop = FALSE]
  prioritize_panel_gene_annotation(out, genes)
}

#' Count variants in a VCF that overlap requested gene symbols (full file scan).
#' @noRd
count_vcf_variants_by_genes <- function(vcf_path, genes, pass_only = FALSE) {
  genes <- parse_gene_filter(genes)
  if (length(genes) == 0L || !file.exists(vcf_path)) {
    return(list(total = 0L, matched_genes = character(), missing_genes = character(),
                per_gene = setNames(integer(), character())))
  }
  genes_norm <- toupper(genes)
  per_gene <- setNames(rep(0L, length(genes_norm)), genes_norm)

  con <- if (grepl("\\.gz$", vcf_path, ignore.case = TRUE)) gzfile(vcf_path, "rt") else file(vcf_path, "rt")
  on.exit(close(con), add = TRUE)

  total <- 0L
  buffer_n <- if (exists("VCF_LINE_BUFFER", inherits = TRUE)) VCF_LINE_BUFFER else 50000L
  repeat {
    lines <- readLines(con, n = buffer_n, warn = FALSE)
    if (length(lines) == 0L) break
    for (line in lines) {
      if (grepl("^#", line)) next
      parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(parts) < 8L) next
      if (isTRUE(pass_only) && !(parts[[7L]] %in% c("PASS", "."))) next
      info_genes <- toupper(extract_all_annotation_genes(parts[[8L]]))
      if (length(info_genes) == 0L) next
      hit <- intersect(info_genes, genes_norm)
      if (length(hit) == 0L) next
      total <- total + 1L
      for (g in hit) per_gene[[g]] <- per_gene[[g]] + 1L
    }
  }

  matched_genes <- genes[vapply(genes, function(g) per_gene[[toupper(g)]] > 0L, logical(1))]
  list(
    total = total,
    matched_genes = matched_genes,
    missing_genes = setdiff(genes, matched_genes),
    per_gene = per_gene
  )
}

format_gene_filter_label <- function(genes) {
  genes <- parse_gene_filter(genes)
  if (length(genes) == 0L) return("(none)")
  if (length(genes) <= 5L) paste(genes, collapse = ", ")
  paste(paste(genes[seq_len(5L)], collapse = ", "), sprintf("(+%d more)", length(genes) - 5L))
}

pathogenic_tier_count <- function(classification_counts) {
  if (is.null(classification_counts) || length(classification_counts) == 0L) return(0L)
  sum(as.integer(classification_counts[["Pathogenic"]] %||% 0L),
      as.integer(classification_counts[["Likely Pathogenic"]] %||% 0L))
}
