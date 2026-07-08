#!/usr/bin/env Rscript
# Import GWAS Catalog associations as SUPPLEMENTARY evidence only.
#
# Policy (ACMGamp professional germline interpretation):
#   - GWAS rows are written to data/reference/gwas_supplementary_evidence.tsv
#   - They are NEVER copied into ps4_case_control_enrichment.tsv
#   - Automated PS4 is NOT awarded from GWAS data
#   - Curators review GWAS notes in the evidence report; Group A manual PS4 if appropriate
#
# Usage:
#   Rscript scripts/import_gwas_supplementary.R path/to/gwas_catalog.tsv
#   Rscript scripts/import_gwas_supplementary.R path/to/gwas.tsv --gene-panel data/gene_panels/pm1_hotspot_genes.csv
#
# Expected GWAS input columns (flexible names):
#   chrom/CHR/chr, pos/BP/start, rsid/SNPS, gene/GENE, trait/DISEASE/TRAIT,
#   pvalue/P/LOG10P, odds_ratio/OR/BETA, maf/FRQ

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  cat(
    "Usage: Rscript scripts/import_gwas_supplementary.R <gwas_catalog.tsv> [--gene-panel genes.csv]\n",
    "Output: data/reference/gwas_supplementary_evidence.tsv (supplementary only; no auto-PS4)\n",
    sep = ""
  )
  quit(status = 1)
}

script_dir <- dirname(normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]),
  winslash = "/", mustWork = FALSE
))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
setwd(project_root)
source("global_cli.R")

input_path <- args[[1L]]
gene_panel_path <- NA_character_
if (length(args) >= 3L && args[[2L]] == "--gene-panel") {
  gene_panel_path <- args[[3L]]
}

if (!file.exists(input_path)) stop("GWAS input not found: ", input_path)

pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) NA_character_ else hit[[1L]]
}

raw <- utils::read.delim(input_path, stringsAsFactors = FALSE, comment.char = "#", check.names = FALSE)
if (nrow(raw) == 0L) stop("GWAS input is empty.")

chrom_col <- pick_col(raw, c("chrom", "CHR", "chr", "Chromosome"))
pos_col <- pick_col(raw, c("pos", "BP", "position", "start", "Position"))
rsid_col <- pick_col(raw, c("rsid", "SNPS", "rsID", "snp"))
gene_col <- pick_col(raw, c("gene", "GENE", "MAPPED_GENE", "Reported Gene"))
trait_col <- pick_col(raw, c("trait", "DISEASE", "DISEASE/TRAIT", "Trait", "Phenotype"))
p_col <- pick_col(raw, c("pvalue", "P", "P-VALUE", "p_value"))
or_col <- pick_col(raw, c("odds_ratio", "OR", "OR/BETA", "oddsratio"))
maf_col <- pick_col(raw, c("maf", "FRQ", "freq", "MAF"))

if (is.na(chrom_col) || is.na(pos_col)) {
  stop("GWAS file must include chromosome and position columns.")
}

out <- data.frame(
  chrom = as.character(raw[[chrom_col]]),
  pos = suppressWarnings(as.integer(raw[[pos_col]])),
  ref = ".",
  alt = ".",
  rsid = if (!is.na(rsid_col)) as.character(raw[[rsid_col]]) else NA_character_,
  gene = if (!is.na(gene_col)) as.character(raw[[gene_col]]) else NA_character_,
  trait = if (!is.na(trait_col)) as.character(raw[[trait_col]]) else NA_character_,
  pvalue = if (!is.na(p_col)) suppressWarnings(as.numeric(raw[[p_col]])) else NA_real_,
  odds_ratio = if (!is.na(or_col)) suppressWarnings(as.numeric(raw[[or_col]])) else NA_real_,
  maf = if (!is.na(maf_col)) suppressWarnings(as.numeric(raw[[maf_col]])) else NA_real_,
  source = "GWAS_Catalog_import",
  stringsAsFactors = FALSE
)

out$chrom <- sub("^chr", "chr", out$chrom, ignore.case = TRUE)
out$chrom <- ifelse(grepl("^chr", out$chrom, ignore.case = TRUE), out$chrom, paste0("chr", out$chrom))
out$gene <- toupper(trimws(sub(";.*$", "", out$gene)))
out$rsid <- vapply(out$rsid, normalize_rsid, character(1L))
out <- out[!is.na(out$pos), , drop = FALSE]

if (!is.na(gene_panel_path) && file.exists(gene_panel_path)) {
  panel <- toupper(trimws(utils::read.csv(gene_panel_path, stringsAsFactors = FALSE)$gene))
  panel <- panel[nzchar(panel)]
  out <- out[out$gene %in% panel, , drop = FALSE]
  cat("Filtered to gene panel (", length(panel), " genes): ", nrow(out), " rows\n", sep = "")
}

out$variant_key <- variant_key_chr_pos_ref_alt(out$chrom, out$pos, out$ref, out$alt)
out <- out[!duplicated(out$variant_key), , drop = FALSE]

dest <- GWAS_SUPPLEMENTARY_DB_PATH
dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

write.table(
  out[, c("chrom", "pos", "ref", "alt", "rsid", "gene", "trait", "pvalue", "odds_ratio", "maf", "source")],
  file = dest,
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("Wrote ", nrow(out), " GWAS supplementary rows to:\n  ", dest, "\n", sep = "")
cat("These rows do NOT trigger automated PS4.\n")
cat("Add variant-specific case-control data to ps4_case_control_enrichment.tsv for PS4 automation.\n")
