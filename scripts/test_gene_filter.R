#!/usr/bin/env Rscript
setwd(normalizePath(file.path(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])), ".."), winslash = "/"))
source("global_cli.R")

vcf <- normalizePath(file.path("..", "HMC-1.final.vcf"), winslash = "/", mustWork = TRUE)
genes <- c("BCR-ABL", "MYC", "CD63", "REXO5", "ENO1")

cat("Expanded genes:", paste(parse_gene_filter(genes), collapse = ", "), "\n")
cnt <- count_vcf_variants_by_genes(vcf, genes)
print(cnt)

n <- 0L
stream_vcf_chunks(vcf, chunk_size = 5000L, processor = function(df, id) {
  f <- filter_variants_by_genes(df, genes)
  n <<- n + nrow(f)
}, progress_fn = NULL)
cat("Stream filter matched rows:", n, "\n")

lines <- readLines(vcf)
myc_line <- lines[grep("\\|MYC\\|", lines)][1]
parts <- strsplit(myc_line, "\t", fixed = TRUE)[[1]]
row <- parse_variant_from_vcf_fields(
  parts[[1]], parts[[2]], parts[[4]], strsplit(parts[[5]], ",")[[1]][1],
  suppressWarnings(as.numeric(parts[[6]])), parts[[7]], parts[[8]]
)
cat("MYC variant primary gene:", row$gene, "| all_genes:", row$all_genes, "\n")
f1 <- filter_variants_by_genes(row, "MYC")
cat("MYC filter keeps row:", nrow(f1), "| display gene:", f1$gene, "\n")

if (cnt$total != n) stop("Count mismatch: full scan ", cnt$total, " vs stream ", n)
if (n < 20L) stop("Expected at least 20 gene-matched variants, got ", n)
cat("OK\n")
