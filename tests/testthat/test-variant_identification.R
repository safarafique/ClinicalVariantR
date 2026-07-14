test_that("normalize_chrom adds chr prefix consistently", {
  expect_equal(ACMGamp:::normalize_chrom("17"), "chr17")
  expect_equal(ACMGamp:::normalize_chrom("chr17"), "chr17")
  expect_equal(ACMGamp:::normalize_chrom("CHR17"), "chr17")
  expect_equal(ACMGamp:::normalize_chrom("X"), "chrX")
  expect_equal(ACMGamp:::normalize_chrom("MT"), "chrM")
  expect_equal(ACMGamp:::normalize_chrom("chrM"), "chrM")
})

test_that("variant_key_chr_pos_ref_alt normalizes chr for matching", {
  key_with <- ACMGamp:::variant_key_chr_pos_ref_alt("chr17", 41223094, "T", "G")
  key_without <- ACMGamp:::variant_key_chr_pos_ref_alt("17", 41223094, "T", "G")
  expect_equal(key_with, key_without)
  expect_equal(key_with, "chr17:41223094:T:G")
})

test_that("variant_key_from_parts aliases normalized key builder", {
  expect_equal(
    ACMGamp:::variant_key_from_parts("X", 154030906, "T", "C"),
    "chrX:154030906:T:C"
  )
})

test_that("variant_coords_equal matches chr and no-chr coordinates", {
  expect_true(ACMGamp:::variant_coords_equal("chr17", 1, "A", "G", "17", 1, "A", "G"))
  expect_false(ACMGamp:::variant_coords_equal("chr17", 1, "A", "G", "chr17", 2, "A", "G"))
})

test_that("split_vcf_alt_alleles expands multi-allelic ALT fields", {
  expect_equal(ACMGamp:::split_vcf_alt_alleles("A,G,T"), c("A", "G", "T"))
  expect_equal(ACMGamp:::split_vcf_alt_alleles("A"), "A")
  expect_equal(ACMGamp:::split_vcf_alt_alleles("."), ".")
})

test_that("parse_vcf_line expands multi-allelic records", {
  line <- paste(
    c("chr1", "100", ".", "A", "G,T", "60", "PASS", "AF=0.1"),
    collapse = "\t"
  )
  parsed <- ACMGamp:::parse_vcf_line(line, NULL)
  expect_s3_class(parsed, "data.frame")
  expect_equal(nrow(parsed), 2L)
  expect_equal(parsed$alt, c("G", "T"))
  expect_equal(parsed$variant_id, c("chr1:100:A:G", "chr1:100:A:T"))
})

test_that("find_reference_hits matches PS4 DB with or without chr prefix", {
  db <- data.frame(
    chrom = "chrX",
    pos = 154030906L,
    ref = "T",
    alt = "C",
    gene = "MECP2",
    rsid = "rs1273236261",
    hgvs_p = "p.Thr320Ala",
    case_af = 0.18,
    control_af = 0.000001,
    disease = "Rett syndrome",
    source = "literature_cohort",
    ps4_confirmed = TRUE,
    stringsAsFactors = FALSE
  )
  db$chrom <- ACMGamp:::normalize_chrom(db$chrom)
  db$variant_key <- ACMGamp:::variant_key_chr_pos_ref_alt(db$chrom, db$pos, db$ref, db$alt)
  db$hgvs_p_norm <- "Thr320Ala"
  db$source_norm <- "literature_cohort"

  row_chr <- data.frame(
    chrom = "chrX", pos = 154030906L, ref = "T", alt = "C",
    gene = "MECP2", rsids = "rs1273236261", hgvs_p = "p.Thr320Ala",
    stringsAsFactors = FALSE
  )
  row_no_chr <- row_chr
  row_no_chr$chrom <- "X"

  hits_chr <- ACMGamp:::find_reference_hits(row_chr, db)
  hits_no_chr <- ACMGamp:::find_reference_hits(row_no_chr, db)

  expect_equal(nrow(hits_chr), 1L)
  expect_equal(nrow(hits_no_chr), 1L)
  expect_equal(hits_chr$gene[[1]], "MECP2")
})

test_that("find_reference_hits falls back to rsID when genomic key misses", {
  db <- data.frame(
    variant_key = "chr1:1:A:G",
    chrom = "chr1",
    pos = 1L,
    ref = "A",
    alt = "G",
    rsid = "rs9990001",
    gene = "GENE1",
    hgvs_p_norm = NA_character_,
    stringsAsFactors = FALSE
  )
  row <- data.frame(
    chrom = "chr2", pos = 2L, ref = "C", alt = "T",
    gene = "GENE1", rsids = "rs9990001", hgvs_p = NA_character_,
    stringsAsFactors = FALSE
  )
  hits <- ACMGamp:::find_reference_hits(row, db)
  expect_equal(nrow(hits), 1L)
  expect_equal(hits$rsid[[1]], "rs9990001")
})

test_that("parse_protein_change identifies PS1/PM5 protein coordinates", {
  change <- ACMGamp:::parse_protein_change("p.Thr320Ala", NA_character_, NA_character_)
  expect_equal(change$protein_position, 320L)
  expect_equal(change$ref_aa, "THR")
  expect_equal(change$alt_aa, "ALA")

  from_aa_field <- ACMGamp:::parse_protein_change(
    NA_character_, "Thr/Ala", "320/320"
  )
  expect_equal(from_aa_field$protein_position, 320L)
  expect_equal(from_aa_field$ref_aa, "THR")
  expect_equal(from_aa_field$alt_aa, "ALA")
})

test_that("annotate_variants joins reference rows across chr conventions", {
  variants <- data.frame(
    chrom = "17",
    pos = 41223094L,
    ref = "T",
    alt = "G",
    stringsAsFactors = FALSE
  )
  refs <- list(
    gnomad = data.table::data.table(
      chrom = "chr17", pos = 41223094L, ref = "T", alt = "G",
      AF = 0.001, popmax_AF = 0.002
    ),
    clinvar = data.table::data.table(
      chrom = "chr17", pos = 41223094L, ref = "T", alt = "G",
      clinical_significance = "Pathogenic", review_status = "reviewed"
    ),
    revel = data.table::data.table(
      chrom = "chr17", pos = 41223094L, ref = "T", alt = "G",
      REVEL = 0.91
    )
  )
  for (ref_name in names(refs)) {
    data.table::set(refs[[ref_name]], j = "chrom", value = ACMGamp:::normalize_chrom(refs[[ref_name]]$chrom))
    data.table::set(
      refs[[ref_name]],
      j = "variant_key",
      value = ACMGamp:::variant_key_chr_pos_ref_alt(
        refs[[ref_name]]$chrom,
        refs[[ref_name]]$pos,
        refs[[ref_name]]$ref,
        refs[[ref_name]]$alt
      )
    )
  }

  annotated <- ACMGamp:::annotate_variants(variants, refs)
  expect_equal(annotated$gnomad_af[[1]], 0.001)
  expect_equal(annotated$revel_score[[1]], 0.91)
  expect_equal(annotated$clinvar_classification[[1]], "Pathogenic")
})

test_that("match_variant_rows selects normalized coordinates from parsed table", {
  parsed <- data.frame(
    chrom = c("chr1", "chr1"),
    pos = c(100L, 100L),
    ref = "A",
    alt = c("G", "T"),
    stringsAsFactors = FALSE
  )
  hit <- ACMGamp:::match_variant_rows(parsed, "1", 100, "A", "T")
  expect_equal(nrow(hit), 1L)
  expect_equal(hit$alt[[1]], "T")
})
