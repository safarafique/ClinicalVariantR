#' Open VCF connection (plain or gzip).
open_vcf_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "r") else file(path, "r")
}

#' GATK and many callers use FILTER="." for passing variants; VCF spec also allows PASS.
vcf_filter_is_pass <- function(filter) {
  is.na(filter) || filter %in% c("PASS", ".")
}

extract_info_field <- function(info_string, key) {
  if (is.na(info_string) || info_string == ".") {
    return(if (key %in% c("AF", "REVEL")) NA_real_ else NA_character_)
  }
  pattern <- paste0("(^|;)", key, "=([^;]+)")
  m <- regexpr(pattern, info_string, perl = TRUE)
  if (m[1] == -1) {
    return(if (key %in% c("AF", "REVEL")) NA_real_ else NA_character_)
  }
  val <- sub(pattern, "\\2", regmatches(info_string, m))
  if (key %in% c("AF", "REVEL")) {
    num <- suppressWarnings(as.numeric(val))
    return(if (length(num) == 0 || is.na(num)) NA_real_ else num)
  }
  as.character(val)
}

#' Extract first VEP CSQ consequence token without retaining full CSQ string.
extract_vep_consequence <- function(csq_value) {
  if (is.na(csq_value) || csq_value == ".") return(NA_character_)
  first <- strsplit(csq_value, ",")[[1]][1]
  if (is.na(first)) return(NA_character_)
  parts <- strsplit(first, "|", fixed = TRUE)[[1]]
  if (length(parts) >= 2) parts[2] else first
}

#' Extract GENE from VEP CSQ (SYMBOL field) or INFO/GENE.
extract_gene_from_info <- function(info_string) {
  gene <- extract_info_field(info_string, "GENE")
  if (!is.na(gene) && nzchar(gene)) return(gene)
  csq <- extract_info_field(info_string, "CSQ")
  if (is.na(csq)) return(NA_character_)
  first <- strsplit(csq, ",")[[1]][1]
  parts <- strsplit(first, "|", fixed = TRUE)[[1]]
  if (length(parts) >= 4) parts[4] else NA_character_
}

#' Parse one VCF data line into a variant row (CSQ or ANN unified parser).
parse_vcf_line <- function(line, header_cols) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(fields) < 8) return(NULL)

  chrom <- fields[1]
  pos <- fields[2]
  ref <- fields[4]
  alt_raw <- fields[5]
  alt <- strsplit(alt_raw, ",")[[1]][1]
  qual <- suppressWarnings(as.numeric(fields[6]))
  filter <- if (length(fields) >= 7) fields[7] else "."
  info <- if (length(fields) >= 8) fields[8] else "."

  row <- parse_variant_from_vcf_fields(chrom, pos, ref, alt, qual, filter, info)

  if (length(fields) >= 9L && !is.null(header_cols) && length(header_cols) >= 9L) {
    format_str <- fields[9]
    sample_names <- header_cols[10:length(header_cols)]
    sample_fields <- if (length(fields) >= 10L) fields[10:length(fields)] else character()
    if (length(sample_names) > 0L && nzchar(format_str) && format_str != ".") {
      gts <- parse_vcf_genotypes(format_str, sample_fields, sample_names)
      row$sample_genotypes <- serialize_sample_genotypes(gts)
    } else {
      row$sample_genotypes <- "{}"
    }
  } else {
    row$sample_genotypes <- "{}"
  }

  row
}

#' Count all variant rows in a VCF (streaming, no memory load).
count_vcf_variants <- function(vcf_path, pass_only = FALSE, min_qual = 0) {
  con <- open_vcf_connection(vcf_path)
  on.exit(close(con), add = TRUE)

  header_cols <- NULL
  total <- 0L
  pass_n <- 0L

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0) break
    if (grepl("^#CHROM\t", line)) {
      header_cols <- strsplit(sub("^#", "", line), "\t")[[1]]
      next
    }
    if (grepl("^#", line)) next

    total <- total + 1L
    if (pass_only || min_qual > 0) {
      row <- parse_vcf_line(line, header_cols)
      if (is.null(row)) next
      if (pass_only && !vcf_filter_is_pass(row$filter)) next
      if (min_qual > 0 && (is.na(row$qual) || row$qual < min_qual)) next
      pass_n <- pass_n + 1L
    }
  }

  if (pass_only || min_qual > 0) {
    list(total = total, analyzed = pass_n, skipped = total - pass_n)
  } else {
    list(total = total, analyzed = total, skipped = 0L)
  }
}

#' Stream VCF variants in chunks; call processor(chunk_df, chunk_id) per batch.
stream_vcf_chunks <- function(
    vcf_path,
    chunk_size = 10000L,
    pass_only = FALSE,
    min_qual = 0,
    max_variants = Inf,
    processor,
    progress_fn = NULL) {

  con <- open_vcf_connection(vcf_path)
  on.exit(close(con), add = TRUE)

  header_cols <- NULL
  batch <- list()
  chunk_id <- 0L
  read_total <- 0L
  kept_total <- 0L
  skipped_total <- 0L

  flush_batch <- function() {
    if (length(batch) == 0) return(invisible(NULL))
    chunk_id <<- chunk_id + 1L
    chunk_df <- do.call(rbind, batch)
    rownames(chunk_df) <- NULL
    processor(chunk_df, chunk_id)
    batch <<- list()
    invisible(NULL)
  }

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0) break

    if (grepl("^#CHROM\t", line)) {
      header_cols <- strsplit(sub("^#", "", line), "\t")[[1]]
      next
    }
    if (grepl("^#", line)) next

    read_total <- read_total + 1L

    if (kept_total >= max_variants) {
      skipped_total <- skipped_total + 1L
      next
    }

    row <- parse_vcf_line(line, header_cols)
    if (is.null(row)) {
      skipped_total <- skipped_total + 1L
      next
    }
    if (pass_only && !vcf_filter_is_pass(row$filter)) {
      skipped_total <- skipped_total + 1L
      next
    }
    if (min_qual > 0 && (is.na(row$qual) || row$qual < min_qual)) {
      skipped_total <- skipped_total + 1L
      next
    }

    row$qual <- NULL
    row$filter <- NULL
    batch[[length(batch) + 1L]] <- row
    kept_total <- kept_total + 1L

    if (!is.null(progress_fn) && read_total %% 50000L == 0L) {
      progress_fn(read_total, kept_total, skipped_total)
    }

    if (length(batch) >= chunk_size) flush_batch()
  }

  flush_batch()

  list(
    rows_read = read_total,
    rows_analyzed = kept_total,
    rows_skipped = skipped_total,
    chunks = chunk_id
  )
}

bcftools_available <- function() {
  nzchar(Sys.which("bcftools"))
}

#' Run bcftools query on Ubuntu/WSL/Linux for faster VCF field extraction.
bcftools_stream_chunks <- function(
    vcf_path,
    chunk_size = 10000L,
    pass_only = FALSE,
    min_qual = 0,
    max_variants = Inf,
    processor,
    progress_fn = NULL) {

  if (!bcftools_available()) {
    stop("bcftools not found on PATH.")
  }

  query_fmt <- paste(
    "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO\\n",
    sep = ""
  )

  bcftools_cmd <- Sys.which("bcftools")
  if (pass_only || min_qual > 0) {
    view_args <- c("view")
    if (pass_only) view_args <- c(view_args, "-i", shQuote('FILTER="PASS" || FILTER="."'))
    if (min_qual > 0) view_args <- c(view_args, "-i", sprintf("QUAL>=%s", min_qual))
    view_args <- c(view_args, shQuote(vcf_path))
    query_args <- c("query", "-f", shQuote(query_fmt))

    cmd <- paste(
      shQuote(bcftools_cmd), paste(view_args, collapse = " "),
      "|",
      shQuote(bcftools_cmd), paste(query_args, collapse = " "),
      sep = " "
    )
  } else {
    cmd <- paste(
      shQuote(bcftools_cmd), "query", "-f", shQuote(query_fmt), shQuote(vcf_path),
      sep = " "
    )
  }

  pipe <- pipe(cmd, "r")
  on.exit(close(pipe), add = TRUE)

  batch <- list()
  chunk_id <- 0L
  kept_total <- 0L
  read_total <- 0L

  flush_batch <- function() {
    if (length(batch) == 0) return(invisible(NULL))
    chunk_id <<- chunk_id + 1L
    chunk_df <- do.call(rbind, batch)
    rownames(chunk_df) <- NULL
    processor(chunk_df, chunk_id)
    batch <<- list()
    invisible(NULL)
  }

  repeat {
  line <- readLines(pipe, n = 1L, warn = FALSE)
    if (length(line) == 0) break
    read_total <- read_total + 1L
    if (kept_total >= max_variants) next

    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < 7) next

    row <- parse_variant_from_vcf_fields(
      chrom = parts[1],
      pos = parts[2],
      ref = parts[3],
      alt = strsplit(parts[4], ",")[[1]][1],
      qual = suppressWarnings(as.numeric(parts[5])),
      filter = parts[6],
      info = if (length(parts) >= 7) parts[7] else "."
    )

    batch[[length(batch) + 1L]] <- row
    kept_total <- kept_total + 1L

    if (!is.null(progress_fn) && read_total %% 50000L == 0L) {
      progress_fn(read_total, kept_total, 0L)
    }
    if (length(batch) >= chunk_size) flush_batch()
  }
  flush_batch()

  list(rows_read = read_total, rows_analyzed = kept_total, rows_skipped = 0L, chunks = chunk_id)
}

#' Analyze complete VCF — all rows, chunked, results written to CSV on disk.
analyze_complete_vcf <- function(
    vcf_path,
    mode = c("full", "rapid"),
    output_csv = NULL,
    pass_only = FALSE,
    min_qual = 0,
    chunk_size = 10000L,
    use_bcftools = TRUE,
    refs = NULL,
    manual_inputs = list(),
    manual_by_variant = list(),
    clinical_context = NULL,
    pedigree_context = NULL,
    session_id = NA_character_,
    profile_id = DEFAULT_PROFILE_ID,
    write_audit = TRUE,
    gene_filter = character(),
    progress_fn = NULL) {

  mode <- match.arg(mode)
  if (is.null(output_csv)) {
    output_csv <- file.path("logs", paste0("report_", session_id, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
  }
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

  run_metadata <- build_run_metadata(
    vcf_path = vcf_path,
    profile_id = profile_id,
    session_id = session_id,
    mode = mode
  )
  metadata_path <- sub("\\.csv$", ".metadata.json", output_csv)
  write_run_metadata_json(run_metadata, metadata_path)

  first_write <- TRUE
  preview_rows <- list()
  preview_limit <- 1000L
  preview_count <- 0L
  classification_counts <- list()
  rows_classified <- 0L
  rows_gene_skipped <- 0L
  gene_filter <- parse_gene_filter(gene_filter)

  process_chunk <- function(chunk_df, chunk_id) {
    if (length(gene_filter) > 0L) {
      n_before <- nrow(chunk_df)
      chunk_df <- filter_variants_by_genes(chunk_df, gene_filter)
      rows_gene_skipped <<- rows_gene_skipped + (n_before - nrow(chunk_df))
      if (nrow(chunk_df) == 0L) return(invisible(NULL))
    }

    if (!is.null(refs)) {
      chunk_df <- dedupe_variants_by_key(annotate_variants(chunk_df, refs))
    }

    report <- run_acmg_pro_chunk(
      variants_df = chunk_df,
      mode = mode,
      manual_inputs = manual_inputs,
      manual_by_variant = manual_by_variant,
      clinical_context = clinical_context,
      pedigree_context = pedigree_context,
      refs = NULL,
      session_id = session_id,
      profile_id = profile_id,
      run_metadata = run_metadata,
      write_audit = write_audit
    )
    if (nrow(report) > 0) {
      rows_classified <<- rows_classified + nrow(report)
      for (cl in unique(report$classification)) {
        classification_counts[[cl]] <<- (classification_counts[[cl]] %||% 0L) + sum(report$classification == cl)
      }
      utils::write.table(
        report, file = output_csv, sep = ",",
        row.names = FALSE, col.names = first_write, append = !first_write
      )
      first_write <<- FALSE

      display_report <- report
      if (nrow(display_report) > 0 && preview_count < preview_limit) {
        n_take <- min(nrow(display_report), preview_limit - preview_count)
        preview_rows[[length(preview_rows) + 1L]] <<- display_report[seq_len(n_take), , drop = FALSE]
        preview_count <<- preview_count + n_take
      }
    }

    if (!is.null(progress_fn)) {
      progress_fn(detail = sprintf("Chunk %d classified (%d variants)", chunk_id, nrow(chunk_df)))
    }
    invisible(NULL)
  }

  stream_fun <- if (isTRUE(use_bcftools) && bcftools_available()) {
    bcftools_stream_chunks
  } else {
    stream_vcf_chunks
  }

  stats <- tryCatch(
    stream_fun(
      vcf_path = vcf_path,
      chunk_size = chunk_size,
      pass_only = pass_only,
      min_qual = min_qual,
      max_variants = Inf,
      processor = process_chunk,
      progress_fn = progress_fn
    ),
    error = function(e) {
      if (identical(stream_fun, bcftools_stream_chunks)) {
        message("bcftools failed, falling back to R streaming: ", conditionMessage(e))
        stream_vcf_chunks(
          vcf_path = vcf_path,
          chunk_size = chunk_size,
          pass_only = pass_only,
          min_qual = min_qual,
          max_variants = Inf,
          processor = process_chunk,
          progress_fn = progress_fn
        )
      } else {
        stop(e)
      }
    }
  )

  engine <- if (isTRUE(use_bcftools) && bcftools_available()) {
    paste0("bcftools+stream+", ACMG_PRO_ENGINE)
  } else {
    paste0("R-stream+", ACMG_PRO_ENGINE)
  }
  if (!is.null(attr(stats, "fallback"))) engine <- "R-stream (bcftools fallback)"

  preview_df <- if (length(preview_rows) > 0) {
    df <- do.call(rbind, preview_rows)
    rownames(df) <- NULL
    df
  } else {
    empty_report()
  }

  list(
    output_csv = output_csv,
    metadata_path = metadata_path,
    run_metadata = run_metadata,
    preview = preview_df,
    stats = stats,
    classification_counts = classification_counts,
    engine = engine,
    rows_analyzed = stats$rows_analyzed,
    rows_skipped = stats$rows_skipped,
    rows_classified = rows_classified,
    rows_displayed = nrow(preview_df),
    rows_gene_skipped = rows_gene_skipped,
    gene_filter = gene_filter
  )
}
