#' Open VCF connection (plain or gzip).
open_vcf_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "r") else file(path, "r")
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

#' Parse one VCF data line into a variant row.
parse_vcf_line <- function(line, header_cols) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(fields) < 5) return(NULL)

  col <- function(name, fallback) {
    idx <- match(name, header_cols)
    if (is.na(idx) || idx > length(fields)) fallback else fields[idx]
  }

  alt_raw <- col("ALT", ".")
  alt <- strsplit(alt_raw, ",")[[1]][1]
  info <- col("INFO", ".")
  chrom <- fields[1]
  pos <- suppressWarnings(as.integer(fields[2]))
  ref <- fields[4]

  consequence <- extract_info_field(info, "CSQ")
  if (!is.na(consequence) && grepl("\\|", consequence, fixed = TRUE)) {
    consequence <- extract_vep_consequence(consequence)
  }

  data.frame(
    variant_id = paste(chrom, pos, ref, alt, sep = ":"),
    chrom = chrom,
    pos = pos,
    ref = ref,
    alt = alt,
    gene = extract_gene_from_info(info),
    consequence = consequence,
    AF = extract_info_field(info, "AF"),
    REVEL = extract_info_field(info, "REVEL"),
    ClinVar = NA_character_,
    qual = suppressWarnings(as.numeric(col("QUAL", "."))),
    filter = col("FILTER", "."),
    stringsAsFactors = FALSE
  )
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
      if (pass_only && row$filter != "PASS") next
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
    if (pass_only && row$filter != "PASS") {
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
    "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t",
    "%INFO/AF\\t%INFO/REVEL\\t%INFO/CSQ\\n",
    sep = ""
  )

  bcftools_cmd <- Sys.which("bcftools")
  if (pass_only || min_qual > 0) {
    view_args <- c("view")
    if (pass_only) view_args <- c(view_args, "-f", "PASS")
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
    if (length(parts) < 9) next

    alt <- strsplit(parts[4], ",")[[1]][1]
    csq <- if (parts[9] == ".") NA_character_ else extract_vep_consequence(parts[9])

    row <- data.frame(
      variant_id = paste(parts[1], parts[2], parts[3], alt, sep = ":"),
      chrom = parts[1],
      pos = as.integer(parts[2]),
      ref = parts[3],
      alt = alt,
      gene = if (parts[9] == ".") NA_character_ else extract_gene_from_info(paste0("CSQ=", parts[9])),
      consequence = csq,
      AF = suppressWarnings(as.numeric(parts[7])),
      REVEL = suppressWarnings(as.numeric(parts[8])),
      ClinVar = NA_character_,
      stringsAsFactors = FALSE
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
    clinical_context = NULL,
    pedigree_context = NULL,
    session_id = NA_character_,
    filter_pathogenic_only = FALSE,
    progress_fn = NULL) {

  mode <- match.arg(mode)
  if (is.null(output_csv)) {
    output_csv <- file.path("logs", paste0("report_", session_id, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
  }
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

  first_write <- TRUE
  preview_rows <- list()
  preview_limit <- 1000L
  preview_count <- 0L
  classification_counts <- list()

  process_chunk <- function(chunk_df, chunk_id) {
    if (!is.null(refs)) {
      chunk_df <- annotate_variants(chunk_df, refs)
    }

    result <- run_pipeline(
      variants_df = chunk_df,
      mode = mode,
      manual_inputs = manual_inputs,
      clinical_context = clinical_context,
      pedigree_context = pedigree_context,
      session_id = session_id
    )

    report <- result$report
    if (filter_pathogenic_only && nrow(report) > 0) {
      report <- report[report$classification %in% c("Pathogenic", "Likely Pathogenic"), , drop = FALSE]
    }

    if (nrow(report) > 0) {
      for (cl in unique(report$classification)) {
        classification_counts[[cl]] <<- (classification_counts[[cl]] %||% 0L) + sum(report$classification == cl)
      }
      utils::write.table(
        report, file = output_csv, sep = ",",
        row.names = FALSE, col.names = first_write, append = !first_write
      )
      first_write <<- FALSE

      if (preview_count < preview_limit) {
        n_take <- min(nrow(report), preview_limit - preview_count)
        preview_rows[[length(preview_rows) + 1L]] <<- report[seq_len(n_take), , drop = FALSE]
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

  engine <- if (isTRUE(use_bcftools) && bcftools_available()) "bcftools+stream" else "R-stream"
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
    preview = preview_df,
    stats = stats,
    classification_counts = classification_counts,
    engine = engine,
    rows_analyzed = stats$rows_analyzed,
    rows_skipped = stats$rows_skipped
  )
}
