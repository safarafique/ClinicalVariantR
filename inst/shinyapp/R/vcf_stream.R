#' Open VCF connection (plain or gzip).
#' @noRd
open_vcf_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "r") else file(path, "r")
}

#' Buffered line reads for large VCF / pipe streams (much faster than n = 1).
#' @noRd
VCF_LINE_BUFFER <- 50000L

#' GATK and many callers use FILTER="." for passing variants; VCF spec also allows PASS.
#' @noRd
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
    num <- scalar_num(val)
    return(if (length(num) == 0 || is.na(num)) NA_real_ else num)
  }
  as.character(val)
}

#' Extract first VEP CSQ consequence token without retaining full CSQ string.
#' @noRd
extract_vep_consequence <- function(csq_value) {
  if (is.na(csq_value) || csq_value == ".") return(NA_character_)
  first <- strsplit(csq_value, ",")[[1]][1]
  if (is.na(first)) return(NA_character_)
  parts <- strsplit(first, "|", fixed = TRUE)[[1]]
  if (length(parts) >= 2) parts[2] else first
}

#' Extract GENE from VEP CSQ (SYMBOL field) or INFO/GENE.
#' @noRd
extract_gene_from_info <- function(info_string) {
  gene <- extract_info_field(info_string, "GENE")
  if (!is.na(gene) && nzchar(gene)) return(gene)
  csq <- extract_info_field(info_string, "CSQ")
  if (is.na(csq)) return(NA_character_)
  first <- strsplit(csq, ",")[[1]][1]
  parts <- strsplit(first, "|", fixed = TRUE)[[1]]
  if (length(parts) >= 4) parts[4] else NA_character_
}

#' Bind list of 1-row (or multi-row) data.frames efficiently.
#' @noRd
rbind_parsed_rows <- function(batch) {
  if (length(batch) == 0L) return(NULL)
  df <- data.table::rbindlist(batch, use.names = TRUE, fill = TRUE)
  data.table::setDF(df)
  rownames(df) <- NULL
  df
}

#' Parse one VCF data line into variant row(s); multi-allelic records expand to one row per ALT.
#' @noRd
parse_vcf_line <- function(line, header_cols) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(fields) < 8) return(NULL)

  chrom <- fields[1]
  pos <- fields[2]
  ref <- fields[4]
  alt_raw <- fields[5]
  qual <- scalar_num(fields[6])
  filter <- if (length(fields) >= 7) fields[7] else "."
  info <- if (length(fields) >= 8) fields[8] else "."

  sample_genotypes <- "{}"
  if (length(fields) >= 9L && !is.null(header_cols) && length(header_cols) >= 9L) {
    format_str <- fields[9]
    sample_names <- header_cols[10:length(header_cols)]
    sample_fields <- if (length(fields) >= 10L) fields[10:length(fields)] else character()
    if (length(sample_names) > 0L && nzchar(format_str) && format_str != ".") {
      gts <- parse_vcf_genotypes(format_str, sample_fields, sample_names)
      sample_genotypes <- serialize_sample_genotypes(gts)
    }
  }

  alts <- split_vcf_alt_alleles(alt_raw)
  rows <- lapply(alts, function(alt) {
    row <- parse_variant_from_vcf_fields(chrom, pos, ref, alt, qual, filter, info)
    row$sample_genotypes <- sample_genotypes
    row
  })

  if (length(rows) == 0L) return(NULL)
  if (length(rows) == 1L) return(rows[[1L]])
  rbind_parsed_rows(rows)
}

append_parsed_vcf_rows <- function(batch, parsed_rows) {
  if (is.null(parsed_rows)) return(batch)
  if (is.data.frame(parsed_rows) && nrow(parsed_rows) > 1L) {
    for (i in seq_len(nrow(parsed_rows))) {
      batch[[length(batch) + 1L]] <- parsed_rows[i, , drop = FALSE]
    }
    return(batch)
  }
  batch[[length(batch) + 1L]] <- parsed_rows
  batch
}

#' Count all variant rows in a VCF (streaming, no memory load).
#' @noRd
count_vcf_variants <- function(vcf_path, pass_only = FALSE, min_qual = 0) {
  con <- open_vcf_connection(vcf_path)
  on.exit(close(con), add = TRUE)

  header_cols <- NULL
  total <- 0L
  pass_n <- 0L

  repeat {
    lines <- readLines(con, n = VCF_LINE_BUFFER, warn = FALSE)
    if (length(lines) == 0L) break

    for (line in lines) {
      if (grepl("^#CHROM\t", line)) {
        header_cols <- strsplit(sub("^#", "", line), "\t")[[1]]
        next
      }
      if (grepl("^#", line)) next

      total <- total + 1L
      if (pass_only || min_qual > 0) {
        rows_df <- parse_vcf_line(line, header_cols)
        if (is.null(rows_df)) next
        n_rows <- if (is.data.frame(rows_df)) nrow(rows_df) else 1L
        for (i in seq_len(n_rows)) {
          row <- if (is.data.frame(rows_df)) rows_df[i, , drop = FALSE] else rows_df
          if (pass_only && !vcf_filter_is_pass(row$filter)) next
          if (min_qual > 0 && (is.na(row$qual) || row$qual < min_qual)) next
          pass_n <- pass_n + 1L
        }
      }
    }
  }

  if (pass_only || min_qual > 0) {
    list(total = total, analyzed = pass_n, skipped = total - pass_n)
  } else {
    list(total = total, analyzed = total, skipped = 0L)
  }
}

#' Stream VCF variants in chunks; call processor(chunk_df, chunk_id) per batch.
#' @noRd
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
  stream_state <- new.env(parent = emptyenv())
  stream_state$batch <- list()
  stream_state$chunk_id <- 0L
  read_total <- 0L
  kept_total <- 0L
  skipped_total <- 0L

  flush_batch <- function() {
    if (length(stream_state$batch) == 0) return(invisible(NULL))
    stream_state$chunk_id <- stream_state$chunk_id + 1L
    chunk_df <- rbind_parsed_rows(stream_state$batch)
    processor(chunk_df, stream_state$chunk_id)
    stream_state$batch <- list()
    invisible(NULL)
  }

  repeat {
    lines <- readLines(con, n = VCF_LINE_BUFFER, warn = FALSE)
    if (length(lines) == 0L) break

    for (line in lines) {
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

      rows_df <- parse_vcf_line(line, header_cols)
      if (is.null(rows_df)) {
        skipped_total <- skipped_total + 1L
        next
      }
      n_rows <- if (is.data.frame(rows_df)) nrow(rows_df) else 1L
      for (i in seq_len(n_rows)) {
        if (kept_total >= max_variants) {
          skipped_total <- skipped_total + (n_rows - i + 1L)
          break
        }
        row <- if (is.data.frame(rows_df)) rows_df[i, , drop = FALSE] else rows_df
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
        stream_state$batch <- append_parsed_vcf_rows(stream_state$batch, row)
        kept_total <- kept_total + 1L
      }

      if (!is.null(progress_fn) && read_total %% 50000L == 0L) {
        progress_fn(read_total, kept_total, skipped_total)
      }

      if (length(stream_state$batch) >= chunk_size) flush_batch()
    }
  }

  flush_batch()

  list(
    rows_read = read_total,
    rows_analyzed = kept_total,
    rows_skipped = skipped_total,
    chunks = stream_state$chunk_id
  )
}

bcftools_available <- function() {
  nzchar(Sys.which("bcftools"))
}

#' Build a single bcftools -i expression (only one -i/-e is allowed).
#' @noRd
bcftools_filter_expression <- function(pass_only = FALSE, min_qual = 0) {
  parts <- character()
  if (isTRUE(pass_only)) {
    parts <- c(parts, '(FILTER="PASS" || FILTER=".")')
  }
  if (isTRUE(min_qual > 0)) {
    parts <- c(parts, sprintf("QUAL>=%s", min_qual))
  }
  if (length(parts) == 0L) return(NULL)
  paste(parts, collapse = " && ")
}

#' Run bcftools query on Ubuntu/WSL/Linux for faster VCF field extraction.
#' Uses a single `bcftools query` (no view|query pipe) to avoid stdin/type errors.
#' @noRd
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

  # bcftools format escapes: pass literal \t and \n through the shell.
  query_fmt <- "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO\\n"
  bcftools_cmd <- Sys.which("bcftools")
  filt <- bcftools_filter_expression(pass_only = pass_only, min_qual = min_qual)

  cmd_parts <- c(shQuote(bcftools_cmd), "query", "-f", shQuote(query_fmt))
  if (!is.null(filt)) {
    cmd_parts <- c(cmd_parts, "-i", shQuote(filt))
  }
  cmd_parts <- c(cmd_parts, shQuote(vcf_path))
  cmd <- paste(cmd_parts, collapse = " ")

  err_file <- tempfile("bcftools_err_")
  on.exit(unlink(err_file), add = TRUE)
  # Capture stderr so silent pipe failures become hard errors and trigger fallback.
  pipe_cmd <- paste(cmd, "2>", shQuote(err_file))
  pipe <- pipe(pipe_cmd, "r")
  on.exit(close(pipe), add = TRUE)

  stream_state <- new.env(parent = emptyenv())
  stream_state$batch <- list()
  stream_state$chunk_id <- 0L
  kept_total <- 0L
  read_total <- 0L

  flush_batch <- function() {
    if (length(stream_state$batch) == 0) return(invisible(NULL))
    stream_state$chunk_id <- stream_state$chunk_id + 1L
    chunk_df <- rbind_parsed_rows(stream_state$batch)
    processor(chunk_df, stream_state$chunk_id)
    stream_state$batch <- list()
    invisible(NULL)
  }

  repeat {
    lines <- readLines(pipe, n = VCF_LINE_BUFFER, warn = FALSE)
    if (length(lines) == 0L) break

    for (line in lines) {
      # Ignore accidental stderr leakage into the pipe
      if (!nzchar(line) || startsWith(line, "Error:") || startsWith(line, "Failed ")) next

      read_total <- read_total + 1L
      if (kept_total >= max_variants) next

      parts <- strsplit(line, "\t")[[1]]
      if (length(parts) < 7) next

      alts <- split_vcf_alt_alleles(parts[4])
      for (alt in alts) {
        if (kept_total >= max_variants) break

        row <- parse_variant_from_vcf_fields(
          chrom = parts[1],
          pos = parts[2],
          ref = parts[3],
          alt = alt,
          qual = scalar_num(parts[5]),
          filter = parts[6],
          info = if (length(parts) >= 7) parts[7] else "."
        )

        stream_state$batch[[length(stream_state$batch) + 1L]] <- row
        kept_total <- kept_total + 1L
      }

      if (!is.null(progress_fn) && read_total %% 50000L == 0L) {
        progress_fn(read_total, kept_total, 0L)
      }
      if (length(stream_state$batch) >= chunk_size) flush_batch()
    }
  }
  flush_batch()

  err_lines <- if (file.exists(err_file)) {
    tryCatch(readLines(err_file, warn = FALSE), error = function(e) character())
  } else {
    character()
  }
  err_text <- paste(err_lines, collapse = "\n")
  if (nzchar(err_text) &&
      grepl("unknown file type|only one -i|Failed|Error:", err_text, ignore.case = TRUE)) {
    # Surface bcftools failures so run_vcf_stream_with_fallback can use R streaming.
    stop(err_text, call. = FALSE)
  }
  if (read_total == 0L) {
    if (nzchar(err_text)) {
      message(sprintf(
        "bcftools query returned 0 rows with current filters (%s).",
        err_text
      ))
    } else {
      message("bcftools query returned 0 rows with current filters.")
    }
  }

  list(rows_read = read_total, rows_analyzed = kept_total, rows_skipped = 0L, chunks = stream_state$chunk_id)
}

new_complete_analysis_state <- function() {
  state <- new.env(parent = emptyenv())
  state$first_write <- TRUE
  state$preview_rows <- list()
  state$preview_count <- 0L
  state$classification_counts <- list()
  state$rows_classified <- 0L
  state$rows_gene_skipped <- 0L
  state$rows_seen <- 0L
  state$total_variants <- NA_integer_
  state
}

#' Instant variant total when a bcftools index exists.
#' Never does a full-file scan (that delayed analysis start).
#' Optional env CLINICALVARIANTR_PRECOUNT=1 enables slow full count.
#' @noRd
count_vcf_variants_fast <- function(vcf_path, pass_only = FALSE, min_qual = 0) {
  if (!isTRUE(pass_only) && !(isTRUE(min_qual > 0)) && bcftools_available()) {
    has_index <- file.exists(paste0(vcf_path, ".tbi")) || file.exists(paste0(vcf_path, ".csi"))
    if (isTRUE(has_index)) {
      out <- tryCatch(
        system2("bcftools", c("index", "-n", vcf_path), stdout = TRUE, stderr = FALSE),
        error = function(e) character()
      )
      n <- suppressWarnings(as.integer(trimws(out[1])))
      if (!is.na(n) && n >= 0L) {
        return(list(total = n, analyzed = n, skipped = 0L, method = "bcftools-index"))
      }
    }
  }

  # Opt-in full scan only (slow on large VCFs).
  if (identical(Sys.getenv("CLINICALVARIANTR_PRECOUNT", unset = ""), "1")) {
    counts <- count_vcf_variants(vcf_path, pass_only = pass_only, min_qual = min_qual)
    counts$method <- "stream"
    return(counts)
  }

  list(total = NA_integer_, analyzed = NA_integer_, skipped = 0L, method = "deferred")
}

#' Quick estimate from a small sample + file size (plain VCF only; skip for .gz).
#' @noRd
estimate_vcf_variants_quick <- function(vcf_path, sample_lines = 20000L) {
  if (grepl("\\.gz$", vcf_path, ignore.case = TRUE)) return(NA_integer_)
  sz <- suppressWarnings(as.numeric(file.info(vcf_path)$size))
  if (is.na(sz) || sz < 1) return(NA_integer_)

  con <- tryCatch(open_vcf_connection(vcf_path), error = function(e) NULL)
  if (is.null(con)) return(NA_integer_)
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  sample_n <- 0L
  bytes_approx <- 0L
  in_variants <- FALSE
  target <- as.integer(sample_lines)

  repeat {
    lines <- readLines(con, n = min(5000L, max(1L, target - sample_n)), warn = FALSE)
    if (length(lines) == 0L) break
    for (line in lines) {
      bytes_approx <- bytes_approx + nchar(line, type = "bytes") + 1L
      if (!in_variants) {
        if (startsWith(line, "#CHROM")) {
          in_variants <- TRUE
          next
        }
        if (startsWith(line, "#")) next
        in_variants <- TRUE
      }
      if (!startsWith(line, "#")) sample_n <- sample_n + 1L
      if (sample_n >= target) break
    }
    if (sample_n >= target) break
    if (bytes_approx > 8e6) break
  }

  if (sample_n < 50L || bytes_approx < 1000L) return(NA_integer_)
  as.integer(max(sample_n, round(as.numeric(sample_n) * (sz / bytes_approx))))
}

#' Push percent-complete to Shiny progress (value 0..1 + detail text).
#' @noRd
emit_analysis_progress <- function(progress_fn, done, total, extra = NULL) {
  if (is.null(progress_fn)) return(invisible(NULL))
  done <- as.integer(done %||% 0L)
  total <- as.integer(total %||% NA_integer_)
  if (!is.na(total) && total > 0L) {
    pct <- max(0, min(0.99, done / total))
    pct_i <- as.integer(round(100 * done / total))
    if (pct_i > 99L && done < total) pct_i <- 99L
    detail <- sprintf(
      "%d%% complete — %s / %s variants",
      pct_i,
      format(done, big.mark = ","),
      format(total, big.mark = ",")
    )
    if (!is.null(extra) && nzchar(extra)) detail <- paste0(detail, " | ", extra)
    progress_fn(
      value = pct,
      message = sprintf("Running... %d%%", pct_i),
      detail = detail
    )
  } else {
    # No pre-count: soft progress so analysis can start immediately.
    soft <- min(0.95, 1 - 1 / (1 + done / 50000))
    detail <- sprintf("%s variants processed", format(done, big.mark = ","))
    if (!is.null(extra) && nzchar(extra)) detail <- paste0(detail, " | ", extra)
    progress_fn(
      value = soft,
      message = sprintf("Running... %s processed", format(done, big.mark = ",")),
      detail = detail
    )
  }
  invisible(NULL)
}

record_analysis_report <- function(report, state, output_csv, preview_limit = 1000L) {
  if (nrow(report) == 0L) return(invisible(NULL))

  state$rows_classified <- state$rows_classified + nrow(report)
  for (cl in unique(report$classification)) {
    state$classification_counts[[cl]] <-
      (state$classification_counts[[cl]] %||% 0L) + sum(report$classification == cl)
  }

  data.table::fwrite(
    report,
    file = output_csv,
    sep = ",",
    append = !state$first_write,
    col.names = state$first_write
  )
  state$first_write <- FALSE

  if (state$preview_count < preview_limit) {
    n_take <- min(nrow(report), preview_limit - state$preview_count)
    state$preview_rows[[length(state$preview_rows) + 1L]] <- report[seq_len(n_take), , drop = FALSE]
    state$preview_count <- state$preview_count + n_take
  }
  invisible(NULL)
}

empty_report <- function() {
  if (exists("REPORT_COLUMNS", inherits = TRUE)) {
    cols <- get("REPORT_COLUMNS", inherits = TRUE)
  } else {
    cols <- character()
  }
  df <- as.data.frame(matrix(ncol = max(length(cols), 0L), nrow = 0))
  if (length(cols) > 0L) names(df) <- cols
  df
}

analysis_preview_df <- function(state) {
  if (length(state$preview_rows) == 0L) return(empty_report())
  df <- rbind_parsed_rows(state$preview_rows)
  df
}

select_vcf_streamer <- function(use_bcftools = TRUE) {
  if (isTRUE(use_bcftools) && bcftools_available()) bcftools_stream_chunks else stream_vcf_chunks
}

run_vcf_stream_with_fallback <- function(stream_fun, vcf_path, chunk_size, pass_only,
                                         min_qual, processor, progress_fn) {
  tryCatch(
    stream_fun(
      vcf_path = vcf_path,
      chunk_size = chunk_size,
      pass_only = pass_only,
      min_qual = min_qual,
      max_variants = Inf,
      processor = processor,
      progress_fn = progress_fn
    ),
    error = function(e) {
      if (!identical(stream_fun, bcftools_stream_chunks)) stop(e)
      message("bcftools failed, falling back to R streaming: ", conditionMessage(e))
      out <- stream_vcf_chunks(
        vcf_path = vcf_path,
        chunk_size = chunk_size,
        pass_only = pass_only,
        min_qual = min_qual,
        max_variants = Inf,
        processor = processor,
        progress_fn = progress_fn
      )
      attr(out, "fallback") <- TRUE
      out
    }
  )
}

# Parallel helpers live in R/parallel_pipeline.R (sourced before this file).
# detect_system_ram_gb / resolve_* / score_chunk_jobs / performance_tuning_ui

#' Analyze complete VCF - all rows, chunked, results written to CSV on disk.
#' Parallelizes ACMG scoring across up to parallel_chunks batches (default 5).
#' @noRd
analyze_complete_vcf <- function(
    vcf_path,
    mode = c("full", "rapid"),
    output_csv = NULL,
    pass_only = FALSE,
    min_qual = 0,
    chunk_size = NULL,
    parallel_chunks = 5L,
    use_bcftools = TRUE,
    refs = NULL,
    manual_inputs = list(),
    manual_by_variant = list(),
    clinical_context = NULL,
    pedigree_context = NULL,
    session_id = NA_character_,
    profile_id = DEFAULT_PROFILE_ID,
    write_audit = FALSE,
    gene_filter = character(),
    progress_fn = NULL) {

  mode <- match.arg(mode)
  parallel_chunks <- resolve_parallel_chunks(parallel_chunks)
  auto_chunk <- is.null(chunk_size)
  if (is.null(output_csv)) {
    output_csv <- file.path("logs", paste0("report_", session_id, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
  }
  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  failed_chunks_log <- sub("\\.csv$", ".failed_chunks.csv", output_csv)

  run_metadata <- build_run_metadata(
    vcf_path = vcf_path,
    profile_id = profile_id,
    session_id = session_id,
    mode = mode
  )
  metadata_path <- sub("\\.csv$", ".metadata.json", output_csv)
  write_run_metadata_json(run_metadata, metadata_path)

  analysis_state <- new_complete_analysis_state()
  analysis_state$n_chunks_ok <- 0L
  analysis_state$n_chunks_failed <- 0L
  gene_filter <- parse_gene_filter(gene_filter)
  chunk_buffer <- new.env(parent = emptyenv())
  chunk_buffer$jobs <- list()

  if (!is.null(progress_fn)) {
    progress_fn(value = 0.01, message = "Running...", detail = NULL)
  }
  counts <- tryCatch(
    count_vcf_variants_fast(vcf_path, pass_only = pass_only, min_qual = min_qual),
    error = function(e) list(total = NA_integer_, analyzed = NA_integer_, skipped = 0L, method = "failed")
  )
  total_variants <- suppressWarnings(as.integer(counts$analyzed %||% NA_integer_))
  if (is.na(total_variants) || total_variants < 1L) {
    total_variants <- suppressWarnings(as.integer(counts$total %||% NA_integer_))
  }
  if ((is.na(total_variants) || total_variants < 1L) &&
      !grepl("\\.gz$", vcf_path, ignore.case = TRUE)) {
    est <- tryCatch(estimate_vcf_variants_quick(vcf_path), error = function(e) NA_integer_)
    if (!is.na(est) && est > 0L) {
      total_variants <- est
      counts$method <- "estimate"
    }
  }
  analysis_state$total_variants <- total_variants

  chunk_size <- resolve_chunk_size(
    chunk_size = if (auto_chunk) NULL else chunk_size,
    parallel_chunks = parallel_chunks,
    auto = auto_chunk,
    n_variants = total_variants
  )
  message(sprintf(
    "ClinicalVariantR performance: parallel_chunks=%d, chunk_size=%s, cores=%d (RAM ~%s GB)",
    parallel_chunks,
    format(chunk_size, big.mark = ","),
    detect_usable_cores(),
    {
      ram <- detect_system_ram_gb()
      if (is.na(ram)) "unknown" else as.character(round(ram, 1))
    }
  ))

  if (!is.null(progress_fn)) {
    if (!is.na(total_variants) && total_variants > 0L) {
      progress_fn(
        value = 0.02,
        message = "Running...",
        detail = sprintf("0%% — 0 / %s variants", format(total_variants, big.mark = ","))
      )
    } else {
      progress_fn(value = 0.02, message = "Running...", detail = NULL)
    }
  }

  flush_chunk_buffer <- function() {
    jobs <- chunk_buffer$jobs
    if (length(jobs) == 0L) return(invisible(NULL))
    chunk_buffer$jobs <- list()

    scored <- score_chunk_jobs(
      jobs = jobs,
      mode = mode,
      manual_inputs = manual_inputs,
      manual_by_variant = manual_by_variant,
      clinical_context = clinical_context,
      pedigree_context = pedigree_context,
      refs = refs,
      session_id = session_id,
      profile_id = profile_id,
      run_metadata = run_metadata,
      workers = parallel_chunks,
      failed_log = failed_chunks_log
    )
    analysis_state$n_chunks_ok <- analysis_state$n_chunks_ok + scored$n_ok
    analysis_state$n_chunks_failed <- analysis_state$n_chunks_failed + scored$n_failed

    for (res in scored$results) {
      if (!isTRUE(res$ok)) next
      report <- res$report
      if (is.null(report) || !is.data.frame(report) || nrow(report) == 0L) next
      if (isTRUE(write_audit)) {
        audit_batch <- build_audit_entries_from_report(report, session_id = session_id)
        append_audit_log(audit_batch)
      }
      record_analysis_report(report, analysis_state, output_csv)
    }

    emit_analysis_progress(
      progress_fn,
      done = analysis_state$rows_seen,
      total = analysis_state$total_variants,
      extra = sprintf(
        "wave ok=%d fail=%d (%d-way parallel)",
        scored$n_ok, scored$n_failed, parallel_chunks
      )
    )
    invisible(NULL)
  }

  process_chunk <- function(chunk_df, chunk_id) {
    n_before <- nrow(chunk_df)
    analysis_state$rows_seen <- analysis_state$rows_seen + n_before

    if (length(gene_filter) > 0L) {
      chunk_df <- filter_variants_by_genes(chunk_df, gene_filter)
      analysis_state$rows_gene_skipped <- analysis_state$rows_gene_skipped + (n_before - nrow(chunk_df))
      if (nrow(chunk_df) == 0L) {
        emit_analysis_progress(
          progress_fn,
          done = analysis_state$rows_seen,
          total = analysis_state$total_variants,
          extra = sprintf("chunk %d (no panel genes in batch)", chunk_id)
        )
        return(invisible(NULL))
      }
    }

    chunk_buffer$jobs[[length(chunk_buffer$jobs) + 1L]] <- list(df = chunk_df, id = chunk_id)
    if (length(chunk_buffer$jobs) >= parallel_chunks) flush_chunk_buffer()
    invisible(NULL)
  }

  # Stream may call progress_fn(read, kept, skipped); map that to %.
  stream_progress_fn <- function(read_total = NULL, kept_total = NULL, skipped_total = NULL,
                                 value = NULL, detail = NULL, message = NULL, ...) {
    if (!is.null(value) || !is.null(detail) || !is.null(message)) {
      if (!is.null(progress_fn)) {
        progress_fn(value = value, detail = detail, message = message, ...)
      }
      return(invisible(NULL))
    }
    if (!is.null(read_total) && !is.na(analysis_state$total_variants) && analysis_state$total_variants > 0L) {
      done <- max(analysis_state$rows_seen, as.integer(read_total))
      emit_analysis_progress(
        progress_fn,
        done = done,
        total = analysis_state$total_variants,
        extra = "streaming"
      )
    }
    invisible(NULL)
  }

  stream_fun <- select_vcf_streamer(use_bcftools)

  stats <- run_vcf_stream_with_fallback(
    stream_fun = stream_fun,
    vcf_path = vcf_path,
    chunk_size = chunk_size,
    pass_only = pass_only,
    min_qual = min_qual,
    processor = process_chunk,
    progress_fn = stream_progress_fn
  )
  flush_chunk_buffer()

  if (!is.null(progress_fn)) {
    progress_fn(
      value = 1,
      message = "Analysis complete — 100%",
      detail = sprintf(
        "100%% complete — %s variants processed",
        format(analysis_state$rows_seen %||% stats$rows_analyzed %||% 0L, big.mark = ",")
      )
    )
  }

  engine <- if (isTRUE(use_bcftools) && bcftools_available()) {
    paste0("bcftools+stream+", ACMG_PRO_ENGINE)
  } else {
    paste0("R-stream+", ACMG_PRO_ENGINE)
  }
  if (!is.null(attr(stats, "fallback"))) engine <- "R-stream (bcftools fallback)"
  if (parallel_chunks > 1L) {
    engine <- paste0(engine, "+parallel", parallel_chunks)
  }

  preview_df <- analysis_preview_df(analysis_state)

  list(
    output_csv = output_csv,
    metadata_path = metadata_path,
    run_metadata = run_metadata,
    preview = preview_df,
    stats = stats,
    classification_counts = analysis_state$classification_counts,
    engine = engine,
    rows_analyzed = stats$rows_analyzed,
    rows_skipped = stats$rows_skipped,
    rows_classified = analysis_state$rows_classified,
    rows_displayed = nrow(preview_df),
    rows_gene_skipped = analysis_state$rows_gene_skipped,
    gene_filter = gene_filter,
    parallel_chunks = parallel_chunks,
    chunk_size = chunk_size,
    total_variants = analysis_state$total_variants,
    n_chunks_ok = analysis_state$n_chunks_ok %||% 0L,
    n_chunks_failed = analysis_state$n_chunks_failed %||% 0L,
    failed_chunks_log = if (isTRUE((analysis_state$n_chunks_failed %||% 0L) > 0L)) failed_chunks_log else NULL
  )
}
