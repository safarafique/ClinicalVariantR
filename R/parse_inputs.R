#' Parse uploaded VCF into a variant data.frame (small files / legacy).
#' For complete large VCF analysis use analyze_complete_vcf() instead.
parse_vcf_upload <- function(vcf_path, max_rows = Inf) {
  parse_vcf_fallback(vcf_path, max_rows = max_rows)
}

parse_vcf_fallback <- function(vcf_path, max_rows = Inf) {
  con <- open_vcf_connection(vcf_path)
  on.exit(close(con), add = TRUE)

  header_cols <- NULL
  batch <- list()

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0) break
    if (grepl("^#CHROM\t", line)) {
      header_cols <- strsplit(sub("^#", "", line), "\t")[[1]]
      next
    }
    if (grepl("^#", line)) next
    rows_df <- parse_vcf_line(line, header_cols)
    if (!is.null(rows_df)) {
      if (is.data.frame(rows_df)) {
        for (i in seq_len(nrow(rows_df))) {
          row <- rows_df[i, , drop = FALSE]
          row$qual <- NULL
          row$filter <- NULL
          batch[[length(batch) + 1L]] <- row
        }
      } else {
        rows_df$qual <- NULL
        rows_df$filter <- NULL
        batch[[length(batch) + 1L]] <- rows_df
      }
    }
    if (length(batch) >= max_rows) break
  }

  if (length(batch) == 0) {
    return(data.frame(
      variant_id = character(), chrom = character(), pos = integer(),
      ref = character(), alt = character(), gene = character(),
      consequence = character(), stringsAsFactors = FALSE
    ))
  }

  df <- do.call(rbind, batch)
  rownames(df) <- NULL
  if (!"sample_genotypes" %in% names(df)) df$sample_genotypes <- "{}"
  df
}

#' Build a lightweight VCF preview (header + first N variant rows).
preview_vcf <- function(vcf_path, max_rows = 50L, max_header_lines = 100000L) {
  if (!file.exists(vcf_path)) {
    stop("VCF file not found.")
  }

  con <- open_vcf_connection(vcf_path)
  on.exit(close(con), add = TRUE)

  header_cols <- NULL
  preview_lines <- character()
  total_variants <- 0L
  has_more <- FALSE
  line_count <- 0L

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0) break
    line_count <- line_count + 1L

    if (is.null(header_cols) && line_count > max_header_lines) {
      stop(sprintf(
        "Invalid VCF: no #CHROM header found within the first %s lines.",
        format(max_header_lines, big.mark = ",")
      ))
    }

    if (grepl("^#CHROM\t", line)) {
      header_cols <- strsplit(sub("^#", "", line), "\t")[[1]]
      next
    }
    if (grepl("^#", line)) next

    total_variants <- total_variants + 1L
    if (length(preview_lines) < max_rows) {
      preview_lines[length(preview_lines) + 1L] <- line
    } else {
      has_more <- TRUE
      break
    }
  }

  if (is.null(header_cols)) {
    stop("Invalid VCF: missing #CHROM header line.")
  }

  total_display <- if (has_more) paste0(max_rows, "+") else as.character(total_variants)

  preview_df <- if (length(preview_lines) == 0) {
    data.frame(
      chrom = character(), pos = character(), ref = character(), alt = character(),
      qual = character(), filter = character(), info = character(),
      stringsAsFactors = FALSE
    )
  } else {
    col_idx <- function(name, fallback) {
      idx <- match(name, header_cols)
      if (is.na(idx)) fallback else idx
    }
    mat <- strsplit(preview_lines, "\t")
    data.frame(
      chrom = vapply(mat, function(x) x[1], FUN.VALUE = character(1)),
      pos = vapply(mat, function(x) x[2], FUN.VALUE = character(1)),
      ref = vapply(mat, function(x) x[4], FUN.VALUE = character(1)),
      alt = vapply(mat, function(x) {
        alleles <- strsplit(x[5], ",")[[1]]
        if (length(alleles) == 0) "." else alleles[1]
      }, FUN.VALUE = character(1)),
      qual = vapply(mat, function(x) x[col_idx("QUAL", 6)], FUN.VALUE = character(1)),
      filter = vapply(mat, function(x) x[col_idx("FILTER", 7)], FUN.VALUE = character(1)),
      info = vapply(mat, function(x) x[col_idx("INFO", 8)], FUN.VALUE = character(1)),
      stringsAsFactors = FALSE
    )
  }

  info_idx <- match("INFO", header_cols)
  if (!is.na(info_idx) && nrow(preview_df) > 0) {
    mat <- strsplit(preview_lines, "\t")
    parsed <- lapply(mat, function(x) {
      if (length(x) < info_idx) return(NULL)
      parse_variant_from_vcf_fields(
        chrom = x[1], pos = x[2], ref = x[4],
        alt = strsplit(x[5], ",")[[1]][1],
        qual = suppressWarnings(as.numeric(x[6])),
        filter = if (length(x) >= 7) x[7] else ".",
        info = x[info_idx]
      )
    })
    parsed <- Filter(Negate(is.null), parsed)
    if (length(parsed) > 0) {
      p_df <- do.call(rbind, parsed)
      rownames(p_df) <- NULL
      preview_df$gene <- p_df$gene
      preview_df$all_genes <- p_df$all_genes
      preview_df$consequence <- p_df$consequence
      preview_df$annotation_source <- p_df$annotation_source
      preview_df$population_af <- p_df$population_af
    }
  }

  list(
    preview = preview_df[, intersect(
      c("chrom", "pos", "ref", "alt", "gene", "all_genes", "consequence", "annotation_source",
        "population_af", "qual", "filter"),
      names(preview_df)
    ), drop = FALSE],
    total_variants = total_variants,
    total_display = total_display,
    has_more = has_more,
    preview_rows = nrow(preview_df),
    file_size_mb = round(file.info(vcf_path)$size / 1024^2, 2)
  )
}

#' Parse clinical logs CSV for Group A manual criteria context.
parse_clinical_logs <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  required <- c("sample_id", "phenotype", "cml_phase", "tki_response")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Clinical logs missing columns: ", paste(missing, collapse = ", "))
  }
  optional <- c("hpo_terms", "hpo_ids", "disease", "clinical_summary", "omim_id")
  for (col in optional) {
    if (!col %in% names(df)) df[[col]] <- NA_character_
  }
  as.data.frame(df)
}

#' Parse pedigree CSV for segregation / de novo context.
parse_pedigree <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  required <- c("sample_id", "relation", "affected_status")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Pedigree data missing columns: ", paste(missing, collapse = ", "))
  }
  as.data.frame(df)
}
