render_vcf_preview_table <- function(preview_df) {
  DT::datatable(
    preview_df,
    options = list(pageLength = 10, scrollX = TRUE, dom = "tip"),
    rownames = FALSE,
    class = "cell-border stripe hover compact"
  )
}

render_validation_table <- function(checks_df) {
  if (is.null(checks_df) || nrow(checks_df) == 0) {
    return(DT::datatable(
      data.frame(Message = "No validation results."),
      rownames = FALSE, options = list(dom = "t")
    ))
  }

  DT::datatable(
    checks_df,
    options = list(pageLength = 15, scrollX = TRUE, dom = "tip"),
    rownames = FALSE,
    class = "cell-border stripe hover compact"
  ) |>
    DT::formatStyle(
      "status",
      backgroundColor = DT::styleEqual(
        c("PASS", "WARN", "FAIL"),
        c("#d4edda", "#fff3cd", "#f8d7da")
      ),
      fontWeight = "bold"
    )
}

readiness_indicator_ui <- function(validation, label = "Analysis readiness", processing = FALSE) {
  if (is.null(validation) && isTRUE(processing)) {
    return(div(
      class = "readiness-badge readiness-neutral",
      tags$span(class = "readiness-dot readiness-dot-grey"),
      tags$span(class = "readiness-text", icon("spinner"), " Validating uploaded file...")
    ))
  }
  if (is.null(validation)) {
    return(div(
      class = "readiness-badge readiness-neutral",
      tags$span(class = "readiness-dot readiness-dot-grey"),
      tags$span(class = "readiness-text", icon("hourglass-half"), " Upload files to check readiness")
    ))
  }

  ready <- isTRUE(validation$can_analyze)
  badge_class <- if (ready) "readiness-ready" else "readiness-not-ready"
  dot_class <- if (ready) "readiness-dot-green" else "readiness-dot-red"
  icon_name <- if (ready) "check-circle" else "times-circle"
  status_text <- if (ready) {
    "Ready for analysis - click Run Analysis"
  } else {
    "Not ready for analysis"
  }

  missing_block <- NULL
  missing_items <- validation$missing_items %||% character()
  if (!ready && length(missing_items) > 0) {
    missing_block <- tags$div(
      class = "readiness-missing mt-2",
      tags$strong("Missing:"),
      tags$ul(class = "mb-0 ps-3", lapply(missing_items, tags$li))
    )
  } else if (!ready) {
    missing_block <- tags$div(
      class = "readiness-missing mt-2",
      tags$strong(validation$summary %||% "Fix validation failures before analysis.")
    )
  }

  div(
    class = paste("readiness-badge", badge_class),
    div(
      class = "d-flex align-items-center gap-2",
      tags$span(class = paste("readiness-dot", dot_class)),
      div(
        div(class = "readiness-label", label),
        div(class = "readiness-text", icon(icon_name), status_text)
      )
    ),
    missing_block
  )
}

validation_summary_ui <- function(validation, file_label = "VCF") {
  if (is.null(validation)) {
    return(div(class = "alert alert-secondary", icon("clipboard-check"), " Upload a file to run requirement checks."))
  }
  readiness_indicator_ui(validation, label = paste(file_label, "requirement check"))
}

render_results_table <- function(report_df, selection = "single") {
  display_cols <- intersect(c(EVIDENCE_DETAIL_COLUMNS, REPORT_COLUMNS), names(report_df))
  dat <- report_df[, unique(display_cols), drop = FALSE]

  tbl <- DT::datatable(
    dat,
    options = list(
      pageLength = 15, scrollX = TRUE, dom = "Bfrtip",
      selection = selection
    ),
    rownames = FALSE,
    class = "cell-border stripe hover compact"
  )

  if ("classification" %in% names(dat)) {
    tbl <- tbl |>
      DT::formatStyle(
        "classification",
        backgroundColor = DT::styleEqual(
          c("Pathogenic", "Likely Pathogenic", "VUS", "Benign", "Likely Benign"),
          c("#dc3545", "#fd7e14", "#ffc107", "#28a745", "#6cbf6c")
        ),
        color = DT::styleEqual(
          c("Pathogenic", "Likely Pathogenic", "VUS", "Benign", "Likely Benign"),
          c("#ffffff", "#ffffff", "#212529", "#ffffff", "#ffffff")
        ),
        fontWeight = "bold"
      )
  }

  round_cols <- intersect(c("gnomad_af", "revel_score", "cadd_score", "spliceai_max", "alphamissense_score"), names(dat))
  if (length(round_cols) > 0) {
    tbl <- tbl |> DT::formatRound(columns = round_cols, digits = 4)
  }

  tbl
}

render_evidence_detail_table <- function(evidence_df) {
  if (is.null(evidence_df) || nrow(evidence_df) == 0) {
    return(DT::datatable(
      data.frame(Message = "Select a variant to view criterion-level evidence."),
      rownames = FALSE, options = list(dom = "t")
    ))
  }
  show <- evidence_df[, intersect(
    c("criterion", "status", "strength", "observed_value", "threshold", "reason", "description"),
    names(evidence_df)
  ), drop = FALSE]
  DT::datatable(
    show,
    options = list(pageLength = 18, scrollX = TRUE, dom = "tip"),
    rownames = FALSE,
    class = "cell-border stripe hover compact"
  ) |>
    DT::formatStyle(
      "status",
      backgroundColor = DT::styleEqual(
        c("Triggered", "Not triggered", "Review required"),
        c("#fff3cd", "#f8f9fa", "#e7f1ff")
      ),
      fontWeight = DT::styleEqual("Triggered", "bold")
    )
}

filter_report_by_classifications <- function(report_df, selected) {
  if (is.null(report_df) || nrow(report_df) == 0) return(report_df)
  if (is.null(selected) || length(selected) == 0) {
    return(report_df[0, , drop = FALSE])
  }
  selected <- if (length(selected) == 1L) selected else selected
  report_df[report_df$classification %in% selected, , drop = FALSE]
}

#' Map DT summary-table row selection to ACMG class.
resolve_classification_from_row_selection <- function(row_selected, summary_df) {
  if (is.null(row_selected) || length(row_selected) == 0) return(NULL)
  row_num <- as.integer(row_selected[1])
  if (row_num < 1L) return(NULL)

  if (!is.null(summary_df) && nrow(summary_df) > 0 && "classification" %in% names(summary_df)) {
    if (row_num <= nrow(summary_df)) {
      return(as.character(summary_df$classification[row_num]))
    }
    return(NULL)
  }

  if (row_num <= length(ACMG_CLASSIFICATIONS)) {
    return(ACMG_CLASSIFICATIONS[row_num])
  }
  NULL
}

load_full_analysis_report <- function(csv_path, preview_fallback = NULL) {
  if (!is.null(csv_path) && nzchar(csv_path) && file.exists(csv_path)) {
    df <- tryCatch(
      utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    if (!is.null(df) && nrow(df) > 0) return(df)
  }
  preview_fallback
}

build_classification_summary_df <- function(classification_counts, full_report = NULL) {
  counts <- setNames(rep(0L, length(ACMG_CLASSIFICATIONS)), ACMG_CLASSIFICATIONS)
  if (!is.null(classification_counts) && length(classification_counts) > 0) {
    for (nm in names(classification_counts)) {
      if (nm %in% names(counts)) counts[[nm]] <- as.integer(classification_counts[[nm]])
    }
  } else if (!is.null(full_report) && nrow(full_report) > 0 && "classification" %in% names(full_report)) {
    tab <- table(full_report$classification)
    for (nm in names(tab)) {
      if (nm %in% names(counts)) counts[[nm]] <- as.integer(tab[[nm]])
    }
  }
  data.frame(
    classification = ACMG_CLASSIFICATIONS,
    variant_count = as.integer(unlist(counts[ACMG_CLASSIFICATIONS], use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}

classification_colors <- function() {
  list(
    bg = c("#dc3545", "#fd7e14", "#ffc107", "#28a745", "#6cbf6c"),
    fg = c("#ffffff", "#ffffff", "#212529", "#ffffff", "#ffffff")
  )
}

render_classification_summary_table <- function(summary_df, table_id = "classification_summary") {
  cols <- classification_colors()
  DT::datatable(
    summary_df,
    options = list(
      pageLength = 10, dom = "t", scrollX = TRUE,
      selection = "single",
      ordering = FALSE
    ),
    rownames = FALSE,
    class = "cell-border stripe hover compact classification-summary-table",
    selection = "single"
  ) |>
    DT::formatStyle(
      "classification",
      backgroundColor = DT::styleEqual(ACMG_CLASSIFICATIONS, cols$bg),
      color = DT::styleEqual(ACMG_CLASSIFICATIONS, cols$fg),
      fontWeight = "bold"
    ) |>
    DT::formatStyle(
      columns = colnames(summary_df),
      cursor = "pointer"
    )
}

render_variant_detail_table <- function(report_df, selection = "single") {
  if (is.null(report_df) || nrow(report_df) == 0) {
    return(DT::datatable(
      data.frame(Message = "No variants in this classification."),
      rownames = FALSE, options = list(dom = "t")
    ))
  }
  display_cols <- intersect(VARIANT_DETAIL_COLUMNS, names(report_df))
  dat <- report_df[, unique(display_cols), drop = FALSE]

  tbl <- DT::datatable(
    dat,
    options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip", selection = selection),
    rownames = FALSE,
    class = "cell-border stripe hover compact"
  )

  if ("classification" %in% names(dat)) {
    cols <- classification_colors()
    tbl <- tbl |>
      DT::formatStyle(
        "classification",
        backgroundColor = DT::styleEqual(ACMG_CLASSIFICATIONS, cols$bg),
        color = DT::styleEqual(ACMG_CLASSIFICATIONS, cols$fg),
        fontWeight = "bold"
      )
  }

  round_cols <- intersect(c("gnomad_af", "revel_score", "cadd_score"), names(dat))
  if (length(round_cols) > 0) {
    tbl <- tbl |> DT::formatRound(columns = round_cols, digits = 4)
  }
  tbl
}

format_analysis_notification <- function(result) {
  classified <- result$rows_classified %||% result$rows_analyzed %||% 0L
  displayed <- result$rows_displayed %||% nrow(result$preview %||% empty_report())
  if (classified == 0L) {
    return(paste0(
      "No variants loaded. Preview counts rows in the VCF; analysis applies filters first. ",
      "Uncheck passing-filter only if your VCF uses FILTER='.' (common for GATK/VEP)."
    ))
  }
  msg <- sprintf(
    "Analysis complete: %s variant(s) classified (%s).",
    format(classified, big.mark = ","),
    result$engine
  )
  paste0(msg, " Click a classification in the summary table to view variant details.")
}

category_selection_hint_ui <- function(selected_category, report_df) {
  if (is.null(selected_category) || !nzchar(selected_category)) {
    return(div(
      class = "alert alert-info mb-2 py-2",
      icon("hand-pointer"),
      " Click a classification row in the summary table to open variant details."
    ))
  }
  n <- if (is.null(report_df)) 0L else nrow(report_df)
  div(
    class = "alert alert-light border mb-2 py-2",
    tags$b(selected_category), " - ",
    format(n, big.mark = ","), " variant(s). ",
    tags$span(class = "text-muted", "Select a variant row below to view all 18 automated criteria.")
  )
}

render_reproducibility_card <- function(metadata) {
  if (is.null(metadata)) {
    return(div(class = "alert alert-secondary", "Run analysis to view reproducibility metadata."))
  }
  tags$div(
    class = "repro-card",
    tags$table(
      class = "table table-sm",
      tags$tbody(
        tags$tr(tags$td(tags$b("App version")), tags$td(metadata$app_version %||% "")),
        tags$tr(tags$td(tags$b("Engine")), tags$td(metadata$engine %||% "")),
        tags$tr(tags$td(tags$b("ACMG guideline")), tags$td(metadata$acmg_guideline_version %||% "")),
        tags$tr(tags$td(tags$b("Disease profile")), tags$td(metadata$disease_profile_name %||% metadata$disease_profile %||% "")),
        tags$tr(tags$td(tags$b("Input VCF")), tags$td(metadata$input_vcf %||% "")),
        tags$tr(tags$td(tags$b("VCF checksum")), tags$td(tags$code(metadata$input_vcf_checksum %||% "NA"))),
        tags$tr(tags$td(tags$b("R version")), tags$td(metadata$r_version %||% "")),
        tags$tr(tags$td(tags$b("Run timestamp")), tags$td(metadata$run_timestamp %||% "")),
        tags$tr(tags$td(tags$b("Thresholds")), tags$td(tags$small(metadata$thresholds_snapshot %||% "")))
      )
    )
  )
}

profile_choices <- function() {
  tbl <- load_disease_profiles_table()
  if (nrow(tbl) == 0) return(c("General germline" = "general_germline"))
  stats::setNames(tbl$profile_id, tbl$profile_name)
}

disease_profile_help_ui <- function() {
  helpText(
    class = "text-muted small",
    tags$strong("Disease profile"),
    " adjusts population-frequency cutoffs (mainly BS1). ",
    "BA1 (AF > 5%) and PM2 (rare/absent) are the same across all presets in this config. ",
    "Stricter profiles (e.g. hematologic, neurological) use a lower BS1 threshold (0.5% vs 1%)."
  )
}

pdf_truncate_cell <- function(x, max_chars = 120L) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(val) {
    if (!nzchar(val)) return("")
    if (nchar(val) <= max_chars) return(val)
    paste0(substr(val, 1L, max_chars - 3L), "...")
  }, FUN.VALUE = character(1L))
}

pdf_format_numeric_cols <- function(df) {
  num_cols <- intersect(
    c("gnomad_af", "revel_score", "cadd_score", "spliceai_max", "alphamissense_score", "confidence_score"),
    names(df)
  )
  for (col in num_cols) {
    df[[col]] <- vapply(df[[col]], function(v) {
      if (is.null(v) || length(v) == 0L || is.na(v)) return("")
      format(round(as.numeric(v), 4), scientific = FALSE, trim = TRUE)
    }, FUN.VALUE = character(1L))
  }
  df
}

pdf_column_widths <- function(col_names) {
  weights <- vapply(col_names, function(nm) {
    switch(nm,
      variant_id = 1.6,
      gene = 0.7,
      classification = 0.9,
      confidence_score = 0.55,
      confidence_label = 0.65,
      criteria_met = 0.9,
      gnomad_af = 0.55,
      disease_profile = 0.8,
      evidence_summary = 2.4,
      1
    )
  }, numeric(1L))
  weights / sum(weights)
}

pdf_classification_color <- function(classification) {
  switch(as.character(classification),
    "Pathogenic" = "#b02a37",
    "Likely Pathogenic" = "#ca6510",
    "VUS" = "#997404",
    "Benign" = "#1e7e34",
    "Likely Benign" = "#2f8f44",
    "#212529"
  )
}

pdf_export_settings <- function(n_rows) {
  if (n_rows <= 500L) {
    list(
      rows_per_page = 22L,
      export_cols = PDF_EXPORT_COLUMNS,
      evidence_chars = 160L,
      cell_chars = 80L
    )
  } else if (n_rows <= 5000L) {
    list(
      rows_per_page = 40L,
      export_cols = c(
        "variant_id", "gene", "classification", "criteria_met",
        "gnomad_af", "confidence_score"
      ),
      evidence_chars = 0L,
      cell_chars = 60L
    )
  } else {
    list(
      rows_per_page = 50L,
      export_cols = c("variant_id", "gene", "classification", "criteria_met", "gnomad_af"),
      evidence_chars = 0L,
      cell_chars = 50L
    )
  }
}

pdf_format_chunk <- function(chunk, evidence_chars = 160L, cell_chars = 80L) {
  if (is.null(chunk) || nrow(chunk) == 0L) return(chunk)
  chunk <- pdf_format_numeric_cols(chunk)
  if ("evidence_summary" %in% names(chunk) && evidence_chars > 0L) {
    chunk$evidence_summary <- pdf_truncate_cell(chunk$evidence_summary, evidence_chars)
  }
  skip_cols <- c("evidence_summary", "gnomad_af", "confidence_score")
  for (col in setdiff(names(chunk), skip_cols)) {
    chunk[[col]] <- pdf_truncate_cell(chunk[[col]], cell_chars)
  }
  chunk
}

pdf_draw_page_header <- function(
    title, subtitle, generated_at, page, n_pages, start, end, total) {
  graphics::text(0.5, 0.97, title, cex = 1.35, font = 2)
  if (nzchar(subtitle)) {
    graphics::text(0.02, 0.935, subtitle, cex = 0.85, adj = c(0, 0.5), col = "#495057")
  }
  footer <- sprintf(
    "Generated %s | Page %d of %d | Variants %s-%s of %s",
    generated_at, page, n_pages,
    format(start, big.mark = ","),
    format(end, big.mark = ","),
    format(total, big.mark = ",")
  )
  graphics::text(0.98, 0.935, footer, cex = 0.75, adj = c(1, 0.5), col = "#6c757d")
}

pdf_draw_table_page <- function(df, top, bottom, left, right) {
  ncols <- ncol(df)
  nrows <- nrow(df)
  if (ncols == 0L || nrows == 0L) return(invisible(NULL))

  widths <- pdf_column_widths(names(df))
  x_edges <- left + c(0, cumsum(widths) * (right - left))
  header_h <- 0.035
  body_top <- top - header_h
  row_h <- (body_top - bottom) / max(nrows, 1L)

  graphics::rect(
    left, body_top, right, top,
    col = "#e9ecef", border = "#adb5bd"
  )

  col_labels <- gsub("_", " ", names(df), fixed = TRUE)
  for (j in seq_len(ncols)) {
    graphics::text(
      (x_edges[j] + x_edges[j + 1L]) / 2,
      top - header_h / 2,
      col_labels[j],
      cex = 0.55,
      font = 2
    )
  }

  for (i in seq_len(nrows)) {
    y_top <- body_top - (i - 1L) * row_h
    y_mid <- body_top - (i - 0.5) * row_h
    y_bot <- body_top - i * row_h
    if (i %% 2L == 0L) {
      graphics::rect(left, y_bot, right, y_top, col = "#f8f9fa", border = NA)
    }
    for (j in seq_len(ncols)) {
      val <- as.character(df[i, j])
      if (is.na(val)) val <- ""
      col_nm <- names(df)[j]
      graphics::text(
        (x_edges[j] + x_edges[j + 1L]) / 2,
        y_mid,
        val,
        cex = 0.48,
        font = if (col_nm == "classification") 2L else 1L,
        col = if (col_nm == "classification") pdf_classification_color(val) else "#212529"
      )
    }
    graphics::segments(left, y_bot, right, y_bot, col = "#dee2e6", lwd = 0.5)
  }

  graphics::rect(left, bottom, right, top, border = "#adb5bd")
  invisible(NULL)
}

resolve_explorer_export_df <- function(
    report_df,
    scope = c("all", "filtered"),
    filtered_row_indices = NULL) {

  scope <- match.arg(scope)
  if (is.null(report_df) || nrow(report_df) == 0L) {
    stop("No results available to export.")
  }
  if (scope == "filtered") {
    if (is.null(filtered_row_indices) || length(filtered_row_indices) == 0L) {
      stop("No rows match the current table search/filter.")
    }
    idx <- as.integer(filtered_row_indices) + 1L
    idx <- unique(idx[idx >= 1L & idx <= nrow(report_df)])
    if (length(idx) == 0L) {
      stop("No rows match the current table search/filter.")
    }
    return(report_df[idx, , drop = FALSE])
  }
  report_df
}

export_evidence_summary_pdf <- function(
    summary_df,
    file,
    title = "ClinicalVariantR Evidence Report - Summary",
    subtitle = "",
    total_variants = NA_integer_) {

  if (is.null(summary_df) || nrow(summary_df) == 0L) {
    stop("No summary data available to export.")
  }

  generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M UTC")
  grDevices::pdf(file, width = 11.69, height = 8.27, paper = "a4r", onefile = TRUE, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))

  pdf_draw_page_header(
    title, subtitle, generated_at, 1L, 1L, 1L,
    sum(summary_df$variant_count, na.rm = TRUE),
    if (is.na(total_variants)) sum(summary_df$variant_count, na.rm = TRUE) else total_variants
  )

  graphics::text(0.02, 0.88, "Classification summary", cex = 1.0, font = 2, adj = c(0, 1))

  tbl <- summary_df
  tbl$variant_count <- format(tbl$variant_count, big.mark = ",")
  nrows <- nrow(tbl)
  top <- 0.82
  bottom <- 0.45
  left <- 0.15
  right <- 0.55
  row_h <- (top - bottom) / max(nrows, 1L)

  graphics::rect(left, bottom, right, top, border = "#adb5bd")
  graphics::rect(left, top - row_h, right, top, col = "#e9ecef", border = "#adb5bd")
  graphics::text((left + right) / 2, top - row_h / 2, "Classification", cex = 0.75, font = 2)
  graphics::text(right + 0.12, top - row_h / 2, "Variant count", cex = 0.75, font = 2)

  for (i in seq_len(nrows)) {
    y_mid <- top - row_h * (i - 0.5) - row_h
    cls <- as.character(tbl$classification[i])
    if (i %% 2L == 0L) {
      graphics::rect(left, y_mid - row_h / 2, right + 0.2, y_mid + row_h / 2, col = "#f8f9fa", border = NA)
    }
    graphics::text(left + 0.02, y_mid, cls, adj = c(0, 0.5), cex = 0.7, font = 2,
                   col = pdf_classification_color(cls))
    graphics::text(right + 0.12, y_mid, tbl$variant_count[i], adj = c(0.5, 0.5), cex = 0.7)
    graphics::segments(left, y_mid - row_h / 2, right + 0.2, y_mid - row_h / 2, col = "#dee2e6")
  }

  if (!is.na(total_variants)) {
    graphics::text(
      0.02, 0.38,
      sprintf("Total variants analyzed: %s", format(total_variants, big.mark = ",")),
      adj = c(0, 0.5), cex = 0.85, font = 2
    )
  }
  graphics::text(
    0.02, 0.32,
    "Use CSV export for the complete variant table. Detailed PDF includes paginated variant rows.",
    adj = c(0, 0.5), cex = 0.75, col = "#6c757d"
  )

  invisible(file)
}

#' Write the Evidence Explorer results table to a paginated landscape PDF.
export_evidence_report_pdf <- function(
    report_df,
    file,
    title = "ClinicalVariantR Evidence Report",
    subtitle = "",
    progress_callback = NULL) {

  if (is.null(report_df) || nrow(report_df) == 0L) {
    stop("No results available to export.")
  }

  settings <- pdf_export_settings(nrow(report_df))
  cols <- intersect(settings$export_cols, names(report_df))
  if (length(cols) == 0L) {
    cols <- intersect(c(EVIDENCE_DETAIL_COLUMNS, REPORT_COLUMNS), names(report_df))
  }

  rows_per_page <- settings$rows_per_page
  n_total <- nrow(report_df)
  n_pages <- max(1L, as.integer(ceiling(n_total / rows_per_page)))
  generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M UTC")

  grDevices::pdf(file, width = 11.69, height = 8.27, paper = "a4r", onefile = TRUE, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (page in seq_len(n_pages)) {
    if (is.function(progress_callback)) {
      progress_callback(page / n_pages)
    }
    start <- (page - 1L) * rows_per_page + 1L
    end <- min(page * rows_per_page, n_total)
    chunk <- report_df[start:end, cols, drop = FALSE]
    chunk <- pdf_format_chunk(
      chunk,
      evidence_chars = settings$evidence_chars,
      cell_chars = settings$cell_chars
    )

    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))

    pdf_draw_page_header(title, subtitle, generated_at, page, n_pages, start, end, n_total)
    pdf_draw_table_page(chunk, top = 0.88, bottom = 0.05, left = 0.02, right = 0.98)
  }

  invisible(file)
}

render_reference_readiness_banner <- function() {
  chk <- check_reference_readiness()
  if (isTRUE(chk$ready)) {
    return(div(
      class = "alert alert-success py-2 mt-3",
      icon("database"),
      " Reference databases loaded (gnomAD, ClinVar, REVEL)."
    ))
  }
  items <- vapply(chk$status, function(x) {
    sprintf("%s: %s (%s rows)", basename(x$path), x$label, format(x$rows, big.mark = ","))
  }, character(1))
  div(
    class = "alert alert-warning py-2 mt-3",
    tags$b("Reference data incomplete. "),
    "Predictions rely heavily on VCF annotation until full gnomAD/ClinVar/REVEL files are installed.",
    tags$br(),
    tags$small(paste(items, collapse = " | "))
  )
}

explorer_source_label <- function(source_id) {
  switch(source_id,
    a = "Group A - Clinical prediction",
    b = "Group B - Automated prediction",
    c = "Group C - Gene panel prediction",
    "Evidence Explorer"
  )
}
