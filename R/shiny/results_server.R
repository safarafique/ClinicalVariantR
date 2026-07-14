#' Classification summary, status, results tables, evidence, repro, and downloads.
register_results_server <- function(ctx) {
  input <- ctx$input
  output <- ctx$output
  report_a_data <- ctx$report_a_data
  report_a_full <- ctx$report_a_full
  report_b_data <- ctx$report_b_data
  report_b_full <- ctx$report_b_full
  report_c_data <- ctx$report_c_data
  report_c_full <- ctx$report_c_full
  selected_category_a <- ctx$selected_category_a
  selected_category_b <- ctx$selected_category_b
  selected_category_c <- ctx$selected_category_c
  report_a_csv <- ctx$report_a_csv
  report_b_csv <- ctx$report_b_csv
  report_c_csv <- ctx$report_c_csv
  analysis_stats_a <- ctx$analysis_stats_a
  analysis_stats_b <- ctx$analysis_stats_b
  analysis_stats_c <- ctx$analysis_stats_c
  run_metadata_a <- ctx$run_metadata_a
  run_metadata_b <- ctx$run_metadata_b
  run_metadata_c <- ctx$run_metadata_c
  analysis_run_a <- ctx$analysis_run_a
  analysis_run_b <- ctx$analysis_run_b
  analysis_run_c <- ctx$analysis_run_c
  get_selected_row <- ctx$get_selected_row
  selected_variant_card <- ctx$selected_variant_card

  select_category_from_summary <- function(table_input, full_rv, selected_rv, data_rv, stats_rv) {
    sel <- table_input
    full <- full_rv()
    if (is.null(sel) || length(sel) == 0 || is.null(full)) return(invisible(NULL))

    stats <- stats_rv()
    summary_df <- build_classification_summary_df(
      if (!is.null(stats)) stats$classification_counts else NULL,
      full
    )
    category <- resolve_classification_from_row_selection(sel, summary_df)
    if (is.null(category) || !nzchar(category)) return(invisible(NULL))

    selected_rv(category)
    filtered <- filter_report_by_classifications(full, category)
    data_rv(filtered)
    if (!is.null(stats)) {
      stats_rv(modifyList(stats, list(rows_displayed = nrow(filtered))))
    }
    invisible(category)
  }

  observeEvent(input$classification_summary_a_rows_selected, {
    select_category_from_summary(
      input$classification_summary_a_rows_selected,
      report_a_full, selected_category_a, report_a_data, analysis_stats_a
    )
  }, ignoreNULL = TRUE)

  observeEvent(input$classification_summary_b_rows_selected, {
    select_category_from_summary(
      input$classification_summary_b_rows_selected,
      report_b_full, selected_category_b, report_b_data, analysis_stats_b
    )
  }, ignoreNULL = TRUE)

  observeEvent(input$classification_summary_c_rows_selected, {
    select_category_from_summary(
      input$classification_summary_c_rows_selected,
      report_c_full, selected_category_c, report_c_data, analysis_stats_c
    )
  }, ignoreNULL = TRUE)

  output$category_hint_a <- renderUI({
    req(analysis_run_a())
    category_selection_hint_ui(selected_category_a(), report_a_data())
  })

  output$category_hint_b <- renderUI({
    req(analysis_run_b())
    category_selection_hint_ui(selected_category_b(), report_b_data())
  })

  output$category_hint_c <- renderUI({
    req(analysis_run_c())
    category_selection_hint_ui(selected_category_c(), report_c_data())
  })

  output$classification_summary_a <- DT::renderDT({
    req(analysis_run_a())
    stats <- analysis_stats_a()
    summary_df <- build_classification_summary_df(
      stats$classification_counts,
      report_a_full()
    )
    render_classification_summary_table(summary_df)
  })

  output$classification_summary_b <- DT::renderDT({
    req(analysis_run_b())
    stats <- analysis_stats_b()
    summary_df <- build_classification_summary_df(
      stats$classification_counts,
      report_b_full()
    )
    render_classification_summary_table(summary_df)
  })

  output$classification_summary_c <- DT::renderDT({
    req(analysis_run_c())
    stats <- analysis_stats_c()
    summary_df <- build_classification_summary_df(
      stats$classification_counts,
      report_c_full()
    )
    render_classification_summary_table(summary_df)
  })

  output$status_b <- renderUI({
    req(analysis_run_b())
    stats <- analysis_stats_b()
    df <- report_b_data()
    if (is.null(stats) || stats$rows_analyzed == 0) {
      return(div(
        class = "alert alert-warning",
        "No variants loaded. ",
        "Uncheck ", tags$strong("passing-filter only"), " if your VCF uses FILTER='.' (common for GATK/VEP). ",
        "Also confirm the uploaded file is a non-empty VCF."
      ))
    }
    counts <- stats$classification_counts
    path_n <- pathogenic_tier_count(counts)
    count_text <- if (length(counts) > 0) {
      paste(names(counts), unlist(counts), sep = ": ", collapse = " | ")
    } else {
      "See download for classifications"
    }
    classified <- stats$rows_classified %||% stats$rows_analyzed
    div(
      class = "alert alert-success",
      tags$b(format(classified, big.mark = ",")), " variant(s) classified across ",
      tags$b("5 ACMG categories"), ". ",
      tags$strong(sprintf("Pathogenic tier (P + LP): %s", format(path_n, big.mark = ","))), " | ",
      if (is.null(selected_category_b())) {
        tags$em("Click a row in the summary table to drill down.")
      } else {
        tagList(
          "Viewing ", tags$b(selected_category_b()), " (",
          format(nrow(df %||% empty_report()), big.mark = ","), " variant(s))."
        )
      },
      if (stats$rows_skipped > 0) tags$span(class = "text-muted", sprintf("(%s skipped by filters) ", format(stats$rows_skipped, big.mark = ","))),
      tags$br(),
      "Engine: ", tags$code(stats$engine), " | ",
      count_text,
      tags$br(),
      tags$em("Expand ", tags$strong("Variant evidence detail"), " below a selected variant for all 18 automated criteria.")
    )
  })

  output$status_c <- renderUI({
    req(analysis_run_c())
    stats <- analysis_stats_c()
    df <- report_c_data()
    genes <- stats$gene_filter %||% character()

    if (is.null(stats) || (stats$rows_classified %||% 0L) == 0L) {
      return(div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        " No variants found for ", tags$b(format_gene_filter_label(genes)), " in this VCF. ",
        "Check gene symbols match VEP/SnpEff SYMBOL fields (case-insensitive)."
      ))
    }

    counts <- stats$classification_counts
    count_text <- if (length(counts) > 0) {
      paste(names(counts), unlist(counts), sep = ": ", collapse = " | ")
    } else {
      "See download for classifications"
    }
    classified <- stats$rows_classified %||% 0L
    path_n <- pathogenic_tier_count(counts)

    tagList(
      if (path_n > 0L) {
        div(
          class = "alert alert-danger py-2",
          icon("triangle-exclamation"),
          tags$b(path_n), " variant(s) classified as ",
          tags$b("Pathogenic"), " or ", tags$b("Likely Pathogenic"), " in this gene panel."
        )
      },
      div(
        class = "alert alert-success",
        tags$b(format(classified, big.mark = ",")), " variant(s) in ",
        tags$b(format_gene_filter_label(genes)), " classified. ",
        if (is.null(selected_category_c())) {
          tags$em("Click a row in the summary table to drill down.")
        } else {
          tagList(
            "Viewing ", tags$b(selected_category_c()), " (",
            format(nrow(df %||% empty_report()), big.mark = ","), " variant(s))."
          )
        },
        if ((stats$rows_gene_skipped %||% 0L) > 0L) {
          tagList(
            tags$br(),
            tags$span(
              class = "text-muted",
              sprintf(
                "%s other variant(s) in VCF skipped (not in gene list).",
                format(stats$rows_gene_skipped, big.mark = ",")
              )
            )
          )
        },
        tags$br(),
        "Engine: ", tags$code(stats$engine), " | ",
        count_text,
        tags$br(),
        tags$em("Expand ", tags$strong("Variant evidence detail"), " for ACMG criterion breakdown.")
      )
    )
  })

  output$status_a <- renderUI({
    req(analysis_run_a())
    stats <- analysis_stats_a()
    df <- report_a_data()
    if (is.null(stats) || stats$rows_analyzed == 0) {
      return(div(
        class = "alert alert-warning",
        "No variants loaded. ",
        "If ", tags$strong("passing-filter only"), " is enabled, GATK/VEP VCFs often use FILTER='.' instead of PASS — this is now accepted. ",
        "Otherwise check that the uploaded VCF contains variant rows."
      ))
    }
    counts <- stats$classification_counts
    count_text <- if (length(counts) > 0) {
      paste(names(counts), unlist(counts), sep = ": ", collapse = " | ")
    } else {
      "See download for classifications"
    }
    classified <- stats$rows_classified %||% stats$rows_analyzed
    div(
      class = "alert alert-info",
      tags$b(format(classified, big.mark = ",")), " variant(s) classified across ",
      tags$b("5 ACMG categories"), ". ",
      if (is.null(selected_category_a())) {
        tags$em("Click a row in the summary table to drill down.")
      } else {
        tagList(
          "Viewing ", tags$b(selected_category_a()), " (",
          format(nrow(df %||% empty_report()), big.mark = ","), " variant(s))."
        )
      },
      if (stats$rows_skipped > 0) tags$span(class = "text-muted", sprintf(" (%s skipped by filters) ", format(stats$rows_skipped, big.mark = ","))),
      tags$br(),
      "Engine: ", tags$code(stats$engine), " | ",
      count_text,
      tags$br(),
      tags$em("Expand ", tags$strong("Variant evidence detail"), " below a selected variant for all 18 automated criteria.")
    )
  })

  output$results_a <- DT::renderDT({
    req(analysis_run_a())
    if (is.null(selected_category_a())) {
      return(DT::datatable(
        data.frame(Message = "Select a classification in the summary table above."),
        rownames = FALSE, options = list(dom = "t")
      ))
    }
    render_variant_detail_table(report_a_data(), selection = "single")
  })

  output$results_b <- DT::renderDT({
    req(analysis_run_b())
    if (is.null(selected_category_b())) {
      return(DT::datatable(
        data.frame(Message = "Select a classification in the summary table above."),
        rownames = FALSE, options = list(dom = "t")
      ))
    }
    render_variant_detail_table(report_b_data(), selection = "single")
  })

  output$results_c <- DT::renderDT({
    req(analysis_run_c())
    if (is.null(selected_category_c())) {
      return(DT::datatable(
        data.frame(Message = "Select a classification in the summary table above."),
        rownames = FALSE, options = list(dom = "t")
      ))
    }
    render_variant_detail_table(report_c_data(), selection = "single")
  })

  output$selected_variant_a <- renderUI({
    row <- get_selected_row("results_a", report_a_data())
    selected_variant_card(row)
  })

  output$selected_variant_b <- renderUI({
    row <- get_selected_row("results_b", report_b_data())
    selected_variant_card(row)
  })

  output$selected_variant_c <- renderUI({
    row <- get_selected_row("results_c", report_c_data())
    selected_variant_card(row)
  })

  output$evidence_detail_a <- DT::renderDT({
    row <- get_selected_row("results_a", report_a_data())
    if (is.null(row)) return(render_evidence_detail_table(NULL))
    render_evidence_detail_table(parse_evidence_json(row$evidence_json[1]))
  })

  output$evidence_detail_b <- DT::renderDT({
    row <- get_selected_row("results_b", report_b_data())
    if (is.null(row)) return(render_evidence_detail_table(NULL))
    render_evidence_detail_table(parse_evidence_json(row$evidence_json[1]))
  })

  output$evidence_detail_c <- DT::renderDT({
    row <- get_selected_row("results_c", report_c_data())
    if (is.null(row)) return(render_evidence_detail_table(NULL))
    render_evidence_detail_table(parse_evidence_json(row$evidence_json[1]))
  })

  output$repro_a <- renderUI({
    render_reproducibility_card(run_metadata_a())
  })

  output$repro_b <- renderUI({
    render_reproducibility_card(run_metadata_b())
  })

  output$repro_c <- renderUI({
    render_reproducibility_card(run_metadata_c())
  })

  output$download_a <- downloadHandler(
    filename = function() paste0("ACMGamp_Prediction_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      path <- report_a_csv()
      if (!is.null(path) && file.exists(path)) {
        file.copy(path, file)
      } else {
        req(report_a_data())
        write.csv(report_a_data(), file, row.names = FALSE)
      }
    }
  )

  output$download_b <- downloadHandler(
    filename = function() paste0("ACMGamp_Prediction_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      path <- report_b_csv()
      if (!is.null(path) && file.exists(path)) {
        file.copy(path, file)
      } else {
        req(report_b_data())
        write.csv(report_b_data(), file, row.names = FALSE)
      }
    }
  )

  output$download_c <- downloadHandler(
    filename = function() {
      genes <- parse_gene_filter(isolate(input$genes_c))
      tag <- if (length(genes) > 0L) paste(genes[1:min(3, length(genes))], collapse = "_") else "genes"
      paste0("Gene_Panel_", tag, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      path <- report_c_csv()
      if (!is.null(path) && file.exists(path)) {
        file.copy(path, file)
      } else {
        req(report_c_data())
        write.csv(report_c_data(), file, row.names = FALSE)
      }
    }
  )

  invisible(ctx)
}

