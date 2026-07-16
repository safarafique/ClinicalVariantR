#' Evidence explorer tab - full report browse, PDF/CSV export.
register_explorer_server <- function(ctx) {
  input <- ctx$input
  output <- ctx$output
  report_a_data <- ctx$report_a_data
  report_a_full <- ctx$report_a_full
  report_b_data <- ctx$report_b_data
  report_b_full <- ctx$report_b_full
  report_c_data <- ctx$report_c_data
  report_c_full <- ctx$report_c_full
  report_a_csv <- ctx$report_a_csv
  report_b_csv <- ctx$report_b_csv
  report_c_csv <- ctx$report_c_csv
  analysis_stats_a <- ctx$analysis_stats_a
  analysis_stats_b <- ctx$analysis_stats_b
  analysis_stats_c <- ctx$analysis_stats_c
  get_selected_row <- ctx$get_selected_row
  selected_variant_card <- ctx$selected_variant_card

  explorer_report <- reactive({
    if (identical(input$explorer_source, "a")) {
      report_a_full() %||% report_a_data()
    } else if (identical(input$explorer_source, "c")) {
      report_c_full() %||% report_c_data()
    } else {
      report_b_full() %||% report_b_data()
    }
  })

  explorer_stats <- reactive({
    if (identical(input$explorer_source, "a")) {
      analysis_stats_a()
    } else if (identical(input$explorer_source, "c")) {
      analysis_stats_c()
    } else {
      analysis_stats_b()
    }
  })

  explorer_csv_path <- reactive({
    if (identical(input$explorer_source, "a")) {
      report_a_csv()
    } else if (identical(input$explorer_source, "c")) {
      report_c_csv()
    } else {
      report_b_csv()
    }
  })

  output$explorer_status <- renderUI({
    df <- explorer_report()
    if (is.null(df) || nrow(df) == 0) {
      return(div(class = "alert alert-secondary", "No results available for the selected pipeline."))
    }
    n <- nrow(df)
    settings <- pdf_export_settings(n)
    n_pages <- ceiling(n / settings$rows_per_page)
    div(
      class = "alert alert-light border mb-0",
      tags$b(format(n, big.mark = ",")), " variant(s) loaded.",
      tags$br(),
      tags$small(
        class = "text-muted",
        "Full PDF: ~", format(n_pages, big.mark = ","), " page(s). ",
        "Use ", tags$b("Current search filter"), " to export a subset, or ",
        tags$b("Save CSV"), " for the complete table."
      )
    )
  })

  output$explorer_results <- DT::renderDT({
    df <- explorer_report()
    if (is.null(df) || nrow(df) == 0) {
      return(DT::datatable(data.frame(Message = "Run Group A, B, or C analysis first."), rownames = FALSE, options = list(dom = "t")))
    }
    render_results_table(df, selection = "single")
  })

  output$explorer_selected_variant <- renderUI({
    row <- get_selected_row("explorer_results", explorer_report())
    selected_variant_card(row)
  })

  output$explorer_evidence <- DT::renderDT({
    row <- get_selected_row("explorer_results", explorer_report())
    if (is.null(row)) return(render_evidence_detail_table(NULL))
    render_evidence_detail_table(parse_evidence_json(row$evidence_json[1]))
  })

  output$download_explorer_pdf <- downloadHandler(
    filename = function() {
      src <- switch(isolate(input$explorer_source),
        a = "GroupA", b = "GroupB", c = "GroupC", "Report"
      )
      scope <- switch(isolate(input$explorer_pdf_scope) %||% "all",
        all = "Full",
        filtered = "Filtered",
        summary = "Summary",
        "Report"
      )
      paste0("ClinicalVariantR_Evidence_Report_", src, "_", scope, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      scope <- isolate(input$explorer_pdf_scope) %||% "all"
      subtitle <- explorer_source_label(isolate(input$explorer_source))

      if (identical(scope, "summary")) {
        full <- explorer_report()
        req(!is.null(full), nrow(full) > 0L)
        stats <- explorer_stats()
        summary_df <- build_classification_summary_df(
          stats$classification_counts %||% NULL,
          full_report = full
        )
        export_evidence_summary_pdf(
          summary_df,
          file,
          title = "ClinicalVariantR Evidence Report - Summary",
          subtitle = subtitle,
          total_variants = nrow(full)
        )
        return(invisible(NULL))
      }

      df_full <- explorer_report()
      req(!is.null(df_full), nrow(df_full) > 0L)
      filtered_idx <- isolate(input$explorer_results_rows_all)
      df <- resolve_explorer_export_df(
        df_full,
        scope = if (identical(scope, "filtered")) "filtered" else "all",
        filtered_row_indices = filtered_idx
      )

      withProgress(
        message = "Generating PDF report...",
        value = 0, min = 0, max = 1,
        {
          export_evidence_report_pdf(
            df,
            file,
            title = "ClinicalVariantR Evidence Report",
            subtitle = sprintf("%s | %s variants", subtitle, format(nrow(df), big.mark = ",")),
            progress_callback = function(p) {
              setProgress(value = p, detail = sprintf("%d%% complete", round(100 * p)))
            }
          )
        }
      )
    }
  )

  output$download_explorer_csv <- downloadHandler(
    filename = function() {
      src <- switch(isolate(input$explorer_source),
        a = "GroupA", b = "GroupB", c = "GroupC", "Report"
      )
      scope <- switch(isolate(input$explorer_pdf_scope) %||% "all",
        all = "Full",
        filtered = "Filtered",
        summary = "Summary",
        "Report"
      )
      paste0("ClinicalVariantR_Evidence_Report_", src, "_", scope, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      scope <- isolate(input$explorer_pdf_scope) %||% "all"

      if (identical(scope, "summary")) {
        full <- explorer_report()
        req(!is.null(full), nrow(full) > 0L)
        stats <- explorer_stats()
        summary_df <- build_classification_summary_df(
          stats$classification_counts %||% NULL,
          full_report = full
        )
        write.csv(summary_df, file, row.names = FALSE)
        return(invisible(NULL))
      }

      if (identical(scope, "all")) {
        csv_path <- explorer_csv_path()
        if (!is.null(csv_path) && nzchar(csv_path) && file.exists(csv_path)) {
          file.copy(csv_path, file, overwrite = TRUE)
          return(invisible(NULL))
        }
      }

      df_full <- explorer_report()
      req(!is.null(df_full), nrow(df_full) > 0L)
      filtered_idx <- isolate(input$explorer_results_rows_all)
      df <- resolve_explorer_export_df(
        df_full,
        scope = if (identical(scope, "filtered")) "filtered" else "all",
        filtered_row_indices = filtered_idx
      )
      write.csv(df, file, row.names = FALSE)
    }
  )

  invisible(ctx)
}

