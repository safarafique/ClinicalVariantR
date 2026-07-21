#' Group B - automated VCF-only prediction tab.
register_group_b_server <- function(ctx) {
  input <- ctx$input
  output <- ctx$output
  session <- ctx$session
  authorized <- ctx$authorized
  vcf_preview_b <- ctx$vcf_preview_b
  vcf_validation_b <- ctx$vcf_validation_b
  vcf_path_b <- ctx$vcf_path_b
  report_b_data <- ctx$report_b_data
  report_b_full <- ctx$report_b_full
  selected_category_b <- ctx$selected_category_b
  analysis_run_b <- ctx$analysis_run_b
  report_b_csv <- ctx$report_b_csv
  run_metadata_b <- ctx$run_metadata_b
  analysis_stats_b <- ctx$analysis_stats_b
  load_vcf_on_upload <- ctx$load_vcf_on_upload
  defer_secure_upload <- ctx$defer_secure_upload
  wait_for_upload_datapath <- ctx$wait_for_upload_datapath
  run_complete_analysis <- ctx$run_complete_analysis
  preview_status_ui <- ctx$preview_status_ui
  results_placeholder_ui <- ctx$results_placeholder_ui
  auth_required_msg <- ctx$auth_required_msg

  run_b_running <- reactiveVal(FALSE)
  run_b_feedback <- reactiveVal(list(
    type = "info",
    message = "Upload a VCF, then click Run Analysis."
  ))

  can_run_b <- reactive({
    v <- vcf_validation_b()
    !is.null(v) && isTRUE(v$can_analyze)
  })

  # Keep static button label/state in sync (do NOT recreate button in renderUI; clicks get lost).
  observe({
    ready <- can_run_b()
    running <- isTRUE(run_b_running())
    shiny::updateActionButton(
      session,
      "run_b",
      label = if (running) {
        "Running analysis..."
      } else if (ready) {
        "Run Analysis"
      } else {
        "Run Analysis (not ready)"
      },
      icon = icon(if (running) "spinner" else if (ready) "play" else "ban")
    )
    shinyjs_available <- requireNamespace("shinyjs", quietly = TRUE)
    if (shinyjs_available) {
      if (!ready || running) shinyjs::disable("run_b") else shinyjs::enable("run_b")
    }
  })

  output$run_b_state <- renderUI({
    if (isTRUE(run_b_running())) {
      return(div(
        class = "run-state-card mt-2",
        tags$div(class = "small fw-semibold mb-1", icon("spinner"), " Analysis running..."),
        tags$div(
          class = "progress",
          tags$div(
            class = "progress-bar progress-bar-striped progress-bar-animated bg-success",
            role = "progressbar",
            style = "width: 100%"
          )
        ),
        tags$p(class = "text-muted small mt-2 mb-0", "Please wait - classifying variants now.")
      ))
    }
    feedback <- run_b_feedback()
    if (is.null(feedback) || is.null(feedback$message) || !nzchar(feedback$message)) return(NULL)
    alert_class <- switch(
      as.character(feedback$type %||% "info"),
      success = "alert alert-success mt-2 py-2",
      error = "alert alert-danger mt-2 py-2",
      "alert alert-secondary mt-2 py-2"
    )
    div(class = alert_class, feedback$message)
  })

  output$engine_status_b <- renderText(ctx$engine_status_text())

  observeEvent(input$vcf_b, {
    if (is.null(input$vcf_b)) {
      vcf_preview_b(NULL)
      vcf_validation_b(NULL)
      vcf_path_b(NULL)
      run_b_running(FALSE)
      run_b_feedback(list(type = "info", message = "Upload a VCF file to run analysis."))
      return()
    }
    if (!isTRUE(authorized())) auth_required_msg()
    report_b_data(NULL)
    report_b_full(NULL)
    selected_category_b(NULL)
    analysis_run_b(FALSE)
    run_b_running(FALSE)
    run_b_feedback(list(type = "info", message = "File ready. Click Run Analysis to start."))
    load_vcf_on_upload(input$vcf_b, vcf_preview_b, vcf_validation_b, mode = "rapid", session = session)
    if (!is.null(input$vcf_b$datapath) && file.exists(input$vcf_b$datapath)) {
      vcf_path_b(input$vcf_b$datapath)
    }
    defer_secure_upload(input$vcf_b, "group_b_vcf")
  }, ignoreNULL = TRUE)

  output$readiness_b <- renderUI({
    processing <- !is.null(input$vcf_b) && is.null(vcf_validation_b())
    readiness_indicator_ui(vcf_validation_b(), label = "VCF readiness", processing = processing)
  })

  output$validation_status_b <- renderUI({
    processing <- !is.null(input$vcf_b) && is.null(vcf_validation_b())
    if (processing) {
      return(div(class = "alert alert-info", icon("spinner"), " Validating uploaded file..."))
    }
    validation_summary_ui(vcf_validation_b(), "VCF")
  })

  output$validation_b <- DT::renderDT({
    v <- vcf_validation_b()
    if (is.null(v)) return(render_validation_table(data.frame()))
    checks <- v$checks
    if (is.null(checks) || !is.data.frame(checks) || nrow(checks) == 0L) {
      checks <- data.frame(
        category = "Summary",
        requirement = "VCF validation",
        status = if (isTRUE(v$can_analyze)) "PASS" else "FAIL",
        detail = v$summary %||% "No detailed checks recorded.",
        stringsAsFactors = FALSE
      )
    }
    render_validation_table(checks)
  })

  output$preview_status_b <- renderUI({
    preview_status_ui(vcf_preview_b(), if (!is.null(input$vcf_b)) input$vcf_b$name else NULL)
  })

  output$preview_b <- DT::renderDT({
    preview <- vcf_preview_b()
    if (is.null(preview) || !is.null(preview$error) || is.null(preview$preview)) {
      return(DT::datatable(data.frame(Message = "No preview available."), rownames = FALSE, options = list(dom = "t")))
    }
    render_vcf_preview_table(preview$preview)
  })

  output$results_placeholder_b <- renderUI({
    results_placeholder_ui(analysis_run_b())
  })

  # Critical: ignoreInit so startup does not fire a fake run; ignoreNULL so clicks always register.
  observeEvent(input$run_b, {
    message(sprintf("[Group B] Run clicked at %s", format(Sys.time(), "%H:%M:%S")))
    run_b_feedback(list(type = "info", message = "Run clicked - checking input..."))

    if (isTRUE(run_b_running())) {
      run_b_feedback(list(type = "error", message = "Analysis already running. Please wait."))
      return()
    }

    is_auth <- isTRUE(isolate(authorized()))
    ready <- isTRUE(isolate(can_run_b()))
    vcf_path <- isolate(vcf_path_b())
    if (is.null(vcf_path) || !nzchar(as.character(vcf_path %||% ""))) {
      vcf_path <- isolate(input$vcf_b$datapath)
    }

    if (!is_auth) {
      run_b_feedback(list(type = "error", message = "Not started: authentication blocked (set AUTH_ENABLED = FALSE)."))
      showNotification("Authentication blocked analysis.", type = "error", duration = 8)
      return()
    }
    if (!ready) {
      v <- isolate(vcf_validation_b())
      detail <- if (is.null(v)) {
        "No validation result. Re-upload the VCF."
      } else {
        v$summary %||% "Validation did not pass."
      }
      run_b_feedback(list(type = "error", message = paste("Not started:", detail)))
      showNotification(paste("Not ready:", detail), type = "error", duration = 8)
      return()
    }
    if (is.null(vcf_path) || !nzchar(as.character(vcf_path))) {
      run_b_feedback(list(type = "error", message = "Not started: VCF path missing. Re-upload the file."))
      showNotification("VCF path missing. Re-upload the file.", type = "error", duration = 8)
      return()
    }
    if (!isTRUE(wait_for_upload_datapath(vcf_path, max_wait_sec = 3))) {
      run_b_feedback(list(
        type = "error",
        message = paste("Not started: file not found on server:", basename(vcf_path))
      ))
      showNotification("Uploaded file not found. Re-upload and try again.", type = "error", duration = 8)
      return()
    }
    # Ensure complete_vcf is treated as TRUE if user never touched the checkbox after UI load.
    if (!isTRUE(isolate(input$complete_vcf_b))) {
      updateCheckboxInput(session, "complete_vcf_b", value = TRUE)
    }

    run_b_running(TRUE)
    run_b_feedback(list(type = "info", message = "Running Group B analysis..."))
    showNotification("Starting Group B analysis...", type = "message", duration = 4)

    tryCatch({
      withProgress(message = "Running complete Group B analysis...", value = 0, {
        incProgress(0.05, detail = "Starting streaming pipeline")
        out <- run_complete_analysis(vcf_path, mode = "rapid", suffix = "b")
        result <- out$result

        report_b_full(load_full_analysis_report(result$output_csv, result$preview))
        selected_category_b(NULL)
        report_b_data(NULL)
        report_b_csv(result$output_csv)
        run_metadata_b(result$run_metadata)
        analysis_stats_b(list(
          rows_analyzed = result$rows_analyzed,
          rows_skipped = result$rows_skipped,
          rows_classified = result$rows_classified %||% result$rows_analyzed,
          rows_displayed = 0L,
          engine = result$engine,
          classification_counts = result$classification_counts,
          metadata_path = result$metadata_path
        ))
        analysis_run_b(TRUE)
        if (is.function(ctx$load_audit)) ctx$load_audit()
        n <- result$rows_analyzed %||% 0L
        run_b_feedback(list(
          type = "success",
          message = sprintf("Analysis complete: %s variant(s) processed.", format(n, big.mark = ","))
        ))
        showNotification(format_analysis_notification(result), type = "message", duration = 10)
        play_analysis_complete_sound(session)
      })
    }, error = function(e) {
      report_b_data(NULL)
      report_b_full(NULL)
      selected_category_b(NULL)
      report_b_csv(NULL)
      run_metadata_b(NULL)
      analysis_stats_b(NULL)
      analysis_run_b(FALSE)
      run_b_feedback(list(type = "error", message = paste("Analysis failed:", conditionMessage(e))))
      showNotification(paste("Analysis failed:", conditionMessage(e)), type = "error", duration = NULL)
      message("[Group B] Analysis error: ", conditionMessage(e))
    }, finally = {
      run_b_running(FALSE)
    })
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  invisible(ctx)
}
