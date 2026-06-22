# CML Clinical Variant Interpreter — Server Logic

server <- function(input, output, session) {
  session_id <- paste0("SES-", substr(as.character(session$token), 1, 8))

  ensure_audit_log()
  refs <- reactiveVal(NULL)

  vcf_preview_a <- reactiveVal(NULL)
  vcf_preview_b <- reactiveVal(NULL)
  vcf_validation_a <- reactiveVal(NULL)
  vcf_validation_b <- reactiveVal(NULL)
  group_a_validation <- reactiveVal(NULL)
  report_a_data <- reactiveVal(NULL)
  report_b_data <- reactiveVal(NULL)
  report_a_csv <- reactiveVal(NULL)
  report_b_csv <- reactiveVal(NULL)
  analysis_stats_a <- reactiveVal(NULL)
  analysis_stats_b <- reactiveVal(NULL)
  analysis_run_a <- reactiveVal(FALSE)
  analysis_run_b <- reactiveVal(FALSE)

  observe({
    tryCatch({
      refs(load_reference_data())
    }, error = function(e) {
      showNotification(paste("Reference data load warning:", e$message), type = "warning")
    })
  })

  output$session_label <- renderText(paste("Session:", session_id))

  observeEvent(input$select_group_a, {
    bslib::nav_select("main_nav", "group_a")
  })
  observeEvent(input$select_group_b, {
    bslib::nav_select("main_nav", "group_b")
  })

  preview_status_ui <- function(preview, file_name) {
    if (is.null(preview)) {
      return(div(class = "alert alert-secondary", icon("upload"), " Upload a VCF file to preview variants."))
    }
    if (!is.null(preview$error)) {
      return(div(class = "alert alert-danger", icon("exclamation-triangle"), " ", preview$error))
    }
    div(
      class = "alert alert-light border",
      tags$b(file_name), " — ",
      preview$file_size_mb, " MB | ",
      "Total variants in file: ", format(preview$total_variants, big.mark = ","), " | ",
      "Preview: first ", preview$preview_rows, " row(s)."
    )
  }

  results_placeholder_ui <- function(run_flag) {
    if (!isTRUE(run_flag)) {
      div(
        class = "alert alert-info",
        icon("info-circle"),
        " Analysis has not been run. Complete requirement checks, review preview, then click Run Analysis."
      )
    }
  }

  refresh_group_a_validation <- function() {
    if (is.null(input$vcf_a)) {
      vcf_validation_a(NULL)
      group_a_validation(NULL)
      return()
    }

    v_val <- tryCatch(
      validate_vcf(input$vcf_a$datapath, mode = "full"),
      error = function(e) list(valid = FALSE, can_analyze = FALSE, checks = data.frame(), summary = conditionMessage(e))
    )
    vcf_validation_a(v_val)

    if (!is.null(input$clinical_a) && !is.null(input$pedigree_a)) {
      g_val <- tryCatch(
        validate_group_a_inputs(
          input$vcf_a$datapath,
          input$clinical_a$datapath,
          input$pedigree_a$datapath
        ),
        error = function(e) list(can_analyze = FALSE, checks = data.frame(), summary = conditionMessage(e))
      )
      group_a_validation(g_val)
    } else {
      missing <- c()
      if (is.null(input$clinical_a)) missing <- c(missing, "Clinical Logs CSV file not uploaded")
      if (is.null(input$pedigree_a)) missing <- c(missing, "Pedigree CSV file not uploaded")
      group_a_validation(list(
        can_analyze = FALSE,
        readiness = "NOT_READY",
        checks = v_val$checks,
        missing_items = unique(c(v_val$missing_items %||% character(), missing)),
        summary = if (length(missing) > 0) {
          paste0("Not ready. Missing: ", paste(unique(c(v_val$missing_items %||% character(), missing)), collapse = "; "))
        } else {
          v_val$summary
        }
      ))
    }
  }

  load_vcf_on_upload <- function(file_info, preview_rv, validation_rv, mode) {
    preview_rv(NULL)
    validation_rv(NULL)
    if (is.null(file_info)) return(invisible(NULL))

    tryCatch({
      validation_rv(validate_vcf(file_info$datapath, mode = mode))
      preview_rv(preview_vcf(file_info$datapath))
      showNotification("VCF uploaded. Review requirement checks and preview before analysis.", type = "message", duration = 4)
    }, error = function(e) {
      preview_rv(list(error = conditionMessage(e)))
      validation_rv(list(valid = FALSE, can_analyze = FALSE, checks = data.frame(), summary = conditionMessage(e)))
      showNotification(paste("VCF error:", e$message), type = "error")
    })
    invisible(NULL)
  }

  can_run_a <- reactive({
    g <- group_a_validation()
    !is.null(g) && isTRUE(g$can_analyze)
  })

  can_run_b <- reactive({
    v <- vcf_validation_b()
    !is.null(v) && isTRUE(v$can_analyze)
  })

  output$run_a_ui <- renderUI({
    ready <- can_run_a()
    actionButton(
      "run_a",
      if (ready) "Run Analysis" else "Run Analysis (not ready — see red indicator)",
      class = if (ready) "btn-primary w-100 mt-2" else "btn-danger w-100 mt-2",
      icon = icon(if (ready) "play" else "ban"),
      disabled = !ready
    )
  })

  output$run_b_ui <- renderUI({
    ready <- can_run_b()
    actionButton(
      "run_b",
      if (ready) "Run Analysis" else "Run Analysis (not ready — see red indicator)",
      class = if (ready) "btn-success w-100 mt-2" else "btn-danger w-100 mt-2",
      icon = icon(if (ready) "play" else "ban"),
      disabled = !ready
    )
  })

  # --- Group A ---
  observeEvent(input$vcf_a, {
    report_a_data(NULL)
    analysis_run_a(FALSE)
    load_vcf_on_upload(input$vcf_a, vcf_preview_a, vcf_validation_a, mode = "full")
    refresh_group_a_validation()
  }, ignoreNULL = FALSE)

  observeEvent(input$clinical_a, {
    report_a_data(NULL)
    analysis_run_a(FALSE)
    refresh_group_a_validation()
  }, ignoreNULL = FALSE)

  observeEvent(input$pedigree_a, {
    report_a_data(NULL)
    analysis_run_a(FALSE)
    refresh_group_a_validation()
  }, ignoreNULL = FALSE)

  output$readiness_a <- renderUI({
    readiness_indicator_ui(group_a_validation(), label = "Group A readiness")
  })

  output$validation_status_a <- renderUI({
    validation_summary_ui(group_a_validation(), "Group A inputs")
  })

  output$validation_a <- DT::renderDT({
    g <- group_a_validation()
    checks <- if (is.null(g)) data.frame() else g$checks
    render_validation_table(checks)
  })

  output$preview_status_a <- renderUI({
    preview_status_ui(vcf_preview_a(), if (!is.null(input$vcf_a)) input$vcf_a$name else NULL)
  })

  output$preview_a <- DT::renderDT({
    preview <- vcf_preview_a()
    if (is.null(preview) || !is.null(preview$error) || is.null(preview$preview)) {
      return(DT::datatable(data.frame(Message = "No preview available."), rownames = FALSE, options = list(dom = "t")))
    }
    render_vcf_preview_table(preview$preview)
  })

  output$results_placeholder_a <- renderUI({
    results_placeholder_ui(analysis_run_a())
  })

  manual_inputs <- reactive({
    list(
      PS3_functional = input$manual_ps3,
      PP4_phenotype = input$manual_pp4,
      PS4_case_control = input$manual_ps4,
      PS2_de_novo = input$manual_ps2,
      PM6_de_novo = input$manual_pm6,
      PP1_segregation = input$manual_pp1,
      PP2_missense_mechanism = input$manual_pp2
    )
  })

  output$engine_status_a <- renderText({
    if (bcftools_available()) {
      "bcftools detected — recommended for large VEP VCFs on Ubuntu/WSL."
    } else {
      "bcftools not found — using R streaming (install on Ubuntu: bash scripts/ubuntu_setup.sh)."
    }
  })

  output$engine_status_b <- renderText({
    if (bcftools_available()) {
      "bcftools detected — recommended for large VEP VCFs on Ubuntu/WSL."
    } else {
      "bcftools not found — using R streaming (install on Ubuntu: bash scripts/ubuntu_setup.sh)."
    }
  })

  run_complete_analysis <- function(vcf_path, mode, suffix) {
    pass_only <- isTRUE(input[[paste0("pass_only_", suffix)]])
    min_qual <- as.numeric(input[[paste0("min_qual_", suffix)]]) %||% 0
    use_bcftools <- isTRUE(input[[paste0("use_bcftools_", suffix)]])
    chunk_size <- as.integer(input[[paste0("chunk_size_", suffix)]]) %||% 10000L
    complete <- isTRUE(input[[paste0("complete_vcf_", suffix)]])
    filter_path <- isTRUE(input[[paste0("filter_pathogenic_", suffix)]])

    if (!complete) {
      stop("Enable 'Analyze entire VCF' for complete analysis.")
    }

    clinical <- NULL
    pedigree <- NULL
    manual <- list()
    if (mode == "full") {
      clinical <- parse_clinical_logs(input$clinical_a$datapath)
      pedigree <- parse_pedigree(input$pedigree_a$datapath)
      manual <- manual_inputs()
    }

    r <- refs()
    progress_count <- 0L

    result <- analyze_complete_vcf(
      vcf_path = vcf_path,
      mode = mode,
      pass_only = pass_only,
      min_qual = min_qual,
      chunk_size = chunk_size,
      use_bcftools = use_bcftools,
      refs = r,
      manual_inputs = manual,
      clinical_context = clinical,
      pedigree_context = pedigree,
      session_id = session_id,
      filter_pathogenic_only = filter_path && mode == "rapid",
      progress_fn = function(detail = NULL, ...) {
        progress_count <<- progress_count + 1L
        incProgress(0.1, detail = detail)
      }
    )

    list(result = result, filter_pathogenic = filter_path && mode == "rapid")
  }

  observeEvent(input$run_a, {
    req(can_run_a())
    req(input$vcf_a, input$clinical_a, input$pedigree_a)

    tryCatch({
      withProgress(message = "Running complete Group A analysis...", value = 0, {
        incProgress(0.05, detail = "Starting streaming pipeline")
        out <- run_complete_analysis(input$vcf_a$datapath, mode = "full", suffix = "a")
        result <- out$result

        report_a_data(result$preview)
        report_a_csv(result$output_csv)
        analysis_stats_a(list(
          rows_analyzed = result$rows_analyzed,
          rows_skipped = result$rows_skipped,
          engine = result$engine,
          classification_counts = result$classification_counts
        ))
        analysis_run_a(TRUE)
        load_audit()
        showNotification(
          sprintf("Complete analysis: %s variants classified (%s). Download for full report.",
                  format(result$rows_analyzed, big.mark = ","), result$engine),
          type = "message", duration = 8
        )
      })
    }, error = function(e) {
      report_a_data(NULL)
      report_a_csv(NULL)
      analysis_stats_a(NULL)
      analysis_run_a(FALSE)
      showNotification(paste("Analysis failed:", e$message), type = "error", duration = NULL)
    })
  })

  # --- Group B ---
  observeEvent(input$vcf_b, {
    report_b_data(NULL)
    analysis_run_b(FALSE)
    load_vcf_on_upload(input$vcf_b, vcf_preview_b, vcf_validation_b, mode = "rapid")
  }, ignoreNULL = FALSE)

  output$readiness_b <- renderUI({
    readiness_indicator_ui(vcf_validation_b(), label = "VCF readiness")
  })

  output$validation_status_b <- renderUI({
    validation_summary_ui(vcf_validation_b(), "VCF")
  })

  output$validation_b <- DT::renderDT({
    v <- vcf_validation_b()
    render_validation_table(if (is.null(v)) data.frame() else v$checks)
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

  observeEvent(input$run_b, {
    req(can_run_b())
    req(input$vcf_b)

    tryCatch({
      withProgress(message = "Running complete Group B analysis...", value = 0, {
        incProgress(0.05, detail = "Starting streaming pipeline")
        out <- run_complete_analysis(input$vcf_b$datapath, mode = "rapid", suffix = "b")
        result <- out$result

        report_b_data(result$preview)
        report_b_csv(result$output_csv)
        analysis_stats_b(list(
          rows_analyzed = result$rows_analyzed,
          rows_skipped = result$rows_skipped,
          engine = result$engine,
          classification_counts = result$classification_counts
        ))
        analysis_run_b(TRUE)
        load_audit()
        showNotification(
          sprintf("Complete analysis: %s variants classified (%s). Download for full report.",
                  format(result$rows_analyzed, big.mark = ","), result$engine),
          type = "message", duration = 8
        )
      })
    }, error = function(e) {
      report_b_data(NULL)
      report_b_csv(NULL)
      analysis_stats_b(NULL)
      analysis_run_b(FALSE)
      showNotification(paste("Analysis failed:", e$message), type = "error", duration = NULL)
    })
  })

  output$status_a <- renderUI({
    req(analysis_run_a())
    stats <- analysis_stats_a()
    df <- report_a_data()
    if (is.null(stats) || stats$rows_analyzed == 0) {
      return(div(class = "alert alert-warning", "No variants classified."))
    }
    counts <- stats$classification_counts
    count_text <- if (length(counts) > 0) {
      paste(names(counts), unlist(counts), sep = ": ", collapse = " | ")
    } else {
      "See download for classifications"
    }
    div(
      class = "alert alert-info",
      tags$b(format(stats$rows_analyzed, big.mark = ",")), " variant(s) analyzed (complete VCF). ",
      if (stats$rows_skipped > 0) tags$span(class = "text-muted", sprintf("(%s skipped by filters) ", format(stats$rows_skipped, big.mark = ","))),
      tags$br(),
      "Engine: ", tags$code(stats$engine), " | ",
      count_text,
      tags$br(),
      tags$em("Table shows first 1,000 rows. Download CSV for the full report.")
    )
  })

  output$status_b <- renderUI({
    req(analysis_run_b())
    stats <- analysis_stats_b()
    if (is.null(stats) || stats$rows_analyzed == 0) {
      return(div(class = "alert alert-warning", "No variants matched filter criteria."))
    }
    counts <- stats$classification_counts
    count_text <- if (length(counts) > 0) {
      paste(names(counts), unlist(counts), sep = ": ", collapse = " | ")
    } else {
      "See download for classifications"
    }
    div(
      class = "alert alert-success",
      tags$b(format(stats$rows_analyzed, big.mark = ",")), " variant(s) analyzed (complete VCF). ",
      if (stats$rows_skipped > 0) tags$span(class = "text-muted", sprintf("(%s skipped by filters) ", format(stats$rows_skipped, big.mark = ","))),
      tags$br(),
      "Engine: ", tags$code(stats$engine), " | ",
      count_text,
      tags$br(),
      tags$em("Table shows first 1,000 rows. Download CSV for the full report.")
    )
  })

  output$results_a <- DT::renderDT({
    req(analysis_run_a(), report_a_data())
    render_results_table(report_a_data())
  })

  output$results_b <- DT::renderDT({
    req(analysis_run_b(), report_b_data())
    render_results_table(report_b_data())
  })

  output$download_a <- downloadHandler(
    filename = function() paste0("Final_Clinical_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
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
    filename = function() paste0("Final_Clinical_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
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

  audit_data <- reactiveVal(data.frame())

  load_audit <- function() {
    path <- AUDIT_LOG_PATH
    audit_data(if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else data.frame())
  }

  observe({ load_audit() })
  observeEvent(input$refresh_audit, { load_audit() })

  output$audit_table <- DT::renderDT({
    DT::datatable(
      audit_data(),
      options = list(pageLength = 20, scrollX = TRUE, order = list(list(1, "desc"))),
      rownames = FALSE
    )
  })
}
