#' Group A — full clinical + pedigree analysis tab.
register_group_a_server <- function(ctx) {
  input <- ctx$input
  output <- ctx$output
  session <- ctx$session
  session_id <- ctx$session_id
  authorized <- ctx$authorized
  refs <- ctx$refs
  vcf_preview_a <- ctx$vcf_preview_a
  vcf_validation_a <- ctx$vcf_validation_a
  group_a_validation <- ctx$group_a_validation
  report_a_data <- ctx$report_a_data
  report_a_full <- ctx$report_a_full
  selected_category_a <- ctx$selected_category_a
  report_a_csv <- ctx$report_a_csv
  analysis_stats_a <- ctx$analysis_stats_a
  run_metadata_a <- ctx$run_metadata_a
  analysis_run_a <- ctx$analysis_run_a
  manual_evidence_a <- ctx$manual_evidence_a
  analysis_context_a <- ctx$analysis_context_a
  load_vcf_on_upload <- ctx$load_vcf_on_upload
  defer_secure_upload <- ctx$defer_secure_upload
  auth_required_msg <- ctx$auth_required_msg
  preview_status_ui <- ctx$preview_status_ui
  results_placeholder_ui <- ctx$results_placeholder_ui
  get_selected_row <- ctx$get_selected_row
  run_complete_analysis <- ctx$run_complete_analysis

  refresh_group_a_validation <- function() {
    if (is.null(input$vcf_a)) {
      vcf_validation_a(NULL)
      group_a_validation(NULL)
      return()
    }

    v_val <- vcf_validation_a()
    if (is.null(v_val)) {
      v_val <- tryCatch(
        validate_vcf(input$vcf_a$datapath, mode = "full"),
        error = function(e) list(valid = FALSE, can_analyze = FALSE, checks = data.frame(), summary = conditionMessage(e))
      )
      vcf_validation_a(v_val)
    }

    if (!is.null(input$clinical_a) && !is.null(input$pedigree_a)) {
      g_val <- tryCatch(
        validate_group_a_inputs(
          input$vcf_a$datapath,
          input$clinical_a$datapath,
          input$pedigree_a$datapath,
          vcf_val = v_val
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
  ctx$refresh_group_a_validation <- refresh_group_a_validation

  can_run_a <- reactive({
    g <- group_a_validation()
    !is.null(g) && isTRUE(g$can_analyze)
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

  observeEvent(input$vcf_a, {
    if (is.null(input$vcf_a)) {
      vcf_preview_a(NULL)
      vcf_validation_a(NULL)
      group_a_validation(NULL)
      return()
    }
    if (!isTRUE(authorized())) {
      auth_required_msg()
    }
    report_a_data(NULL)
    report_a_full(NULL)
    selected_category_a(NULL)
    analysis_run_a(FALSE)
    load_vcf_on_upload(input$vcf_a, vcf_preview_a, vcf_validation_a, mode = "full", session = session)
    refresh_group_a_validation()
    defer_secure_upload(input$vcf_a, "group_a_vcf")
  })

  observeEvent(input$clinical_a, {
    if (is.null(input$clinical_a)) return()
    if (!isTRUE(authorized())) {
      auth_required_msg()
      return()
    }
    defer_secure_upload(input$clinical_a, "group_a_clinical")
    report_a_data(NULL)
    analysis_run_a(FALSE)
    refresh_group_a_validation()
  })

  observeEvent(input$pedigree_a, {
    if (is.null(input$pedigree_a)) return()
    if (!isTRUE(authorized())) {
      auth_required_msg()
      return()
    }
    defer_secure_upload(input$pedigree_a, "group_a_pedigree")
    report_a_data(NULL)
    analysis_run_a(FALSE)
    refresh_group_a_validation()
  })

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

  output$reference_readiness_home <- renderUI({
    render_reference_readiness_banner()
  })

  sync_curation_inputs <- function(variant_id) {
    cur <- manual_inputs_for_variant(manual_evidence_a(), variant_id)
    updateCheckboxInput(session, "cur_ps3", value = isTRUE(cur$PS3_functional))
    updateCheckboxInput(session, "cur_pp4", value = isTRUE(cur$PP4_phenotype))
    updateCheckboxInput(session, "cur_ps4", value = isTRUE(cur$PS4_case_control))
    updateCheckboxInput(session, "cur_ps2", value = isTRUE(cur$PS2_de_novo))
    updateCheckboxInput(session, "cur_pm6", value = isTRUE(cur$PM6_de_novo))
    updateCheckboxInput(session, "cur_pp1", value = isTRUE(cur$PP1_segregation))
    updateCheckboxInput(session, "cur_pp2", value = isTRUE(cur$PP2_missense_mechanism))
  }

  observeEvent(input$results_a_rows_selected, {
    row <- get_selected_row("results_a", report_a_data())
    if (!is.null(row)) sync_curation_inputs(row$variant_id[1])
  }, ignoreNULL = FALSE)

  output$curation_status_a <- renderUI({
    NULL
  })

  observeEvent(input$apply_curation_a, {
    row <- get_selected_row("results_a", report_a_data())
    a_ctx <- analysis_context_a()
    if (is.null(row) || is.null(a_ctx)) {
      showNotification("Select a variant and run analysis first.", type = "warning")
      return()
    }
    manual <- list(
      PS3_functional = isTRUE(input$cur_ps3),
      PP4_phenotype = isTRUE(input$cur_pp4),
      PS4_case_control = isTRUE(input$cur_ps4),
      PS2_de_novo = isTRUE(input$cur_ps2),
      PM6_de_novo = isTRUE(input$cur_pm6),
      PP1_segregation = isTRUE(input$cur_pp1),
      PP2_missense_mechanism = isTRUE(input$cur_pp2)
    )
    vid <- row$variant_id[1]
    store <- manual_evidence_a()
    store[[vid]] <- manual
    manual_evidence_a(store)

    tryCatch({
      vrow <- parse_single_variant_from_vcf(
        a_ctx$vcf_path, row$chrom[1], row$pos[1], row$ref[1], row$alt[1],
        pass_only = isTRUE(a_ctx$pass_only)
      )
      if (is.null(vrow)) stop("Variant not found in VCF (check FILTER/pass settings).")
      scored <- rescore_variant_with_manual(
        vrow,
        manual_inputs = manual,
        clinical_context = a_ctx$clinical,
        pedigree_context = a_ctx$pedigree,
        profile_id = a_ctx$profile_id,
        refs = refs(),
        evidence_scope = "full"
      )
      new_report <- rescore_variant_to_report_row(
        scored, mode = "full", session_id = session_id, run_metadata = run_metadata_a()
      )
      full <- patch_report_row(report_a_full(), new_report)
      report_a_full(full)
      write_report_csv(full, report_a_csv())
      cat <- selected_category_a()
      report_a_data(full[full$classification == cat, , drop = FALSE])
      showNotification(paste("Reclassified", vid, "→", new_report$classification[1]), type = "message")
    }, error = function(e) {
      showNotification(paste("Curation failed:", conditionMessage(e)), type = "error")
    })
  })

  output$engine_status_a <- renderText(ctx$engine_status_text())

  observeEvent(input$run_a, {
    req(authorized())
    req(can_run_a())
    req(input$vcf_a, input$clinical_a, input$pedigree_a)

    tryCatch({
      withProgress(message = "Running complete Group A analysis...", value = 0, {
        incProgress(0.05, detail = "Starting streaming pipeline")
        out <- run_complete_analysis(input$vcf_a$datapath, mode = "full", suffix = "a")
        result <- out$result

        report_a_full(load_full_analysis_report(result$output_csv, result$preview))
        selected_category_a(NULL)
        report_a_data(NULL)
        report_a_csv(result$output_csv)
        run_metadata_a(result$run_metadata)
        analysis_stats_a(list(
          rows_analyzed = result$rows_analyzed,
          rows_skipped = result$rows_skipped,
          rows_classified = result$rows_classified %||% result$rows_analyzed,
          rows_displayed = result$rows_displayed %||% nrow(result$preview),
          engine = result$engine,
          classification_counts = result$classification_counts,
          metadata_path = result$metadata_path
        ))
        analysis_run_a(TRUE)
        analysis_context_a(list(
          vcf_path = input$vcf_a$datapath,
          clinical = parse_clinical_logs(input$clinical_a$datapath),
          pedigree = parse_pedigree(input$pedigree_a$datapath),
          profile_id = input$profile_a %||% DEFAULT_PROFILE_ID,
          pass_only = isTRUE(input$pass_only_a)
        ))
        ctx$load_audit()
        showNotification(
          format_analysis_notification(result),
          type = "message", duration = 10
        )
      })
    }, error = function(e) {
      report_a_data(NULL)
      report_a_full(NULL)
      selected_category_a(NULL)
      report_a_csv(NULL)
      run_metadata_a(NULL)
      analysis_stats_a(NULL)
      analysis_run_a(FALSE)
      showNotification(paste("Analysis failed:", e$message), type = "error", duration = NULL)
    })
  })

  invisible(ctx)
}
