#' Shared UI helpers and lightweight outputs used across pipeline tabs.
register_shared_server <- function(ctx) {
  input <- ctx$input
  output <- ctx$output

  report_a_full <- ctx$report_a_full
  report_b_full <- ctx$report_b_full
  report_c_full <- ctx$report_c_full

  register_expert_worklist_downloads(output, "a", report_a_full, ctx$session)
  register_expert_worklist_downloads(output, "b", report_b_full, ctx$session)
  register_expert_worklist_downloads(output, "c", report_c_full, ctx$session)

  output$expert_review_stats_a <- renderUI({
    render_expert_worklist_stats_ui(report_a_full())
  })
  output$expert_review_stats_b <- renderUI({
    render_expert_worklist_stats_ui(report_b_full())
  })
  output$expert_review_stats_c <- renderUI({
    render_expert_worklist_stats_ui(report_c_full())
  })

  expert_handoff_panel <- function(report_full) {
    df <- report_full()
    if (is.null(df) || nrow(df) == 0L) return(NULL)
    expert_csv_handoff_ui()
  }

  output$expert_handoff_a <- renderUI({ expert_handoff_panel(report_a_full) })
  output$expert_handoff_b <- renderUI({ expert_handoff_panel(report_b_full) })
  output$expert_handoff_c <- renderUI({ expert_handoff_panel(report_c_full) })

  output$session_label <- renderText(paste("Session:", ctx$session_id))

  # Keep idle browser tabs connected for at least SESSION_IDLE_TIMEOUT_SEC (default 1 hour).
  if (!is.null(ctx$session$allowReconnect)) {
    ctx$session$allowReconnect(TRUE)
  }
  observe({
    invalidateLater(5 * 60 * 1000) # heartbeat every 5 minutes
    invisible(TRUE)
  })

  output$auth_label <- renderText({
    if (!isTRUE(AUTH_ENABLED)) return("Auth: off")
    if (ctx$authorized()) paste("User:", ctx$auth_user()) else "Not signed in"
  })

  observeEvent(input$select_group_a, {
    bslib::nav_select("main_nav", "group_a")
  })
  observeEvent(input$select_group_b, {
    bslib::nav_select("main_nav", "group_b")
  })
  observeEvent(input$select_group_c, {
    bslib::nav_select("main_nav", "group_c")
  })

  observe({
    tryCatch({
      ctx$refs(load_reference_data())
    }, error = function(e) {
      showNotification(paste("Reference data load warning:", e$message), type = "warning")
    })
  })

  ctx$preview_status_ui <- function(preview, file_name) {
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
      "Total variants in file: ", preview$total_display %||% format(preview$total_variants, big.mark = ","), " | ",
      "Preview: first ", preview$preview_rows, " row(s)."
    )
  }

  ctx$selected_variant_card <- function(row) {
    if (is.null(row) || nrow(row) == 0) {
      return(div(class = "alert alert-light border", "No variant selected."))
    }
    strength <- row$evidence_strength[1] %||% row$confidence_label[1] %||% "—"
    div(
      class = "alert alert-light border",
      tags$b(row$variant_id[1]), " | Gene: ", row$gene[1], " | ",
      "Prediction: ", tags$span(class = "badge bg-secondary", row$classification[1]), " | ",
      "Evidence strength: ", tags$span(class = "badge bg-info text-dark", strength),
      tags$br(),
      tags$small(class = "text-muted", "In silico: ", row$prediction_scores[1]),
      if ("prediction_limitations" %in% names(row) && nzchar(row$prediction_limitations[1])) {
        tagList(
          tags$br(),
          tags$small(class = "text-warning", "Limitations: ", row$prediction_limitations[1])
        )
      }
    )
  }

  ctx$get_selected_row <- function(table_id, report_df) {
    sel <- input[[paste0(table_id, "_rows_selected")]]
    if (is.null(sel) || length(sel) == 0 || is.null(report_df)) return(NULL)
    idx <- as.integer(sel[1])
    if (idx < 1L || idx > nrow(report_df)) return(NULL)
    report_df[idx, , drop = FALSE]
  }

  ctx$results_placeholder_ui <- function(run_flag) {
    if (!isTRUE(run_flag)) {
      div(
        class = "alert alert-info",
        icon("info-circle"),
        " Analysis has not been run. Complete requirement checks, review preview, then click Run Analysis."
      )
    }
  }

  ctx$engine_status_text <- function() {
    engine <- ACMG_PRO_ENGINE
    bcftools_msg <- if (bcftools_available()) {
      "bcftools detected — recommended for large VCFs on Ubuntu/WSL."
    } else {
      "bcftools not found — using R streaming (install: bash scripts/ubuntu_setup.sh)."
    }
    paste(engine, "|", bcftools_msg)
  }

  invisible(ctx)
}
