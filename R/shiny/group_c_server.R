#' Group C - gene panel analysis tab.

register_group_c_server <- function(ctx) {

  input <- ctx$input

  output <- ctx$output

  authorized <- ctx$authorized

  vcf_preview_c <- ctx$vcf_preview_c

  vcf_validation_c <- ctx$vcf_validation_c

  gene_scan_c <- ctx$gene_scan_c

  report_c_data <- ctx$report_c_data

  report_c_full <- ctx$report_c_full

  selected_category_c <- ctx$selected_category_c

  report_c_csv <- ctx$report_c_csv

  analysis_stats_c <- ctx$analysis_stats_c

  run_metadata_c <- ctx$run_metadata_c

  analysis_run_c <- ctx$analysis_run_c

  load_vcf_on_upload <- ctx$load_vcf_on_upload

  defer_secure_upload <- ctx$defer_secure_upload

  auth_required_msg <- ctx$auth_required_msg

  preview_status_ui <- ctx$preview_status_ui

  results_placeholder_ui <- ctx$results_placeholder_ui

  run_complete_analysis <- ctx$run_complete_analysis



  can_run_c <- reactive({

    v <- vcf_validation_c()

    genes <- parse_gene_filter(input$genes_c)

    !is.null(v) && isTRUE(v$can_analyze) && length(genes) > 0L

  })



  output$run_c_ui <- renderUI({

    ready <- can_run_c()

    genes_ok <- length(parse_gene_filter(input$genes_c)) > 0L

    label <- if (ready) {

      "Run Gene Panel Analysis"

    } else if (!genes_ok) {

      "Enter at least one gene symbol"

    } else {

      "Run Analysis (not ready - see red indicator)"

    }

    actionButton(

      "run_c",

      label,

      class = if (ready) "btn-warning w-100 mt-2" else "btn-danger w-100 mt-2",

      icon = icon(if (ready) "play" else "ban"),

      disabled = !ready

    )

  })



  observeEvent(input$vcf_c, {

    if (is.null(input$vcf_c)) {

      vcf_preview_c(NULL)

      vcf_validation_c(NULL)

      gene_scan_c(NULL)

      return()

    }

    if (!isTRUE(authorized())) {

      auth_required_msg()

    }

    report_c_data(NULL)

    report_c_full(NULL)

    selected_category_c(NULL)

    analysis_run_c(FALSE)

    gene_scan_c(NULL)

    load_vcf_on_upload(input$vcf_c, vcf_preview_c, vcf_validation_c, mode = "rapid", session = ctx$session)

    defer_secure_upload(input$vcf_c, "group_c_vcf")

  })



  observeEvent(list(input$vcf_c, input$genes_c), {

    report_c_data(NULL)

    selected_category_c(NULL)

    analysis_run_c(FALSE)

    gene_scan_c(NULL)

    genes <- parse_gene_filter(isolate(input$genes_c))

    vcf_info <- isolate(input$vcf_c)

    if (length(genes) == 0L || is.null(vcf_info) || is.null(vcf_info$datapath)) return()

    tryCatch({

      gene_scan_c(count_vcf_variants_by_genes(vcf_info$datapath, genes))

    }, error = function(e) {

      gene_scan_c(list(error = conditionMessage(e)))

    })

  }, ignoreNULL = FALSE)



  output$readiness_c <- renderUI({

    readiness_indicator_ui(vcf_validation_c(), label = "Gene panel readiness")

  })



  output$engine_status_c <- renderText(ctx$engine_status_text())



  output$validation_status_c <- renderUI({

    validation_summary_ui(vcf_validation_c(), "VCF")

  })



  output$validation_c <- DT::renderDT({

    v <- vcf_validation_c()

    render_validation_table(if (is.null(v)) data.frame() else v$checks)

  })



  output$preview_status_c <- renderUI({

    preview_status_ui(vcf_preview_c(), if (!is.null(input$vcf_c)) input$vcf_c$name else NULL)

  })



  output$gene_preview_hint_c <- renderUI({

    preview <- vcf_preview_c()

    genes <- parse_gene_filter(input$genes_c)

    if (is.null(preview) || !is.null(preview$error) || is.null(preview$preview) || length(genes) == 0L) {

      return(NULL)

    }

    matched <- filter_variants_by_genes(preview$preview, genes)

    scan <- gene_scan_c()

    full_total <- if (!is.null(scan) && is.null(scan$error)) scan$total else NA_integer_

    missing <- if (!is.null(scan) && is.null(scan$error)) scan$missing_genes else character()



    div(

      class = "alert alert-warning py-2 mb-2",

      icon("filter"),

      " Preview (first ", preview$preview_rows, " rows): ", tags$b(nrow(matched)),

      " match ", tags$b(format_gene_filter_label(genes)), ".",

      if (!is.na(full_total)) {

        tagList(

          tags$br(),

          "Full VCF scan: ", tags$b(full_total), " variant row(s) annotated to ",

          tags$b(format_gene_filter_label(genes)), "."

        )

      },

      if (length(missing) > 0L) {

        tagList(

          tags$br(),

          icon("info-circle"),

          " Not found in this VCF: ", tags$b(paste(missing, collapse = ", ")), ". ",

          "Use HGNC symbols (e.g. ", tags$b("ABL1"), " and ", tags$b("BCR"), ", not ", tags$b("BCR-ABL"), ")."

        )

      },

      if (nrow(matched) == 0L && !is.na(full_total) && full_total > 0L) {

        tagList(

          tags$br(),

          icon("lightbulb"),

          " Genes are present later in the file - click ", tags$b("Run Gene Panel Analysis"), " to classify them."

        )

      }

    )

  })



  output$preview_c <- DT::renderDT({

    preview <- vcf_preview_c()

    if (is.null(preview) || !is.null(preview$error) || is.null(preview$preview)) {

      return(DT::datatable(data.frame(Message = "No preview available."), rownames = FALSE, options = list(dom = "t")))

    }

    genes <- parse_gene_filter(input$genes_c)

    df <- if (length(genes) > 0L) filter_variants_by_genes(preview$preview, genes) else preview$preview

    if (nrow(df) == 0L) {

      scan <- gene_scan_c()

      msg <- if (!is.null(scan) && is.null(scan$error) && scan$total > 0L) {

        paste0(

          "No matches in the first ", preview$preview_rows,

          " preview rows, but ", scan$total,

          " variant(s) match in the full VCF. Run Gene Panel Analysis to classify them."

        )

      } else if (!is.null(scan) && is.null(scan$error) && length(scan$missing_genes) > 0L) {

        paste0(

          "No variants match in preview. Not found in VCF: ",

          paste(scan$missing_genes, collapse = ", "),

          ". Check HGNC gene symbols."

        )

      } else {

        "No variants in preview match the entered gene symbols."

      }

      return(DT::datatable(

        data.frame(Message = msg),

        rownames = FALSE, options = list(dom = "t")

      ))

    }

    render_vcf_preview_table(df)

  })



  output$results_placeholder_c <- renderUI({

    results_placeholder_ui(analysis_run_c())

  })



  observeEvent(input$run_c, {

    req(authorized())

    req(can_run_c())

    req(input$vcf_c)



    genes <- parse_gene_filter(input$genes_c)



    tryCatch({

      withProgress(message = paste("Running Group C gene panel:", format_gene_filter_label(genes)), value = 0, {

        incProgress(0.05, detail = "Filtering by gene symbols")

        out <- run_complete_analysis(

          input$vcf_c$datapath,

          mode = "rapid",

          suffix = "c",

          gene_filter = genes

        )

        result <- out$result



        report_c_full(load_full_analysis_report(result$output_csv, result$preview))

        selected_category_c(NULL)

        report_c_data(NULL)

        report_c_csv(result$output_csv)

        run_metadata_c(result$run_metadata)

        analysis_stats_c(list(

          rows_analyzed = result$rows_analyzed,

          rows_skipped = result$rows_skipped,

          rows_classified = result$rows_classified %||% 0L,

          rows_gene_skipped = result$rows_gene_skipped %||% 0L,

          gene_filter = genes,

          rows_displayed = 0L,

          engine = result$engine,

          classification_counts = result$classification_counts,

          metadata_path = result$metadata_path

        ))

        analysis_run_c(TRUE)

        ctx$load_audit()



        classified <- result$rows_classified %||% 0L

        if (classified == 0L) {

          showNotification(

            paste0("No variants found for: ", format_gene_filter_label(genes)),

            type = "warning", duration = 10

          )

        } else {

          path_n <- pathogenic_tier_count(result$classification_counts)

          msg <- sprintf(

            "Gene panel complete: %s variant(s) in %s.",

            format(classified, big.mark = ","),

            format_gene_filter_label(genes)

          )

          if (path_n > 0L) {

            msg <- paste0(msg, sprintf(" %d pathogenic / likely pathogenic.", path_n))

          }

          showNotification(msg, type = if (path_n > 0L) "warning" else "message", duration = 10)

        }

        play_analysis_complete_sound(session)

      })

    }, error = function(e) {

      report_c_data(NULL)

      report_c_full(NULL)

      selected_category_c(NULL)

      report_c_csv(NULL)

      run_metadata_c(NULL)

      analysis_stats_c(NULL)

      analysis_run_c(FALSE)

      showNotification(paste("Gene panel analysis failed:", e$message), type = "error", duration = NULL)

    })

  })



  invisible(ctx)

}


