#' VCF upload, validation, and deferred secure storage.
register_upload_server <- function(ctx) {
  session <- ctx$session
  input <- ctx$input
  authorized <- ctx$authorized
  auth_user <- ctx$auth_user

  ctx$auth_required_msg <- function() {
    showNotification(
      "Sign in first. Files selected before login are processed after sign-in.",
      type = "warning",
      duration = 8
    )
  }

  ctx$wait_for_upload_datapath <- function(datapath, max_wait_sec = 5) {
    if (is.null(datapath) || !nzchar(datapath)) return(FALSE)
    if (file.exists(datapath)) return(TRUE)
    deadline <- Sys.time() + max_wait_sec
    while (!file.exists(datapath) && Sys.time() < deadline) {
      Sys.sleep(0.1)
    }
    file.exists(datapath)
  }

  ctx$upload_path_unavailable_validation <- function(file_name = NULL) {
    list(
      valid = FALSE,
      can_analyze = FALSE,
      checks = data.frame(
        category = "File",
        requirement = "Server upload path",
        status = "FAIL",
        detail = "File not found on server. Re-select the file after sign-in.",
        stringsAsFactors = FALSE
      ),
      summary = if (nzchar(file_name %||% "")) {
        sprintf("Upload of '%s' incomplete — re-select the file.", file_name)
      } else {
        "Upload incomplete — re-select the file."
      },
      missing_items = "File not available on server"
    )
  }

  ctx$load_vcf_on_upload <- function(file_info, preview_rv, validation_rv, mode, session = NULL) {
    wait_for_upload_datapath <- ctx$wait_for_upload_datapath
    upload_path_unavailable_validation <- ctx$upload_path_unavailable_validation

    if (is.null(file_info)) {
      preview_rv(NULL)
      validation_rv(NULL)
      return(invisible(NULL))
    }

    datapath <- file_info$datapath
    if (is.null(datapath) || !wait_for_upload_datapath(datapath)) {
      preview_rv(NULL)
      validation_rv(upload_path_unavailable_validation(file_info$name))
      return(invisible(NULL))
    }

    preview_rv(NULL)

    tryCatch({
      validation_result <- validate_vcf(datapath, mode = mode)
      validation_rv(validation_result)
      tryCatch({
        preview_rv(preview_vcf(datapath))
      }, error = function(e) {
        preview_rv(list(error = conditionMessage(e)))
      })

      showNotification(
        sprintf("VCF loaded: %s", file_info$name %||% basename(datapath)),
        type = "message",
        duration = 4
      )
      invisible(validation_result)
    }, error = function(e) {
      preview_rv(list(error = conditionMessage(e)))
      validation_rv(list(
        valid = FALSE, can_analyze = FALSE,
        checks = data.frame(), summary = conditionMessage(e)
      ))
      showNotification(paste("VCF error:", conditionMessage(e)), type = "error", duration = NULL)
    })
    invisible(NULL)
  }

  ctx$defer_secure_upload <- function(fileinfo, label) {
    if (is.null(fileinfo)) return(invisible(NULL))
    if (!isTRUE(isolate(authorized()))) return(invisible(NULL))
    user <- isolate(auth_user())
    fileinfo_capture <- fileinfo
    label_capture <- label

    session$onFlushed(function() {
      rec <- secure_store_shiny_upload(fileinfo_capture, label_capture)
      if (!is.null(rec)) {
        append_access_audit(
          user,
          paste0("secure_upload:", label_capture),
          paste(rec$original_name, rec$method, sep = " | ")
        )
      }
      invisible(rec)
    }, once = TRUE)
    invisible(NULL)
  }

  ctx$reprocess_pending_uploads <- function() {
    if (!isTRUE(authorized())) return(invisible(NULL))

    load_vcf_on_upload <- ctx$load_vcf_on_upload
    defer_secure_upload <- ctx$defer_secure_upload

    if (!is.null(input$vcf_a)) {
      ctx$report_a_data(NULL)
      ctx$report_a_full(NULL)
      ctx$selected_category_a(NULL)
      ctx$analysis_run_a(FALSE)
      load_vcf_on_upload(input$vcf_a, ctx$vcf_preview_a, ctx$vcf_validation_a, mode = "full", session = session)
      if (!is.null(ctx$refresh_group_a_validation)) ctx$refresh_group_a_validation()
      defer_secure_upload(input$vcf_a, "group_a_vcf")
    }
    if (!is.null(input$clinical_a)) {
      defer_secure_upload(input$clinical_a, "group_a_clinical")
      if (!is.null(ctx$refresh_group_a_validation)) ctx$refresh_group_a_validation()
    }
    if (!is.null(input$pedigree_a)) {
      defer_secure_upload(input$pedigree_a, "group_a_pedigree")
      if (!is.null(ctx$refresh_group_a_validation)) ctx$refresh_group_a_validation()
    }
    if (!is.null(input$vcf_b)) {
      ctx$report_b_data(NULL)
      ctx$report_b_full(NULL)
      ctx$selected_category_b(NULL)
      ctx$analysis_run_b(FALSE)
      load_vcf_on_upload(input$vcf_b, ctx$vcf_preview_b, ctx$vcf_validation_b, mode = "rapid", session = session)
      if (!is.null(input$vcf_b$datapath) && file.exists(input$vcf_b$datapath)) {
        ctx$vcf_path_b(input$vcf_b$datapath)
      }
      defer_secure_upload(input$vcf_b, "group_b_vcf")
    }
    if (!is.null(input$vcf_c)) {
      ctx$report_c_data(NULL)
      ctx$report_c_full(NULL)
      ctx$selected_category_c(NULL)
      ctx$analysis_run_c(FALSE)
      ctx$gene_scan_c(NULL)
      load_vcf_on_upload(input$vcf_c, ctx$vcf_preview_c, ctx$vcf_validation_c, mode = "rapid", session = session)
      defer_secure_upload(input$vcf_c, "group_c_vcf")
    }
    invisible(NULL)
  }

  invisible(ctx)
}
