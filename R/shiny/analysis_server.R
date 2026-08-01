#' Streaming VCF analysis runner (shared by Groups A/B/C).
register_analysis_server <- function(ctx) {
  input <- ctx$input
  session_id <- ctx$session_id
  refs <- ctx$refs
  manual_evidence_a <- ctx$manual_evidence_a

  ctx$run_complete_analysis <- function(vcf_path, mode, suffix, gene_filter = character()) {
    pass_only <- isTRUE(input[[paste0("pass_only_", suffix)]])
    write_audit <- !isTRUE(input[[paste0("skip_audit_", suffix)]])

    # Always 5-way parallel; chunk size auto-scales with RAM unless user disables auto.
    parallel_chunks <- 5L
    auto_chunk <- input[[paste0("auto_chunk_", suffix)]]
    use_auto <- if (is.null(auto_chunk)) TRUE else isTRUE(auto_chunk)
    chunk_raw <- input[[paste0("chunk_size_", suffix)]]
    chunk_manual <- scalar_int(chunk_raw)
    chunk_size <- resolve_chunk_size(
      chunk_size = if (use_auto) NULL else chunk_manual,
      parallel_chunks = parallel_chunks,
      auto = use_auto
    )

    complete_val <- input[[paste0("complete_vcf_", suffix)]]
    # Default to TRUE when checkbox has not been rendered yet / is NULL
    complete <- if (is.null(complete_val)) TRUE else isTRUE(complete_val)

    if (!complete) {
      stop("Enable 'Analyze entire VCF' for complete analysis.")
    }

    clinical <- NULL
    pedigree <- NULL
    if (mode == "full") {
      clinical <- parse_clinical_logs(input$clinical_a$datapath)
      pedigree <- parse_pedigree(input$pedigree_a$datapath)
    }

    r <- refs()
    profile_id <- input[[paste0("profile_", suffix)]] %||% DEFAULT_PROFILE_ID
    manual_map <- if (mode == "full") manual_evidence_a() else list()

    min_qual_raw <- input[[paste0("min_qual_", suffix)]]
    min_qual <- scalar_num(min_qual_raw)
    if (is.na(min_qual)) min_qual <- 0

    showNotification(
      sprintf(
        "Running with 5 parallel chunks; chunk size %s (%s)",
        format(chunk_size, big.mark = ","),
        if (use_auto) "auto from RAM" else "manual"
      ),
      type = "message",
      duration = 5
    )

    result <- analyze_complete_vcf(
      vcf_path = vcf_path,
      mode = mode,
      pass_only = pass_only,
      min_qual = min_qual,
      chunk_size = chunk_size,
      parallel_chunks = parallel_chunks,
      use_bcftools = isTRUE(input[[paste0("use_bcftools_", suffix)]]),
      refs = r,
      manual_by_variant = manual_map,
      clinical_context = clinical,
      pedigree_context = pedigree,
      session_id = session_id,
      profile_id = profile_id,
      write_audit = write_audit,
      gene_filter = gene_filter,
      progress_fn = function(value = NULL, detail = NULL, message = NULL, ...) {
        # Absolute % progress (0..1). Avoids overflowing the Shiny progress bar.
        if (!is.null(value) && is.finite(value)) {
          shiny::setProgress(
            value = max(0, min(1, as.numeric(value))),
            message = message,
            detail = detail
          )
        } else {
          shiny::setProgress(message = message, detail = detail)
        }
      }
    )

    list(result = result)
  }

  invisible(ctx)
}
