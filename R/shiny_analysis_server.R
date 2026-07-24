#' Streaming VCF analysis runner (shared by Groups A/B/C).
#' @noRd
register_analysis_server <- function(ctx) {
  input <- ctx$input
  session_id <- ctx$session_id
  refs <- ctx$refs
  manual_evidence_a <- ctx$manual_evidence_a

  ctx$run_complete_analysis <- function(vcf_path, mode, suffix, gene_filter = character()) {
    pass_only <- isTRUE(input[[paste0("pass_only_", suffix)]])
    write_audit <- !isTRUE(input[[paste0("skip_audit_", suffix)]])
    chunk_raw <- input[[paste0("chunk_size_", suffix)]]
    chunk_size <- scalar_int(chunk_raw)
    if (is.na(chunk_size) || chunk_size < 100L) chunk_size <- 10000L

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
    progress_state <- new.env(parent = emptyenv())
    progress_state$count <- 0L
    profile_id <- input[[paste0("profile_", suffix)]] %||% DEFAULT_PROFILE_ID
    manual_map <- if (mode == "full") manual_evidence_a() else list()

    min_qual_raw <- input[[paste0("min_qual_", suffix)]]
    min_qual <- scalar_num(min_qual_raw)
    if (is.na(min_qual)) min_qual <- 0

    result <- analyze_complete_vcf(
      vcf_path = vcf_path,
      mode = mode,
      pass_only = pass_only,
      min_qual = min_qual,
      chunk_size = chunk_size,
      use_bcftools = isTRUE(input[[paste0("use_bcftools_", suffix)]]),
      refs = r,
      manual_by_variant = manual_map,
      clinical_context = clinical,
      pedigree_context = pedigree,
      session_id = session_id,
      profile_id = profile_id,
      write_audit = write_audit,
      gene_filter = gene_filter,
      progress_fn = function(detail = NULL, ...) {
        progress_state$count <- progress_state$count + 1L
        # Cap progress increments so Shiny progress UI does not overflow
        step <- if (progress_state$count <= 8L) 0.1 else 0.02
        incProgress(step, detail = detail)
      }
    )

    list(result = result)
  }

  invisible(ctx)
}
