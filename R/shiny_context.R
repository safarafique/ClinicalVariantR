#' Create shared Shiny session context (reactive values + handles).
#' @noRd
init_shiny_context <- function(input, output, session) {
  session_id <- paste0("SES-", substr(as.character(session$token), 1, 8))

  ctx <- new.env(parent = emptyenv())
  ctx$input <- input
  ctx$output <- output
  ctx$session <- session
  ctx$session_id <- session_id

  ctx$authorized <- reactiveVal(!isTRUE(AUTH_ENABLED))
  ctx$auth_user <- reactiveVal("system")
  ctx$refs <- reactiveVal(NULL)

  ctx$vcf_preview_a <- reactiveVal(NULL)
  ctx$vcf_preview_b <- reactiveVal(NULL)
  ctx$vcf_preview_c <- reactiveVal(NULL)
  ctx$gene_scan_c <- reactiveVal(NULL)
  ctx$vcf_validation_a <- reactiveVal(NULL)
  ctx$vcf_validation_b <- reactiveVal(NULL)
  ctx$vcf_path_b <- reactiveVal(NULL)
  ctx$vcf_validation_c <- reactiveVal(NULL)
  ctx$group_a_validation <- reactiveVal(NULL)

  ctx$report_a_data <- reactiveVal(NULL)
  ctx$report_a_full <- reactiveVal(NULL)
  ctx$report_b_data <- reactiveVal(NULL)
  ctx$report_b_full <- reactiveVal(NULL)
  ctx$report_c_data <- reactiveVal(NULL)
  ctx$report_c_full <- reactiveVal(NULL)

  ctx$selected_category_a <- reactiveVal(NULL)
  ctx$selected_category_b <- reactiveVal(NULL)
  ctx$selected_category_c <- reactiveVal(NULL)

  ctx$report_a_csv <- reactiveVal(NULL)
  ctx$report_b_csv <- reactiveVal(NULL)
  ctx$report_c_csv <- reactiveVal(NULL)

  ctx$analysis_stats_a <- reactiveVal(NULL)
  ctx$analysis_stats_b <- reactiveVal(NULL)
  ctx$analysis_stats_c <- reactiveVal(NULL)

  ctx$run_metadata_a <- reactiveVal(NULL)
  ctx$run_metadata_b <- reactiveVal(NULL)
  ctx$run_metadata_c <- reactiveVal(NULL)

  ctx$analysis_run_a <- reactiveVal(FALSE)
  ctx$analysis_run_b <- reactiveVal(FALSE)
  ctx$analysis_run_c <- reactiveVal(FALSE)

  ctx$manual_evidence_a <- reactiveVal(list())
  ctx$analysis_context_a <- reactiveVal(NULL)
  ctx$audit_data <- reactiveVal(data.frame())

  ensure_audit_log()
  ctx
}
