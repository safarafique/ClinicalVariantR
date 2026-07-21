#' @title ClinicalVariantR: ACMG/AMP germline variant interpretation
#'
#' @description
#' ClinicalVariantR classifies germline sequence variants under ACMG/AMP 2015
#' guidelines from VEP-, SnpEff-, or ANNOVAR-annotated VCF files. The package
#' provides an interactive Shiny application with three workflows (full clinical,
#' automated rapid, and gene-panel) plus streaming analysis for large call sets.
#'
#' @details
#' Launch the app with:
#' \preformatted{
#' app <- ClinicalVariantR()
#' if (interactive()) {
#'   shiny::runApp(app)
#' }
#' }
#'
#' See the package vignette for workflows, input requirements, and Bioconductor
#' packaging notes.
#'
#' @seealso \code{\link{ClinicalVariantR}}, \code{\link{ClinicalVariantRApp}}
#'
#' @import shiny
#' @importFrom bslib bs_theme page_navbar
#' @importFrom VariantAnnotation scanVcfHeader
#' @importFrom methods is
#' @importFrom digest digest
#' @importFrom openssl aes_gcm_encrypt aes_gcm_decrypt rand_bytes
#' @importFrom stats setNames
#' @importFrom utils packageVersion installed.packages modifyList read.csv read.delim write.csv
#'
#' @keywords internal
"_PACKAGE"

# Shiny app symbols assigned in inst/shinyapp/global.R and cross-file helpers.
utils::globalVariables(c(
    "ACMG_CLASSIFICATIONS",
    "ACMG_GUIDELINE_VERSION",
    "ACMG_PRO_ENGINE",
    "APP_VERSION",
    "AUDIT_LOG_PATH",
    "EVIDENCE_DETAIL_COLUMNS",
    "PDF_EXPORT_COLUMNS",
    "REFERENCE_PATHS",
    "REPORT_COLUMNS",
    "VARIANT_DETAIL_COLUMNS"
))
