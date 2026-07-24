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
#' Imports are selective (\code{@importFrom}) so shiny / DT / jsonlite do not
#' clash on \code{renderDataTable}, \code{dataTableOutput}, or \code{validate}.
#' Shiny UI/server code lives under \code{inst/shinyapp/} and is not loaded into
#' the package namespace.
#'
#' @importFrom shiny shinyAppDir
#' @importFrom bslib bs_theme page_navbar
#' @importFrom DT datatable renderDT DTOutput formatStyle styleEqual formatRound
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom data.table fread as.data.table fifelse set
#' @importFrom readr read_csv
#' @importFrom VariantAnnotation scanVcfHeader
#' @importFrom methods is
#' @importFrom digest digest
#' @importFrom openssl aes_gcm_encrypt aes_gcm_decrypt rand_bytes
#' @importFrom stats setNames
#' @importFrom utils packageVersion installed.packages modifyList read.csv read.delim write.csv
#'
#' @keywords internal
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        "ClinicalVariantR ", as.character(utils::packageVersion("ClinicalVariantR")),
        "\nLaunch: shiny::runApp(ClinicalVariantR())",
        "\nInstall once with BiocManager (dependencies = TRUE), then library(ClinicalVariantR)."
    )
}

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
    "VARIANT_DETAIL_COLUMNS",
    "session"
))
