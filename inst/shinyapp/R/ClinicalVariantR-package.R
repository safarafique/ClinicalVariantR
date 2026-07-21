#' @title ClinicalVariantR: ACMG/AMP germline variant interpretation
#'
#' @description
#' ClinicalVariantR classifies germline sequence variants under ACMG/AMP 2015 guidelines
#' from VEP-, SnpEff-, or ANNOVAR-annotated VCF files. The package provides an
#' interactive Shiny application with three workflows (full clinical, automated
#' rapid, and gene-panel) plus streaming analysis for large call sets.
#'
#' @details
#' Install from GitHub and launch the app with:
#' \preformatted{
#' remotes::install_github("safarafique/ClinicalVariantR")
#' library(ClinicalVariantR)
#' app <- ClinicalVariantR()
#' if (interactive()) {
#'   shiny::runApp(app, launch.browser = TRUE)
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
