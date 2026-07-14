#' @title ACMGamp: ACMG/AMP germline variant interpretation
#'
#' @description
#' ACMGamp classifies germline sequence variants under ACMG/AMP 2015 guidelines
#' from VEP-, SnpEff-, or ANNOVAR-annotated VCF files. The package provides an
#' interactive Shiny application with three workflows (full clinical, automated
#' rapid, and gene-panel) plus streaming analysis for large call sets.
#'
#' @details
#' Launch the app with:
#' \preformatted{
#' app <- ACMGamp()
#' if (interactive()) {
#'   shiny::runApp(app)
#' }
#' }
#'
#' See the package vignette for workflows, input requirements, and Bioconductor
#' packaging notes.
#'
#' @seealso \code{\link{ACMGamp}}, \code{\link{ACMGampApp}}
#'
#' @keywords internal
"_PACKAGE"
