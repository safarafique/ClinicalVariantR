#' Build the ClinicalVariantR Shiny application
#'
#' Returns a Shiny app object for ACMG/AMP germline variant classification.
#' Following Bioconductor Shiny guidelines, this function **returns** the app
#' and does not launch it. Call \code{shiny::runApp()} yourself.
#'
#' @param ... Reserved for future options (currently ignored).
#'
#' @return A \code{shiny.appobj} created by \code{shiny::shinyApp()}.
#'
#' @examples
#' stopifnot(is.function(ClinicalVariantR))
#' stopifnot(is.function(ClinicalVariantRApp))
#'
#' @seealso \code{\link{ClinicalVariantRApp}}
#' @export
ClinicalVariantR <- function(...) {
    .ClinicalVariantR_shiny_app()
}

#' Alias for \code{\link{ClinicalVariantR}}
#'
#' @inheritParams ClinicalVariantR
#' @inherit ClinicalVariantR return
#' @examples
#' stopifnot(is.function(ClinicalVariantRApp))
#' @export
ClinicalVariantRApp <- function(...) {
    ClinicalVariantR(...)
}

#' Locate Shiny app directory (installed inst/shinyapp or source tree).
#'
#' @keywords internal
#' @noRd
.ClinicalVariantR_app_dir <- function() {
    env_root <- Sys.getenv("CLINICALVARIANTR_APP_ROOT", unset = "")
    pkg_path <- tryCatch(
        find.package("ClinicalVariantR", quiet = TRUE),
        error = function(e) character()
    )
    if (!length(pkg_path)) pkg_path <- ""

    candidates <- c(
        env_root,
        system.file("shinyapp", package = "ClinicalVariantR"),
        if (nzchar(pkg_path)) file.path(pkg_path, "shinyapp") else "",
        system.file(package = "ClinicalVariantR"),
        pkg_path,
        normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    )

    for (root in unique(candidates[nzchar(candidates)])) {
        if (!dir.exists(root)) next
        if (file.exists(file.path(root, "global.R")) &&
            file.exists(file.path(root, "ui.R")) &&
            file.exists(file.path(root, "server.R"))) {
            return(normalizePath(root, winslash = "/", mustWork = FALSE))
        }
    }
    NULL
}

#' Internal: construct shinyApp from installed shinyapp/ or source tree.
#'
#' @keywords internal
#' @noRd
.ClinicalVariantR_shiny_app <- function() {
    app_dir <- .ClinicalVariantR_app_dir()
    if (is.null(app_dir)) {
        stop(
            "Unable to locate ClinicalVariantR Shiny sources (global.R / ui.R / server.R). ",
            "Install the package or set CLINICALVARIANTR_APP_ROOT / working directory ",
            "to the package root (or inst/shinyapp).",
            call. = FALSE
        )
    }

    Sys.setenv(CLINICALVARIANTR_APP_ROOT = app_dir)
    shiny::shinyAppDir(app_dir)
}

#' Read a VCF header via VariantAnnotation (keeps the Import wired for checks).
#'
#' @param path Character path to a VCF / VCF.gz file.
#' @return A \code{VCFHeader} object, or \code{NULL} on failure.
#' @keywords internal
#' @noRd
.ClinicalVariantR_scan_vcf_header <- function(path) {
    if (!is.character(path) || length(path) != 1L || !nzchar(path) || !file.exists(path)) {
        return(NULL)
    }
    hdr <- try(VariantAnnotation::scanVcfHeader(path), silent = TRUE)
    if (inherits(hdr, "try-error") || !methods::is(hdr, "VCFHeader")) {
        return(NULL)
    }
    hdr
}
