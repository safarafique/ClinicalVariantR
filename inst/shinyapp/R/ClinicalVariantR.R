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
#' # To launch interactively:
#' # shiny::runApp(ClinicalVariantR())
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
#' # To launch interactively:
#' # shiny::runApp(ClinicalVariantRApp())
#' @export
ClinicalVariantRApp <- function(...) {
    ClinicalVariantR(...)
}

#' Locate Shiny app directory (installed inst/shinyapp or source tree).
#'
#' @keywords internal
#' @noRd
.ClinicalVariantR_app_dir <- function() {
    candidates <- c(
        system.file("shinyapp", package = "ClinicalVariantR"),
        system.file(package = "ClinicalVariantR"),
        normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    )

    # When running from a source checkout of this repo, also try typical roots.
    pkg_path <- tryCatch(
        find.package("ClinicalVariantR", quiet = TRUE),
        error = function(e) character()
    )
    if (length(pkg_path) == 1L && nzchar(pkg_path)) {
        candidates <- c(
            file.path(pkg_path, "shinyapp"),
            pkg_path,
            candidates
        )
    }

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
            "Reinstall with a package that includes inst/shinyapp/, or set the working ",
            "directory to the package source root and run shiny::runApp('.').",
            call. = FALSE
        )
    }

    Sys.setenv(CLINICALVARIANTR_APP_ROOT = app_dir)
    shiny::shinyAppDir(app_dir)
}
