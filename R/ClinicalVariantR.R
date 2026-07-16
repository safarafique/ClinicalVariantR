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
        normalizePath(getwd(), winslash = "/", mustWork = FALSE),
        # Common local checkout paths (development fallback)
        "E:/ACGM/ClinicalVariantR/inst/shinyapp",
        "E:/ACGM/ClinicalVariantR",
        "e:/ACGM/ClinicalVariantR/inst/shinyapp",
        "e:/ACGM/ClinicalVariantR"
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
            "Unable to locate ClinicalVariantR Shiny sources (global.R / ui.R / server.R).\n",
            "The package was likely not reinstalled (it was still loaded).\n",
            "Fix:\n",
            "  1) Session -> Restart R\n",
            "  2) remotes::install_local('E:/ACGM/ClinicalVariantR', force = TRUE, upgrade = 'never')\n",
            "  3) library(ClinicalVariantR); shiny::runApp(ClinicalVariantR(), launch.browser = TRUE)\n",
            "Or launch from the source tree without reinstalling:\n",
            "  shiny::runApp('E:/ACGM/ClinicalVariantR', launch.browser = TRUE)",
            call. = FALSE
        )
    }

    Sys.setenv(CLINICALVARIANTR_APP_ROOT = app_dir)
    shiny::shinyAppDir(app_dir)
}
