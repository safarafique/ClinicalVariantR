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
#' is.function(ClinicalVariantR)
#' is.function(ClinicalVariantRApp)
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
#' is.function(ClinicalVariantRApp)
#' @export
ClinicalVariantRApp <- function(...) {
    ClinicalVariantR(...)
}

#' Internal: construct shinyApp from existing ui/server sources
#'
#' @keywords internal
#' @noRd
.ClinicalVariantR_shiny_app <- function() {
    pkg_root <- system.file(package = "ClinicalVariantR")
    # During development before install, fall back to source tree.
    if (!nzchar(pkg_root) || !dir.exists(file.path(pkg_root, "R"))) {
        pkg_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    }

    # Prefer installed layout: UI/server live as package sources; for the
    # transitional Shiny-app repo, source app entry files from package root.
    global_file <- file.path(pkg_root, "global.R")
    ui_file <- file.path(pkg_root, "ui.R")
    server_file <- file.path(pkg_root, "server.R")

    if (file.exists(global_file) && file.exists(ui_file) && file.exists(server_file)) {
        env <- new.env(parent = globalenv())
        sys.source(global_file, envir = env)
        sys.source(ui_file, envir = env)
        sys.source(server_file, envir = env)
        return(shiny::shinyApp(ui = env$ui, server = env$server))
    }

    stop(
        "Unable to locate ClinicalVariantR Shiny sources (global.R / ui.R / server.R). ",
        "Install the package or set the working directory to the package root.",
        call. = FALSE
    )
}
