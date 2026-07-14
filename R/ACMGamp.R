#' Build the ACMGamp Shiny application
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
#' \dontrun{
#' app <- ACMGamp()
#' if (interactive()) {
#'     shiny::runApp(app)
#' }
#' }
#'
#' @seealso \code{\link{ACMGampApp}}
#' @export
ACMGamp <- function(...) {
    .ACMGamp_shiny_app()
}

#' Alias for \code{\link{ACMGamp}}
#'
#' @inheritParams ACMGamp
#' @inherit ACMGamp return
#' @examples
#' \dontrun{
#' app <- ACMGampApp()
#' if (interactive()) {
#'     shiny::runApp(app)
#' }
#' }
#' @export
ACMGampApp <- function(...) {
    ACMGamp(...)
}

#' Internal: construct shinyApp from existing ui/server sources
#'
#' @keywords internal
#' @noRd
.ACMGamp_shiny_app <- function() {
    pkg_root <- system.file(package = "ACMGamp")
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
        "Unable to locate ACMGamp Shiny sources (global.R / ui.R / server.R). ",
        "Install the package or set the working directory to the package root.",
        call. = FALSE
    )
}
