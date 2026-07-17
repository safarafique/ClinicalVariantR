# ClinicalVariantR — package-root launcher
#
# Prefer these (no package-R/ auto-load warning):
#   shiny::runApp("E:/ACGM/ClinicalVariantR/inst/shinyapp", launch.browser = TRUE)
#   library(ClinicalVariantR); shiny::runApp(ClinicalVariantR(), launch.browser = TRUE)
#
# This file lets shiny::runApp("E:/ACGM/ClinicalVariantR") still work by
# delegating to inst/shinyapp (which is not an R package directory).

options(shiny.autoload.r = FALSE)

.app_dir <- file.path("inst", "shinyapp")
if (!file.exists(file.path(.app_dir, "server.R")) ||
    !file.exists(file.path(.app_dir, "ui.R")) ||
    !file.exists(file.path(.app_dir, "global.R"))) {
  stop(
    "Cannot find inst/shinyapp/{global,ui,server}.R.\n",
    "Run shiny::runApp() from the ClinicalVariantR package root, or use:\n",
    "  shiny::runApp('E:/ACGM/ClinicalVariantR/inst/shinyapp', launch.browser = TRUE)",
    call. = FALSE
  )
}

.app_dir <- normalizePath(.app_dir, winslash = "/", mustWork = TRUE)
Sys.setenv(CLINICALVARIANTR_APP_ROOT = .app_dir)
shiny::shinyAppDir(.app_dir)
