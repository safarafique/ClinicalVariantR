#!/usr/bin/env Rscript
# Install all ClinicalVariantR app dependencies (CRAN + Bioconductor).
# Usage:
#   Rscript scripts/install_app_deps.R
# Or from R:
#   source("R/bootstrap_deps.R"); clinicalvariantr_ensure_dependencies()

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
root <- normalizePath(file.path(script_dir, ".."))
setwd(root)

Sys.unsetenv("CLINICALVARIANTR_NO_AUTO_INSTALL")
boot <- file.path(root, "scripts", "bootstrap_deps.R")
if (!file.exists(boot)) boot <- file.path(root, "R", "bootstrap_deps.R")
source(boot, local = FALSE)
clinicalvariantr_ensure_dependencies(auto_install = TRUE)
cat("OK: app dependencies installed.\n")
cat("Launch with:\n")
cat("  shiny::runApp(\"", root, "\", launch.browser = TRUE)\n", sep = "")
