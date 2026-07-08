#!/usr/bin/env Rscript
# Install R packages needed for CLI scripts (WSL/Linux, no Shiny).
# Usage: Rscript scripts/install_r_cli_deps.R

repos <- "https://cloud.r-project.org"
pkgs <- c("data.table", "readr", "jsonlite")

ensure_user_library <- function() {
  user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
  if (!nzchar(user_lib)) {
    rver <- paste(R.version$major, R.version$minor, sep = ".")
    user_lib <- file.path(
      Sys.getenv("HOME"),
      "R",
      paste0(R.version$platform, "-library"),
      rver
    )
  }
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
  user_lib
}

user_lib <- ensure_user_library()
cat("Using R library:", user_lib, "\n")

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) == 0L) {
  cat("CLI dependencies already installed:", paste(pkgs, collapse = ", "), "\n")
  quit(status = 0)
}

cat("Installing:", paste(missing, collapse = ", "), "\n")
install.packages(missing, repos = repos, lib = user_lib, dependencies = c("Depends", "Imports"))

still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0L) {
  stop("Failed to install: ", paste(still_missing, collapse = ", "), call. = FALSE)
}

cat("Done. Re-run:\n")
cat("  Rscript scripts/verify_group_a_28.R\n")
cat("  Rscript scripts/verify_group_b_c.R\n")
cat("  Rscript scripts/generate_validation_report.R\n")
