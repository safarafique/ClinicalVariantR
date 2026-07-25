#' Ensure ClinicalVariantR runtime packages are installed (new-user bootstrap).
#'
#' Sourced by global.R before modules load. Missing packages are installed
#' automatically unless CLINICALVARIANTR_NO_AUTO_INSTALL=1.
#' @noRd

.clinicalvariantr_cran_packages <- function() {
  c(
    # Declared Imports (app runtime)
    "shiny",
    "bslib",
    "DT",
    "crosstalk",
    "data.table",
    "digest",
    "readr",
    "jsonlite",
    "openssl",
    # Common DT / shiny transitive deps that break launch if missing
    "htmltools",
    "htmlwidgets",
    "httpuv",
    "promises",
    "later",
    "rlang",
    "jquerylib",
    "sass",
    "fontawesome",
    "yaml",
    "mime",
    "R6",
    "cli",
    "glue",
    "lifecycle",
    "fastmap",
    "cachem",
    "memoise",
    "withr",
    "magrittr",
    "stringr",
    "vctrs",
    "tidyselect",
    "cpp11",
    "bit64",
    "bit"
  )
}

.clinicalvariantr_bioc_packages <- function() {
  c("VariantAnnotation")
}

.clinicalvariantr_auto_install_enabled <- function() {
  !identical(Sys.getenv("CLINICALVARIANTR_NO_AUTO_INSTALL", unset = ""), "1")
}

.clinicalvariantr_ensure_user_library <- function() {
  libs <- .libPaths()
  writable <- vapply(libs, function(p) {
    dir.exists(p) && file.access(p, 2L) == 0L
  }, logical(1))
  if (any(writable)) return(libs[[which(writable)[1L]]])

  user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
  if (!nzchar(user_lib)) {
    rver <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
    user_lib <- file.path(
      path.expand("~"),
      "R",
      paste0(R.version$platform, "-library"),
      rver
    )
  }
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
  user_lib
}

.clinicalvariantr_missing_pkgs <- function(pkgs) {
  pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
}

.clinicalvariantr_cran_repos <- function() {
  repos <- getOption("repos")
  cran <- if (!is.null(repos) && "CRAN" %in% names(repos)) unname(repos[["CRAN"]]) else NA_character_
  if (is.null(repos) || is.na(cran) || !nzchar(cran) || identical(cran, "@CRAN@")) {
    return(c(CRAN = "https://cloud.r-project.org"))
  }
  repos
}

.clinicalvariantr_install_cran <- function(pkgs, lib) {
  if (length(pkgs) == 0L) return(invisible(TRUE))
  message("ClinicalVariantR: installing CRAN package(s): ", paste(pkgs, collapse = ", "))
  utils::install.packages(
    pkgs,
    lib = lib,
    repos = .clinicalvariantr_cran_repos(),
    dependencies = c("Depends", "Imports", "LinkingTo"),
    quiet = FALSE
  )
  invisible(TRUE)
}

.clinicalvariantr_install_bioc <- function(pkgs, lib) {
  if (length(pkgs) == 0L) return(invisible(TRUE))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("ClinicalVariantR: installing BiocManager...")
    utils::install.packages(
      "BiocManager",
      lib = lib,
      repos = .clinicalvariantr_cran_repos(),
      quiet = FALSE
    )
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    warning(
      "BiocManager could not be installed; skipping Bioconductor packages: ",
      paste(pkgs, collapse = ", "),
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  message("ClinicalVariantR: installing Bioconductor package(s): ", paste(pkgs, collapse = ", "))
  BiocManager::install(
    pkgs,
    lib = lib,
    update = FALSE,
    ask = FALSE,
    dependencies = TRUE
  )
  invisible(TRUE)
}

#' Install any missing app dependencies, then verify.
#' @noRd
clinicalvariantr_ensure_dependencies <- function(auto_install = NULL) {
  if (is.null(auto_install)) auto_install <- .clinicalvariantr_auto_install_enabled()

  cran_pkgs <- .clinicalvariantr_cran_packages()
  bioc_pkgs <- .clinicalvariantr_bioc_packages()

  missing_cran <- .clinicalvariantr_missing_pkgs(cran_pkgs)
  missing_bioc <- .clinicalvariantr_missing_pkgs(bioc_pkgs)

  if (length(missing_cran) == 0L && length(missing_bioc) == 0L) {
    return(invisible(TRUE))
  }

  if (!isTRUE(auto_install)) {
    stop(
      "Missing required package(s): ",
      paste(c(missing_cran, missing_bioc), collapse = ", "),
      ".\nSet CLINICALVARIANTR_NO_AUTO_INSTALL unset (default) to auto-install,\n",
      "or install manually:\n  install.packages(c(\"",
      paste(missing_cran, collapse = "\", \""),
      "\"))\n  BiocManager::install(c(\"",
      paste(missing_bioc, collapse = "\", \""),
      "\"))",
      call. = FALSE
    )
  }

  lib <- .clinicalvariantr_ensure_user_library()
  message("ClinicalVariantR: using library ", lib)

  .clinicalvariantr_install_cran(missing_cran, lib)
  .clinicalvariantr_install_bioc(missing_bioc, lib)

  still_cran <- .clinicalvariantr_missing_pkgs(cran_pkgs)
  still_bioc <- .clinicalvariantr_missing_pkgs(bioc_pkgs)
  if (length(still_cran) > 0L || length(still_bioc) > 0L) {
    stop(
      "Failed to install required package(s): ",
      paste(c(still_cran, still_bioc), collapse = ", "),
      ".\nTry in R:\n  source(\"scripts/bootstrap_deps.R\"); clinicalvariantr_ensure_dependencies()\n",
      "or:\n  Rscript scripts/install_app_deps.R",
      call. = FALSE
    )
  }

  message("ClinicalVariantR: all required packages are available.")
  invisible(TRUE)
}
