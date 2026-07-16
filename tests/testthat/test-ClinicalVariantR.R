test_that("ClinicalVariantR constructors exist", {
    expect_true(is.function(ClinicalVariantR))
    expect_true(is.function(ClinicalVariantRApp))
})

test_that("ClinicalVariantR returns a shiny.appobj when sources are available", {
    app_dir <- system.file("shinyapp", package = "ClinicalVariantR")
    has_installed_app <- nzchar(app_dir) &&
        file.exists(file.path(app_dir, "global.R")) &&
        file.exists(file.path(app_dir, "ui.R")) &&
        file.exists(file.path(app_dir, "server.R"))
    has_source_app <- file.exists("global.R") &&
        file.exists("ui.R") &&
        file.exists("server.R")
    skip_if_not(
        has_installed_app || has_source_app,
        "Shiny entry files not in install layout (inst/shinyapp) or working directory"
    )
    app <- ClinicalVariantR()
    expect_s3_class(app, "shiny.appobj")
})
