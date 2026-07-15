test_that("ClinicalVariantR constructors exist", {
    expect_true(is.function(ClinicalVariantR))
    expect_true(is.function(ClinicalVariantRApp))
})

test_that("ClinicalVariantR returns a shiny.appobj when sources are available", {
    skip_if_not(
        file.exists("global.R") && file.exists("ui.R") && file.exists("server.R"),
        "Shiny entry files not in working directory / install layout yet"
    )
    app <- ClinicalVariantR()
    expect_s3_class(app, "shiny.appobj")
})
