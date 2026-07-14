test_that("ACMGamp constructors exist", {
    expect_true(is.function(ACMGamp))
    expect_true(is.function(ACMGampApp))
})

test_that("ACMGamp returns a shiny.appobj when sources are available", {
    skip_if_not(
        file.exists("global.R") && file.exists("ui.R") && file.exists("server.R"),
        "Shiny entry files not in working directory / install layout yet"
    )
    app <- ACMGamp()
    expect_s3_class(app, "shiny.appobj")
})
