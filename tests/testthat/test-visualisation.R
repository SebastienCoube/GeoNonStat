set.seed(123)
obs_locs <- matrix(rnorm(200), ncol = 2)
vecchia_approx <- createVecchia(obs_locs)
pepito <- suppressMessages(
  createPP(vecchia_approx, plot = FALSE)
)

test_that("plotKnotsPP doesn't produce errors", {
  expect_error(
    res <- plotKnotsPP(pepito),
    NA
  )
  set.seed(123)
  mar_var <- rnorm(125)
  expect_error(
    res <- plotKnotsPP(pepito, mar_var_loss = mar_var),
    NA
  )
  expect_error(
    res <- plotKnotsPP(pepito, mar_var_loss = mar_var, show_knots = FALSE),
    NA
  )
})

test_that("plot doesn't produce errors", {
  expect_message(
    res <- plot(pepito, mar_var_loss = FALSE),
    NA
  )
  expect_message(
    res <- plot(pepito, mar_var_loss = TRUE),
    "1.7% of marginal variance"
  )
  expect_message(
    res <- plot(pepito, mar_var_loss = TRUE, separate = TRUE),
    "1.7% of marginal variance"
  )
  
  suppressMessages(mar_var <- varLossPP(pepito))
  expect_error(
    res <- plotKnotsPP(pepito, mar_var_loss = mar_var),
    NA
  )
  expect_error(
    res <- plotKnotsPP(pepito, mar_var_loss = mar_var, show_knots = FALSE),
    NA
  )
})


test_that("getColors produce expected output", {
  expect_error(
    res<- getColors(c(1,2,3)),
          NA
  )
  expect_identical(res, c("#B10026", "#D8806C", "#FFFFB2"))
})


test_that("getColorsCat produce expected output", {
  expect_error(
    res<- getColorsCat(3),
    NA
  )
  expect_identical(res, c("#000000", "#E73F74", "#F1CE63"))
})


test_that("plotPointillistPainting produce expected errors", {
  locs=c(2,3,6)
  field=c(1,2,3)
  if(!is.null(dev.list())) dev.off()
  expect_error(
    plotPointillistPainting(locs, field, add=TRUE),
    "plot.new has not been called yet"
  )
  expect_error(
    plotPointillistPainting(locs, field),
    NA
  )
})

test_that("pointillistColorscale produce expected output", {
  op <- par("mfrow")
  expect_error(
    pointillistColorscale(c(1,2,3)),
    NA
  )
  # Test that par is correctly reset
  testthat::expect_identical(op, par("mfrow"))
})
