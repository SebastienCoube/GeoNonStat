test_that("plot.PP executes without error for a minimal PP object", {
  pp <- structure(
    list(
      knots = matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE),
      matern_range = 1,
      sparse_chol = Matrix::Diagonal(4),
      t_sparse_chol = Matrix::Diagonal(4),
      n_knots = 2,
      vecchia_locs = matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
    ),
    class = "PP"
  )

  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  expect_error(plot.PP(pp, mar_var_loss = TRUE, samples = FALSE), NA)
  dev.off()
})

test_that("plotKnotsPP works when mar_var_loss is provided", {
  pp <- structure(
    list(
      knots = matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE),
      vecchia_locs = matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
    ),
    class = "PP"
  )
  mar_var_loss <- c(10, 20, 30, 40)
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  expect_error(plotKnotsPP(pp, mar_var_loss = mar_var_loss, show_knots = TRUE), NA)
  dev.off()
})

test_that("plotPointillistPainting works for a simple dataset", {
  locs <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  field <- c(0.1, 0.4, 0.7, 0.9)
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  expect_error(plotPointillistPainting(locs, field, add = FALSE, legend = FALSE), NA)
  dev.off()
})

test_that("pointillistColorscale handles non-empty field", {
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  expect_error(pointillistColorscale(c(0.1, 0.5, 0.9)), NA)
  dev.off()
})

test_that("plotEllipses works for 1-column and 3-column log_range", {
  locs <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  expect_error(plotEllipses(locs, matrix(c(log(0.5), log(0.8)), ncol = 1), shrink = 0.1), NA)
  expect_error(plotEllipses(locs, matrix(rnorm(6), ncol = 3), shrink = 0.1), NA)
  dev.off()
})


