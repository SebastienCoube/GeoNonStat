# knotsFromKmeans tests
test_that("knotsFromKmeans returns the expected matrix shape", {
  set.seed(123)
  locs <- matrix(runif(20), ncol = 2)
  centers <- knotsFromKmeans(2, locs)
  expect_equal(dim(centers), c(2, 2))
  expect_true(is.matrix(centers))
})

test_that("knotsFromKmeans rejects too many knots", {
  locs <- matrix(runif(4), ncol = 2)
  expect_error(knotsFromKmeans(2, locs), "49999 knots maximum allowed")
})

# createPP validation tests
test_that("createPP rejects invalid vecchia_approx and matern_range", {
  expect_error(createPP(vecchia_approx = NULL, plot = FALSE, verbose = FALSE), "Argument 'vecchia_approx' must be a list")

  vecchia_approx <- list(
    locs = matrix(runif(8), ncol = 2),
    NNarray = matrix(NA_integer_, nrow = 2, ncol = 1),
    n_locs = 2
  )
  expect_error(createPP(vecchia_approx = vecchia_approx, matern_range = -1, knots = 1, plot = FALSE, verbose = FALSE), "'matern_range' must be positive")
})

test_that("createPP rejects invalid knot matrix dimensions", {
  vecchia_approx <- list(
    locs = matrix(runif(8), ncol = 2),
    NNarray = matrix(NA_integer_, nrow = 2, ncol = 1),
    n_locs = 2
  )
  bad_knots <- matrix(runif(9), ncol = 3)
  expect_error(createPP(vecchia_approx = vecchia_approx, knots = bad_knots, plot = FALSE, verbose = FALSE), "Matrix 'knots' must have the same number of columns as spatial locations")
})

test_that("createPP return expected format", {
  vecchia_approx <- list(
    locs = matrix(runif(8), ncol = 2),
    NNarray = matrix(NA_integer_, nrow = 2, ncol = 1),
    n_locs = 2
  )
  expect_message(
    expect_message(
      pptest <- expect_error(createPP(vecchia_approx = vecchia_approx, 
                            knots = 2, 
                            plot = FALSE, verbose = FALSE),
                   NA),
      "knot placement done by default using"
    ),"Matern range"
  )
  expect_true(is(pptest, "PP"))
  expect_named(pptest, c("knots", "matern_range", "sparse_chol", "t_sparse_chol", "n_knots", "vecchia_locs"))
  expect_identical(pptest$n_knots, 2L)
  expect_true(is.matrix(pptest$knots) && identical(dim(pptest$knots),c(2L, 2L)))
  expect_true(is.numeric(pptest$matern_range) & pptest$matern_range>0)
  expect_true(is(pptest$sparse_chol, "Matrix") & identical(dim(pptest$sparse_chol), c(6L, 6L)))
  expect_identical(pptest$sparse_chol, t(pptest$t_sparse_chol))
  expect_true(is.matrix(pptest$vecchia_locs) && identical(dim(pptest$vecchia_locs),c(4L, 2L)))
})


# varLossPP tests
test_that("varLossPP computes loss for a valid PP object", {
  sparse_chol <- Matrix::sparseMatrix(
    i = c(1, 2, 2),
    j = c(1, 1, 2),
    x = c(1, 0.2, 1),
    triangular = TRUE
  )
  pp <- list(
    knots = matrix(1:2, ncol = 2),
    sparse_chol = sparse_chol
  )
  res <- varLossPP(pp, verbose = FALSE)
  expect_equal(length(res), nrow(sparse_chol))
  expect_true(all(res >= 0))
})

