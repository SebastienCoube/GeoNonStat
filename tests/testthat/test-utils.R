set.seed(123)

make_test_vecchia <- function(n) {
  observed_locs <- cbind(runif(n), runif(n))
  return(createVecchia(observed_locs, m = 10))
}

make_test_PP <- function(vecchia, knots=2) {
  suppressMessages(
    createPP(vecchia_approx = vecchia, knots = knots)
  )
}


vecchia_test <- make_test_vecchia(100)
pepito_test <- make_test_PP(vecchia_test, knots = 10)


test_that("naiveGreedyColoring produce expected results", {
  n <- 5
  i <- c(rep(1, n), 2:n) # 1re ligne + 1re colonne sauf la 1re case
  j <- c(1:n, rep(1, n - 1))
  M <- sparseMatrix(i = i, j = j, x = 1, dims = c(n, n))
  expect_error(
    tt <- naiveGreedyColoring(M),
    NA
  )
  expect_identical(tt, c(1, 2, 2, 2, 2))
  M[2, 3] <- 1
  M[3, 2] <- 1
  expect_identical(
    naiveGreedyColoring(M),
    c(1, 2, 3, 2, 2)
  )
  i <- c(1, 3:8)
  j <- c(2, 9, 6:10)
  x <- 7 * (1:7)
  A <- sparseMatrix(i, j, x = x)
  expect_error(
    naiveGreedyColoring(A),
    "M must be symmetric"
  )
})

test_that("decompressChol produce expected results", {
  # TODO
})

test_that("expmat produce expected output", {
  coords <- c(1, 2, 3, 4, 5, 6)
  expect_error(
    res <- expmat(coords),
    NA
  )
  expect_true(is(res, "array"))
  expect_identical(dim(res), c(3L, 3L))

  expected_mat <- matrix(
    c(47859.7700011562, 55641.9937959158, 62462.1209161336,
      55641.9937959158, 64689.8331149192, 72618.8933675331,
      62462.1209161336, 72618.8933675331, 81519.8962286822),
    nrow=3)
  expect_equal(res, expected_mat)
  expect_error(
    res <- expmat(coords),
    NA
  )
  diag(expected_mat) <- diag(expected_mat) - 0.0001
  expect_equal(res, expected_mat)
})

test_that("symmat produce expected output", {
  # 1x1 matrix
  coords <- c(1)
  expect_error(
    res <- symmat(coords),
    NA
  )
  expect_true(is(res, "matrix"))
  expect_identical(dim(res), c(1L, 1L))
  expect_equal(res, matrix(c(1)))

  # 2x2 matrix
  coords <- c(1, 2, 3)
  expect_error(
    res <- symmat(coords),
    NA
  )
  expect_true(is(res, "matrix"))
  expect_identical(dim(res), c(2L, 2L))
  expect_equal(
    res,
    matrix(c(1, 3, 3, 2),
      nrow = 2
    )
  )

  # 3x3 matrix
  coords <- c(1, 2, 3, 4, 5, 6)
  expect_error(
    res <- symmat(coords),
    NA
  )
  expect_true(is(res, "matrix"))
  expect_identical(dim(res), c(3L, 3L))
  expect_equal(
    res,
    matrix(c(1, 4, 5, 4, 2, 6, 5, 6, 3),
      nrow = 3
    )
  )

  # Vector incompatible
  coords <- c(1, 2, 3, 4, 5)
  expect_error(
    res <- symmat(coords),
    "length of coords incompatible with a symetric matrix"
  )
})


## betaPriorLogDens ##################################
test_that("betaPriorLogDens produce expected output", {
  set.seed(123)
  beta1 <- matrix(rnorm(100), 100, ncol = 1)
  beta3 <- matrix(rnorm(300), 100, ncol = 3)
  expect_error(
    tmp <- betaPriorLogDens(
      beta = beta3,
      n_PP = 120,
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3)
    ),
    "n_PP can't be greater than"
  )
  expect_error(
    tmp <- betaPriorLogDens(
      beta = beta3,
      n_PP = 90,
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3)
    ),
    "log_scale is supposed to be of length 2"
  )

  expect_error(
    tmp <- betaPriorLogDens(
      beta = beta3,
      n_PP = 90,
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3, -2)
    ),
    NA
  )
  expect_true(inherits(tmp, "numeric"))
  # expect_equal(tmp, -2697.276)

  expect_error(
    tmp <- betaPriorLogDens(
      beta = beta1,
      n_PP = 90,
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3)
    ),
    NA
  )
  expect_true(inherits(tmp, "numeric"))
  # expect_equal(tmp, -1020.22585)
})

## betaPriorLogDensDerivative ##################################
test_that("betaPriorLogDensDerivative produce expected output", {
  set.seed(123)
  beta1 <- matrix(rnorm(100), 100, ncol = 1)
  beta3 <- matrix(rnorm(300), 100, ncol = 3)
  expect_error(
    tmp <- betaPriorLogDensDerivative(
      beta = beta3,
      n_PP = 120,
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3, 2)
    ),
    "n_PP can't be greater than"
  )
  expect_error(
    tmp <- betaPriorLogDensDerivative(
      beta = beta3,
      n_PP = 90,
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3)
    ),
    "log_scale is supposed to be of length 2"
  )

  expect_error(
    tmp <- betaPriorLogDensDerivative(
      beta = beta3,
      n_PP = 90,
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3, -2)
    ),
    NA
  )
  expect_true(inherits(tmp, "matrix"))
  expect_identical(dim(tmp), dim(beta3))
  # expect_equal(mean(tmp), 1.00837737)

  expect_error(
    tmp <- betaPriorLogDensDerivative(
      beta = beta1,
      n_PP = 90,
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3)
    ),
    NA
  )
  expect_true(inherits(tmp, "matrix"))
  expect_identical(dim(tmp), dim(beta1))
  # expect_equal(mean(tmp), -2.9948913)
})

test_that("getting correct values using betaPriorLogDens and betaPriorLogDensDerivative", {
  set.seed(123)
  beta <- matrix(rnorm(300), 100)
  res <- matrix(0, 100, 3)
  for (i in seq(100)) {
    for (j in seq(3)) {
      beta1 <- beta
      beta1[i, j] <- beta1[i, j] + .0001
      res[i, j] <- (
        betaPriorLogDens(beta1, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = c(0, 2)) -
          betaPriorLogDens(beta, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = c(0, 2))
      ) * 10000 / betaPriorLogDensDerivative(beta, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = c(0, 2))[i, j]
    }
  }
  expect_equal(mean(res), 1.00004157)
  expect_equal(var(c(res)), 3.722309e-07, tolerance = 1e-7)
})

## xPPCrossprod ##################################
test_that("xPPCrossprod Simple crossprod if no PP, vecchia unused", {
  X <- matrix(rnorm(100), 50)
  Y <- matrix(rnorm(30 * nrow(X)), nrow(X))
  expect_error(
    res1 <- xPPCrossprod(
      X = X,
      PP = NULL,
      Y = Y,
      vecchia_approx = vecchia_approx
    ),
    NA
  )

  expect_true(is(res1, "matrix"))
  expect_identical(dim(res1), c(2L, 30L))
  expect_identical(res1, crossprod(X, Y))

  # Simple corssprod, vecchia unused.
  expect_error(
    res2 <- xPPCrossprod(X = X, PP = NULL, Y = Y, vecchia_approx = NULL),
    NA
  )
  expect_identical(res1, res2)
})

test_that("xPPCrossprod with PP", {
  set.seed(123)
  X <- matrix(rnorm(500), nrow=100)
  Y <- matrix(rnorm(30 * nrow(X)), nrow(X))
  expect_error(
    res2 <- xPPCrossprod(X = X, PP = pepito_test, Y = Y, vecchia_approx = NULL),
    "vecchia_approx should be provided too"
  )
})

test_that("xPPCrossprod with PP and permute obs", {
  set.seed(123)
  X <- matrix(rnorm(500), nrow=100)
  Y <- matrix(rnorm(30 * nrow(X)), nrow(X))
  expect_error(
    res3 <- xPPCrossprod(X = X, PP = pepito_test, Y = Y, 
                         vecchia_approx = vecchia_test, 
                         permutate_PP_to_obs = TRUE),
    NA
  )
  expect_true(is(res3, "matrix"))
  expect_identical(dim(res3), c(15L, 30L))

  # First line is crossprod
  expect_identical(res3[1, ], crossprod(X, Y)[1, ])
})

# ## xPPMultRight ##################################
# test_that("xPPMultRight with PP and permute obs", {
#   set.seed(123)
#   X <- matrix(rnorm(100), nrow = 50)
#   Y <- matrix(rnorm(30 * nrow(X)), nrow(X))
#   expect_error(
#     resmr <- xPPMultRight(
#       X = X,
#       PP = pepito_test,
#       Y = Y,
#       vecchia_approx = vecchia_test, permutate_PP_to_obs = FALSE
#     ),
#     "Y should have"
#   )
#   # No PP
#   expect_error(
#     resmr <- xPPMultRight(
#       X = X,
#       PP = NULL,
#       Y = matrix(rep(c(0.5, 1, 1.1, 0.2), 2), nrow = 2),
#       vecchia_approx = vecchia_approx,
#       permutate_PP_to_obs = FALSE
#     ),
#     NA
#   )
#   expect_true(is(resmr, "matrix"))
#   expect_identical(dim(resmr), c(50L, 4L))
#   # expect_equal(mean(resmr), 0.1153678, tolerance = 1e-5)
# 
#   # No X
#   expect_error(
#     resmr <- xPPMultRight(X = NULL, PP = pepito_test, 
#                           Y = matrix(rnorm(250), nrow = 10), 
#                           vecchia_approx = vecchia_test, 
#                           permutate_PP_to_obs = FALSE),
#     NA
#   )
#   expect_true(is(resmr, "matrix"))
#   expect_identical(dim(resmr), c(100L, 25L))
#   # expect_equal(mean(as.matrix(resmr)), -0.160772, tolerance = 1e-5)
# 
#   # X and PP
#   expect_error(
#     resmr <- xPPMultRight(X = X, PP = pepito_test,
#                           Y = matrix(nrow(X)*5, nrow = nrow(X)), 
#                           vecchia_approx = vecchia_test, 
#                           permutate_PP_to_obs = FALSE),
#     NA
#   )
#   expect_true(is(resmr, "matrix"))
#   expect_identical(dim(resmr), c(100L, 25L))
#   # expect_equal(mean(as.matrix(resmr)), -0.03502011, tolerance = 1e-5)
# })
