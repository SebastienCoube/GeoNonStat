set.seed(123)
locs <- cbind(runif(1000), runif(1000))
vecchia_approx = createVecchia(locs, 10, ncores = 1)

set.seed(123)
PP = createPP(vecchia_approx, plot=FALSE)
X = matrix(rnorm(500), nrow(locs))
Y = matrix(rnorm(30*nrow(X)), nrow(X))
range_X = matrix(1, nrow(locs))

test_that("naive_greedy_coloring produce expected results", {
  n <- 5  
  i <- c(rep(1, n), 2:n)       # 1re ligne + 1re colonne sauf la 1re case
  j <- c(1:n, rep(1, n - 1))
  M <- sparseMatrix(i = i, j = j, x = 1, dims = c(n, n))
  expect_error(
    tt <- naive_greedy_coloring(M), 
    NA
  )
  expect_identical(tt, c(1, 2, 2, 2, 2))
  M[2,3]=1
  M[3,2]=1
  expect_identical(naive_greedy_coloring(M),
                   c(1, 2, 3, 2, 2))
  i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
  A <- sparseMatrix(i, j, x = x)
  expect_error(naive_greedy_coloring(A),
               "M must be symmetric")
})

test_that("decompress_chol produce expected results", {
# TODO  
})

test_that("expmat produce expected output", {
  coords <- c(1,2,3,4,5,6)
  expect_error(
    res <- expmat(coords),
    NA
  )
  expect_true(is(res, "array"))
  expect_identical(dim(res), c(3L,3L))
  # 
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
  # coords <- c(1)
  # expect_error(
  #   res <- symmat(coords), 
  #   NA
  # )
  # expect_true(is(res, "matrix"))
  # expect_identical(dim(res), c(1L,1L))
  # expect_equal(res, matrix(c(1)))
  # 
  # # 2x2 matrix
  # coords <- c(1,2,3)
  # expect_error(
  #   res <- symmat(coords), 
  #   NA
  # )
  # expect_true(is(res, "matrix"))
  # expect_identical(dim(res), c(2L,2L))
  # expect_equal(res, 
  #              matrix(c(1,3,3,2),
  #                     nrow=2))
  # 
  # # 3x3 matrix
  # coords <- c(1,2,3,4,5,6)
  # expect_error(
  #   res <- symmat(coords), 
  #   NA
  # )
  # expect_true(is(res, "matrix"))
  # expect_identical(dim(res), c(3L,3L))
  # expect_equal(res, 
  #              matrix(c(1,4,5,4,2,6,5,6,3),
  #                     nrow=3))
  # 
  # # Vector incompatible
  # coords <- c(1,2,3,4,5)
  # expect_error(
  #   res <- symmat(coords), 
  #   "length of coords incompatible with a symetric matrix"
  # )
})

test_that("compute_sparse_chol produce expected output", {
  set.seed(123)
  X <- data.frame(cbind(runif(vecchia_approx$n_obs), 
                        rnorm(vecchia_approx$n_obs), 
                        rpois(vecchia_approx$n_obs, 5)))
  range_X <- process_covariates(
    X = X,
    vecchia_approx = vecchia_approx)

  range_beta <- matrix(rnorm(4 + PP$n_knots))
  
  expect_error(
    tmp <- compute_sparse_chol(range_beta = range_beta,
                        vecchia_approx = vecchia_approx,
                        range_X = range_X,
                        PP = PP,
                        matern_smoothness = 1.5,
                        compute_derivative = T),
    NA)
  
  expect_true(inherits(tmp, "array"))
  expect_identical(dim(tmp), c(11L, 12L, 1000L))
  expect_equal(mean(tmp), 0.0008364, tolerance = 1e-5)
})


## beta_prior_log_dens ##################################
test_that("beta_prior_log_dens produce expected output", {
  set.seed(123)
  beta1 <- matrix(rnorm(100), 100, ncol=1)
  beta3 <- matrix(rnorm(300), 100, ncol=3)
  expect_error(
    tmp <- beta_prior_log_dens(
      beta = beta3, 
      n_PP = 120, 
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3)
    ),
    "n_PP can't be greater than")
  expect_error(
    tmp <- beta_prior_log_dens(
      beta = beta3, 
      n_PP = 90, 
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3)
    ),
    "log_scale is supposed to be of length 2")
  
  expect_error(
    tmp <- beta_prior_log_dens(
      beta = beta3, 
      n_PP = 90, 
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3, -2)
     ),
    NA)
  expect_equal(tmp, -3012.276)
  
  expect_error(
    tmp <- beta_prior_log_dens(
      beta = beta1, 
      n_PP = 90, 
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3)
    ),
    NA)
  expect_equal(tmp, -1155.22585)
})

## beta_prior_log_dens_derivative ##################################
test_that("beta_prior_log_dens_derivative produce expected output", {
  set.seed(123)
  beta1 <- matrix(rnorm(100), 100, ncol=1)
  beta3 <- matrix(rnorm(300), 100, ncol=3)
  expect_error(
    tmp <- beta_prior_log_dens_derivative(
      beta = beta3, 
      n_PP = 120, 
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3, 2)
    ),
    "n_PP can't be greater than")
  expect_error(
    tmp <- beta_prior_log_dens_derivative(
      beta = beta3, 
      n_PP = 90, 
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3)
    ),
    "log_scale is supposed to be of length 2")
  
  expect_error(
    tmp <- beta_prior_log_dens_derivative(
      beta = beta3, 
      n_PP = 90, 
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3, -2)
    ),
    NA)
  expect_true(inherits(tmp, "matrix"))
  expect_identical(dim(tmp), dim(beta3))
  expect_equal(mean(tmp), 1.00837737)
  
  expect_error(
    tmp <- beta_prior_log_dens_derivative(
      beta = beta1, 
      n_PP = 90, 
      beta0_mean = -5,
      beta0_var = 2,
      log_scale = c(-3)
    ),
    NA)
  expect_true(inherits(tmp, "matrix"))
  expect_identical(dim(tmp), dim(beta1))
  expect_equal(mean(tmp), -2.9948913)
})

test_that("getting correct values using beta_prior_log_dens and beta_prior_log_dens_derivative", {
  set.seed(123)
  beta = matrix(rnorm(300), 100)
  res <-matrix(0, 100, 3)
  for(i in seq(100)) {
    for(j in seq(3)) {
      beta1 = beta 
      beta1[i,j] = beta1[i,j]+ .0001
      res[i,j] <- (
        beta_prior_log_dens(beta1, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = c(0, 2)) -
          beta_prior_log_dens(beta, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = c(0, 2))
      )*10000 / beta_prior_log_dens_derivative(beta, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = c(0, 2))[i,j]
   }
  }
  expect_equal(mean(res), 1.00004157)
  expect_equal(var(c(res)), 3.722309e-07, tolerance = 1e-7)
})

## X_PP_crossprod ##################################
test_that("X_PP_crossprod Simple crossprod if no PP, vecchia unused", {
  X = matrix(rnorm(100), 50)
  Y = matrix(rnorm(30*nrow(X)), nrow(X))
  expect_error(
    res1 <- X_PP_crossprod(X = X, 
                           PP = NULL, 
                           Y = Y, 
                           vecchia_approx = vecchia_approx),
    NA
  )
  
  expect_true(is(res1, "matrix"))
  expect_identical(dim(res1), c(2L, 30L))
  expect_identical(res1,crossprod(X, Y))
  
  # Simple corssprod, vecchia unused. 
  expect_error(
    res2 <- X_PP_crossprod(X = X, PP = NULL, Y = Y, vecchia_approx = NULL),
    NA
  )
  expect_identical(res1,res2)
})

test_that("X_PP_crossprod with PP", {
  set.seed(123)
  X = matrix(rnorm(100), 50)
  Y = matrix(rnorm(30*nrow(X)), nrow(X))
  expect_error(
    res2 <- X_PP_crossprod(X = X, PP = PP, Y = Y, vecchia_approx = vecchia_approx),
    NA
  )
  expect_true(is(res2, "matrix"))
  expect_identical(dim(res2), c(27L, 30L))
  
  # First line is crossprod
  expect_identical(res2[1,],crossprod(X, Y)[1,])
  expect_equal(colMeans(res2)[1:5], 
               c(-1.2353077, -0.2666516, -0.5238175,  0.6992207,  0.2442785),
               tolerance = 1e-5)
})

test_that("X_PP_crossprod with PP and permute obs", {
  set.seed(123)
  X = matrix(rnorm(100), 50)
  Y = matrix(rnorm(30*nrow(X)), nrow(X))
  expect_error(
    res3 <- X_PP_crossprod(X = X, PP = PP, Y = Y, vecchia_approx = vecchia_approx, permutate_PP_to_obs = TRUE),
    NA
  )
  expect_true(is(res3, "matrix"))
  expect_identical(dim(res3), c(27L, 30L))
  
  # First line is crossprod
  expect_identical(res3[1,],crossprod(X, Y)[1,])
  expect_equal(colMeans(res3)[1:5], 
               c(-1.2334375, -0.2906708, -0.6123371, 0.9166211,  0.2427284),
               tolerance = 1e-5)
})

## X_PP_mult_right ##################################
test_that("X_PP_mult_right with PP and permute obs", {
  set.seed(123)
  X = matrix(rnorm(100), 50)
  Y = matrix(rnorm(30*nrow(X)), nrow(X))
  expect_error(
    resmr <- X_PP_mult_right(X = X, PP = PP, Y = Y, vecchia_approx = vecchia_approx, permutate_PP_to_obs = FALSE),
    "Y should have"
  )
  # No PP
  expect_error(
    resmr <- X_PP_mult_right(X = X, 
                             PP = NULL, 
                             Y = matrix(rep(c(0.5, 1, 1.1, 0.2),2), nrow = 2), 
                             vecchia_approx = vecchia_approx, 
                             permutate_PP_to_obs = FALSE),
    NA
  )
  expect_true(is(resmr, "matrix"))
  expect_identical(dim(resmr), c(50L, 4L))
  expect_equal(mean(resmr), 0.1153678, tolerance = 1e-5)
  
  # No X
  expect_error(
    resmr <- X_PP_mult_right(X = NULL, PP = PP, Y = matrix(rnorm(250), nrow = 25), vecchia_approx = vecchia_approx, permutate_PP_to_obs = FALSE),
    NA
  )
  expect_true(is(resmr, "matrix"))
  expect_identical(dim(resmr), c(50L, 10L))
  expect_equal(mean(as.matrix(resmr)), -0.160772, tolerance = 1e-5)
  
  # X and PP
  expect_error(
    resmr <- X_PP_mult_right(X = X, PP = PP, Y =  matrix(rnorm(270), nrow = 27), vecchia_approx = vecchia_approx, permutate_PP_to_obs = FALSE),
    NA
  )
  expect_true(is(resmr, "matrix"))
  expect_identical(dim(resmr), c(50L, 10L))
  expect_equal(mean(as.matrix(resmr)), -0.03502011, tolerance = 1e-5)
})


