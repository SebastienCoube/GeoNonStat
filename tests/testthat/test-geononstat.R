# Tests for functions in GeoNonStat.R
test_that("processTransitionKernels returns expected structure and values", {
  res <- processTransitionKernels(init = -3)
  expect_type(res, "list")
  expect_true(all(c("range_log_scale_sufficient", "noise_log_scale") %in% names(res)))
  expect_equal(res$range_log_scale_sufficient, -3)
  expect_equal(res$noise_log_scale, -3)
})

test_that("generateLocationPartitions returns matrix with correct rows", {
  set.seed(1)
  locs <- matrix(runif(40), ncol = 2)
  parts <- generateLocationPartitions(locs, n = nrow(locs))
  expect_true(is.matrix(parts) || is.data.frame(parts))
  expect_equal(nrow(parts), nrow(locs))
})

test_that("processPPPrior handles NULL and invalid inputs correctly", {
  # both NULL -> returns NULL
  expect_null(processPPPrior(PP = NULL, log_scale_bounds = NULL, parameter_name = "p"))

  # PP NULL but bounds provided -> error
  expect_error(processPPPrior(PP = NULL, log_scale_bounds = c(-3, 2), parameter_name = "p"))

  # PP not of class PP -> error
  fake_pp <- list()
  class(fake_pp) <- "notPP"
  expect_error(processPPPrior(PP = fake_pp, log_scale_bounds = NULL, parameter_name = "p"))

  # PP of class PP and NULL bounds -> returns default bounds and emits message
  good_pp <- list()
  class(good_pp) <- "PP"
  expect_message(out <- processPPPrior(PP = good_pp, log_scale_bounds = NULL, parameter_name = "p"), "automatically set")
  expect_equal(out, c(-5, 6))

  # provided bounds numeric of length 2 are sorted and returned
  class(good_pp) <- "PP"
  expect_equal(processPPPrior(PP = good_pp, log_scale_bounds = c(5, -1), parameter_name = "p"), c(-1, 5))
})

test_that("processHierarchicalModel rejects invalid matern_smoothness early", {
  # only check early validation; other args can be dummies because check occurs first
  expect_error(
    processHierarchicalModel(
      vecchia_approx = list(locs = matrix(runif(4), ncol = 2)),
      noise_PP = NULL,
      noise_log_scale_bounds = NULL,
      range_PP = NULL,
      range_log_scale_bounds = NULL,
      observed_locs = matrix(runif(4), ncol = 2),
      matern_smoothness = 2.0,
      observed_field = rnorm(2),
      covariates = list(X = list(X = matrix(1, 1))),
      anisotropic = FALSE
    ), "only matern_smoothness = 1.5 or matern_smoothness = 0.5"
  )
})

test_that("processHierarchicalModel is deterministic with the same seed and inputs", {
  set.seed(123)
  obs_locs <- cbind(runif(20), runif(20))
  vecchia_approx <- createVecchia(obs_locs, m = 3)
  X <- as.data.frame(cbind(runif(nrow(obs_locs)), rnorm(nrow(obs_locs))))
  range_X <- processCovariates(X, vecchia_approx)
  noise_X <- processCovariates(X, vecchia_approx)
  observed_field <- rnorm(nrow(obs_locs))

  set.seed(123)
  hm1 <- processHierarchicalModel(
    vecchia_approx = vecchia_approx,
    noise_PP = NULL,
    noise_log_scale_bounds = NULL,
    range_PP = NULL,
    range_log_scale_bounds = NULL,
    observed_locs = obs_locs,
    matern_smoothness = 1.5,
    observed_field = observed_field,
    covariates = list(X = range_X, range_X = range_X, noise_X = noise_X),
    anisotropic = TRUE
  )

  set.seed(123)
  hm2 <- processHierarchicalModel(
    vecchia_approx = vecchia_approx,
    noise_PP = NULL,
    noise_log_scale_bounds = NULL,
    range_PP = NULL,
    range_log_scale_bounds = NULL,
    observed_locs = obs_locs,
    matern_smoothness = 1.5,
    observed_field = observed_field,
    covariates = list(X = range_X, range_X = range_X, noise_X = noise_X),
    anisotropic = TRUE
  )

  expect_equal(hm1$noise, hm2$noise)
  expect_equal(hm1$scale, hm2$scale)
  expect_equal(hm1$range, hm2$range)
  expect_equal(hm1$anisotropic, hm2$anisotropic)
  expect_equal(hm1$matern_smoothness, hm2$matern_smoothness)
  expect_equal(coef(hm1$naive_ols), coef(hm2$naive_ols))
  expect_equal(hm1$naive_ols$residuals, hm2$naive_ols$residuals)
})

test_that("processHierarchicalModel is deterministic with PP for range and noise", {
  set.seed(1234)
  obs_locs <- cbind(runif(50), runif(50))
  vecchia_approx <- createVecchia(obs_locs, m = 5)

  # build small PP objects for range and noise with fixed seed
  set.seed(42)
  PP_range <- createPP(vecchia_approx, knots = 4, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)
  set.seed(42)
  PP_noise <- createPP(vecchia_approx, knots = 4, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)

  X <- as.data.frame(cbind(runif(nrow(obs_locs)), rnorm(nrow(obs_locs))))
  range_X <- processCovariates(X, vecchia_approx)
  noise_X <- processCovariates(X, vecchia_approx)
  observed_field <- rnorm(nrow(obs_locs))

  set.seed(999)
  hm1 <- processHierarchicalModel(
    vecchia_approx = vecchia_approx,
    noise_PP = PP_noise,
    noise_log_scale_bounds = NULL,
    range_PP = PP_range,
    range_log_scale_bounds = NULL,
    observed_locs = obs_locs,
    matern_smoothness = 1.5,
    observed_field = observed_field,
    covariates = list(X = range_X, range_X = range_X, noise_X = noise_X),
    anisotropic = TRUE
  )

  set.seed(999)
  hm2 <- processHierarchicalModel(
    vecchia_approx = vecchia_approx,
    noise_PP = PP_noise,
    noise_log_scale_bounds = NULL,
    range_PP = PP_range,
    range_log_scale_bounds = NULL,
    observed_locs = obs_locs,
    matern_smoothness = 1.5,
    observed_field = observed_field,
    covariates = list(X = range_X, range_X = range_X, noise_X = noise_X),
    anisotropic = TRUE
  )

  expect_equal(hm1$noise$PP$n_knots, hm2$noise$PP$n_knots)
  expect_equal(hm1$range$PP$n_knots, hm2$range$PP$n_knots)
  expect_equal(hm1$noise, hm2$noise)
  expect_equal(hm1$scale, hm2$scale)
  expect_equal(hm1$range, hm2$range)
  expect_equal(coef(hm1$naive_ols), coef(hm2$naive_ols))
  expect_equal(hm1$naive_ols$residuals, hm2$naive_ols$residuals)
})

test_that("GeoNonStat object reproducible with anisotropic, range_PP and noise_PP", {
  set.seed(2026)
  nlocs <- 50
  obs_locs <- cbind(runif(nlocs), runif(nlocs))
  vecchia_approx <- createVecchia(obs_locs, m = 4)

  # Create small PP objects reproducibly
  set.seed(11)
  PP_range <- createPP(vecchia_approx, knots = 6, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)
  set.seed(11)
  PP_noise <- createPP(vecchia_approx, knots = 6, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)

  X <- as.data.frame(cbind(obs_locs, rnorm(nlocs)))
  observed_field <- rnorm(nlocs)

  set.seed(7)
  g1 <- GeoNonStat(
    vecchia_approx = vecchia_approx,
    data_list = list(
      observed_field = observed_field,
      X = X,
      range_X = X,
      noise_X = X),
    matern_smoothness = 1.5,
    anisotropic = TRUE,
    n_chains = 1,
    range_PP = PP_range,
    noise_PP = PP_noise
  )

  set.seed(7)
  g2 <- GeoNonStat(
    vecchia_approx = vecchia_approx,
    data_list = list(
      observed_field = observed_field,
      X = X,
      range_X = X,
      noise_X = X),
    matern_smoothness = 1.5,
    anisotropic = TRUE,
    n_chains = 1,
    range_PP = PP_range,
    noise_PP = PP_noise
  )

  # Compare covariates design matrices
  expect_equal(as.matrix(g1$covariates$X$X), as.matrix(g2$covariates$X$X))

  # Compare PP metadata
  expect_equal(g1$hierarchical_model$range$PP$n_knots, g2$hierarchical_model$range$PP$n_knots)
  expect_equal(g1$hierarchical_model$noise$PP$n_knots, g2$hierarchical_model$noise$PP$n_knots)

  # Compare hierarchical model numeric parts and OLS results
  expect_equal(coef(g1$hierarchical_model$naive_ols), coef(g2$hierarchical_model$naive_ols))
  expect_equal(g1$hierarchical_model$naive_ols$residuals, g2$hierarchical_model$naive_ols$residuals)

  # Compare a parameter from the MCMC state
  expect_equal(g1$states$chain_1$params$range_beta, g2$states$chain_1$params$range_beta)
})




make_setup <- function(nlocs = 100, nobs = 200, m = 10) {
  set.seed(123)
  unique_locs <- cbind(runif(nlocs), runif(nlocs))
  # chaque localisation apparait au moins une fois, avec des doublons
  idx <- c(seq_len(nlocs), sample(nlocs, nobs - nlocs, replace = TRUE))
  observed_locs <- unique_locs[idx, ]
  X <- as.data.frame(cbind(runif(nobs), rnorm(nobs), rpois(nobs, 5)))
  vecchia_approx <- createVecchia(observed_locs, m)
  PP <- createPP(vecchia_approx, plot = FALSE)
  list(
    nobs = nobs,
    observed_locs = observed_locs,
    X = X,
    vecchia_approx = vecchia_approx,
    PP = PP
  )
}

s <- make_setup()
nobs <- s$nobs
vecchia_approx <- s$vecchia_approx
PP <- s$PP
X <- s$X
n_unique_locs <- length(vecchia_approx$hctam_scol_1)

test_that("processCovariates works with no PP and no X (intercept only)", {
  res <- processCovariates(
    X = NULL, vecchia_approx, PP = NULL, one_obs_per_locs = FALSE
  )
  
  expect_type(res, "list")
  expect_null(res$arg)
  expect_true(methods::is(res$X, "sparseMatrix"))
  expect_equal(dim(res$X), c(nobs, 1))
  expect_equal(colnames(res$X), "(Intercept)")
  expect_equal(res$n_regressors, 1)
  expect_equal(res$which_locs, 1)
  expect_equal(dim(res$X_locs), c(n_unique_locs, 1))
  expect_equal(colnames(res$X_locs), "(Intercept)")
  expect_equal(dim(res$crossprod_X), c(1, 1))
  expect_equal(as.numeric(as.matrix(res$crossprod_X)), nobs)
})

test_that("processCovariates works with no PP and an X", {
  res <- processCovariates(
    X = X, vecchia_approx, PP = NULL, one_obs_per_locs = FALSE
  )
  
  expect_identical(res$arg, X)
  expect_true(methods::is(res$X, "sparseMatrix"))
  expect_equal(dim(res$X), c(nobs, ncol(X) + 1))
  expect_equal(colnames(res$X)[1], "(Intercept)")
  expect_equal(res$n_regressors, ncol(X) + 1)
  # random covariates vary within a location:
  # only the intercept is constant within locations
  expect_equal(res$which_locs, 1)
  expect_equal(dim(res$X_locs), c(n_unique_locs, 1))
  expect_equal(dim(res$crossprod_X), rep(ncol(X) + 1, 2))
})

test_that("processCovariates works with a PP and an X (crossprod_X includes PP basis)", {
  res <- processCovariates(
    X = X, vecchia_approx, PP = PP, one_obs_per_locs = FALSE
  )
  
  expect_identical(res$arg, X)
  expect_equal(res$n_regressors, ncol(X) + 1)
  expect_equal(dim(res$crossprod_X), rep(ncol(X) + 1 + PP$n_knots, 2))
  # X_locs does not contain the PP columns
  expect_equal(ncol(res$X_locs), length(res$which_locs))
})

test_that("processCovariates works with a PP and no X", {
  res <- processCovariates(
    X = NULL, vecchia_approx, PP = PP, one_obs_per_locs = FALSE
  )
  
  expect_null(res$arg)
  expect_equal(res$n_regressors, 1)
  expect_equal(dim(res$crossprod_X), rep(1 + PP$n_knots, 2))
})

test_that("processCovariates works with one_obs_per_locs = TRUE and covariates constant within locations", {
  locs <- vecchia_approx$observed_locs
  X_loc <- as.data.frame(cbind(locs, locs[, 1]^2 + locs[, 2]^2))
  
  res <- processCovariates(
    X = X_loc, vecchia_approx, PP = NULL, one_obs_per_locs = TRUE
  )
  
  expect_identical(res$arg, X_loc)
  expect_true(is(res$X, "Matrix")) # removed when one_obs_per_locs = TRUE
  expect_equal(res$which_locs, seq_len(ncol(X_loc) + 1))
  expect_equal(dim(res$X_locs), c(n_unique_locs, ncol(X_loc) + 1))
  expect_equal(colnames(res$X_locs)[1], "(Intercept)")
  expect_equal(dim(res$crossprod_X), rep(ncol(X_loc) + 1, 2))
})

test_that("processCovariates works with one_obs_per_locs = TRUE and no X (intercept only)", {
  res <- processCovariates(
    X = NULL, vecchia_approx, PP = NULL, one_obs_per_locs = TRUE
  )
  
  expect_true("list" %in% class(res))
  expect_true(is(res$X_locs, "Matrix"))
  expect_identical(dim(res$X_locs), c(100L, 1L))
})

# Error cases ------------------------------------------------------------------

test_that("processCovariates errors if X varies within a location and one_obs_per_locs = TRUE", {
  expect_error(
    processCovariates(
      X = X, vecchia_approx, PP = NULL, one_obs_per_locs = TRUE
    ),
    "cannot vary within one spatial location"
  )
})

test_that("processCovariates errors if X has the wrong number of rows", {
  X_bad <- X[-1, ]
  expect_error(
    processCovariates(
      X = X_bad, vecchia_approx, PP = NULL, one_obs_per_locs = FALSE
    ),
    "should have the same number of rows"
  )
  expect_error(
    processCovariates(
      X = X_bad, vecchia_approx, PP = NULL, one_obs_per_locs = TRUE
    ),
    "should have the same number of rows"
  )
})

test_that("processCovariates errors if X is not a data.frame", {
  X_matrix <- as.matrix(X)
  expect_error(
    processCovariates(
      X = X_matrix, vecchia_approx, PP = NULL, one_obs_per_locs = FALSE
    ),
    "should be a data.frame or NULL"
  )
  expect_error(
    processCovariates(
      X = X_matrix, vecchia_approx, PP = NULL, one_obs_per_locs = TRUE
    ),
    "should be a data.frame or NULL"
  )
})

test_that("processCovariates errors if covariates are not independent", {
  locs <- vecchia_approx$observed_locs
  X_collinear <- as.data.frame(
    cbind(locs, locs, locs[, 1]^2 + locs[, 2]^2)
  )
  expect_error(
    processCovariates(
      X = X_collinear, vecchia_approx, PP = NULL, one_obs_per_locs = TRUE
    ),
    "does not induce an independent collection of covariates"
  )
  expect_error(
    processCovariates(
      X = X_collinear, vecchia_approx, PP = NULL, one_obs_per_locs = FALSE
    ),
    "does not induce an independent collection of covariates"
  )
})

test_that("processCovariates errors if X has a constant column (collinear with the intercept)", {
  X_const <- data.frame(a = rep(1, nobs))
  expect_error(
    processCovariates(
      X = X_const, vecchia_approx, PP = NULL, one_obs_per_locs = FALSE
    ),
    "does not induce an independent collection of covariates"
  )
})


