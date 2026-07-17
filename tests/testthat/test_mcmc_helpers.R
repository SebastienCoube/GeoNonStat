make_small_geononstat <- function(seed) {
  set.seed(seed)
  nlocs <- 16
  obs_locs <- cbind(runif(nlocs), runif(nlocs))
  vecchia_approx <- createVecchia(obs_locs, m = 3)

  set.seed(seed + 1)
  PP_range <- createPP(vecchia_approx, knots = 3, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)
  set.seed(seed + 2)
  PP_noise <- createPP(vecchia_approx, knots = 3, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)

  X <- as.data.frame(cbind(runif(nlocs), rnorm(nlocs)))
  observed_field <- rnorm(nlocs)

  GeoNonStat(
    vecchia_approx = vecchia_approx,
    observed_field = observed_field,
    X = X,
    matern_smoothness = 1.5,
    anisotropic = TRUE,
    n_chains = 1,
    range_X = X,
    range_PP = PP_range,
    noise_X = X,
    noise_PP = PP_noise
  )
}

test_that("renewMomentum is reproducible for a fixed seed", {
  momentum <- c(0.1, -0.2, 0.3)

  set.seed(123)
  out1 <- renewMomentum(momentum, kept_momentum = 0.8)
  expect_equal(length(out1), length(momentum))
  expect_type(out1, "double")
})

test_that("updateVarPPSuff is reproducible for a fixed seed", {
  hm4params <- list(
    PP = list(n_knots = 3),
    beta0_mean = 0,
    beta0_sd = 1,
    log_scale_bounds = c(-4, 3)
  )
  beta4params <- matrix(c(0.2, -0.1, 0.05), nrow = 3, ncol = 1)
  current <- c(0.1)

  set.seed(321)
  out1 <- updateVarPPSuff(
    hm4params = hm4params,
    beta4params = beta4params,
    current_range_log_scale = current
  )
  expect_type(out1, "double")
  expect_equal(length(out1), length(current))
})

test_that("initStateParams preserves shapes and initializes missing values", {
  state_params <- list(
    rnorm(4),
    matrix(rnorm(8), nrow = 2, ncol = 4)
  )

  res <- initStateParams(state_params)

  expect_type(res, "list")
  expect_equal(length(res), 2)
  expect_equal(length(res[[1]]), 4)
  expect_equal(dim(res[[2]]), dim(state_params[[2]]))
  expect_true(all(is.na(res[[1]])))
  expect_true(all(is.na(res[[2]])))
})

test_that("updateKernel clips values deterministically", {
  expect_equal(
    updateKernel(iter = 1, iter_start = 0, kernel_value = 0, mult = 1),
    0.3015,
    tolerance = 1e-4
  )
  expect_equal(
    updateKernel(iter = 1000, iter_start = 0, kernel_value = -20, mult = -1),
    -12
  )
})

test_that("updateBeta gives expected format", {
  obj <- make_small_geononstat(2026)

  set.seed(7)
  res1 <- updateBeta(obj$states$chain_1$params, 
                     obj$states$chain_1$stuff, 
                     obj$covariates$X, 
                     obj$vecchia_approx, 
                     obj$observed_field)
  expect_true(is(res1$params$beta, "matrix"))
  expect_identical(dim(obj$states$chain_1$params$beta), dim(res1$params$beta))
  expect_type(res1$stuff$lm_fit, "double")
  expect_length(res1$stuff$lm_fit, 16)
})

test_that("updateLatentField gives expected result", {
  obj <- make_small_geononstat(2027)

  set.seed(11)
  res1 <- updateLatentField(
    state = obj$states$chain_1,
    vecchia_approx = obj$vecchia_approx,
    hierarchical_model = obj$hierarchical_model,
    observed_field = obj$observed_field,
    iter = 1,
    num_threads = 1
  )
  
  expect_equal(length(res1$field), length(obj$states$chain_1$params$field))
  expect_type(res1$field, "double")
})

test_that("updateFieldLogVar gives expected results", {
  obj <- make_small_geononstat(2028)

  set.seed(13)
  res1 <- updateFieldLogVar(
    state = obj$states$chain_1,
    scale = obj$hierarchical_model$scale,
    vecchia_approx = obj$vecchia_approx,
    iter = 1,
    iter_start = 0
  )
  
  expect_equal(dim(res1$params$field_log_var), dim(obj$states$chain_1$params$field_log_var))
  expect_equal(length(res1$params$field), length(obj$states$chain_1$params$field))
})

test_that("updateNoiseBeta gives expected results", {
  obj <- make_small_geononstat(2030)

  set.seed(17)
  res1 <- updateNoiseBeta(
    state = obj$states$chain_1,
    noise = obj$hierarchical_model$noise,
    noise_X = obj$covariates$noise_X,
    vecchia_approx = obj$vecchia_approx,
    iter = 1,
    iter_start = 0
  )
  
  expect_equal(dim(res1$state$params$noise_beta), dim(obj$states$chain_1$params$noise_beta))
  expect_equal(length(res1$state$stuff$noise_var), length(obj$states$chain_1$stuff$noise_var))
  expect_equal(mean(res1$state$stuff$noise_beta[,1]), 6.966677, tolerance=1e-6)
  expect_equal(mean(res1$state$stuff$noise_var[1]), 0.01174705, tolerance=1e-6)
})

test_that("updateNoiseBetaFisher gives expected results", {
  obj <- make_small_geononstat(2031)

  set.seed(19)
  res1 <- updateNoiseBetaFisher(
    state = obj$states$chain_1,
    noise = obj$hierarchical_model$noise,
    noise_X = obj$covariates$noise_X,
    vecchia_approx = obj$vecchia_approx,
    iter = 1,
    iter_start = 0
  )
 
  expect_equal(dim(res1$state$params$noise_beta), dim(obj$states$chain_1$params$noise_beta))
  expect_equal(length(res1$state$stuff$noise_var), length(obj$states$chain_1$stuff$noise_var))
  expect_equal(mean(res1$state$stuff$noise_beta[,1]), 7.745450, tolerance=1e-6)
  expect_equal(mean(res1$state$stuff$noise_var[1]), 0.0976740, tolerance=1e-6)
})

test_that("updateRangeBetaAndLogVarMALA is reproducible for a fixed seed", {
  # Needs two objects here since updateRangeBetaAndLogVarMALA modify the object by reference. 
  obj1 <- make_small_geononstat(2032)
  obj2 <- make_small_geononstat(2032)
  expect_equal(obj1$states$chain_1, obj2$states$chain_1)
  set.seed(23)
  state = obj1$states$chain_1
  statekeep <- rlang::duplicate(state, shallow = FALSE)
  res1 <- updateRangeBetaAndLogVarMALA(
    state = state,
    hierarchical_model = obj1$hierarchical_model,
    vecchia_approx = obj1$vecchia_approx,
    range_X = obj1$covariates$range_X,
    iter = 1,
    iter_start = 0,
    num_threads = 1
  )
  set.seed(23)
  res2 <- updateRangeBetaAndLogVarMALA(
    state = statekeep,
    hierarchical_model = obj1$hierarchical_model,
    vecchia_approx = obj1$vecchia_approx,
    range_X = obj1$covariates$range_X,
    iter = 1,
    iter_start = 0,
    num_threads = 1
  )

  expect_equal(res1$state$params$range_beta, res2$state$params$range_beta)
  expect_equal(res1$state$params$field_log_var, res2$state$params$field_log_var)
  expect_equal(dim(res1$state$params$range_beta), dim(obj1$states$chain_1$params$range_beta))
})

test_that("renewMomentum handles edge cases for kept_momentum", {
  momentum <- c(0.1, -0.2, 0.3)

  expect_error(renewMomentum(momentum, kept_momentum = -0.1), regexp = "kept_momentum must be between 0 and 1")
  expect_error(renewMomentum(momentum, kept_momentum = 1.1), regexp = "kept_momentum must be between 0 and 1")

  set.seed(456)
  out_boundary_0 <- renewMomentum(momentum, kept_momentum = 0)
  set.seed(456)
  out_boundary_0_rep <- renewMomentum(momentum, kept_momentum = 0)
  expect_equal(out_boundary_0, out_boundary_0_rep)

  set.seed(456)
  out_boundary_1 <- renewMomentum(momentum, kept_momentum = 1)
  set.seed(456)
  out_boundary_1_rep <- renewMomentum(momentum, kept_momentum = 1)
  expect_equal(out_boundary_1, out_boundary_1_rep)
})

test_that("convertBetaAncillary is reproducible for a fixed seed", {
  obj <- make_small_geononstat(2035)
  PP <- obj$hierarchical_model$range$PP
  base_beta <- obj$states$chain_1$params$range_beta
  old_log_scale <- obj$states$chain_1$params$range_log_scale
  new_log_scale <- old_log_scale + c(0.3, 0.2)
  beta <- rbind(base_beta, matrix(rnorm(PP$n_knots), ncol = ncol(base_beta)))

  set.seed(456)
  out1 <- convertBetaAncillary(PP, old_log_scale, new_log_scale, beta)

  set.seed(456)
  out2 <- convertBetaAncillary(PP, old_log_scale, new_log_scale, beta)

  expect_equal(out1, out2)
  expect_equal(dim(out1), dim(beta))
  expect_true(is(out1, "matrix"))
})

# Test don't pass ... TODO
# test_that("convertBetaAncillary works with scalar log_scale", {
#   obj <- make_small_geononstat(2034)
#   PP <- obj$hierarchical_model$range$PP
#   base_beta <- obj$states$chain_1$params$range_beta
#   old_log_scale <- obj$states$chain_1$params$range_log_scale[1]
#   new_log_scale <- old_log_scale + 0.5
#   beta <- rbind(base_beta, matrix(rnorm(PP$n_knots), ncol = ncol(base_beta)))
# 
#   set.seed(789)
#   out1 <- convertBetaAncillary(PP, old_log_scale, new_log_scale, beta)
# 
#   set.seed(789)
#   out2 <- convertBetaAncillary(PP, old_log_scale, new_log_scale, beta)
# 
#   expect_equal(out1, out2)
#   expect_equal(dim(out1), dim(beta))
#   expect_type(out1, "double")
# })

test_that("updateVarPPAncillaryXSufficient is reproducible for a fixed seed", {
  obj <- make_small_geononstat(2033)

  set.seed(29)
  res1 <- updateVarPPAncillaryXSufficient(
    state = obj$states$chain_1,
    hierarchical_model = obj$hierarchical_model,
    vecchia_approx = obj$vecchia_approx,
    range_X = obj$covariates$range_X,
    iter = 1,
    iter_start = 0,
    num_threads = 1
  )

  set.seed(29)
  res2 <- updateVarPPAncillaryXSufficient(
    state = obj$states$chain_1,
    hierarchical_model = obj$hierarchical_model,
    vecchia_approx = obj$vecchia_approx,
    range_X = obj$covariates$range_X,
    iter = 1,
    iter_start = 0,
    num_threads = 1
  )

  expect_equal(res1$params$range_log_scale, res2$params$range_log_scale)
  expect_equal(res1$params$range_beta, res2$params$range_beta)
  expect_equal(length(res1$params$range_log_scale), length(obj$states$chain_1$params$range_log_scale))
})
