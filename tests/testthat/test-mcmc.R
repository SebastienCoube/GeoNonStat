test_that("runOneChainMcmc returns expected structure and is reproducible for fixed inputs", {
  set.seed(2026)
  nlocs <- 30
  obs_locs <- cbind(runif(nlocs), runif(nlocs))
  vecchia_approx <- createVecchia(obs_locs, m = 3)

  # small reproducible PP objects
  set.seed(11)
  PP_range <- createPP(vecchia_approx, knots = 4, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)
  set.seed(11)
  PP_noise <- createPP(vecchia_approx, knots = 4, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)

  X <- as.data.frame(cbind(runif(nlocs), rnorm(nlocs)))
  observed_field <- rnorm(nlocs)

  # Build GeoNonStat object with 1 chain and small PPs
  set.seed(7)
  g <- GeoNonStat(
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

  state <- g$states$chain_1

  # Run a short chain (no parallelism) and check structure
  set.seed(42)
  res1 <- runOneChainMcmc(
    covariates = g$covariates,
    observed_field = g$observed_field,
    hierarchical_model = g$hierarchical_model,
    vecchia_approx = g$vecchia_approx,
    state = state,
    n_iterations = 3,
    num_threads = 1,
    iter_start = 0
  )

  expect_type(res1, "list")
  expect_true(all(c("state", "params_records") %in% names(res1)))
  expect_length(res1$params_records, 3)

  # Re-run with same inputs and expect same numeric results
  set.seed(42)
  res2 <- runOneChainMcmc(
    covariates = g$covariates,
    observed_field = g$observed_field,
    hierarchical_model = g$hierarchical_model,
    vecchia_approx = g$vecchia_approx,
    state = state,
    n_iterations = 3,
    num_threads = 1,
    iter_start = 0
  )

  # Compare numeric parts of the returned states and records
  expect_equal(res1$params_records[[3]]$field, res2$params_records[[3]]$field)
  expect_equal(res1, res2)
})

test_that("runOneChainMcmc is reproducible after more iterations", {
  set.seed(2026)
  nlocs <- 24
  obs_locs <- cbind(runif(nlocs), runif(nlocs))
  vecchia_approx <- createVecchia(obs_locs, m = 3)

  set.seed(11)
  PP_range <- createPP(vecchia_approx, knots = 3, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)
  set.seed(11)
  PP_noise <- createPP(vecchia_approx, knots = 3, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)

  X <- as.data.frame(cbind(runif(nlocs), rnorm(nlocs)))
  observed_field <- rnorm(nlocs)

  set.seed(7)
  g1 <- GeoNonStat(
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
  set.seed(7)
  g2 <- GeoNonStat(
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

  set.seed(123)
  res1 <- runOneChainMcmc(
    covariates = g1$covariates,
    observed_field = g1$observed_field,
    hierarchical_model = g1$hierarchical_model,
    vecchia_approx = g1$vecchia_approx,
    state = g1$states$chain_1,
    n_iterations = 65,
    num_threads = 1,
    iter_start = 0
  )

  set.seed(123)
  res2 <- runOneChainMcmc(
    covariates = g2$covariates,
    observed_field = g2$observed_field,
    hierarchical_model = g2$hierarchical_model,
    vecchia_approx = g2$vecchia_approx,
    state = g2$states$chain_1,
    n_iterations = 65,
    num_threads = 1,
    iter_start = 0
  )

  expect_equal(res1, res2)
})

test_that("GeoNonStatMcmc is reproducible when run twice with same seed", {
  set.seed(2027)
  nlocs <- 30
  obs_locs <- cbind(runif(nlocs), runif(nlocs))
  vecchia_approx <- createVecchia(obs_locs, m = 3)

  set.seed(13)
  PP_range <- createPP(vecchia_approx, knots = 4, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)
  set.seed(13)
  PP_noise <- createPP(vecchia_approx, knots = 4, plot = FALSE, verbose = FALSE, reorder_knots = FALSE)

  X <- as.data.frame(cbind(runif(nlocs), rnorm(nlocs)))
  observed_field <- rnorm(nlocs)

  set.seed(9)
  g <- GeoNonStat(
    vecchia_approx = vecchia_approx,
    observed_field = observed_field,
    X = X,
    matern_smoothness = 1.5,
    anisotropic = TRUE,
    n_chains = 2,
    range_X = X,
    range_PP = PP_range,
    noise_X = X,
    noise_PP = PP_noise
  )

  # run twice with same seed; function warns for n_iterations < 60
  set.seed(100)
  expect_warning(
    out1 <- GeoNonStatMcmc(g, n_chains_in_parallel = 1, n_threads_per_chain = 1, n_iterations = 5),
    "The algorithm is implemented to update the spatial range parameters"
  )
  
  set.seed(9)
  g2 <- GeoNonStat(
    vecchia_approx = vecchia_approx,
    observed_field = observed_field,
    X = X,
    matern_smoothness = 1.5,
    anisotropic = TRUE,
    n_chains = 2,
    range_X = X,
    range_PP = PP_range,
    noise_X = X,
    noise_PP = PP_noise
  )
  # run twice with same seed; function warns for n_iterations < 60
  set.seed(100)
  out1 <- GeoNonStatMcmc(g, n_chains_in_parallel = 3, 
                         n_threads_per_chain = 1, 
                         n_iterations = 65)
  
  set.seed(100)
  out2 <- GeoNonStatMcmc(g2, n_chains_in_parallel = 3, 
                         n_threads_per_chain = 1, 
                         n_iterations = 65)
  
  # Compare numeric records across runs
  expect_equal(lapply(out1$records, function(x) lapply(x, function(p) p$noise_beta)), 
               lapply(out2$records, function(x) lapply(x, function(p) p$noise_beta)))
  
})
