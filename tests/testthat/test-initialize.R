test_that("process_covariates works", {
  expect_equal(2 * 2, 4)
})


test_that("initialize works", {
  set.seed(100)
  size <- 2000
  locs = cbind(runif(size), runif(size))
  PP = PP(
    observed_locs = locs[seq(size),], # spatial sites
    matern_range = .1,
    knots = 50, # number of knots
    m = 15 # number of NNGP parents
  )
  range_beta = rbind(c(-4, 0, 0), # intercept
                     matrix(.5*rnorm(150), 50)
                     ) %*% diag(c(1, 2,.5))
  NNarray_aniso = GpGp::find_ordered_nn(locs[seq(size),], 10)
  
  # getting the coefficients
  chol_precision = compute_sparse_chol(
      range_beta = range_beta, # range parameters
      NNarray = NNarray_aniso, # Vecchia approx DAG
      locs = locs[seq(size),], # spatial coordinates
      range_X = matrix(1, size), # covariates for the range (just an intercept)
      PP = PP, use_PP = T, # predictive process
      nu = 1.5, # smoothness
      anisotropic = T # anisotropy
    )[[1]]
  # putting coefficients in precision Cholesly
  chol_precision = Matrix::sparseMatrix(
    x = chol_precision[!is.na(NNarray_aniso)], 
    i = row(NNarray_aniso)[!is.na(NNarray_aniso)], 
    j = (NNarray_aniso)[!is.na(NNarray_aniso)], 
    triangular = T
  )
  # sampling the anisotropic process
  seed_vector = rnorm(size)
  aniso_latent_field = as.vector(Matrix::solve(chol_precision, seed_vector)) 
  aniso_observed_field = aniso_latent_field + .8*rnorm(size)
  #plot_pointillist_painting(locs, aniso_latent_field)
  expect_message(
    MCMC_NNGP <- initialize_nngp(
      observed_locs = locs[seq(size),], 
      observed_field = aniso_observed_field, 
      nu = 1.5, 
      n_chains = 5,
      range_PP = TRUE, 
      PP = PP, # use PP for range
      anisotropic = TRUE # Covariance will be anisotropic
    ),
    "range_log_scale_prior was automatically set to an"
  )
  
  expect_true(is(MCMC_NNGP, "list"))
  expect_named(MCMC_NNGP, c('data', 'hierarchical_model', 'vecchia_approx', 'states', 'records', 't_begin', 'seed', 'iterations'))
  
  # data
  expect_named(MCMC_NNGP$data, c('locs', 'observed_field', 'observed_locs', 'covariates'))
  expect_identical(MCMC_NNGP$data$observed_locs, locs[seq(size),])
  expect_identical(MCMC_NNGP$data$observed_field, aniso_observed_field)
  expect_identical(mean(MCMC_NNGP$data$locs), mean(locs[seq(size),]))
  
  expect_named(MCMC_NNGP$data$covariates,c('X', 'range_X', 'scale_X', 'noise_X'))
  
  
})

