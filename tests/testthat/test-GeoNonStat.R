test_that("process_vecchia GeoNonStat works", {
  set.seed(100)
  size <- 2000
  observed_locs = cbind(runif(size), runif(size))
  expect_message(
  res <- process_vecchia(observed_locs, m=10),
    "building DAGs and indices for Vecchia approximation"
  )
  expect_true("list" %in% class(res))
  expect_named(res, c('n_locs', 'n_obs', 'locs', 'locs_match', 
                      'locs_match_matrix', 'hctam_scol', 
                      'hctam_scol_1', 'obs_per_loc', 
                      'NNarray', 'NNarray_non_NA', 
                      'sparse_chol_column_idx', 'sparse_chol_row_idx', 
                      'locs_partition'))
  expect_identical(res$n_locs, 2000L)
  expect_identical(res$n_obs, 2000L)
  expect_true(is.numeric(res$locs))
  expect_identical(dim(res$locs), c(2000L, 2L))
  expect_equal(mean(res$locs), 0.5035152, tolerance = 1e-5)
  expect_identical(length(res$locs_match), 2000L)
  expect_true(all(res$locs_match %in% seq_len(res$n_locs)))
  
  expect_true("Matrix" %in% class(res$locs_match_matrix))
  expect_identical(dim(res$locs_match_matrix), c(2000L, 2000L))
  
  expect_true("list" %in% class(res$hctam_scol))
  expect_true(all(res$hctam_scol %in% seq_len(res$n_locs)))
  expect_true(all(res$hctam_scol_1 %in% seq_len(res$n_locs)))
  
  expect_true(all(unlist(res$hctam_scol) == res$hctam_scol_1))
  expect_true(all(res$obs_per_loc == 1))
  
  expect_identical(dim(res$NNarray), c(11L, 2000L))
  expect_identical(length(res$sparse_chol_column_idx), 21945L)
  expect_identical(length(res$sparse_chol_row_idx), 21945L)
  expect_identical(dim(res$locs_partition), c(2000L, 1L))
})


test_that("initialize class GeoNonStat works", {
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
    myobj <- GeoNonStat(
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
  
  expect_true(is(myobj, "GeoNonStat"))
  expect_named(myobj, c('data', 'hierarchical_model', 'vecchia_approx', 'states', 'records', 'seed', 'iterations'))
  
  # data
  expect_named(myobj$data, c('locs', 'observed_field', 'observed_locs', 'covariates'))
  expect_identical(myobj$data$observed_locs, locs[seq(size),])
  expect_identical(myobj$data$observed_field, aniso_observed_field)
  expect_identical(mean(myobj$data$locs), mean(locs[seq(size),]))
  
  expect_named(myobj$data$covariates,c('X', 'range_X', 'scale_X', 'noise_X'))
})
