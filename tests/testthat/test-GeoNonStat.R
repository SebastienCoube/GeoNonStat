set.seed(100)
size <- 2000
observed_locs = cbind(runif(size), runif(size))

test_that("createVecchia works", {
  expect_error(
    res <- createVecchia(observed_locs, m=10, ncores=1),
    NA
  )
  expect_true(is(res, "list"))
  expect_named(res, c('n_locs', 'n_obs', 'observed_locs', 'locs', 't_locs', 
                      'locs_match', 'locs_match_matrix', 
                      'hctam_scol', 'hctam_scol_1', 
                      'obs_per_loc', 
                      'NNarray', 'NNarray_non_NA', 
                      'sparse_chol_i', 'sparse_chol_p', 'sparse_chol_x_reorder', 
                      'locs_partition'))
  expect_identical(res$n_locs, 2000L)
  expect_identical(res$n_obs, 2000L)
  expect_identical(res$observed_locs, observed_locs)
  expect_identical(dim(res$locs), c(2000L, 2L))
  expect_identical(colMeans(res$locs), colMeans(observed_locs))
  expect_identical(res$t_locs, t(res$locs))
  
  expect_identical(length(res$locs_match), 2000L)
  expect_true(all(res$locs_match %in% seq_len(res$n_locs)))
  
  expect_true("dgCMatrix" %in% class(res$locs_match_matrix))
  expect_identical(dim(res$locs_match_matrix), c(2000L, 2000L))
  
  expect_true("list" %in% class(res$hctam_scol))
  expect_true(all(res$hctam_scol %in% seq_len(res$n_locs)))
  expect_true(all(res$hctam_scol_1 %in% seq_len(res$n_locs)))
  
  expect_true(all(unlist(res$hctam_scol) == res$hctam_scol_1))
  expect_true(all(res$obs_per_loc == 1))
  
  expect_identical(dim(res$NNarray), c(11L, 2000L))
  expect_identical(length(res$sparse_chol_i), 21945L)
  expect_identical(length(res$sparse_chol_x_reorder), 21945L)
  expect_identical(dim(res$locs_partition), c(2000L, 2L))
})


test_that("process_covariates work as expected", {
  nlocs = 50L
  nobs = 100L
  unique_locs = cbind(runif(nlocs), runif(nlocs))
  observed_locs = rbind(unique_locs, unique_locs)
  X = as.data.frame(cbind(runif(nobs), rnorm(nobs), rpois(nobs, 5)))
  vecchia_approx = createVecchia(observed_locs, 12, ncores=1)
  PP = createPP(vecchia_approx, plot=FALSE)
  
  # no PP, an X
  expect_error(
    res <- process_covariates(X = X, 
                             vecchia_approx = vecchia_approx, 
                             PP = NULL, one_obs_per_locs = FALSE),
    NA
  )
  expect_true(is(res, "list"))
  expect_named(res, c('arg', 'X', 'chol_crossprod_X', 'n_regressors', 
                      'which_locs', 'X_locs', 'crossprod_X_locs', 
                      'chol_crossprod_X_locs'))
  
  expect_identical(res$arg, X)
  expect_true(is(res$X, "matrix"))
  expect_true(is(res$chol_crossprod_X, "matrix"))
  expect_identical(res$n_regressors, 4L)
  expect_identical(res$which_locs, 1L)
  expect_true(is(res$X_locs, "matrix"))
  expect_true(is(res$crossprod_X_locs, "matrix"))
  expect_true(is(res$chol_crossprod_X_locs, "matrix"))
  
  expect_identical(dim(res$X), c(nobs, 4L))
  expect_identical(dim(res$chol_crossprod_X), c(4L, 4L))
  expect_identical(dim(res$X_locs), c(nlocs, 1L))
  expect_identical(dim(res$crossprod_X_locs), c(1L, 1L))
  expect_identical(dim(res$chol_crossprod_X_locs), c(1L, 1L))
  
  # TODO tester les valeurs en sortie.
  # PP, no X
  expect_error(
    res <- process_covariates(X = NULL, 
                              vecchia_approx, 
                              PP = PP, 
                              one_obs_per_locs = FALSE),
    NA
  )
  
  # no PP, no X
  expect_error(
    res <- process_covariates(X = NULL, 
                             vecchia_approx, 
                             PP = NULL, 
                             one_obs_per_locs = FALSE),
    NA
  )
  
  # PP and X
  expect_error(
    res <- process_covariates(X = X, 
                              vecchia_approx, 
                              PP = PP, 
                              one_obs_per_locs = FALSE),
    NA
  )

  # one obs of x per loc (no duplicates)
  X = as.data.frame(cbind(observed_locs, observed_locs[,1]^2+ observed_locs[,2]^2))
  expect_error(
  res <- process_covariates(X = X, 
                            vecchia_approx, 
                            PP = NULL, 
                            one_obs_per_locs = TRUE),
  NA)

})

test_that("process_PP_priors works as expected whithout PP", {
  # Whithout PP
  expect_error(
    PPPP <- process_PP_prior(NULL, NULL, "example"),
    NA
  )
  expect_null(PPPP)
  
  expect_error(
    PPPP <- process_PP_prior(NULL, 1, "example"),
    "No PP object was provided for the example parameters, but log - marginal variance bounds were provided"
  )
})

test_that("process_PP_priors works as expected whith PP", {
  # With PP
  suppressMessages({
    myVecchia <- createVecchia(observed_locs, m=10, ncores = 1)
    myPP <- createPP(myVecchia) 
  })
  
  expect_message(
    PPPP <- process_PP_prior(myPP, NULL, "example"),
    "automatically set to \\(-6, 3\\)"
  )
  expect_identical(PPPP, c(-6, 3))
    
  expect_error(
    PPPP <- process_PP_prior(myPP, "toto", "example"),
    "must be a numeric vector of length 2"
  )
  
  expect_error(
    PPPP <- process_PP_prior(myPP, c(1,2,3), "example"),
    "must be a numeric vector of length 2"
  )
  
  expect_error(
    PPPP <- process_PP_prior(myPP, c(1,2), "tatato"),
    NA
  )
  expect_identical(PPPP, c(1,2))
  
})


test_that("process_transition_kernels return expected result", {
  init <- -5
  expect_error(
    res <- process_transition_kernels(init = init),
    NA
  )
  expect_true(is.list(res))
  expect_identical(res, 
               list( range_log_scale_sufficient = init,
                  range_log_scale_ancillary =  init,
                  range_beta_sufficient = c(init, init, init, init),
                  range_beta_ancillary  = c(init, init, init, init),
                  scale_beta_sufficient_mala = c(init, init),
                  scale_beta_ancillary_mala  = c(init, init),
                  scale_log_scale_sufficient = init,
                  scale_log_scale_ancillary =  init,
                  noise_beta_mala = c(init, init),
                  noise_log_scale = init))
})

set.seed(100)
suppressMessages({
  VA <- createVecchia(observed_locs, m=10, ncores=1)
  PP = createPP(VA, plot=FALSE)
})

test_that("initialize class GeoNonStat works", {
  set.seed(123)
  obs_field <- rnorm(nrow(observed_locs))
  expect_message(
    myobj <- GeoNonStat(
      vecchia_approx = VA,
      observed_field = obs_field
    ),
    "Setup done,"
  )
  
  expect_s3_class(myobj, "GeoNonStat")
  expect_named(myobj, c('covariates', 'observed_field', 
                        'hierarchical_model', 'vecchia_approx', 
                        'states', 'records', 'seed', 'checkpoints'))
  expect_true(is(myobj$covariates, "list"))
  expect_true(is(myobj$observed_field, "numeric"))
  expect_true(is(myobj$hierarchical_model, "list"))
  expect_true(is(myobj$vecchia_approx, "list"))
  expect_true(is(myobj$states, "list"))
  expect_true(is(myobj$records, "list"))
  expect_true(is(myobj$seed, "numeric"))
  expect_true(is(myobj$checkpoints, "matrix"))
  
  # covariates
  expect_named(myobj$covariates,c('X', 'range_X', 'scale_X', 'noise_X'))
  
  # observed_filed
  expect_identical(myobj$observed_field, obs_field)
  
  # hierarchical_model
  expect_named(myobj$hierarchical_model,c('noise_PP', 'noise_log_scale_bounds', 
                                          'scale_PP', 'scale_log_scale_bounds', 
                                          'range_PP', 'range_log_scale_bounds', 
                                          'anisotropic', 'matern_smoothness', 
                                          'beta_priors', 
                                          'range_beta0_mean', 'range_beta0_var', 'noise_beta0_mean', 'noise_beta0_var', 'scale_beta0_mean', 'scale_beta0_var', 
                                          'naive_ols'))
  
  # data
  expect_identical(myobj$vecchia_approx, VA)
  # TODO finir les tests
})
