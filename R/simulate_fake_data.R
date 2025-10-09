

# ingredients needed to simulate observations from GeoNonStat model
GeoNonStatSimulator = function(vecchia_approx,
                               X = NULL,
                               # Response variable
                               matern_smoothness = 1.5,
                               #Matern smoothness
                               noise_X = NULL,
                               noise_PP = NULL,
                               range_X = NULL,
                               range_PP = NULL,
                               anisotropic = F) {
  covariates = list()
  covariates$X = process_covariates(
    X = X,
    vecchia_approx = vecchia_approx,
    PP = NULL,
    one_obs_per_locs = F
  )
  # fixed effects and PP for range
  covariates$range_X = process_covariates(
    X = range_X,
    vecchia_approx = vecchia_approx,
    PP = range_PP,
    one_obs_per_locs = T
  )
  # fixed effects and PP for noise
  covariates$noise_X = process_covariates(
    X = noise_X,
    vecchia_approx = vecchia_approx,
    PP = noise_PP,
    one_obs_per_locs = F
  )
  return(
    list(
      "covariates" =  covariates,
      "vecchia_approx" = vecchia_approx,
      "noise_PP" = noise_PP,
      "range_PP" = range_PP,
      "matern_smoothness" = matern_smoothness,
      "anisotropic" = anisotropic
    )
  )
}


# From a GeoNonStatSimulator object, creates regression coefficients at the right format
create_coeff_list = function(GNSSimulator) {
  covs <- GNSSimulator$covariates
  res = list(
    "X_coeff" =       matrix(
      0,
      nrow = ncol(covs$X$X)      ,
      dimnames = list(colnames(covs$X$X))
    ),
    "noise_X_coeff" = matrix(
      0,
      nrow = ncol(covs$noise_X$X),
      dimnames = list(colnames(covs$noise_X$X))
    ),
    "range_X_coeff" = matrix(
      0,
      nrow = ncol(covs$range_X$X_locs),
      dimnames = list(colnames(covs$range_X$X_locs), 
                      c("det", "aniso", "aniso")[seq(1 + 2 *GNSSimulator$anisotropic)]),
      ncol = 1 + 2 * GNSSimulator$anisotropic
    ),
    "noise_PP_coeff" = NULL,
    "field_log_var" = 0
  )
  if (!is.null(GNSSimulator$noise_PP))
    res$noise_PP_coeff = matrix(0,
                                nrow = GNSSimulator$noise_PP$n_knots,
                                dimnames = list(paste(
                                  "PP", seq(GNSSimulator$noise_PP$n_knots), sep = ""
                                )))
  if (!is.null(GNSSimulator$range_PP))
    res$range_PP_coeff = matrix(0,
                                nrow = GNSSimulator$range_PP$n_knots, ncol = 1 + 2*GNSSimulator$anisotropic,
                                dimnames = list(paste(
                                  "PP", seq(GNSSimulator$range_PP$n_knots), sep = ""
                                )))
  res
}


# From a GeoNonStatSimulator object and regression coefficients at the right format,
# simulates fake data
simulate = function(GNSSimulator, coeff_list, num_threads = 5) {
  log_range_field = X_PP_mult_right(
    X = GNSSimulator$covariates$range_X$X_locs,
    PP = GNSSimulator$range_PP,
    vecchia_approx = GNSSimulator$vecchia_approx,
    Y = rbind(coeff_list$range_X_coeff, coeff_list$range_PP_coeff),
    permutate_PP_to_obs = F
  )
  
  
  log_noise_field = X_PP_mult_right(
    X = GNSSimulator$covariates$noise_X$X,
    PP = GNSSimulator$noise_PP,
    vecchia_approx = GNSSimulator$vecchia_approx,
    Y = rbind(coeff_list$noise_X_coeff, coeff_list$noise_PP_coeff),
    permutate_PP_to_obs = T
  )
  
  sparse_chol = decompress_chol(
    vecchia_approx = GNSSimulator$vecchia_approx,
    compute_sparse_chol(
      range_beta = rbind(coeff_list$range_X_coeff, coeff_list$range_PP_coeff),
      vecchia_approx = GNSSimulator$vecchia_approx,
      range_X = GNSSimulator$covariates$range_X,
      PP = GNSSimulator$range_PP,
      matern_smoothness = GNSSimulator$matern_smoothness,
      compute_derivative = F,
      num_threads = num_threads
    )
  )
  
  latent_field = as.vector(exp(.5 * coeff_list$field_log_var) * Matrix::solve(sparse_chol, rnorm(nrow(sparse_chol))))
  noise = as.vector(exp(.5 * log_noise_field) * rnorm(length(log_noise_field)))
  fixed_effects = as.vector(GNSSimulator$covariates$X$X %*% coeff_list$X_coeff)
  
  observed_field =
    as.vector(t(GNSSimulator$vecchia_approx$locs_match_matrix) %*% latent_field) +
    noise +
    fixed_effects
  
  res = list(
    "log_range_field" = log_range_field,
    "log_noise_var_field" = log_noise_field,
    "sparse_chol" = sparse_chol,
    "noise" = noise,
    "fixed_effects" = fixed_effects,
    "latent_field" = latent_field,
    "observed_field" = observed_field
  )
}
