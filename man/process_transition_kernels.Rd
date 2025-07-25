% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeoNonStat.R
\name{process_transition_kernels}
\alias{process_transition_kernels}
\title{Initialize transition kernels}
\usage{
process_transition_kernels(init = -4, hm)
}
\arguments{
\item{init}{a numeric value}
}
\value{
a list
}
\description{
Initialize transition kernels
}
\examples{
tk <- process_transition_kernels()
nobs = 10000
observed_locs = cbind(runif(nobs), runif(nobs))
observed_field = rnorm(nobs)
vecchia_approx = createVecchia(observed_locs)
PP = createPP(vecchia_approx)
X = as.data.frame(cbind(rnorm(nobs), runif(nobs)))
covariates = list()
covariates$X = process_covariates(X, vecchia_approx = vecchia_approx)
covariates$X_range = process_covariates(X, vecchia_approx = vecchia_approx)
covariates$X_noise = process_covariates(X, vecchia_approx = vecchia_approx)
covariates$X_scale = process_covariates(X, vecchia_approx = vecchia_approx)

hm = process_hierarchical_model(
  noise_PP = PP, noise_log_scale_bounds = NULL,
  scale_PP = PP, scale_log_scale_bounds = NULL,
  range_PP = PP, range_log_scale_bounds = NULL,
  observed_locs,
  matern_smoothness = 1.5,
  observed_field = observed_field,
  covariates = covariates,
  anisotropic = T) 
process_transition_kernels(hm = hm)
}
