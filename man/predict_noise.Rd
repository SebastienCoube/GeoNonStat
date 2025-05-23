% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict_noise}
\alias{predict_noise}
\title{Predicts the noise variance at unobserved locations}
\usage{
predict_noise(
  mcmc_nngp_list,
  X_noise_pred = NULL,
  burn_in = 0.5,
  predicted_locs = NULL
)
}
\arguments{
\item{mcmc_nngp_list}{A mcmc_nngp_list object generated using mcmc_nngp_initialize and ran using mcmc_nngp_run}

\item{X_noise_pred}{Covariates for the noise observed at predicted_locs}

\item{burn_in}{MCMC burn-in}

\item{predicted_locs}{A set of predicted locations, who can be the same as the locations in mcmc_nngp_list}
}
\value{
a list with predicted samples, summaries and predicted locations TODO QUESTION : useless to output predicted_locs?
}
\description{
Predicts the noise variance at unobserved locations
}
\examples{
\dontrun{TODO}
}
