% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Useful_stuff.R
\name{beta_prior_log_dens}
\alias{beta_prior_log_dens}
\title{Compute prior logarithmic density of a matrix
TODO A reformater.}
\usage{
beta_prior_log_dens(
  beta,
  n_PP,
  beta0_mean,
  beta0_var,
  chol_crossprod_X,
  log_scale
)
}
\arguments{
\item{beta}{TODO}

\item{n_PP}{TODO}

\item{log_scale}{TODO}

\item{beta_mean}{TODO}

\item{beta_precision}{TODO}
}
\value{
a numeric value
}
\description{
Compute prior logarithmic density of a matrix
TODO A reformater.
}
\examples{
beta = matrix(rnorm(10))
bm = matrix(rnorm(5), 5, 1)
bp = diag(exp(rnorm(5)), 5, 5)
ls = rnorm(1)
# TODO a changer beta_prior_log_dens(beta, n_PP = 5, beta_mean = bm, beta_precision = bp, log_scale = ls)
}
