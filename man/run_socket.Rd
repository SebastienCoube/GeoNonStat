% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{run_socket}
\alias{run_socket}
\title{mcmc_nngp_run_nonstationary_socket}
\usage{
run_socket(
  mcmc_nngp_list,
  burn_in = 0.5,
  seed = 1,
  n_cores = NULL,
  thinning = 0.1,
  plot_diags = T,
  plot_PSRF_fields = F,
  debug_outfile = NULL,
  num_threads_per_chain = 1,
  lib.loc = NULL
)
}
\arguments{
\item{mcmc_nngp_list}{a mcmc_nngp_list generated by mcmc_nngp_initialize_nonstationary}

\item{burn_in}{burn-in to compute diagnostics}

\item{seed}{random seed}

\item{n_cores}{}

\item{thinning}{}

\item{plot_diags}{}

\item{plot_PSRF_fields}{}

\item{debug_outfile}{}

\item{num_threads_per_chain}{}

\item{lib.loc}{}
}
\value{
a list
}
\description{
mcmc_nngp_run_nonstationary_socket
}
