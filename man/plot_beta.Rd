% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualisation.R
\name{plot_beta}
\alias{plot_beta}
\title{Title TODO}
\usage{
plot_beta(
  beta_arrays,
  iterations,
  starting_proportion = 0.5,
  varname,
  var_names = NULL
)
}
\arguments{
\item{beta_arrays}{list of arrays TODO}

\item{iterations}{numerical vector TODO}

\item{starting_proportion}{numeric value, default to 0.5 TODO}

\item{varname}{name of variable}

\item{var_names}{? TODO}
}
\value{
several plots
}
\description{
Title TODO
}
\examples{
beta_arrays = list(
  array(rnorm(400), dim = c(2, 2, 100)),
  array(rnorm(400), dim = c(2, 2, 100)),
  array(rnorm(400), dim = c(2, 2, 100))
)
plot_beta(beta_arrays, seq(100), starting_proportion = .5, varname = "Example", var_names = c(1, 2))
beta_arrays = lapply(seq(3), function(i){
   res = array(data = 0, dim = c(10, 3, 100))
   res[,1,] = rnorm(length(res[,1,]))
   res
 })
plot_beta(beta_arrays, seq(100), starting_proportion = .5, varname = "Example", var_names = c(1, 2))
}
