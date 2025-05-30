% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeoNonStat.R
\name{process_vecchia}
\alias{process_vecchia}
\title{Vecchia approximation setup}
\usage{
process_vecchia(observed_locs, m)
}
\arguments{
\item{observed_locs}{a matrix of spatial coordinates where observations are done}

\item{m}{number of nearest neighbors to do Vecchia's approximation}
}
\value{
a list
}
\description{
Vecchia approximation setup
}
\examples{
  set.seed(100)
  size <- 2000
  observed_locs = cbind(runif(size), runif(size))
  res <- process_vecchia(observed_locs, m=10) 
}
