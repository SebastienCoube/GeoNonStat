% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PP.R
\name{summary.PP}
\alias{summary.PP}
\title{Summary of a 'PP' object}
\usage{
\method{summary}{PP}(object, ...)
}
\arguments{
\item{object}{an object of class \code{PP}}

\item{...}{additional arguments}
}
\description{
Summary of a 'PP' object
}
\examples{
observed_locs = cbind(runif(1000), runif(1000))
observed_locs = observed_locs[ceiling(nrow(observed_locs)*runif(3000)),]
pepito = createPP(observed_locs)
summary(pepito)
}
