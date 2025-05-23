% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualisation.R
\name{plot_pointillist_painting}
\alias{plot_pointillist_painting}
\title{Plots a spatial variable like a pointillist painting using R base's points. Stupid, but handy.}
\usage{
plot_pointillist_painting(locs, field, cex = 1, main = NULL, add = FALSE)
}
\arguments{
\item{locs}{spatial locations}

\item{field}{interest variable}

\item{cex}{shrinks or inflates the points}

\item{main}{main title}

\item{add}{logical (default to FALSE)}
}
\value{
a plot component (a `plot` if `add == FALSE`, `points` if `add ==TRUE`)
}
\description{
Plots a spatial variable like a pointillist painting using R base's points. Stupid, but handy.
}
\examples{
plot_pointillist_painting(locs=c(2,3,6), field=c(1,2,3))
}
