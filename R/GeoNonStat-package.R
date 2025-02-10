#' @title GeoNonStat: Hierarchical Nonstationary Nearest Neighbor Gaussian Process
#' @description This package provides functions for hierarchical nonstationary nearest neighbor Gaussian processes.
#'
#' @docType package
#' @name GeoNonStat
#'
#' @import Rcpp
#' @import abind
#' @import Matrix
#' @import GpGp
#' @import parallel
#' @import fields
#' @import expm
#' @import ellipse
#' @import BH
#' @importFrom Rcpp evalCpp
#' @importFrom methods as
#' @importFrom stats rnorm runif
#' @importFrom graphics plot
#' @useDynLib GeoNonStat, .registration = TRUE
"_PACKAGE"