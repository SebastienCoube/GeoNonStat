#' @title GeoNonStat: Hierarchical Nonstationary Nearest Neighbor Gaussian Process
#' @description This package provides functions for hierarchical nonstationary nearest neighbor Gaussian processes.
#'
#' @docType package
#' @name GeoNonStat
#'
#' @import Rcpp
#' @importFrom abind abind
#' @importFrom Matrix Diagonal crossprod solve sparseMatrix t tcrossprod chol 
#' @importFrom GpGp find_ordered_nn order_maxmin vecchia_Linv matern_isotropic fast_Gp_sim
#' @import parallel
#' @import fields
#' @importFrom ellipse ellipse
#' @import BH
#' @importFrom expm expm
#' @importFrom Rcpp evalCpp
#' @importFrom methods as
#' @importFrom stats rnorm runif
#' @importFrom graphics plot
#' @useDynLib GeoNonStat, .registration = TRUE
"_PACKAGE"