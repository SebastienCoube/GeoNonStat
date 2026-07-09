#' @title GeoNonStat: Hierarchical Nonstationary Nearest Neighbor Gaussian Process
#' @description 
#' Hierarchical Nonstationary Nearest Neighbor Gaussian Process. 
#' This package allows to model observations with point coordinates in a 2D, 
#' geographic space. 
#' Explanatory variable impact on the interest variable is modeled using a 
#' linear model.  
#' Spatial auto-correlation of the data is captured using a possibly 
#' non-stationary Nearest Neighbor Gaussian process. 
#' Independent observation noise is captured using a possibly heteroskedastic 
#' Gaussian noise.   
#' The complexity levels of the Gaussian Process and the observation noise are 
#' specified by the user through an interpretable parametrization.  
#' @docType package
#' @name GeoNonStat
#'
#' @import future
#' @import future.apply
#' @importFrom parallelly supportsMulticore
#' @importFrom methods new
#' @importFrom utils capture.output
#' @import Rcpp
#' @importFrom abind abind
#' @importFrom Matrix Diagonal crossprod solve sparseMatrix t tcrossprod chol
#' @importFrom GpGp find_ordered_nn order_maxmin vecchia_Linv matern_isotropic fast_Gp_sim
#' @import fields
#' @importFrom ellipse ellipse
#' @import BH
#' @importFrom expm expm
#' @importFrom methods as
#' @importFrom stats rnorm runif ecdf predict
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom FNN get.knn
#' @importFrom coda effectiveSize
#' @importFrom stats dist dnorm kmeans lm median model.matrix quantile sd var vcov
#' @importFrom graphics plot abline barplot hist layout legend lines par points text matplot
#' @useDynLib GeoNonStat, .registration = TRUE
"_PACKAGE"
