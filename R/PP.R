
# S3 class PP 
#' Create an object of class PP, a low rank Predictive Process used as a prior to describe 
#' spatial variations of covariance parameters
#' @param vecchia_approx obtained by `createVecchia()`
#' @param matern_range either a positive number used to set the Matérn range of the PP or NULL (default), 
#' if NULL, a default PP range is guessed from the spatial locations.
#' @param knots either (i) a matrix of spatial knots used for the PP, 
#' (ii) or a positive integer used to set the number of knots, the knots being then found using kmeans, 
#' (iii) or NULL (default), a default number and placement of knots being guessed from the spatial locations.
#' @param seed integer. Used as a seed to reproduce results.
#' @param plot logical, should diagnostic plots of PP be produced ? Default to TRUE.
#' @returns a list
#' @export
#'
#' @examples
#' vecchia_approx = createVecchia(cbind(runif(10000), runif(10000)), 10, ncores=1)
#' # automatic 
#' pepito = createPP(vecchia_approx)
#' # choosing manually Matérn range, too small wrt number of knots
#' pepito = createPP(vecchia_approx, matern_range = .1)
#' # choosing manually Matérn range, way too small wrt number of knots
#' pepito = createPP(vecchia_approx, matern_range = .01)
#' # choosing manually number of knots in order to adjust to Matérn range
#' pepito = createPP(vecchia_approx, knots = 1000, matern_range = .1)
#' # choosing manually number of knots, but picking too few for default Matérn range
#' pepito = createPP(vecchia_approx, knots = 20)
#' # choosing manually Matérn range in order to adjust to the number of knots
#' pepito = createPP(vecchia_approx, knots = 20, matern_range = .5)
#' # inputing an user-specified grid of knots
#' pepito = createPP(vecchia_approx, knots = as.matrix(expand.grid(seq(-.05, 1.05, .05), seq(-.05, 1.05, .05))))
#' # inputing an user-specified grid of knots in order to adjust to small Matérn range
#' pepito = createPP(vecchia_approx, knots = as.matrix(expand.grid(seq(-.05, 1.05, .05), seq(-.05, 1.05, .05))), matern_range = .1)
#' pepito = createPP(vecchia_approx, knots = as.matrix(expand.grid(seq(-.05, 1.05, .05), seq(-.05, 1.05, .05))), matern_range = .05)
#' pepito = createPP(vecchia_approx, knots = as.matrix(expand.grid(seq(-.05, 1.05, .025), seq(-.05, 1.05, .025))), matern_range = .05)
createPP = function(vecchia_approx, matern_range = NULL, knots = NULL, seed=1234, plot=TRUE){
  # TODO : est-ce que le nombre de knots peut être plus grand que la dimention des locs de vecchia_approx ?
  # Generate knots
  if(is.null(knots)){
    knots = min(100, vecchia_approx$n_locs-1)
    message(paste("number of knots set to", knots))
  }
  if(!is.matrix(knots)){
    knots = min(knots, nrow(vecchia_approx$observed_locs)-1)
    knots = generate_knots_from_kmeans(knots, vecchia_approx$locs)
    message("knot placement done by default using k-means")
  }
  
  # matern range
  if(is.null(matern_range)){
    matern_range = max(dist(knots))/5
    message(paste("Matérn range set to ", signif(matern_range, 3))," the fifth of the space pseudo-diameter")
  }
  
  # knots order
  if(is.null(rownames(knots))) rownames(knots) <- seq_len(nrow(knots))
  knots = knots[GpGp::order_maxmin(knots),]
  
  # NNarray
  NNarray = rbind(
    GpGp::find_ordered_nn(knots, nrow(vecchia_approx$NNarray)-1), 
    cbind(nrow(knots) + seq(vecchia_approx$n_locs), 
          FNN::get.knnx(query = vecchia_approx$locs, 
                        data = knots, 
                        k = nrow(vecchia_approx$NNarray)-1)$nn.index)
  )
  
  # Cholesky matrix
  combined_locs <- rbind(knots, vecchia_approx$locs)
  Linv_vals <- GpGp::vecchia_Linv(
    covparms = c(1, matern_range, 1e-6),
    covfun_name = "matern15_isotropic",
    locs = combined_locs, 
    NNarray = NNarray
  )
  notnaNNarray <- !is.na(NNarray)
  sparse_chol = Matrix::sparseMatrix(
    i = row(NNarray)[notnaNNarray], 
    j = NNarray[notnaNNarray], 
    x = Linv_vals[notnaNNarray], 
    triangular = TRUE
  )
  
  res = structure(
    list(
      "knots" = knots,
      "matern_range" = matern_range,
      "sparse_chol" = sparse_chol,
      "n_knots" = nrow(knots),
      "vecchia_locs" = vecchia_approx$locs
    ), class = "PP"
  )
  
  if(plot) {
    plot.PP(res, mar_var_loss=TRUE)
  } else {
    varloss <- var_loss_percentage.PP(res)
  }
  
  return(res)
}

#' Generate spatial knots using k-means clustering
#'
#' Selects a set of spatial knots using k-means clustering from the observed locations.
#' If the number of locations is large, a subsample of up to 10,000 points is used,
#' and a small random perturbation is added to break ties and improve cluster separation.
#'
#' @param knots_number Integer. The number of knots (i.e., clusters) to generate.
#' @param locs A matrix of spatial coordinates (typically with two columns for 2D locations).
#' 
#' @return A matrix of size \code{knots_number} × ncol(\code{locs}) containing the spatial knot coordinates (cluster centers).

#' @examples
#' locs <- cbind(runif(5000), runif(5000))
#' knots <- generate_knots_from_kmeans(100, locs)
#' plot(locs, col = "grey", pch = 16, cex = 0.5)
#' points(knots, col = "red", pch = 19)
generate_knots_from_kmeans <- function(knots_number, locs) {
  n_sample <- min(nrow(locs), 10000)
  sampled_locs <- locs[seq(n_sample), ]
  noise <- matrix(rnorm(2 * n_sample, 0, max(dist(sampled_locs)) / 20), ncol = 2)
  centers <- kmeans(sampled_locs + noise, knots_number,
                    algorithm = "Hartigan-Wong", iter.max = 50)$centers
  return(centers)
}

#' Summary of a 'PP' object 
#'
#' @param object an object of class \code{PP}
#' @param ... additional arguments
#' @export
#' @examples
#' vecchia_approx = createVecchia(cbind(runif(100), runif(100)), 10, ncores=1)
#' pepito = createPP(vecchia_approx, plot=FALSE)
#' summary(pepito)
summary.PP <- function(object, ...) {
  cat("Object of class \'PP\' with", 
      object$n_knots, "knots,", 
      "based on", dim(object$vecchia_locs)[1], "locations,",
      "matérn range =", object$matern_range)
  var_loss_percentage.PP(object)
  return(invisible(NULL))
}


#' @title Compute the percentage of marginal variance who is lost because of the use of a PP
#' @param x an object of class \code{PP}
#' @examples
#' vecchia <- createVecchia(cbind(runif(1000), runif(1000)), ncores=1)
#' pepito = createPP(vecchia, plot=FALSE)
#' var_loss_percentage.PP(pepito)
var_loss_percentage.PP = function(x) {
  PP_mar_var = apply(
    Matrix::solve(x$sparse_chol, 
                  Matrix::diag(nrow =  nrow(x$sparse_chol), ncol = nrow(x$knots))), 
    1, 
    function(x) sum(x^2)
  )
  # max(0) because of tiny numerical errors
  PP_mar_var <- (pmax(0, 1.000001 - PP_mar_var)/1.000001) * 100
  mean_mar_var <- mean(PP_mar_var)
  msg <- if (mean_mar_var > 10) {
    "quite a bit of loss, and may be fixed by adding more knots or increasing the Matérn range."
  } else if (mean_mar_var > 3) {
    "fairly good, but it might be improved by adding more knots or increasing the Matérn range."
  } else {
    "great !"
  }
  message(round(mean_mar_var, 1), 
          "% of marginal variance on average is lost with the use of a PP.\nThis is ", msg)
  
  return(PP_mar_var)
}

