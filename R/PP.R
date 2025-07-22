
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
#' vecchia_approx = createVecchia(observed_locs  = cbind(runif(10000), runif(10000)), 10)
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
  # Generate knots
  if(is.null(knots)){
    knots = min(100, vecchia_approx$n_locs-1)
    message(paste("number of knots set to", knots))
  }
  if(!is.matrix(knots)){
    knots = min(knots, nrow(vecchia_approx$observed_locs)-1)
    knots = max(knots, nrow(vecchia_approx$NNarray))
    knots = generate_knots_from_kmeans(knots, vecchia_approx$locs)
    message("knot placement done by default using k-means")
  }
  
  # matern range
  if(is.null(matern_range)){
    matern_range = max(dist(knots))/5
    message(paste("Matérn range set to ", signif(matern_range, 3))," the fifth of the space pseudo-diameter")
  }
  
  # knots order
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
  
  if(plot) plot.PP(res, mar_var_loss=TRUE)
  
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
#' observed_locs = cbind(runif(1000), runif(1000))
#' observed_locs = observed_locs[ceiling(nrow(observed_locs)*runif(3000)),]
#' pepito = createPP(observed_locs)
#' summary(pepito)
summary.PP <- function(object, ...) {
  mar_var_loss = var_loss_percentage.PP(object)
  message_loss = (paste(
    round(mean(mar_var_loss), 1), 
    "percent of marginal variance on average lost with the use of a PP, which is", 
    c(
      "great !", 
      "fairly good, but might be improved by adding more knots or increasing the Matérn range.", 
      "quite a bit of loss, but may be fixed by adding more knots or increasing the Matérn range.")
    [1 + (mean(mar_var_loss)>3) + (mean(mar_var_loss)>10)]))
  cat("Obect of class \'PP\' with ", 
      object$n_knots, 
      "knots, Matérn range =", 
      object$matern_range, 
      "and", message_loss)
}


#' @title Compute the percentage of marginal variance who is lost because of the use of a PP
#' @param x an object of class \code{PP}
#' @examples
#' observed_locs = cbind(runif(1000), runif(1000))
#' observed_locs = observed_locs[ceiling(nrow(observed_locs)*runif(3000)),]
#' pepito = createPP(observed_locs)
#' var_loss_percentage.PP(pepito)
var_loss_percentage.PP = function(x, ...) {
  PP_mar_var = apply(Matrix::solve(x$sparse_chol, Matrix::diag(
    nrow =  nrow(x$sparse_chol), ncol = nrow(x$knots)
  )), 1, function(x)
    sum(x^2))
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
  # max(0) because of tiny numerical errors
}



#' Multiply the concatenation of a matrix of covariates and a PP by a matrix
#' (X|PP) %*% Y
#'
#' @param X a matrix of covariates who will be multiplied by the first rows of Y
#' @param PP either a PP whose basis will be multiplied by the last columns of Y
#' @param locs_idx either a vector of integers who dispatch the PP basis to the covariates, or NULL
#' @param Y the matrix who multiplies the covariates and the PP
#'
#' @export
#' @returns a matrix
#' 
#' @examples
#' locs = cbind(runif(10000), runif(10000))
#' locs = rbind(locs, locs)
#' vecchia_approx = createVecchia(locs, 12)
#' PP = createPP(vecchia_approx)
#' covariate_coefficients = c(4, 1, 1, .5)
#' knots_coeffs = rnorm(PP$n_knots)
#' par(mfrow = c(1,2))
#' # multiplying X alone
#' X = cbind(1, locs, rnorm(nrow(locs)))
#' res1 <- X_PP_mult_right(X = X, Y = covariate_coefficients, vecchia_approx = vecchia_approx)
#' res1
#' plot_pointillist_painting(locs, res1, main = "covariates only, \n one covariate for each observation")
#' # multiplying PP alone
#' res2 <- X_PP_mult_right(PP = PP, Y = knots_coeffs, vecchia_approx = vecchia_approx)
#' plot_pointillist_painting(vecchia_approx$locs, res2, main = "PP only")
#' # multiplying PP and matrix of covariate, one obs for each location
#' X_by_loc = cbind(1, vecchia_approx$locs, rnorm(vecchia_approx$n_locs))
#' res3 <- X_PP_mult_right(PP = PP, X = X_by_loc, Y = c(covariate_coefficients, knots_coeffs), vecchia_approx = vecchia_approx)
#' plot_pointillist_painting(vecchia_approx$locs, res3, main = "PP + covariates, \n one covariate for each location")
#' # multiplying PP and matrix of covariates with an index
#' X_by_obs = cbind(1, vecchia_approx$observed_locs, rnorm(vecchia_approx$n_obs))
#' res4 <- X_PP_mult_right(X = X_by_obs, PP = PP, 
#'                        Y = c(covariate_coefficients, knots_coeffs), 
#'                        vecchia_approx = vecchia_approx
#'                        )
#' plot_pointillist_painting(vecchia_approx$observed_locs, res4, main = "PP + covariates,\n  one covariate for each observation")
#' 
X_PP_mult_right = function(X = NULL, PP = NULL, vecchia_approx, Y, permutate_PP_to_obs = F)
{
  if(is.null(X) & is.null(PP)) stop("X and PP can't be both NULL")
  # Sanity checks
  Y = as.matrix(Y)
  expected_rows <- 0
  if(!is.null(X)) expected_rows <- expected_rows + ncol(X)
  if(!is.null(PP)) expected_rows <- expected_rows + nrow(PP$knots)
  if(nrow(Y) != expected_rows) {
    stop("Y should have ", expected_rows, " rows it has ", nrow(Y))
  }
  if(!permutate_PP_to_obs) locs_idx = seq(vecchia_approx$n_locs)
  if(permutate_PP_to_obs) locs_idx = vecchia_approx$locs_match
  res = matrix(0, length(locs_idx), ncol(Y))
  # Multiply X and Y
  if(!is.null(X)) res = res + X  %*% Y[seq(ncol(X)),]
  if(!is.null(PP)) {
    if(!is.null(X)) Y =  Y[-seq(ncol(X)),, drop = F] 
    V = matrix(0, nrow(PP$sparse_chol), ncol(Y))
    V[seq(nrow(Y)),] = Y
    res = res + as.matrix(Matrix::solve(PP$sparse_chol, V, triangular = T))[-seq(nrow(PP$knots)),,drop =F][locs_idx,,drop =F]
  }
  if(ncol(res)==3)colnames(res) = c("det", "an", "an")
  if(ncol(res)==1)colnames(res) = "det"
  res
}

#' Do the cross-product of the concatenation of a matrix of covariates and a PP, and a matrix
#' t(X|PP) %*% Y
#'
#' @param X a matrix of covariates who will be multiplied by the first rows of Y
#' @param PP either a PP whose basis will be multiplied by the last columns of Y
#' @param locs_idx either a vector of integers who dispatch the PP basis to the covariates, or NULL
#' @param Y the matrix who multiplies the covariates and the PP
#'
#'
#' @returns a matrix
#' @export
#'
#' @examples
#' locs = cbind(runif(5000), runif(5000))
#' locs = rbind(locs, locs)
#' vecchia_approx = createVecchia(locs)
#' PP = createPP(vecchia_approx)
#' 
#' # just surrogate of crossprod
#' X = matrix(rnorm(100000), 10000)
#' Y = matrix(rnorm(30*nrow(X)), nrow(X))
#' res1 <- X_PP_crossprod(X = X, PP = NULL, Y = Y, vecchia_approx = vecchia_approx)
#' crossprod(X, Y) - res1
#' 
#' # crossprod + PP with observations of X on the locs
#' X = matrix(rnorm(20*vecchia_approx$n_locs), vecchia_approx$n_locs)
#' Y = matrix(rnorm(30*nrow(X)), nrow(X))
#' res2 <- X_PP_crossprod(X = X, PP = PP, Y = Y, vecchia_approx = vecchia_approx, permutate_PP_to_obs = F)
#' hist(as.vector(res2 - crossprod(
#'  cbind(X, Matrix::solve(PP$sparse_chol, diag(1, nrow(PP$sparse_chol), PP$n_knots))[-seq(PP$n_knots),]), Y
#' )))
#' 
#' # crossprod + PP with PP dispatched to observations of X
#' X = matrix(rnorm(20*vecchia_approx$n_obs), vecchia_approx$n_obs)
#' Y = matrix(rnorm(30*nrow(X)), nrow(X))
#' res2 <- X_PP_crossprod(X = X, PP = PP, Y = Y, vecchia_approx = vecchia_approx, permutate_PP_to_obs = T)
#' hist(as.vector(res2 - crossprod(
#'  cbind(X, Matrix::solve(PP$sparse_chol, diag(1, nrow(PP$sparse_chol), PP$n_knots))[-seq(PP$n_knots),][vecchia_approx$locs_match,]), Y
#' )))

 X_PP_crossprod = function(X, PP = NULL, Y, vecchia_approx, permutate_PP_to_obs = F)
{
  Y = as.matrix(Y)
  res = crossprod(x = X, y = Y)
  if(!is.null(PP))
  {
    if(permutate_PP_to_obs) Y = vecchia_approx$locs_match_matrix %*% Y
    res = 
      rbind(
        res, 
        Matrix::solve(
          Matrix::t(PP$sparse_chol), 
          rbind(matrix(0, nrow(PP$knots), ncol(Y)), Y)
        )[1:PP$n_knots,,drop=F]
      )
  }
  as.matrix(res)
}

 
