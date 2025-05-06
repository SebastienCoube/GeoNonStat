
# S3 class PP 
#' Create an object of class PP, a low rank Predictive Process used as a prior to describe 
#' spatial variations of covariance parameters
#' @param observed_locs a matrix of spatial coordinates where observations are done
#' @param matern_range either a positive number used to set the Matérn range of the PP or NULL, 
#' if NULL, a default PP range is guessed from the spatial locations.
#' @param knots either (i) a matrix of spatial knots used for the PP, 
#' (ii) or a positive integer used to set the number of knots, the knots being then found using kmeans, 
#' (iii) or NULL, a default number and placement of knots being guessed from the spatial locations.
#' @param m number of nearest neighbors to do Vecchia's approximation who goes into PP
#' @param seed integer used as a seed
#' @returns a list
#' @export
#'
#' @examples
#' observed_locs = cbind(runif(1000), runif(1000))
#' observed_locs = observed_locs[ceiling(nrow(observed_locs)*runif(3000)),]
#' # automatic 
#' pepito = PP(observed_locs)
#' # choosing manually Matérn range, too small wrt number of knots
#' pepito = PP(observed_locs, matern_range = .1)
#' # choosing manually Matérn range, way too small wrt number of knots
#' pepito = PP(observed_locs, matern_range = .01)
#' # choosing manually number of knots in order to adjust to Matérn range
#' pepito = PP(observed_locs, knots = 200, matern_range = .1)
#' # choosing manually number of knots, but picking too few for default Matérn range
#' pepito = PP(observed_locs, knots = 10)
#' # choosing manually Matérn range in order to adjust to the number of knots
#' pepito = PP(observed_locs, knots = 10, matern_range = .5)
#' # inputing an user-specified grid of knots
#' pepito = PP(observed_locs, knots = as.matrix(expand.grid(seq(-.05, 1.05, .1), seq(-.05, 1.05, .1))))
#' # inputing an user-specified grid of knots in order to adjust to small Matérn range
#' pepito = PP(observed_locs, knots = as.matrix(expand.grid(seq(-.05, 1.05, .1), seq(-.05, 1.05, .1))), matern_range = .1)
PP = function(observed_locs, matern_range = NULL, knots = NULL,  m = 10, seed=1234){
  # Suppress duplicates in observed points.
  locs_ = observed_locs[! duplicated(observed_locs),]
  # Random sample. 
  locs_ = locs_[order(runif(nrow(locs_))),]
  # Reordonate first 100000 points to help spatial coverage of vecchia. 
  locs_[seq(min(nrow(locs_), 10000)),] <- locs_[GpGp::order_maxmin(locs_[seq(min(nrow(locs_), 10000)),]),]
  # Corresponding locations
  idx = match(split(observed_locs, row(observed_locs)), split(locs_, row(locs_)))
  
  if(is.null(knots)){
    knots = min(100, nrow(observed_locs)-1)
    message(paste("number of knots set by default to", knots))
    }
  if(!is.matrix(knots)){
    knots = min(knots, nrow(observed_locs)-1)
    knots = kmeans(locs_[seq(min(nrow(locs_), 100000)),] + rnorm(2*min(nrow(locs_), 100000), 0, max(dist(locs_[seq(min(nrow(locs_), 10000)),]))/30), 
                   knots, algorithm = "Hartigan-Wong", iter.max = 50)$centers
    message("knot placement done by default using k-means")
  }
  if(is.null(matern_range)){
    matern_range = max(dist(knots))/5
    message(paste("Matérn range set by default to the fifth of the space pseudo-diameter, that is to ", signif(matern_range, 3)))
  }
  knots = knots[GpGp::order_maxmin(knots),]
  
  # Get the m closed neighbors (ordered)
  NNarray = GpGp::find_ordered_nn(rbind(knots, locs_), m)
  # Lower triangular matrix of inversed cholesky factor 
  # Add a seed here to always get the same results from vecchia
  set.seed(seed)
  sparse_chol = Matrix::sparseMatrix(
    i = row(NNarray)[!is.na(NNarray)], 
    j = NNarray[!is.na(NNarray)], 
    x = GpGp::vecchia_Linv(covparms = c(1, matern_range, .0001), covfun_name = "matern15_isotropic", locs = rbind(knots, locs_), NNarray = NNarray)[!is.na(NNarray)], 
    triangular = T
  )
  
  res = structure(
    list(
      "knots" = knots,
      "unique_reordered_locs" = locs_,
      "idx" = idx,
      "m" = m,
      "matern_range" = matern_range,
      "sparse_chol" = sparse_chol,
      "NNarray" = NNarray,
      "n_knots" = nrow(knots)
    ), class = "PP"
  )
  
  plot_knots.PP(res)
  mar_var_loss = var_loss_percentage.PP(res)
  hist(mar_var_loss, xlab = "percentage of lost variance", main = "histogram of lost marginal variance \n between the PP and the full GP")
  message(paste(
    round(mean(mar_var_loss), 1), 
    "percent of marginal variance on average is lost with the use of a PP.\n This is", 
    c(
      "great !", 
      "fairly good, but it might be improved by adding more knots or increasing the Matérn range.", 
      "quite a bit of loss, and may be fixed by adding more knots or increasing the Matérn range.")
    [1 + (mean(mar_var_loss)>3) + (mean(mar_var_loss)>10)]))
  
  return(res)
}

#' Summary of a 'PP' object 
#'
#' @param object an object of class \code{object}
#' @param ... additional arguments
#' @export
#' @examples
#' observed_locs = cbind(runif(1000), runif(1000))
#' observed_locs = observed_locs[ceiling(nrow(observed_locs)*runif(3000)),]
#' pepito = PP(observed_locs)
#' summary(pepito)
summary.PP <- function(object, ...) {
  mar_var_loss = var_loss_percentage.PP(res)
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

#' @title Plot the knots and the spatial locations of a PP
#' @examples
#' observed_locs = cbind(runif(1000), runif(1000))
#' observed_locs = observed_locs[ceiling(nrow(observed_locs)*runif(3000)),]
#' pepito = PP(observed_locs)
#' plot_knots.PP(pepito)
plot_knots.PP = function(x, ...) {
  nx <- nrow(x$unique_reordered_locs)
  ny <- nrow(x$knots)
  plot(rbind(x$unique_reordered_locs, x$knots), 
       cex = c(rep(.2, nx), rep(1, ny)),
       col = c(rep(8,  nx), rep(2, ny)), 
       pch = c(rep(16,  nx), rep(16, ny)), 
       xlab = "1st spatial coordinate",
       ylab = "2nd spatial coordinate",
       main = "Knot placement of PP"
       )
}

#' @title Compute the percentage of marginal variance who is lost because of the use of a PP
#' @examples
#' observed_locs = cbind(runif(1000), runif(1000))
#' observed_locs = observed_locs[ceiling(nrow(observed_locs)*runif(3000)),]
#' pepito = PP(observed_locs)
#' var_loss_percentage.PP(pepito)
var_loss_percentage.PP = function(x, ...) {
  PP_mar_var = apply(Matrix::solve(x$sparse_chol, Matrix::diag(
    nrow =  nrow(x$sparse_chol), ncol = nrow(x$knots)
  )), 1, function(x)
    sum(x^2))
  return((1.0001 - PP_mar_var) * 100 / 1.0001)
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
#' par(mfrow = c(1,2))
#' pepito = PP(locs)
#' # multiplying PP alone
#' res <- X_PP_mult_right(PP = pepito, Y = rnorm(pepito$n_knots))
#' GeoNonStat::plot_pointillist_painting(locs, res)
#' # multiplying PP and matrix of covariates
#' X = cbind(1, locs, rnorm(nrow(locs)))
#' res <- X_PP_mult_right(PP = pepito, X = X, Y = c(4, 1, 2, .2, rnorm(pepito$n_knots)))
#' GeoNonStat::plot_pointillist_painting(locs, res)
X_PP_mult_right = function(X = NULL, PP = NULL, use_PP = FALSE, locs_idx = NULL, Y)
{
  if(is.null(X) & is.null(PP)) stop("X and PP can't be both NULL")
  if(is.null(locs_idx)) if(!is.null(X)) locs_idx = seq(nrow(X))
  if(is.null(locs_idx)) if(!is.null(PP)) locs_idx = seq(length(PP$idx))
  Y = as.matrix(Y)
  res = matrix(0, length(locs_idx), ncol(Y))
  if(!is.null(X)) res = res + X  %*% Y[seq(ncol(X)),]
  if(!is.null(PP)) {
    if(!is.null(X)) Y =  Y[-seq(ncol(X)),, drop = F] 
    V = matrix(0, nrow(PP$sparse_chol), ncol(Y))
    V[seq(nrow(Y)),] = Y
    res = res + as.matrix(Matrix::solve(PP$sparse_chol, V, triangular = T))[-seq(nrow(PP$knots)),,drop =F][PP$idx[locs_idx],,drop =F]
  }
  colnames(res) = c("det", "an", "an")[seq(ncol(res))]
  res
}


#' Comparison between PP and NNGP
#' Allows to see if a Predictive Process has enough knots.
#' Plots two samples, one from a Predictive Process, and one from the Nearest Neighbor Gaussian Process the PP is obtained from.
#' If there are not enough knots, the PP should be over-smoothed with respect to the NNGP.
#' @param PP a Predictive Process object (produced by `PP()`)
#'
#' @returns a plot
#' @export
#'
#' @examples
#' obs_locs <- matrix(rnorm(50), ncol=2)
#' pepito <- PP(obs_locs)
#' compare_PP_NNGP(pepito, 1)
compare_PP_NNGP = function(PP, cex = .3) {
  op <- par("mfrow")
  seed_vector =  rnorm(PP$n_knots + nrow(PP$unique_reordered_locs))
  par(mfrow = c(1, 2))
  plot_pointillist_painting(
    PP$unique_reordered_locs[PP$idx, ],
    X_PP_mult_right(PP = PP, use_PP = T, Y = seed_vector[seq(PP$n_knots)]),
    cex = cex,
    main = "NNGP into PP"
  )
  points(PP$knots, pch = 3, cex = cex)
  plot_pointillist_painting(
    rbind(PP$knots, PP$unique_reordered_locs),
    as.vector(Matrix::solve(PP$sparse_chol, seed_vector)),
    cex = cex,
    main = "NNGP"
  )
  par(mfrow = op)
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
#' locs = cbind(runif(10000), runif(10000))
#' par(mfrow = c(1,2))
#' # comparing several PP approximations and testing PP mult
#' pepito = PP(locs)
#' res <- X_PP_crossprod(X = matrix(1, 10000), PP = pepito, Y = matrix(rnorm(10000)))

X_PP_crossprod = function(X, PP = NULL, Y, locs_idx = NULL)
{
  if(is.null(locs_idx))locs_idx = seq(nrow(X))
  Y = as.matrix(Y)
  res = crossprod(x = X, y = Y)
  if(!is.null(PP))
  {
    res = 
      rbind(
        res, 
        Matrix::solve(
          Matrix::t(PP$sparse_chol), 
          rbind(
            matrix(0, nrow(PP$knots), ncol(Y)), 
            (
              Matrix::sparseMatrix(x = 1, i = PP$idx, j = seq(length(PP$idx))) %*% # matrix for redunant locations and reordering in PP
                Matrix::sparseMatrix(i = locs_idx, j = seq(nrow(Y)), dims = c(length(PP$idx), nrow(Y)))  # matrix for redunant locations and reordering between Y and PP
            )%*%
              Y
          ))[1:PP$n_knots,,drop=F]
      )
  }
  as.matrix(res)
}

### # simulate locs
### locs = cbind(runif(10000), runif(10000))
### par(mfrow = c(1,2))
### # comparing several PP approximations and testing PP mult
### range= .1
### n_PP = 50
### PP = get_PP(locs, c(1, range, 1.5, 0), n_PP = n_PP, m = 15)
### GeoNonStat::plot_pointillist_painting(locs, field = X_PP_mult_right(PP = PP, use_PP = T, Y = rnorm(n_PP)))
### points(PP$knots, pch = 16, cex = .5)
### n_PP = 100
### PP = get_PP(locs, c(1, range, 1.5, 0), n_PP = n_PP, m = 15)
### GeoNonStat::plot_pointillist_painting(locs, field = X_PP_mult_right(PP = PP, use_PP = T, Y = rnorm(n_PP)))
### points(PP$knots, pch = 16, cex = .5)
### # ploting one PP basis
### GeoNonStat::plot_pointillist_painting(locs, 
###                                   Matrix::solve(PP$sparse_chol, 
###                                                 Matrix::sparseMatrix(
###                                                   i = seq(nrow(PP$knots)), 
###                                                   j = seq(nrow(PP$knots)), 
###                                                   x = 1, 
###                                                   dims = c(nrow(PP$sparse_chol), nrow(PP$knots))
###                                                 ))[-seq(nrow(PP$knots)),][PP$idx,5]
### )
### # testing PP crossprod
### X = matrix(rnorm(10*nrow(PP$unique_reordered_locs)), ncol = 10)
### Y = matrix(rnorm(nrow(X)*3), ncol=3)
### X_PP_crossprod(X = X, PP = PP, use_PP = T, Y = Y)- 
###   Matrix::t(Matrix::crossprod(
###     Y, 
###     cbind(
###       X,
###       Matrix::solve(PP$sparse_chol, 
###                     Matrix::sparseMatrix(
###                       i = seq(nrow(PP$knots)), 
###                       j = seq(nrow(PP$knots)), 
###                       x = 1, 
###                       dims = c(nrow(PP$sparse_chol), nrow(PP$knots))
###                     ))[-seq(nrow(PP$knots)),][PP$idx,]
###     )
###   ))
 
