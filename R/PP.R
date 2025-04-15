
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
#' set.seed(100)
#' locs = cbind(runif(1000), runif(1000))
#' PP = get_PP(
#'   observed_locs = locs, # spatial sites
#'   matern_range = .1,
#'   n_PP = 50, # number of knots
#'   m = 15 # number of NNGP parents
#' )

PP = function(observed_locs, matern_range = NULL, knots = NULL,  m = 10, seed=1234){
  # Suppress duplicates in observed points.
  locs_ = observed_locs[! duplicated(observed_locs),]
  # Random sample. 
  locs_ = locs_[order(runif(nrow(locs_))),]
  # Reordonate first 100000 points to help spatial coverage of vecchia. 
  locs_[seq(min(nrow(locs_), 10000)),] <- locs_[GpGp::order_maxmin(locs_[seq(min(nrow(locs_), 100000)),]),]
  # Corresponding locations
  idx = match(split(observed_locs, row(observed_locs)), split(locs_, row(locs_)))
  
  if(is.null(knots)){
    knots = 100
    message("number of knots set by default to 50")
    }
  if(!is.matrix(knots)){
    knots = kmeans(locs_[seq(min(nrow(locs_), 10000)),] + rnorm(20000, 0, max(dist(locs_[seq(min(nrow(locs_), 10000)),]))/15), 
                   knots, algorithm = "Hartigan-Wong", iter.max = 50)$centers
    knots = knots[GpGp::order_maxmin(knots),]
    message("knot placement done by default using k-means")
  }
  if(is.null(matern_range)){
    matern_range = max(dist(knots))/5
    message(paste("Matérn range set by default to the fifth of the space pseudo-diameter, that is to ", signif(matern_range, 3)))
  }
  
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
  message(paste(round(mean(mar_var_loss), 1), "percent of marginal variance on average is lost with the use of a PP.\n This is", c("great !", "fairly good, but it might be improved by adding more knots or increasing the Matérn range.", "quite a bit of loss, and may be fixed by adding more knots or increasing the Matérn range.")[1 + (mean(mar_var_loss)>3) + (mean(mar_var_loss)>10)]))
  
  return(res)
}




plot_knots.PP = function(x, ...) {
  plot(rbind(x$unique_reordered_locs, x$knots), 
       cex = c(rep(.2, nrow(x$unique_reordered_locs)), rep(1, nrow(x$knots))),
       col = c(rep(8,  nrow(x$unique_reordered_locs)), rep(2, nrow(x$knots))), 
       pch = c(rep(16,  nrow(x$unique_reordered_locs)), rep(16, nrow(x$knots))), 
       xlab = "first spatial coordinate",
       ylab = "second spatial coordinate",
       main = "Knot placement"
       )
}



var_loss_percentage.PP = function(x, ...) {
  PP_mar_var = apply(Matrix::solve(x$sparse_chol, Matrix::diag(nrow=  nrow(x$sparse_chol), ncol = nrow(x$knots))), 1, function(x)sum(x^2))
  return((1.0001 - PP_mar_var)*100/1.0001)
}


set.seed(100)
locs = cbind(runif(10000), runif(10000))
x = PP(
  observed_locs = locs, knots = 100, seed = 4
)
plot_knots.PP(x)
var_loss_percentage.PP(x)






#' Title TODO
#'
#' @param X TODO
#' @param PP TODO
#' @param use_PP Logical, use PP ? Default to FALSE
#' @param locs_idx TODO
#' @param Y TODO
#'
#' @returns a matrix
#' @export
#'
#' @examples
#' locs = cbind(runif(100), runif(100))
#' par(mfrow = c(1,2))
#' # comparing several PP approximations and testing PP mult
#' range= .1
#' n_PP = 50
#' PP = get_PP(locs, c(1, range, 1.5, 0), n_PP = n_PP, m = 15)
#' res <- X_PP_mult_right(PP = PP, use_PP = TRUE, Y = rnorm(n_PP))
X_PP_mult_right = function(X = NULL, PP = NULL, use_PP = FALSE, locs_idx = NULL, Y)
{
  if(is.null(X) & is.null(PP)) stop("X and PP can't be both NULL")
  if(is.null(locs_idx)) if(!is.null(X)) locs_idx = seq(nrow(X))
  if(is.null(locs_idx)) if(!is.null(PP)) locs_idx = seq(length(PP$idx))
  Y = as.matrix(Y)
  res = matrix(0, length(locs_idx), ncol(Y))
  if(!is.null(X)) res = res + X  %*% Y[seq(ncol(X)),]
  if(use_PP) {
    if(!is.null(X)) Y =  Y[-seq(ncol(X)),, drop = F] 
    V = matrix(0, nrow(PP$sparse_chol), ncol(Y))
    V[seq(nrow(Y)),] = Y
    res = res + as.matrix(Matrix::solve(PP$sparse_chol, V, triangular = T))[-seq(nrow(PP$knots)),,drop =F][PP$idx[locs_idx],,drop =F]
  }
  colnames(res) = c("det", "an", "an")[seq(ncol(res))]
  res
}


#' Title TODO
#'
#' @param X TODO
#' @param PP TODO
#' @param use_PP logical, use PP ? Default to FALSE
#' @param Y TODO
#' @param locs_idx TODO 
#'
#' @returns a matrix
#' @export
#'
#' @examples
#' set.seed(123)
#' locs = cbind(runif(100), runif(100))
#' PP = get_PP(locs, c(1, .1, 1.5, 0), n_PP = 50, m = 15)
#' X = matrix(rnorm(10*nrow(PP$unique_reordered_locs)), ncol = 10)
#' Y = matrix(rnorm(nrow(X)*3), ncol=3)
#' res <- X_PP_crossprod(X = X, PP = PP, use_PP = T, Y = Y)
#' res <- X_PP_crossprod(X = X, Y = Y)
X_PP_crossprod = function(X, PP = NULL, use_PP = F,  Y, locs_idx = NULL)
{
  if(is.null(locs_idx))locs_idx = seq(nrow(X))
  Y = as.matrix(Y)
  res = crossprod(x = X, y = Y)
  if(use_PP)
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
          ))[1:PP$n_PP,,drop=F]
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
 
