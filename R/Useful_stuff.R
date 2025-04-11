#' Title TODO
#'
#' @param derivatives TODO
#' @param left_vector TODO
#' @param right_vector TOSO
#' @param NNarray TODO
#'
#' @returns a matrix
#' @export
#'
#' @examples
#' \dontrun{TODO}
derivative_sandwiches = function(
    derivatives, 
    left_vector, 
    right_vector, 
    NNarray
)
{
  M <- matrix(0, length(left_vector), length(derivatives))
  for( i in seq(length(derivatives))) {
    M[,i] <- derivative_sandwich(derivatives[[i]], left_vector, right_vector, NNarray)
  }
  # changing basis between det/aniso and canonical
  if(ncol(M)==3) {
    M <- M %*% t(matrix(
    c(1/sqrt(2), 1/sqrt(2),  0, 
      1/sqrt(2), -1/sqrt(2), 0,
      0,       0,        1), 3))*sqrt(2)
  }
  return(M)
}

#' Title TODO
#'
#' @param sparse_chol_and_grad TODO
#' @param NNarray TODO
#'
#' @returns a matrix
#' @export
#'
#' @examples
#' \dontrun{TODO}
log_determinant_derivatives = function(sparse_chol_and_grad, NNarray)
{
  M <- matrix(0, nrow(NNarray), length(sparse_chol_and_grad[[2]]))
  for( i in seq(length(sparse_chol_and_grad[[2]]))) {
    M[,i] <- log_determinant_derivative(derivative = sparse_chol_and_grad[[2]][[i]], 
                                       compressed_sparse_chol = sparse_chol_and_grad[[1]], 
                                       NNarray = NNarray)
  }
  if(ncol(M)==3) {
    M <- M %*% matrix(
    c(1/sqrt(2), 1/sqrt(2),  0, 
      1/sqrt(2), -1/sqrt(2), 0,
      0,       0,        1), 3)*sqrt(2)
  }
  return(M)
}


#' Exponential of a square matrix, adding a small numeric value on the diagonal
#'
#' @param coords a numeric vector
#' @param eps a numeric value, default to .0001
#'
#' @returns a square matrix
#' @export
#'
#' @examples
#' expmat(c(1,2,3,4,5,6))
expmat = function(coords, eps = .0001)
{
  res = expm::expm(symmat(coords))
  if(eps != 0) diag(res) <- diag(res) + eps
  return(res)
}


#' Create symetric matrix from coordinates
#'
#' @param coords a numeric vector, of length 1, 3 or 6
#'
#' @returns a matrix
#' @export
#'
#' @examples
#' symmat(c(1,2,3,4,5,6))
symmat = function(coords)
{
  # Trouver la taille n de la matrice symétrique n x n
  # Longueur du vecteur doit être égale à n + n*(n-1)/2 = n*(n+1)/2
  n <- (sqrt(8 * length(coords) + 1) - 1) / 2
  if (n != floor(n)) {
    stop("length of coords incompatible with a symetric matrix.")
  }
  mat <- matrix(0, n, n)
  
  diag_indices <- which(row(mat) == col(mat))
  lower_indices <- which(lower.tri(mat))

  # 1st fill diagonal elements
  mat[diag_indices] <- coords[1:n]
  # Then lower and upper mat
  mat[lower_indices] <- coords[(n + 1):length(coords)]
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  
  return(mat)
}

#' Title TODO
#'
#' @param beta TODO
#' @param PP predictive process obtained through get_PP
#' @param use_PP should the PP be used ? Default to FALSE
#' @param X TODO
#' @param locs_idx match between PP basis function and locs.
#'
#' @returns a numeric vector
#' @export
#'
#' @examples
#' locs = cbind(runif(100), runif(100))
#' n_PP = 50
#' PP = get_PP(locs, c(1, .1, 1.5, 0), n_PP = n_PP, m = 15)
#' X = matrix(rnorm(10*nrow(PP$unique_reordered_locs)), ncol = 10)
#' res <- variance_field(beta = rnorm(n_PP), PP = PP, use_PP = TRUE, X = X)
#' res <- variance_field(beta = rnorm(n_PP), X = X)
variance_field = function(beta,
                          PP = NULL,
                          use_PP = F,
                          X,
                          locs_idx = NULL)
{
  as.vector(exp(
    X_PP_mult_right(
      X = X,
      PP = PP,
      use_PP = use_PP,
      locs_idx = locs_idx,
      Y = beta
    )
  ))
}


#' Computes a Vecchia sparse Cholesky factor and its derivatives
#' 
#' The function uses range_beta one one hand, and range_X and PP on the other hand, to compute nonstationary range parameters. 
#' The Vecchia approximation of the sparse Cholesky factor of the precision is then computed. 
#' Derivatives can be computed too. 
#' Warning (for developing users) : the sparse cholesky factor is computed using the radius/anisotropy parametrization for range_beta.  
#' However, for efficiency of implementation, the derivatives are given along the canonical parametrization.  
#' Re-parametrization is done automatically in GeoNonStat::derivative_sandwiches
#' 
#' @param range_beta parameter for the range.
#' If the covariance is anisotropic, it must have 3 columns. It the covariance is isotropic, it must have 1 column. 
#' The first coefficients are multiplied with range_X 
#' The last coefficients are multiplied with the spatial basis functions of PP
#' @param NNarray Vecchia parents array provided by GpGp::find_ordered_nn
#' @param locs matrix of spatial sites
#' @param range_X covariates for range
#' @param PP predictive process obtained through get_PP
#' @param use_PP should the PP be used ? (redundant when using the function "by hand", but handy when automating)
#' @param compute_derivative logical, indicates if derivatives of Vecchia factors are to be computed
#' @param nu Matern smoothness Default to 1.5. Can be 0.5 or 1.5.
#' @param anisotropic Logical, default to FALSE. TODO
#' @param sphere Logical, default to FALSE. TODO
#' @param num_threads numerical, number of treads to use. Default to 1.
#' @param locs_idx match between the duplicated locations used to buile the PP basis function and the non-redundant locs used to compute the sparse chol
#'
#' @returns a list
#' @export
#'
#' @examples
#' locs = cbind(seq(100)/10, 0)
#' NNarray = GpGp::find_ordered_nn(locs, 10)
#' res <- compute_sparse_chol(
#'           range_beta = matrix(.5/sqrt(2),1,1), 
#'           NNarray = NNarray, 
#'           locs = locs,
#'           use_PP = F, 
#'           num_threads = 1, 
#'           anisotropic = F,
#'           range_X = matrix(1, nrow(locs), 1), 
#'           nu = 1.5
#'         )
#' \dontrun{
#' 
#' set.seed(1)
#' observed_locs =  cbind(runif(10000), runif(10000))  # creating spatial locations 
#' observed_locs = rbind(observed_locs, observed_locs) 
#' observed_locs = observed_locs[order(runif(nrow(observed_locs))),]# duplicating observed_locs
#' unique_locs = observed_locs[!duplicated(observed_locs),]# getting unique spatial locations
#' hctam_scol_1 =  
#'   match(
#'     split(unique_locs, row(unique_locs)), 
#'     split(observed_locs, row(observed_locs)), 
#'   )  # match between unique observed_locs and duplicated observed_locs (reverse of locs_match)
#' NNarray = GpGp::find_ordered_nn(unique_locs, 10)  # Vecchia Nearest Neighbor Array
#' range_X = cbind(1, unique_locs) # Covariates for the range
#' PP = GeoNonStat::get_PP(observed_locs = observed_locs, matern_range = .1, lonlat = F, n_PP = 20, m = 10) # Predictive Process is defined on duplicated observed_locs
#' 
#' 
#' # sampling white noise to reuse 
#' v = rnorm(nrow(unique_locs))
#' 
#' # anisotropic case 
#' range_beta = matrix(rnorm(23*3), 23, 3) # regression coefficients for the range
#' # Note : range_beta has 3 col because aniso
#' range_beta[1,1] = -4
#' sparse_chol = compute_sparse_chol(
#'   range_beta =  range_beta, NNarray = NNarray, 
#'   locs = unique_locs, range_X = range_X, 
#'   PP = PP, use_PP = T, compute_derivative = T, 
#'   nu = 1.5, anisotropic = T,# Note : anisotropic is T
#'   sphere = F, num_threads = 1, locs_idx = hctam_scol_1)
#' # plotting a sample generated from sparse chol
#' GeoNonStat::plot_pointillist_painting(
#'   unique_locs, 
#'   Matrix::solve(Matrix::sparseMatrix(
#'     i = row(NNarray)[!is.na(NNarray)], 
#'     j = (NNarray[!is.na(NNarray)]), 
#'     x = (sparse_chol[[1]][!is.na(NNarray)]), 
#'     triangular = T
#'   ), v)
#' )
#' 
#' 
#' # Showing that anisotropic case comprises isotropic case 
#' range_beta[,-1] = 0
#' sparse_chol = compute_sparse_chol(
#'   range_beta =  range_beta, NNarray = NNarray, 
#'   locs = unique_locs, range_X = range_X, 
#'   PP = PP, use_PP = T, compute_derivative = T, 
#'   nu = 1.5, anisotropic = T,# Note : anisotropic is T
#'   sphere = F, num_threads = 1, locs_idx = hctam_scol_1)
#' # plotting a sample generated from sparse chol. It is locally isotropic, and the same as the next !
#' GeoNonStat::plot_pointillist_painting(
#'   unique_locs, 
#'   Matrix::solve(Matrix::sparseMatrix(
#'     i = row(NNarray)[!is.na(NNarray)], 
#'     j = (NNarray[!is.na(NNarray)]), 
#'     x = (sparse_chol[[1]][!is.na(NNarray)]), 
#'     triangular = T
#'   ), v)
#' )
#' 
#' 
#' # isotropic case 
#' range_beta = range_beta[,1,drop = F] # regression coefficients for the range
#' # Note : range_beta has 1 col because iso
#' range_beta[1,1] = -4
#' sparse_chol = compute_sparse_chol(
#'   range_beta =  range_beta, NNarray = NNarray, 
#'   locs = unique_locs, range_X = range_X, 
#'   PP = PP, use_PP = T, compute_derivative = T, 
#'   nu = 1.5, anisotropic = F, 
#'   sphere = F, num_threads = 1, locs_idx = hctam_scol_1)
#' # plotting a sample generated from sparse chol
#' GeoNonStat::plot_pointillist_painting(
#'   unique_locs, 
#'   Matrix::solve(Matrix::sparseMatrix(
#'     i = row(NNarray)[!is.na(NNarray)], 
#'     j = (NNarray[!is.na(NNarray)]), 
#'     x = (sparse_chol[[1]][!is.na(NNarray)]), 
#'     triangular = T
#'   ), v)
#' )
#' 
#' }
compute_sparse_chol = function(range_beta, 
                               NNarray, 
                               locs, 
                               range_X = NULL, 
                               PP = NULL, 
                               use_PP = F, 
                               compute_derivative = T, 
                               nu = 1.5, 
                               anisotropic = F,
                               sphere = F,
                               num_threads = 1,
                               locs_idx = NULL)
{
  if (!nu%in%c(.5, 1.5)) stop("nu must be equal to 0.5 or 1.5")
  # converting to canonical basis
  if(ncol(range_beta)==3) {
    range_beta = range_beta %*% matrix(
    c(1/sqrt(2), 1/sqrt(2),  0, 
      1/sqrt(2), -1/sqrt(2), 0,
      0,       0,        1), 3)*sqrt(2)
  }
  if(ncol(range_beta)==1) {
    range_beta = range_beta# / sqrt(2)
  }
  
  log_range = as.matrix(
    X_PP_mult_right(
      X = range_X, 
      PP = PP, 
      Y = range_beta,  
      use_PP = use_PP, 
      locs_idx = locs_idx))
  #GeoNonStat::plot_ellipses(locs, log_range)
  # exp locally isotropic
  covfun_name <- NULL
  if((!anisotropic) & (nu==0.5)) covfun_name <- "nonstationary_exponential_isotropic" 
  # matern locally isotropic
  if((!anisotropic) & (nu==1.5)) covfun_name <- "nonstationary_matern_isotropic"  
  # exp locally anisotropic
  if(( anisotropic) & (nu==0.5)) covfun_name <- "nonstationary_exponential_anisotropic"
  # matern locally anisotropic
  if(( anisotropic) & (nu==1.5)) covfun_name <- "nonstationary_matern_anisotropic"
  
  res <- nonstat_vecchia_Linv(num_threads=num_threads,
                              log_range = log_range*2, 
                              covfun_name = covfun_name  , 
                              sphere = sphere, 
                              locs = locs, 
                              NNarray = NNarray, 
                              compute_derivative = compute_derivative)
  
  res[[2]] = lapply(res[[2]], function(x)x*2)
  return(res)
}

### # checking equivalence of parametrizations ######
### 
### locs = cbind(seq(100)/10, 0)
### NNarray = GpGp::find_ordered_nn(locs, 10)
### M = 
###   Matrix::tcrossprod(
###     Matrix::solve(
###       Matrix::sparseMatrix(
###         i = row(NNarray)[!is.na(NNarray)],
###         j = (NNarray)[!is.na(NNarray)],
###         x= GeoNonStat::compute_sparse_chol(
###           range_beta = matrix(.5/sqrt(2),1,1), 
###           NNarray = NNarray, 
###           locs = locs,
###           use_PP = F, 
###           num_threads = 1, 
###           anisotropic = F,
###           range_X = matrix(1, nrow(locs), 1), nu = 1.5
###         )[[1]][!is.na(NNarray)],
###       )
###     )
###   ) 
### plot(locs[,1], M[,1])
### 
### locs = cbind(seq(100)/10, 0)
### NNarray = GpGp::find_ordered_nn(locs, 10)
### M = 
###   Matrix::tcrossprod(
###     Matrix::solve(
###       Matrix::sparseMatrix(
###         i = row(NNarray)[!is.na(NNarray)],
###         j = (NNarray)[!is.na(NNarray)],
###         x= GeoNonStat::compute_sparse_chol(
###           range_beta = matrix(c(.5,0,0),1), 
###           NNarray = NNarray, 
###           locs = locs,
###           use_PP = F, 
###           num_threads = 1, 
###           anisotropic = T,
###           range_X = matrix(1, nrow(locs), 1), nu = 1.5
###         )[[1]][!is.na(NNarray)],
###       )
###     )
###   ) 
### points(locs[,1], M[,1], pch=3)


#' Get prediction points (PP) for useful for vecchia 
#'
#' @param observed_locs a matrix of observed locations. TODO
#' @param matern_range a numeric vector. TODO
#' @param lonlat logical, default to FALSE. TODO
#' @param n_PP a numerical value. Number of prediction points, default to 20.
#' @param m  a numerical value. TODO
#' @param seed a numerical value. The random seed, default to 123. Uses mainly to reproduce results in tests. 
#'
#' @returns a list
#' @export
#'
#' @examples
#' obs_locs <- matrix(rnorm(20), ncol=2)
#' get_PP(obs_locs, matern_range=c(1, 1.1, 1.5, 0), n_PP=4)
get_PP = function(observed_locs, matern_range, lonlat = F, n_PP = 20, m = 10, seed=1234)
{
  # Suppress duplicates in observed points.
  locs_ = observed_locs[! duplicated(observed_locs),]
  # Random sample. 
  locs_ = locs_[order(runif(nrow(locs_))),]
  # Reordonate first 100000 points to help spatial coverage of vecchia. 
  locs_[seq(min(nrow(locs_), 100000)),] <- locs_[GpGp::order_maxmin(locs_[seq(min(nrow(locs_), 100000)),]),]
  # Corresponding locations
  idx = match(split(observed_locs, row(observed_locs)), split(locs_, row(locs_)))
  
  knots = kmeans(locs_[seq(min(nrow(locs_), 100000)),], 
                 n_PP, algorithm = "Hartigan-Wong", iter.max = 50)$centers
  knots = knots[GpGp::order_maxmin(knots),]
  
  # Get the m closed neighbors (ordered)
  NNarray = GpGp::find_ordered_nn(rbind(knots, locs_), m, lonlat = lonlat)
  # Lower triangular matrix of inversed cholesky factor 
  # Add a seed here to always get the same results from vecchia
  set.seed(seed)
  sparse_chol = Matrix::sparseMatrix(
    i = row(NNarray)[!is.na(NNarray)], 
    j = NNarray[!is.na(NNarray)], 
    x = GpGp::vecchia_Linv(covparms = c(1, matern_range, .0001), covfun_name = "matern15_isotropic", locs = rbind(knots, locs_), NNarray = NNarray)[!is.na(NNarray)], 
    triangular = T
  )
  return(
    list(
      "knots" = knots,
      "unique_reordered_locs" = locs_,
      "idx" = idx,
      "lonlat" = lonlat,
      "m" = m,
      "matern_range" = matern_range,
      "sparse_chol" = sparse_chol,
      "NNarray" = NNarray,
      "n_PP" = n_PP
    )
  )
}

#' Compute prior logarithmic density of a matrix
#' TODO A reformater.
#'
#' @param beta TODO
#' @param n_PP TODO
#' @param beta_mean TODO
#' @param beta_precision TODO
#' @param log_scale TODO
#'
#' @returns a numeric value
#' @export
#'
#' @examples
#' beta = matrix(rnorm(10))
#' bm = matrix(rnorm(5), 5, 1)
#' bp = diag(exp(rnorm(5)), 5, 5)
#' ls = rnorm(1)
#' # TODO a changer beta_prior_log_dens(beta, n_PP = 5, beta_mean = bm, beta_precision = bp, log_scale = ls)
beta_prior_log_dens = function(beta, 
                               n_PP, 
                               beta0_mean,
                               beta0_var,
                               chol_crossprod_X, 
                               log_scale){
  PP_prior = 0
  if(n_PP>0) 
  {
    scale_mat = GeoNonStat::expmat(-log_scale)
    PP_prior = (
      # PP coefficients follow N(0, scale_mat)
      +.5 * n_PP * determinant(scale_mat, logarithm = T)$mod # determinant is changed by log scale
      -sum(.5 * c(beta[-seq(nrow(beta)-n_PP),,drop = F] %*% scale_mat) * beta[-seq(nrow(beta)-n_PP),,drop = F])
    )
  }
  chol_crossprod_X_ = chol_crossprod_X[seq(nrow(beta)-n_PP),seq(nrow(beta)-n_PP),drop = F]
  chol_crossprod_X_[1,] = chol_crossprod_X_[1,]/chol_crossprod_X_[1,1]
  M = solve(chol_crossprod_X_)
  
  return(
    PP_prior + 
      sum(
        -.5 * (M %*% beta[seq(nrow(beta)-n_PP),] - c(beta0_mean, rep(0, length(beta[seq(nrow(beta)-n_PP),])-1)))^2/c(beta0_var, rep(10, length(beta[seq(nrow(beta)-n_PP),])-1))
      )  
  )
}



#' Compute gradient of logarithmic density prior
#'
#' @param beta a numerical matrix of spatial coordinates. 
#' @param n_PP number of prediction points.
#' @param beta_mean TODO
#' @param beta_precision TODO
#' @param log_scale TODO
#'
#' @returns an array
#' @export
#'
#' @examples
#' beta = matrix(rnorm(10))
#' bm = matrix(rnorm(5), 5, 1)
#' bp = diag(exp(rnorm(5)), 5, 5)
#' ls = rnorm(1)
#' # TODO q corriger beta_prior_log_dens_derivative(beta, n_PP = 5, beta_mean = bm, beta_precision = bp, log_scale = ls)
beta_prior_log_dens_derivative = 
  function(beta, n_PP, 
           beta0_mean,
           beta0_var,
           chol_crossprod_X, 
           log_scale){
    chol_crossprod_X_ = chol_crossprod_X[seq(nrow(beta)-n_PP),seq(nrow(beta)-n_PP),drop = F]
    chol_crossprod_X_[1,] = chol_crossprod_X_[1,]/chol_crossprod_X_[1,1]
    M = solve(chol_crossprod_X_)
    res = t(M) %*% matrix(
      -(M%*%beta[seq(nrow(beta)-n_PP),,drop=F] - c(beta0_mean, rep(0, length(beta[seq(nrow(beta)-n_PP),,drop=F])-1)))/c(beta0_var, rep(10, length(beta[seq(nrow(beta)-n_PP),,drop=F])-1)),
      ncol = ncol(beta)
    ) 
    if(n_PP>0) 
    {
      scale_mat = GeoNonStat::expmat(-log_scale)
      res = rbind(res, 
                  -beta[-seq(nrow(beta)-n_PP),,drop = F] %*% scale_mat
      )
    }
    res
}



## beta = matrix(rnorm(10))
## n_PP = 5
## beta_mean = matrix(rnorm(5), 5, 1)
## beta_precision = diag(exp(rnorm(5)), 5, 5)
## log_scale = rnorm(1)
## 
## beta_prior_log_dens(beta = beta, n_PP = 5, beta_mean = beta_mean, beta_precision = beta_precision, log_scale = log_scale)
## beta_ = beta; beta_[10]=beta_[10]+.0001
## beta_prior_log_dens_derivative(beta = beta, n_PP = 5, beta_mean = beta_mean, beta_precision = beta_precision, log_scale = log_scale)
## 10000*
##   (
##   beta_prior_log_dens(beta = beta_, n_PP = 5, beta_mean = beta_mean, beta_precision = beta_precision, log_scale = log_scale)-
##   beta_prior_log_dens(beta = beta, n_PP = 5, beta_mean = beta_mean, beta_precision = beta_precision, log_scale = log_scale)
##   )


#PP$idx : match between the non redundant locations of PP and the redundant observed locations
#locs_idx : match between the redundant observed locations and those of X

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

#' Title
#'
#' @param coords a numeric vector of length 1 or 6
#' @param eps a numeric value, default to 0.00001
#'
#' @returns an 3 dimensional array
#' @export
#'
#' @examples
#' derivative_chol_expmat(c(1,2,3, 3,2,4))
derivative_chol_expmat = function(coords, eps=0.00001)
{
  dimres = 1
  if(length(coords)==6) dimres = 3
  res = array(data = 0, dim = c(dimres, dimres, length(coords)))
  chol_expmat = chol(expmat(coords))
  for(i in seq(length(coords)))
  {
    coords_ = coords
    coords_[i] = coords_[i] + eps
    res[,,i] = 100000 * (chol(expmat(coords_)) - chol_expmat)
  }
  res
}

#' Title TODO
#'
#' @param field TODO
#' @param coordsTODO
#'
#' @returns an array
#' @export
#'
#' @examples
#' field <- array(c(1,2,3), dim=c(1,3))
#' res <- derivative_field_wrt_scale(field, c(1,2,3,3,2,4))
derivative_field_wrt_scale = function(field, coords)
{
  d_chol_expmat = derivative_chol_expmat(coords)
  white_field = field %*% solve(chol(expmat(coords)))
  res = array(0, dim = c(dim(field), length(coords)))
  for(i in seq(length(coords)))
  {
    res[,,i] = white_field %*% d_chol_expmat[,,i]
  }
  res
}
