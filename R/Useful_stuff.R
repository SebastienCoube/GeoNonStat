

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



#' Computes a Vecchia sparse Cholesky factor and its derivatives
#' 
#' @param range_beta parameter for the range.
#' If the covariance is anisotropic, it must have 3 columns. It the covariance is isotropic, it must have 1 column. 
#' The first coefficients are multiplied with range_X 
#' The last coefficients are multiplied with the spatial basis functions of PP
#' @param vecchia_approx TODO
#' @param range_X covariates for range
#' @param PP predictive process obtained through `createPP()`
#' @param compute_derivative logical, indicates if derivatives of Vecchia factors are to be computed
#' @param smoothness Matern smoothness Default to 1.5. Can be 0.5 or 1.5.
#' @param anisotropic Logical, default to FALSE. TODO
#' @param num_threads numerical, number of treads to use. Default to 1.
#' @param locs_idx match between the duplicated locations used to buile the PP basis function and the non-redundant locs used to compute the sparse chol
#'
#' @returns a list
#' @export
#'
#' @examples
#' locs = cbind(runif(1000), runif(1000))
#' vecchia_approx = createVecchia(locs)
#' Rcpp::sourceCpp("src/vecchia.cpp")
#' range_X = matrix(1, nrow(locs))
#' PP = createPP(vecchia_approx)
#' 
#' 
#' 
#' # test equivalence of parametrizations for locally isotropic Matérn covariance, smoothness= 1.5
#' range_beta = matrix(c(-2))
#' GpGpcov = tcrossprod(solve(
#'   Matrix::sparseMatrix(
#'     i = row(t(vecchia_approx$NNarray))[!is.na(t(vecchia_approx$NNarray))], 
#'     j = t(vecchia_approx$NNarray)[!is.na(t(vecchia_approx$NNarray))], 
#'     x=  
#'   GpGp::vecchia_Linv(
#'   c(1, exp(range_beta[1]), 1.5, .000001), covfun_name = "matern_isotropic", 
#'   (vecchia_approx$locs), NNarray = t(vecchia_approx$NNarray))[!is.na(t(vecchia_approx$NNarray))], 
#'   triangular = T
#'   )
#'   ))
#' image(GpGpcov[vecchia_approx$locs_match, vecchia_approx$locs_match])
#' mycov = tcrossprod(solve(
#'   Matrix::sparseMatrix(
#'     i = col((vecchia_approx$NNarray))[!is.na((vecchia_approx$NNarray))], 
#'     j = (vecchia_approx$NNarray)[!is.na((vecchia_approx$NNarray))], 
#'     x=  
#'       compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = NULL, matern_smoothness = 1.5, compute_derivative = F)[!is.na((vecchia_approx$NNarray))], 
#'     triangular = T
#'   )
#' ))
#' image(mycov[vecchia_approx$locs_match, vecchia_approx$locs_match])
#' hist(mycov- GpGpcov)
#' # test equivalence of parametrizations for locally isotropic Matérn covariance, smoothness = 0.5 aka exponential
#' GpGpcov = tcrossprod(solve(
#'   Matrix::sparseMatrix(
#'     i = row(t(vecchia_approx$NNarray))[!is.na(t(vecchia_approx$NNarray))], 
#'     j = t(vecchia_approx$NNarray)[!is.na(t(vecchia_approx$NNarray))], 
#'     x=  
#'   GpGp::vecchia_Linv(
#'   c(1, exp(range_beta[1]), .5, .000001), covfun_name = "matern_isotropic", 
#'   (vecchia_approx$locs), NNarray = t(vecchia_approx$NNarray))[!is.na(t(vecchia_approx$NNarray))], 
#'   triangular = T
#'   )
#'   ))
#' image(GpGpcov[vecchia_approx$locs_match, vecchia_approx$locs_match])
#' mycov = tcrossprod(solve(
#'   Matrix::sparseMatrix(
#'     i = col((vecchia_approx$NNarray))[!is.na((vecchia_approx$NNarray))], 
#'     j = (vecchia_approx$NNarray)[!is.na((vecchia_approx$NNarray))], 
#'     x=  
#'       compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = NULL, matern_smoothness = .5, compute_derivative = F)[!is.na((vecchia_approx$NNarray))], 
#'     triangular = T
#'   )
#' ))
#' image(mycov[vecchia_approx$locs_match, vecchia_approx$locs_match])
#' hist(mycov- GpGpcov)
#' range_beta = matrix(rnorm(1 + PP$n_knots))
#' compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = PP, matern_smoothness = 1.5, compute_derivative = F)
#' range_beta = matrix(rnorm(1 + PP$n_knots))
#' compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = PP, matern_smoothness = 1.5, compute_derivative = F)
#' 
#' 
#' # test equivalence of parametrizations for locally isotropic and locally aniso
#' range_beta = matrix(c(-2))
#' mycov = tcrossprod(solve(
#'   Matrix::sparseMatrix(
#'     i = col((vecchia_approx$NNarray))[!is.na((vecchia_approx$NNarray))], 
#'     j = (vecchia_approx$NNarray)[!is.na((vecchia_approx$NNarray))], 
#'     x=  
#'       compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = NULL, matern_smoothness = 1.5, compute_derivative = F)[!is.na((vecchia_approx$NNarray))], 
#'     triangular = T
#'   )
#' ))
#' image(mycov[vecchia_approx$locs_match, vecchia_approx$locs_match])
#' 
#' range_beta = matrix(c(-2, 0,0), 1)
#' mycov_aniso = tcrossprod(solve(
#'   Matrix::sparseMatrix(
#'     i = col((vecchia_approx$NNarray))[!is.na((vecchia_approx$NNarray))], 
#'     j = (vecchia_approx$NNarray)[!is.na((vecchia_approx$NNarray))], 
#'     x=  
#'       compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = NULL, matern_smoothness = 1.5, compute_derivative = F)[!is.na((vecchia_approx$NNarray))], 
#'     triangular = T
#'   )
#' ))
#' image(mycov_aniso[vecchia_approx$locs_match, vecchia_approx$locs_match])
#' hist(mycov - mycov_aniso)
#' # test with a PP
#' range_beta = matrix(rnorm(1 + PP$n_knots))
#' compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = PP, matern_smoothness = 1.5, compute_derivative = T)
#' range_beta = matrix(rnorm(3*(1 + PP$n_knots)), ncol = 3)
#' compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = PP, matern_smoothness = 1.5, compute_derivative = T)

compute_sparse_chol = function(range_beta, 
                               vecchia_approx,
                               range_X = NULL, 
                               PP = NULL, 
                               matern_smoothness = 1.5, 
                               compute_derivative = T, 
                               num_threads = 1)
{
  if (!matern_smoothness %in% c(.5, 1.5)) stop("matern_smoothness must be equal to 0.5 or 1.5")
  # converting to canonical basis
  if(ncol(range_beta)==3) {
    range_beta = range_beta %*% matrix(
    c(1/sqrt(2), 1/sqrt(2),  0, 
      1/sqrt(2), -1/sqrt(2), 0,
      0,       0,        1), 3)*sqrt(2)*2
  }
  if(ncol(range_beta)==1) {
    range_beta = range_beta * 2
  }
  
  log_range = as.matrix(
    X_PP_mult_right(
      vecchia_approx=  vecchia_approx, 
      X = range_X, 
      PP = PP, 
      Y = range_beta, 
      permutate_PP_to_obs = F))

  res <- vecchia(num_threads=num_threads,
                 log_range = t(log_range), 
                 locs = vecchia_approx$t_locs, 
                 NNarray = vecchia_approx$NNarray, 
                 compute_derivative = compute_derivative,
                 smoothness = matern_smoothness)
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
###           range_X = matrix(1, nrow(locs), 1), smoothness = 1.5
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
###           range_X = matrix(1, nrow(locs), 1), smoothness = 1.5
###         )[[1]][!is.na(NNarray)],
###       )
###     )
###   ) 
### points(locs[,1], M[,1], pch=3)


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