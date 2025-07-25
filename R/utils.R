naive_greedy_coloring = function(M)
{
  #number of nodes
  n_obs = nrow(M)
  #deducting degrees
  degrees = as.vector(rep(1, n_obs)%*%M)
  #getting adjacent nodes of a given node
  idx = split(M@i+1, rep(seq_along(diff(M@p)),diff(M@p)))
  #creating a color * node matrix of incompatibilities
  incompatibilities = matrix(0, n_obs+1, max(degrees))
  cols = rep(0, n_obs)
  
  for(i in seq(n_obs))
  {
    cols[i] = match(0, incompatibilities[i,])
    incompatibilities[idx[[i]],cols[i]] = 1
  }
  return(cols)
}


decompress_chol = function(vecchia_approx, compressed_sparse_chol){
  Matrix::sparseMatrix(
    i = vecchia_approx$sparse_chol_i, 
    p = vecchia_approx$sparse_chol_p, 
    x = compressed_sparse_chol[,1,][vecchia_approx$sparse_chol_x_reorder], 
    triangular = T
  ) 
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



#' Computes a Vecchia sparse Cholesky factor and its derivatives
#' 
#' @param range_beta parameter for the range.
#' If the covariance is anisotropic, it must have 3 columns. It the covariance is isotropic, it must have 1 column. 
#' The first coefficients are multiplied with range_X 
#' The last coefficients are multiplied with the spatial basis functions of PP
#' @param vecchia_approx TODO
#' @param range_X covariates for range, treated using process_covariates
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
                               range_X, 
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
      X = range_X$X_locs, 
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

#### checking equivalence of parametrizations ##
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
#' locs = cbind(runif(10000), runif(10000))
#' X_range = cbind(locs, locs^2)
#' V = createVecchia(locs)
#' PP = createPP(V)
#' GNS = GeoNonStat(vecchia_approx = V, observed_field = rnorm(nrow(locs)), range_PP = PP, range_X = as.data.frame(X_range))
#' 
#' GNS$covariates$range_X$crossprod_X_locs[seq()]
#' GNS$hierarchical_model$range_beta0_mean 
# 
beta_prior_log_dens = function(beta, 
                               X, 
                               beta0_mean,
                               beta0_var,
                               log_scale){
  PP_prior = 0
  if(n_PP>0) 
  {
    scale_mat = expmat(-log_scale)
    PP_prior = (
      # PP coefficients follow N(0, scale_mat)
      +.5 * n_PP * determinant(scale_mat, logarithm = T)$mod # determinant is changed by log scale
      -sum(.5 * c(beta[-seq(nrow(beta)-n_PP),,drop = F] %*% scale_mat) * beta[-seq(nrow(beta)-n_PP),,drop = F])
    )
  }
  
  return(
    PP_prior +  sum(
      c(
      - .5 * (beta[1,1] - beta0_mean)^2 / beta0_var,
      - .5 * (beta[seq(nrow(beta) - n_PP),][-1])^2 * .01
      )
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
#'
#' @examples
 beta = matrix(rnorm(300), 100)
 beta1 = beta 
 #derived_idx=  cbind(1,1)
 #derived_idx=  cbind(1,2)
 #derived_idx=  cbind(2,1)
 # derived_idx=  cbind(10,3)
 # beta1[derived_idx] = beta1[derived_idx]+ .0001
 # (
 #  beta_prior_log_dens(beta1, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = rep(0, 6)) -
 #   beta_prior_log_dens(beta, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = rep(0, 6))
 # )*10000
 # beta_prior_log_dens_derivative(beta, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = rep(0, 6))[derived_idx]
 

 
 beta_prior_log_dens_derivative = 
  function(beta, n_PP, 
           beta0_mean,
           beta0_var,
           log_scale){
    res = beta[seq(nrow(beta)-n_PP),]
    res[1] = -(res[1] - beta0_mean)/beta0_var
    res[-1] = -(res[-1])/100
    if(n_PP>0) 
    {
      scale_mat = GeoNonStat::expmat(-log_scale)
      res = rbind(res, 
                  beta[-seq(nrow(beta)-n_PP),,drop = F] %*% scale_mat
      )
    }
    res
}

#PP$idx : match between the non redundant locations of PP and the redundant observed locations
#locs_idx : match between the redundant observed locations and those of X


#' Title
#'
#' @param coords a numeric vector of length 1 or 6
#' @param eps a numeric value, default to 0.00001
#'
#' @returns an 3 dimensional array
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




#' Multiply the concatenation of a matrix of covariates and a PP by a matrix
#' (X|PP) %*% Y
#'
#' @param X a matrix of covariates who will be multiplied by the first rows of Y
#' @param PP either a PP whose basis will be multiplied by the last columns of Y
#' @param locs_idx either a vector of integers who dispatch the PP basis to the covariates, or NULL
#' @param Y the matrix who multiplies the covariates and the PP
#'
#' @returns a matrix
#' 
#' @examples
#' locs = cbind(runif(1000), runif(1000))
#' locs = rbind(locs, locs)
#' vecchia_approx = createVecchia(locs, 12, ncores=1)
#' PP = createPP(vecchia_approx, plot=FALSE)
#' covariate_coefficients = c(4, 1, 1, .5)
#' knots_coeffs = rnorm(PP$n_knots)
#' X = cbind(1, vecchia_approx$locs, rnorm(nrow(vecchia_approx$locs)))
#' 
#' # multiplying X alone
#' res1 <- X_PP_mult_right(X = X, Y = covariate_coefficients, vecchia_approx = vecchia_approx)
#' plot_pointillist_painting(locs, res1, main = "covariates only, \n one covariate for each observation")
#' res1 <- X_PP_mult_right(X = X, Y = covariate_coefficients, vecchia_approx = vecchia_approx)
#' plot_pointillist_painting(locs, res1, main = "covariates only, \n one covariate for each observation")
#' 
#' # multiplying PP alone
#' res2 <- X_PP_mult_right(PP = PP, Y = knots_coeffs, vecchia_approx = vecchia_approx)
#' plot_pointillist_painting(vecchia_approx$locs, res2, main = "PP only")
#' 
#' # multiplying PP and matrix of covariate, one obs for each location
#' X_by_loc = cbind(1, vecchia_approx$locs, rnorm(vecchia_approx$n_locs))
#' res3 <- X_PP_mult_right(PP = PP, X = X_by_loc, Y = c(covariate_coefficients, knots_coeffs), vecchia_approx = vecchia_approx)
#' plot_pointillist_painting(vecchia_approx$locs, res3, main = "PP + covariates, \n one covariate for each location")
#' 
#' # multiplying PP and matrix of covariates with an index
#' X_by_obs = cbind(1, vecchia_approx$observed_locs, rnorm(vecchia_approx$n_obs))
#' res4 <- X_PP_mult_right(X = X_by_obs, PP = PP, 
#'                        Y = c(covariate_coefficients, knots_coeffs), 
#'                        vecchia_approx = vecchia_approx
#'                        )
#' plot_pointillist_painting(vecchia_approx$observed_locs, res4, main = "PP + covariates,\n  one covariate for each observation")
#' 
#' # multiplying 
#' X_by_obs = cbind(1, vecchia_approx$observed_locs, rnorm(vecchia_approx$n_obs))
#' res5 <- X_PP_mult_right(X = NULL, PP = PP, 
#'                        Y = diag(1, PP$n_knots), 
#'                        vecchia_approx = vecchia_approx,
#'                        permutate_PP_to_obs = T
#'                        )
X_PP_mult_right = function(X = NULL, PP = NULL, vecchia_approx, Y, permutate_PP_to_obs = F)
{
  if(is.null(X) & is.null(PP)) stop("X and PP can't be both NULL")
  # Sanity checks
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  expected_rows <- 0
  if(!is.null(X)) expected_rows <- expected_rows + ncol(X)
  if(!is.null(PP)) expected_rows <- expected_rows + nrow(PP$knots)
  if(nrow(Y) != expected_rows) {
    stop("Y should have ", expected_rows, " rows it has ", nrow(Y))
  }
  
  if(permutate_PP_to_obs) {
    locs_idx = vecchia_approx$locs_match
  } else {
    locs_idx = seq(vecchia_approx$n_locs)
  }
  res = matrix(0, length(locs_idx), ncol(Y))
  
  # Multiply X and Y
  xrow_offset <- 0
  if(!is.null(X)) {
    xrow_offset <- ncol(X)
    res = res + X  %*% Y[seq_len(xrow_offset), , drop=FALSE]
  } 
  if(!is.null(PP)) {
    # remove X rows from Y if needed
    if(xrow_offset>0) Y =  Y[-seq_len(xrow_offset), , drop = FALSE] 
    V = matrix(0, nrow(PP$sparse_chol), ncol(Y))
    V[seq(nrow(Y)),] = Y
    solved <- Matrix::solve(PP$sparse_chol, V, triangular = TRUE)
    PP_result <- solved[-seq_len(nrow(PP$knots)), , drop = FALSE]
    res <- res + PP_result[locs_idx, , drop = FALSE]
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
#'
#' @examples
#' set.seed(123)
#' locs = cbind(runif(50), runif(50))
#' vecchia_approx = createVecchia(locs, ncores=1)
#' PP = createPP(vecchia_approx, plot=FALSE)
#' X = matrix(rnorm(100), 50)
#' Y = matrix(rnorm(30*nrow(X)), nrow(X))
#' 
#' # just surrogate of crossprod
#' res1 <- X_PP_crossprod(X = X, Y = Y)
#' identical(crossprod(X, Y) , res1)
#' 
#' # crossprod + PP with observations of X on the locs
#' res2 <- X_PP_crossprod(X = X, PP = PP, Y = Y, vecchia_approx = vecchia_approx, permutate_PP_to_obs = F)
#' 
#' # crossprod + PP with observations of X on the obs
#' res3 <- X_PP_crossprod(X = X, PP = PP, Y = Y, vecchia_approx = vecchia_approx, permutate_PP_to_obs = T)
X_PP_crossprod = function(X, PP = NULL, Y, vecchia_approx=NULL, permutate_PP_to_obs = F)
{
  # TODO : est-ce qu'il ne faudrait pas déplacer cette fonction dans un fichier utils ou usefull_stuff ?
  # Si on le fait le faut déplacer les tests avec. 
  if(nrow(X) != nrow(Y)) {
    stop("X and Y should have the same number of rows")
  }
  if(permutate_PP_to_obs & is.null(vecchia_approx)) {
    stop("To permutate PP to observed values vecchia_approx needs to be provided.")
  }
  if(!is.null(PP)){
    if(nrow(X) != nrow(PP$vecchia_locs) ) {
      stop("X should have the same number of rows as locations in vecchia")
    }
    if(vecchia_approx$n_locs != nrow(PP$vecchia_locs) ) {
      stop("vecchia_approx should have the same number of locations as the locations of PP")
    }
  }
  if(!is.matrix(Y)) Y = as.matrix(Y)
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
  if(!is.matrix(res)) res <- as.matrix(res)
  return(res)
}


