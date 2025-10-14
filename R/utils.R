#' Title
#'
#' @param M a SparseMatrix
#'
#' @returns a vector of colors, of length the number of rows of M
#' @examples
#' n <- 5  
#' i <- c(rep(1, n), 2:n)       # 1re ligne + 1re colonne sauf la 1re case
#' j <- c(1:n, rep(1, n - 1))
#' M <- sparseMatrix(i = i, j = j, x = 1, dims = c(n, n))
#' naive_greedy_coloring(M)
#' M[2,3]=1
#' M[3,2]=1
#' naive_greedy_coloring(M)
naive_greedy_coloring = function(M)
{
  #number of nodes
  if(!Matrix::isSymmetric(M, checkDN=FALSE)) 
    stop("M must be symmetric")
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

#' Title
#'
#' @param vecchia_approx  an object created with `createVecchia()`
#' @param compressed_sparse_chol an object created with `compute_sparse_chol()`
#'
#' @returns a sparse triangular matrix
decompress_chol = function(vecchia_approx, compressed_sparse_chol){
  Matrix::sparseMatrix(
    i = vecchia_approx$sparse_chol_i, 
    p = vecchia_approx$sparse_chol_p, 
    x = compressed_sparse_chol[,1,][vecchia_approx$sparse_chol_x_reorder], 
    triangular = T
  ) 
}

#' Computes a Vecchia sparse Cholesky factor and its derivatives
#' 
#' @param range_beta parameter for the range.
#' If the covariance is anisotropic, it must have 3 columns. It the covariance is isotropic, it must have 1 column. 
#' The first coefficients are multiplied with range_X 
#' The last coefficients are multiplied with the spatial basis functions of PP
#' @param vecchia_approx an object created with `createVecchia()`
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
#' compute_sparse_chol(range_beta = range_beta, 
#' vecchia_approx = vecchia_approx, 
#' range_X = range_X, 
#' PP = PP, 
#' matern_smoothness = 1.5, compute_derivative = T)
compute_sparse_chol = function(range_beta, 
                               vecchia_approx,
                               range_X, 
                               PP = NULL, 
                               matern_smoothness = 1.5, 
                               compute_derivative = T, 
                               num_threads = 1)
{
  if (!matern_smoothness %in% c(.5, 1.5)) stop("matern_smoothness must be equal to 0.5 or 1.5")
  if(ncol(range_beta)==3) {
    range_beta = range_beta %*% matrix(
      c(1/sqrt(2), 1/sqrt(2),  0, 
        1/sqrt(2), -1/sqrt(2), 0,
        0,       0,        1), 3)*sqrt(2)*2
  } else if(ncol(range_beta)==1) {
    range_beta = range_beta * 2
  } else {
    stop("range_beta is expected to have 1 (isotropic case) or 3 (anisotropic case) columns")
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
###         x= compute_sparse_chol(
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
###         x= compute_sparse_chol(
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
#'
#' @param beta a numeric matrix, with 1 or 3 columns
#' @param n_PP number of PP
#' @param beta0_mean mean of beta, numeric value
#' @param beta0_var var of beta, numeric value
#' @param log_scale numeric vector of length 1 for a 1 column beta matrix
#' and length 2 for a 1 column beta matrix
#' @description
#' #'   beta follows a Normal distribution, with independent components. 
#'   
#'   -------------------------------------------------------------------------------                    
#'   In the isotropic case, beta has one column. 
#'   It has 3 categories of rows : 
#'   - the first row, corresponding to the Intercept. 
#'   - the next rows, corresponding to the rest of the explnanatory variables if there are any.
#'   - the last rows, corresponding to the PP coefficients if there are any.
#'   
#'                 The mean is a vector     The diagonal of the variance matrix is 
#'                      with shape :                    a vector with shape : 
#'                                                                                                               
#'    (Intercept)      beta_0_mean                           beta_0_var                                              
#'                 |-       0                                   0.01                                                
#'                 |        0                                   0.01                                                
#'            X   -|        0                                   0.01                                                
#'                 |        0                                   0.01                                                
#'                 |_       0                                   0.01                                                
#'                 |-       0                          exp(range_log_scale)                                    
#'                 |        0                          exp(range_log_scale)                                    
#'            PP  -|        0                          exp(range_log_scale)                                    
#'                 |        0                          exp(range_log_scale)                                    
#'                 |        0                          exp(range_log_scale)                                    
#'                 |_       0                          exp(range_log_scale)                                    
#'                    
#'   -------------------------------------------------------------------------------                    
#'                                                                                                            
#'   In the anisotropic case, beta has 3 columns, one for the Range, two for the Anisotropy 
#'   The categories of rows are the same as in the Isotropic case. 
#'   
#'   The mean and the diagonal of the variance matrix are stored under 
#'   the shape of matrices with 3 columns, like beta.
#'                                                                                                            
#'                                     mean                                                        variance                                                                                        
#'                                                                                                                          
#'                      (range)      (aniso)    (aniso)                                                                  
#'                                                                                                                                                   
#'    (Intercept)      beta_0_mean      0           0                     beta_0_var                  0.01                      0.01
#'                 |-       0           0           0                        0.01                     0.01                      0.01                                                                                      
#'                 |        0           0           0                        0.01                     0.01                      0.01                                                                                      
#'            X   -|        0           0           0                        0.01                     0.01                      0.01                                                                                      
#'                 |        0           0           0                        0.01                     0.01                      0.01                                                                                      
#'                 |_       0           0           0                        0.01                     0.01                      0.01                                                                                      
#'                 |-       0           0           0               exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])                                                                          
#'                 |        0           0           0               exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])                                                                          
#'            PP  -|        0           0           0               exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])                                                                          
#'                 |        0           0           0               exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])                                                                          
#'                 |        0           0           0               exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])                                                                          
#'                 |_       0           0           0               exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])                                             
#'                 
#' @returns a numeric value
#' @export
#'
#' @examples
#' beta = matrix(rnorm(300), 100, ncol=3)
#' n_PP = 90 
#' beta0_mean = -5
#' beta0_var = 2
#' log_scale = c(-3,-2)
#' 
#' beta_prior_log_dens(
#'   beta, 
#'   n_PP, 
#'   beta0_mean,
#'   beta0_var,
#'   log_scale
#' )
beta_prior_log_dens = function(beta, 
                               n_PP, 
                               beta0_mean,
                               beta0_var,
                               log_scale){

  nrb <- nrow(beta)
  ncb <- ncol(beta)
  if(!ncb %in% c(1,3)) stop("beta is expected to have 1 or 3 columns")
  if(n_PP > nrb + 2) {
    stop("n_PP can't be greater than nrow(beta) + 2")
  }
  mean_mat <- matrix(0, nrb, ncb) 
  var_mat <- matrix(0.01, nrb, ncb)  
  mean_mat[1,1] <- beta0_mean  # Intercept for range
  var_mat[1,1] <- beta0_var # Intercept for range
  var_mat[-seq(nrb-n_PP), 1] <- exp(log_scale[1])
  if(ncb==3) {
    if(!length(log_scale) == 2) stop("log_scale is supposed to be of length 2")
    var_mat[-seq(nrb-n_PP), c(2,3)] <- exp(log_scale[2])
  }
  return(-0.5 * sum((beta - mean_mat)^2 / var_mat))
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
#' beta_prior_log_dens_derivative(
#'   beta = matrix(rnorm(300), 100, 3), 
#'   n_PP = 90, 
#'   beta0_mean = -5,
#'   beta0_var = 2,
#'   log_scale = c(-3, -2)
#' )
<<<<<<< HEAD
#' 
#' beta = matrix(rnorm(300), 100)
#' n_PP = 90 
#' beta0_mean = -5
#' beta0_var = 2
#' log_scale = c(-3, -5)
#' 
#' beta_prior_log_dens_derivative(
#'   beta, 
#'   n_PP, 
#'   beta0_mean,
#'   beta0_var,
#'   log_scale
#' )
 # TODO Elise : ceci n'est pas un exemple.
 # Par contre passer ça en test pour vérifier 
 # que ça a bonne taille et que ça tourne autour de 1. 
 # TODO : rajouter un exemple pour dire que ça tourne avec beta 1 colonne
 # TODO : rajouter un exemple pour dire que ça tourne et avec n_PP=0
#'res <-matrix(0, 100, 3)
#'for(i in seq(100))
#'{
#'  for(j in seq(3)){
#'    beta1 = beta 
#'    beta1[i,j] = beta1[i,j]+ .000001
#'    res[i,j] <- (
#'      beta_prior_log_dens(beta1, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = c(0, 2)) -
#'        beta_prior_log_dens(beta, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = c(0, 2))
#'    )*1000000
#' }
#'}
#'nalytical_der =  beta_prior_log_dens_derivative(beta, n_PP = 90, beta0_mean = -4, beta0_var = 2, log_scale = c(0, 2))
#'ar(mfrow = c(1,3))
#'oxplot(res- analytical_der)
#'oxplot(res/analytical_der)
#'lot(res,analytical_der)
#'bline(a=0, b=1)
 
=======
>>>>>>> 8fd9265 (tests and sanity checks on beta_prior_log_dens)
beta_prior_log_dens_derivative = 
  function(beta, 
           n_PP, 
           beta0_mean,
           beta0_var,
           log_scale){
<<<<<<< HEAD
    mean_mat <- matrix(0, dim(beta)[1], dim(beta)[2]) 
    var_mat <- matrix(0.1,dim(beta)[1], dim(beta)[2]) 
=======
    nrb <- nrow(beta)
    ncb <- ncol(beta)
    if(!ncb %in% c(1,3)) stop("beta is expected to have 1 or 3 columns")
    if(n_PP > nrb + 2) {
      stop("n_PP can't be greater than nrow(beta) + 2")
    }
    mean_mat <- matrix(0, nrb, ncb) 
    var_mat <- matrix(0.01, nrb, ncb) 
>>>>>>> 8fd9265 (tests and sanity checks on beta_prior_log_dens)
    mean_mat[1,1] <- beta0_mean
    var_mat[1,1] <- beta0_var
    var_mat[-seq(nrb-n_PP), 1] <- exp(log_scale[1])
    if(ncb==3) {
      if(!length(log_scale) == 2) stop("log_scale is supposed to be of length 2")
      var_mat[-seq(nrb-n_PP), c(2,3)]<- exp(log_scale[2])
    }
    return(-(beta - mean_mat)/ var_mat)
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
    # using res[] is important for performance
    res[] = res + X  %*% Y[seq_len(xrow_offset), , drop=FALSE]
  } 
  if(!is.null(PP)) {
    # remove X rows from Y if needed
    if(xrow_offset>0) Y =  Y[-seq_len(xrow_offset), , drop = FALSE] 
    V = matrix(0, nrow(PP$sparse_chol), ncol(Y))
    V[seq(nrow(Y)),] = Y
    solved <- Matrix::solve(PP$sparse_chol, V, triangular = TRUE)
    PP_result <- solved[-seq_len(nrow(PP$knots)), , drop = FALSE]
    # using res[] is important for performance
    res[] <- res + PP_result[locs_idx, , drop = FALSE]
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


