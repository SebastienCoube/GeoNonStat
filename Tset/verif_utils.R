### compute_sparse_chol #######################################################################
# test equivalence of parametrizations for locally isotropic Matérn covariance, smoothness= 1.5
range_beta = matrix(c(-2))
GpGpcov = tcrossprod(solve(
  Matrix::sparseMatrix(
    i = row(t(vecchia_approx$NNarray))[!is.na(t(vecchia_approx$NNarray))], 
    j = t(vecchia_approx$NNarray)[!is.na(t(vecchia_approx$NNarray))], 
    x=  
  GpGp::vecchia_Linv(
  c(1, exp(range_beta[1]), 1.5, .000001), covfun_name = "matern_isotropic", 
  (vecchia_approx$locs), NNarray = t(vecchia_approx$NNarray))[!is.na(t(vecchia_approx$NNarray))], 
  triangular = T
  )
  ))
image(GpGpcov[vecchia_approx$locs_match, vecchia_approx$locs_match])
mycov = tcrossprod(solve(
  Matrix::sparseMatrix(
    i = col((vecchia_approx$NNarray))[!is.na((vecchia_approx$NNarray))], 
    j = (vecchia_approx$NNarray)[!is.na((vecchia_approx$NNarray))], 
    x=  
      compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = NULL, matern_smoothness = 1.5, compute_derivative = F)[!is.na((vecchia_approx$NNarray))], 
    triangular = T
  )
))
image(mycov[vecchia_approx$locs_match, vecchia_approx$locs_match])
hist(mycov- GpGpcov)
# test equivalence of parametrizations for locally isotropic Matérn covariance, smoothness = 0.5 aka exponential
GpGpcov = tcrossprod(solve(
  Matrix::sparseMatrix(
    i = row(t(vecchia_approx$NNarray))[!is.na(t(vecchia_approx$NNarray))], 
    j = t(vecchia_approx$NNarray)[!is.na(t(vecchia_approx$NNarray))], 
    x=  
  GpGp::vecchia_Linv(
  c(1, exp(range_beta[1]), .5, .000001), covfun_name = "matern_isotropic", 
  (vecchia_approx$locs), NNarray = t(vecchia_approx$NNarray))[!is.na(t(vecchia_approx$NNarray))], 
  triangular = T
  )
  ))
image(GpGpcov[vecchia_approx$locs_match, vecchia_approx$locs_match])
mycov = tcrossprod(solve(
  Matrix::sparseMatrix(
    i = col((vecchia_approx$NNarray))[!is.na((vecchia_approx$NNarray))], 
    j = (vecchia_approx$NNarray)[!is.na((vecchia_approx$NNarray))], 
    x=  
      compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = NULL, matern_smoothness = .5, compute_derivative = F)[!is.na((vecchia_approx$NNarray))], 
    triangular = T
  )
))
image(mycov[vecchia_approx$locs_match, vecchia_approx$locs_match])
hist(mycov- GpGpcov)
range_beta = matrix(rnorm(1 + PP$n_knots))
compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = PP, matern_smoothness = 1.5, compute_derivative = F)
range_beta = matrix(rnorm(1 + PP$n_knots))
compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = PP, matern_smoothness = 1.5, compute_derivative = F)


# test equivalence of parametrizations for locally isotropic and locally aniso
range_beta = matrix(c(-2))
mycov = tcrossprod(solve(
  Matrix::sparseMatrix(
    i = col((vecchia_approx$NNarray))[!is.na((vecchia_approx$NNarray))], 
    j = (vecchia_approx$NNarray)[!is.na((vecchia_approx$NNarray))], 
    x=  
      compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = NULL, matern_smoothness = 1.5, compute_derivative = F)[!is.na((vecchia_approx$NNarray))], 
    triangular = T
  )
))
image(mycov[vecchia_approx$locs_match, vecchia_approx$locs_match])

range_beta = matrix(c(-2, 0,0), 1)
mycov_aniso = tcrossprod(solve(
  Matrix::sparseMatrix(
    i = col((vecchia_approx$NNarray))[!is.na((vecchia_approx$NNarray))], 
    j = (vecchia_approx$NNarray)[!is.na((vecchia_approx$NNarray))], 
    x=  
      compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = NULL, matern_smoothness = 1.5, compute_derivative = F)[!is.na((vecchia_approx$NNarray))], 
    triangular = T
  )
))
image(mycov_aniso[vecchia_approx$locs_match, vecchia_approx$locs_match])
hist(mycov - mycov_aniso)


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



## X_PP_mult_right #########################################################
locs = cbind(runif(1000), runif(1000))
locs = rbind(locs, locs)
vecchia_approx = createVecchia(locs, 12, ncores=1)
PP = createPP(vecchia_approx, plot=FALSE)
covariate_coefficients = c(4, 1, 1, .5)
knots_coeffs = rnorm(PP$n_knots)
X = cbind(1, vecchia_approx$locs, rnorm(nrow(vecchia_approx$locs)))

# multiplying X alone
res1 <- X_PP_mult_right(X = X, Y = covariate_coefficients, vecchia_approx = vecchia_approx)
plot_pointillist_painting(locs, res1, main = "covariates only, \n one covariate for each observation")
res1 <- X_PP_mult_right(X = X, Y = covariate_coefficients, vecchia_approx = vecchia_approx)
plot_pointillist_painting(locs, res1, main = "covariates only, \n one covariate for each observation")

# multiplying PP alone
res2 <- X_PP_mult_right(PP = PP, Y = knots_coeffs, vecchia_approx = vecchia_approx)
plot_pointillist_painting(vecchia_approx$locs, res2, main = "PP only")

# multiplying PP and matrix of covariate, one obs for each location
X_by_loc = cbind(1, vecchia_approx$locs, rnorm(vecchia_approx$n_locs))
res3 <- X_PP_mult_right(PP = PP, X = X_by_loc, Y = c(covariate_coefficients, knots_coeffs), vecchia_approx = vecchia_approx)
plot_pointillist_painting(vecchia_approx$locs, res3, main = "PP + covariates, \n one covariate for each location")

# multiplying PP and matrix of covariates with an index
X_by_obs = cbind(1, vecchia_approx$observed_locs, rnorm(vecchia_approx$n_obs))
res4 <- X_PP_mult_right(X = X_by_obs, PP = PP, 
                      Y = c(covariate_coefficients, knots_coeffs), 
                      vecchia_approx = vecchia_approx
                      )
plot_pointillist_painting(vecchia_approx$observed_locs, res4, main = "PP + covariates,\n  one covariate for each observation")

# multiplying 
X_by_obs = cbind(1, vecchia_approx$observed_locs, rnorm(vecchia_approx$n_obs))
res5 <- X_PP_mult_right(X = NULL, PP = PP, 
                      Y = diag(1, PP$n_knots), 
                      vecchia_approx = vecchia_approx,
                      permutate_PP_to_obs = T
                      )

