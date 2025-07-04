% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Useful_stuff.R
\name{compute_sparse_chol}
\alias{compute_sparse_chol}
\title{Computes a Vecchia sparse Cholesky factor and its derivatives}
\usage{
compute_sparse_chol(
  range_beta,
  vecchia_approx,
  range_X = NULL,
  PP = NULL,
  matern_smoothness = 1.5,
  compute_derivative = T,
  num_threads = 1
)
}
\arguments{
\item{range_beta}{parameter for the range.
If the covariance is anisotropic, it must have 3 columns. It the covariance is isotropic, it must have 1 column. 
The first coefficients are multiplied with range_X 
The last coefficients are multiplied with the spatial basis functions of PP}

\item{vecchia_approx}{TODO}

\item{range_X}{covariates for range}

\item{PP}{predictive process obtained through `createPP()`}

\item{compute_derivative}{logical, indicates if derivatives of Vecchia factors are to be computed}

\item{num_threads}{numerical, number of treads to use. Default to 1.}

\item{smoothness}{Matern smoothness Default to 1.5. Can be 0.5 or 1.5.}

\item{anisotropic}{Logical, default to FALSE. TODO}

\item{locs_idx}{match between the duplicated locations used to buile the PP basis function and the non-redundant locs used to compute the sparse chol}
}
\value{
a list
}
\description{
Computes a Vecchia sparse Cholesky factor and its derivatives
}
\examples{
locs = cbind(runif(1000), runif(1000))
vecchia_approx = createVecchia(locs)
Rcpp::sourceCpp("src/vecchia.cpp")
range_X = matrix(1, nrow(locs))
PP = createPP(vecchia_approx)



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
# test with a PP
range_beta = matrix(rnorm(1 + PP$n_knots))
compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = PP, matern_smoothness = 1.5, compute_derivative = T)
range_beta = matrix(rnorm(3*(1 + PP$n_knots)), ncol = 3)
compute_sparse_chol(range_beta = range_beta, vecchia_approx = vecchia_approx, range_X = range_X, PP = PP, matern_smoothness = 1.5, compute_derivative = T)
}
