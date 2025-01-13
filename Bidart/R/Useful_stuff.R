#' @export
derivative_sandwiches = function(
    derivatives, 
    left_vector, 
    right_vector, 
    NNarray
)
{
  M = matrix(0, length(left_vector), length(derivatives))
  for( i in seq(length(derivatives)))M[,i] = Bidart::derivative_sandwich(derivatives[[i]], left_vector, right_vector, NNarray)
  # changing basis between det/aniso and canonical
  if(ncol(M)==3) M = M %*% t(matrix(
    c(1/sqrt(2), 1/sqrt(2),  0, 
      1/sqrt(2), -1/sqrt(2), 0,
      0,       0,        1), 3))*sqrt(2)
  M
}

#' @export
log_determinant_derivatives = function(sparse_chol_and_grad, NNarray)
{
  M = matrix(0, nrow(NNarray), length(sparse_chol_and_grad[[2]]))
  for( i in seq(length(sparse_chol_and_grad[[2]])))M[,i] = Bidart::log_determinant_derivative(derivative = sparse_chol_and_grad[[2]][[i]], compressed_sparse_chol = sparse_chol_and_grad[[1]], NNarray = NNarray)
  if(ncol(M)==3) M = M %*% matrix(
    c(1/sqrt(2), 1/sqrt(2),  0, 
      1/sqrt(2), -1/sqrt(2), 0,
      0,       0,        1), 3)*sqrt(2)
  M
}

#' @export
expmat = function(coords)
{
  res = expm::expm(symmat(coords)) 
  res + diag(.0001,nrow(res), ncol(res))
}

#' @export
symmat = function(coords)
{
  if(length(coords)==1)logm = matrix(coords)
  if(length(coords)==3)
  {
    logm = matrix(coords [c(1, 3, 3, 2)] , 2)
  }
  if(length(coords)==6)
  {
    logm = matrix(0, 3, 3)
    diag(logm) = coords[seq(3)]
    logm[lower.tri(logm)] = coords[-seq(3)]
    logm[upper.tri(logm)] = logm[lower.tri(logm)]
  }
  logm
}

#' @export
variance_field = function(beta, PP = NULL, use_PP = F, X, locs_idx = NULL)
{
  as.vector(exp(X_PP_mult_right(X = X, PP = PP, use_PP = use_PP, locs_idx = locs_idx, Y = beta)))
}


# NB : computes sparse chol wrt determinant/anisotropy basis of range, 
# but gives derivatives wrt canonical

#' Computes a Vecchia sparse Cholesky factor
#' 
#' @param range_beta parameter for the range. If a stationary covariance is used, the parameter is passed to the exponential. Else, 
#' @param NNarray Vecchia parents array provided by GpGp::find_ordered_nn
#' @param locs matrix of spatial sites
#' @param range_X covariates for range
#' @param PP predictive process obtained through get_PP
#' @param use_PP should the PP be used ?
#' @param compute_derivative logical, indicates if derivatives of Vecchia factors are to be computed
#' @param nu Matern smoothness
#' @param locs_idx match between PP basis function and locs

#' @export
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
  if(ncol(range_beta)==3)range_beta = range_beta %*% matrix(
    c(1/sqrt(2), 1/sqrt(2),  0, 
      1/sqrt(2), -1/sqrt(2), 0,
      0,       0,        1), 3)*sqrt(2)
  if(ncol(range_beta)==1)range_beta = range_beta# / sqrt(2)
  
  log_range = as.matrix(
    Bidart::X_PP_mult_right(
      X = range_X, 
      PP = PP, 
      Y = range_beta,  
      use_PP = use_PP, 
      locs_idx = locs_idx))
  #Bidart::plot_ellipses(locs, log_range)
  # exp locally isotropic
  if((!anisotropic)&(nu==0.5))res = Bidart::nonstat_vecchia_Linv(num_threads=num_threads,log_range = log_range*2, covfun_name = "nonstationary_exponential_isotropic"  , sphere = sphere, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
  
  # matern locally isotropic
  if((!anisotropic)&(nu==1.5)) res = Bidart::nonstat_vecchia_Linv(num_threads=num_threads,log_range = log_range*2, covfun_name = "nonstationary_matern_isotropic"  , sphere = sphere, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
  
  # exp locally anisotropic
  if(( anisotropic)&(nu==0.5)) res = Bidart::nonstat_vecchia_Linv(num_threads=num_threads,log_range = log_range*2, covfun_name = "nonstationary_exponential_anisotropic", sphere = sphere, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
  
  # matern locally anisotropic
  if(( anisotropic)&(nu==1.5)) res = Bidart::nonstat_vecchia_Linv(num_threads=num_threads,log_range = log_range*2, covfun_name = "nonstationary_matern_anisotropic", sphere = sphere, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
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
###         x= Bidart::compute_sparse_chol(
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
###         x= Bidart::compute_sparse_chol(
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
#########

#' @export
get_PP = function(observed_locs, matern_range, lonlat = F, n_PP = 20, m = 10)
{
  locs_ = observed_locs[! duplicated(observed_locs),]
  locs_ = locs_[order(runif(nrow(locs_))),]
  locs_[seq(min(nrow(locs_), 100000)),] = locs_[GpGp::order_maxmin(locs_[seq(min(nrow(locs_), 100000)),]),]
  idx = match(split(observed_locs, row(observed_locs)), split(locs_, row(locs_)))
  
  
  knots = kmeans(locs_[seq(min(nrow(locs_), 100000)),], n_PP, algorithm = "Hartigan-Wong", iter.max = 50)$centers
  knots = knots[GpGp::order_maxmin(knots),]
  
  NNarray = GpGp::find_ordered_nn(rbind(knots, locs_), m, lonlat = lonlat)
  sparse_chol = Matrix::sparseMatrix(
    i = row(NNarray)[!is.na(NNarray)], 
    j = NNarray[!is.na(NNarray)], 
    x = GpGp::vecchia_Linv(covparms = c(1, matern_range, .0001), covfun_name = "matern15_isotropic", locs = rbind(knots, locs_), NNarray = NNarray)[!is.na(NNarray)], 
    triangular = T
  )
  return(list("knots" = knots, "unique_reordered_locs" = locs_, "idx" = idx, "lonlat" = lonlat, "m" = m, "matern_range" = matern_range, "sparse_chol" = sparse_chol, "NNarray" = NNarray, "n_PP" = n_PP))
}

#' @export
beta_prior_log_dens = function(beta, n_PP, beta_mean, beta_precision, log_scale)
{
  if(n_PP>0) 
  {
    scale_mat = Bidart::expmat(-log_scale)
    return(
      (
        # PP coefficients follow N(0, scale_mat)
        +.5 * n_PP * determinant(scale_mat, logarithm = T)$mod # determinant is changed by log scale
        -sum(.5 * c(c(t(beta[seq(nrow(beta)-n_PP),, drop = F])-t(beta_mean)) %*% beta_precision) * c(t(beta[seq(nrow(beta)-n_PP),,drop = F])-t(beta_mean)))
        -sum(.5 * c(beta[-seq(nrow(beta)-n_PP),,drop = F] %*% scale_mat) * beta[-seq(nrow(beta)-n_PP),,drop = F])
      )
    )
  }
  return(
    sum(
      c(
        -.5 * (c(t(beta[seq(nrow(beta)-n_PP),,drop=F])-t(beta_mean)) %*% beta_precision) * c(t(beta[seq(nrow(beta)-n_PP),,drop=F])-t(beta_mean))
      )
    )
  )
}


# #' @export
#beta_prior_log_dens_derivative = function(beta, n_PP, beta_mean, beta_precision, log_scale)
#{
#  if(n_PP>0) 
#  {
#    scale_mat = Bidart::expmat(-log_scale)
#    return(
#      matrix(
#        c(
#          -(c(t(beta[seq(nrow(beta)-n_PP),,drop=F])-t(beta_mean)) %*% beta_precision),
#          -c(t(beta[-seq(nrow(beta)-n_PP),,drop = F] %*% scale_mat))
#        ), 
#        nrow(beta)
#      )
#    )
#  }
#  return(
#    matrix(
#      c(
#        -(c(t(beta[seq(nrow(beta)-n_PP),])-t(beta_mean)) %*% beta_precision)
#      ), 
#      nrow(beta)
#    )
#  )
#}

#' @export
beta_prior_log_dens_derivative = function(beta, n_PP, beta_mean, beta_precision, log_scale)
{
  
  res = matrix(
    c(
      -(c(t(beta[seq(nrow(beta)-n_PP),, drop = F])-t(beta_mean)) %*% beta_precision)
    ), 
    nrow(beta)-n_PP
  )
  if(n_PP>0) 
  {
    scale_mat = Bidart::expmat(-log_scale)
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
#' @export
X_PP_mult_right = function(X = NULL, PP = NULL, use_PP = F, locs_idx = NULL, Y)
{
  if(is.null(locs_idx))if(!is.null(X))locs_idx = seq(nrow(X))
  if(is.null(locs_idx))if(!is.null(PP))locs_idx = seq(length(PP$idx))
  Y = as.matrix(Y)
  res = matrix(0, length(locs_idx), ncol(Y))
  if(!is.null(X))res = res + X  %*% Y[seq(ncol(X)),]
  if(use_PP)
  {
    if(!is.null(X))Y =  Y[-seq(ncol(X)),, drop = F] 
    V = matrix(0, nrow(PP$sparse_chol), ncol(Y))
    V[seq(nrow(Y)),] = Y
    res = res + as.matrix(Matrix::solve(PP$sparse_chol, V, triangular = T))[-seq(nrow(PP$knots)),,drop =F][PP$idx[locs_idx],,drop =F]
  }
  colnames(res) = c("det", "an", "an")[seq(ncol(res))]
  res
}

#' @export
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
### Bidart::plot_pointillist_painting(locs, field = X_PP_mult_right(PP = PP, use_PP = T, Y = rnorm(n_PP)))
### points(PP$knots, pch = 16, cex = .5)
### n_PP = 100
### PP = get_PP(locs, c(1, range, 1.5, 0), n_PP = n_PP, m = 15)
### Bidart::plot_pointillist_painting(locs, field = X_PP_mult_right(PP = PP, use_PP = T, Y = rnorm(n_PP)))
### points(PP$knots, pch = 16, cex = .5)
### # ploting one PP basis
### Bidart::plot_pointillist_painting(locs, 
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


#' @export
derivative_chol_expmat = function(coords)
{
  dimres = 1
  if(length(coords)==6)dimres = 3
  res = array(data = 0, dim = c(dimres, dimres, length(coords)))
  chol_expmat = chol(Bidart::expmat(coords))
  for(i in seq(length(coords)))
  {
    coords_ = coords
    coords_[i] = coords_[i] + 0.00001
    res[,,i] = 100000 * (chol(Bidart::expmat(coords_)) - chol_expmat)
  }
  res
}

#' @export
derivative_field_wrt_scale = function(field, coords)
{
  d_chol_expmat = derivative_chol_expmat(coords)
  white_field = field %*% solve(chol(Bidart::expmat(coords)))
  res = array(0, dim = c(dim(field), length(coords)))
  for(i in seq(length(coords)))
  {
    res[,,i] = white_field %*% d_chol_expmat[,,i]
  }
  res
}
