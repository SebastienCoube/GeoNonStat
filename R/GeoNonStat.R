#' Vecchia approximation setup
#'
#' @param observed_locs a matrix of spatial coordinates where observations are done
#' @param m number of nearest neighbors to do Vecchia's approximation
#'
#' @returns a list
#' @examples
#' set.seed(100)
#' size <- 20000
#' observed_locs = cbind(runif(size), runif(size))
#' res <- createVecchia(observed_locs, m=10) 
createVecchia <- function(observed_locs, m = 12){
  #message("building DAGs and indices for Vecchia approximation...")
  # Vecchia approximation ##########################################################################
  # This object gathers the NNarray table used by GpGp package and related objects
  if (!is.matrix(observed_locs)) stop("observed_locs should be a matrix")
  if (ncol(observed_locs)!=2) stop("observed_locs should have 2 columns")
  
  # re-ordering and treating spatial locations 
  locs = observed_locs[!duplicated(observed_locs), ]
  max_locs <- min(nrow(locs), 10000)
  locs_reordering = order(runif(nrow(locs)))
  locs_reordering[seq(max_locs)] = locs_reordering[GpGp::order_maxmin(locs[locs_reordering[seq(max_locs)], ])]
  locs = locs[locs_reordering, ]
  
  # storing numbers
  n_locs = n = nrow(locs)
  n_obs = nrow(observed_locs)
  # matching observed locations with reordered, unrepeated locations
  locs_match = match(split(observed_locs, row(observed_locs)), 
                     split(locs, row(locs)))
  locs_match_matrix = Matrix::sparseMatrix(i = locs_match,j = seq(n_obs), x = 1)
  # doing reversed operation : for a given unrepeated location, tell which observations correspond
  hctam_scol = split(seq(n_obs), locs_match)
  
  #extracting NNarray =  nearest neighbours for Vecchia approximation
  NNarray = t(GpGp::find_ordered_nn(locs, m))
  #computations from createVecchia$NNarray in order to create sparse Cholesky using Matrix::sparseMatrix
  #non_NA indices from createVecchia$NNarray
  sparse_mat = Matrix::sparseMatrix(x= seq(sum(!is.na(NNarray))), i = col(NNarray)[!is.na(NNarray)], j =NNarray[!is.na(NNarray)], triangular = T)
  sparse_chol_x_reorder = (seq(length(NNarray)))[!is.na(NNarray)][match(sparse_mat@x, seq(sum(!is.na(NNarray))),)]
  # sparse_mat@x = (NNarray + 0.0)[sparse_chol_x_reorder]
  # sparse_mat[30,]
  
  # partition of locs for field update
  cl = parallel::makeCluster(5)
  parallel::clusterExport(cl = cl, varlist = c("locs", "n"), envir = environment())
  locs_partition = parallel::parSapply(
    cl = cl, 
    round(seq(n / 10000 + 1, 2 * (n / 10000) + 1, length.out = min(10, ceiling(n / 10000)))), 
    function(i) { kmeans(locs, centers = i, iter.max = 200, algorithm = "Hartigan-Wong")$cluster}
  )
  colnames(locs_partition) = paste(round(seq(n / 10000 + 1, 2 * (n / 10000) + 1, length.out = min(10, ceiling(n / 10000)))), "clust")
  parallel::stopCluster(cl)
  
  return(list(
    n_locs = n_locs,
    n_obs = n_obs,
    observed_locs = observed_locs, 
    locs = locs, 
    t_locs = t(locs), 
    locs_match = locs_match,
    locs_match_matrix = locs_match_matrix,
    hctam_scol = hctam_scol,
    hctam_scol_1 = sapply(hctam_scol, function(x) x[1]),
    obs_per_loc = unlist(sapply(hctam_scol, length)), # count how many observations correspond to one location
    NNarray = NNarray,
    NNarray_non_NA = !is.na(NNarray),
    sparse_chol_i = sparse_mat@i+1, 
    sparse_chol_p = sparse_mat@p, 
    sparse_chol_x_reorder = sparse_chol_x_reorder, 
    locs_partition = locs_partition
  ))
}

#' process_covariates : pre-processes covariates by adding an intercept,
#' creating useful indices, and pre-computing useful matrices and vectors
#'
#' @param X a m
#' @param vecchia_approx TODO
#' @param explicit_PP_basis TODO
#' @param use_PP TODO
#'
#' @returns a list
#'
#' @examples
#' \dontrun{TODO}
#' 
#' 
#' nobs = 10000
#' observed_locs = cbind(runif(nobs), runif(nobs))
#' X = as.data.frame(cbind(runif(nobs), rnorm(nobs), rpois(nobs, 5)))
#' vecchia_approx = createVecchia(observed_locs, 12)
#' PP = createPP(vecchia_approx)
#' 
#' process_covariates(X = X, vecchia_approx)
#' process_covariates(X = NULL, vecchia_approx)
#' process_covariates(X = X, vecchia_approx, PP)
#' process_covariates(X = NULL, vecchia_approx, PP)
 
process_covariates = function(X, vecchia_approx, PP = NULL){
  # covariates in the observed field #
  res = list()
  # creating model matrix
  # extracting a model matrix and storing the original argument
  if(!is.null(X))
  {
    res$arg = X
    res$X = model.matrix(~., X)
  }
  # extracting a model matrix with only intercept and storing a message about the lack of original argument if no X is provided
  if(is.null(X))
  {
    res$arg = "No covariates were provided"
    res$X = matrix(model.matrix(~., as.data.frame(rep(1, vecchia_approx$n_obs)))[,-1], vecchia_approx$n_obs)
  }
  colnames(res$X)[1] = "(Intercept)"
  
  
  X_ = res$X
  if(!is.null(PP)) X_ = cbind(res$X, X_PP_mult_right(PP = PP, vecchia_approx = vecchia_approx, permutate_PP_to_obs = T, Y = diag(1, PP$n_knots)))
  # pre- computing XTX
  crossprod_X = crossprod(X_)
  res$chol_crossprod_X = chol(crossprod_X)
  res$n_regressors = ncol(res$X)
  # identifying  which X do not vary within location
  res$which_locs = c()
  for(i in seq(ncol(res$X))) 
  {
    if(all(duplicated(cbind(vecchia_approx$observed_locs, res$X[,i])) == vecchia_approx$duplicated_locs)) {
      res$which_locs = c(res$which_locs, i)
    }
  }
  res$X_locs = matrix(res$X[vecchia_approx$hctam_scol_1,res$which_locs], ncol = length(res$which_locs))
  colnames(res$X_locs) = colnames(res$X)[res$which_locs]
  X_locs_ = res$X_locs
  if(!is.null(PP))X_locs_ =cbind(X_locs_, X_PP_mult_right(PP = PP, vecchia_approx = vecchia_approx, permutate_PP_to_obs = F, Y = diag(1, PP$n_knots)))
  res$crossprod_X_locs = crossprod(X_locs_)
  res$chol_crossprod_X_locs = chol(res$crossprod_X_locs)
  res
}

#' Safety checks and automatic treatment of PP objects and their marginal variance bounds
#'
#' @param PP a PP object, as given by create_PP
#' @param log_scale_prior a numeric vector of size 2 indicating the lower and upper PP log-variance bounds
#' @param parameter_name a character string indicating the number of the parameter, used for prints
#'
#' @examples 
#' myPP <- createPP(observed_locs = cbind(runif(100), runif(100))) 
#' PPPP <- process_PP_prior(pepito, c(1,2), "example")
process_PP_prior = function(
  PP = NULL, 
  log_scale_bounds = NULL, 
  parameter_name
){
  if(is.null(PP) & is.null(log_scale_bounds)) return(log_scale_bounds)
  if(is.null(PP) & !is.null(log_scale_bounds))stop(paste("No PP object was provided for the", parameter_name, "parameters, but log - marginal variance bounds were provided") )
  if(class(PP)!= "PP")stop(paste("PP who describes the", parameter_name, "parameters must be of class PP"))
  if (!is.null(PP) & is.null(log_scale_bounds)){
    message(paste("The log - marginal variance (log_scale) bounds for the PP who describes the", parameter_name, "parameters were automatically set to (-6, 3)"))
    log_scale_bounds = c(-6, 3)
  }
  if(!is.numeric(log_scale_bounds) | length(log_scale_bounds) != 2)stop(paste("The log - marginal variance (log_scale) bounds for the PP who describes the", parameter_name, "parameters must be a numeric vector of length 2"))
  log_scale_bounds = sort(log_scale_bounds)
  return(log_scale_bounds)
}

#' Initialize hierarchical model
#'
#' @param noise_PP either an object of class PP used to model the field of log-noise parameters, or NULL in which case the model is stationary
#' @param scale_PP either an object of class PP used to model the field of log-scale parameters, or NULL in which case the model is stationary
#' @param range_PP either an object of class PP used to model the field of log-range parameters, or NULL in which case the model is stationary
#' @param noise_log_scale_bounds either a vector containing two numeric values bounding Uniform prior for the log-marginal variance of the noise's PP, or NULL in which case the bounds are set automatically
#' @param scale_log_scale_bounds either a vector containing two numeric values bounding Uniform prior for the log-marginal variance of the scale's PP, or NULL in which case the bounds are set automatically
#' @param range_log_scale_bounds either a vector containing two numeric values bounding Uniform prior for the log-marginal variance of the range's PP, or NULL in which case the bounds are set automatically
#' @param matern_smoothness a Matérn smoothness parameter, either 0.5 (aka ``exponential kernel'') or 1.5
#' @param observed_field a vector of observations. 
#' @param covariates The list of covariates obtained with process_covariates, contianing X, X_noise, X_range, X_scale
#' @param anisotropic Boolean indicating if anisotropic
#'
#' @returns a list
#' @examples
#' 
#' 
#' nobs = 10000
#' observed_locs = cbind(runif(nobs), runif(nobs))
#' observed_field = rnorm(nobs)
#' vecchia_approx = createVecchia(observed_locs)
#' PP = createPP(vecchia_approx)
#' X = as.data.frame(cbind(rnorm(nobs), runif(nobs)))
#' covariates = list()
#' covariates$X = process_covariates(X, vecchia_approx = vecchia_approx)
#' 
#' hm1 = process_hierarchical_model(
#'   noise_PP = PP, noise_log_scale_bounds = NULL,
#'   scale_PP = PP, scale_log_scale_bounds = NULL,
#'   range_PP = PP, range_log_scale_bounds = NULL,
#'   observed_locs,
#'   matern_smoothness = 1.5,
#'   observed_field = observed_field,
#'   covariates = covariates,
#'   anisotropic = T) 
#' 
#' hm2 = process_hierarchical_model(
#'   noise_PP = PP, noise_log_scale_bounds = c(1, 3),
#'   scale_PP = PP, scale_log_scale_bounds = c(1, 3),
#'   range_PP = PP, range_log_scale_bounds = c(1, 3),
#'   observed_locs,
#'   matern_smoothness = 1.5,
#'   observed_field = observed_field,
#'   covariates = covariates,
#'   anisotropic = T) 
#' 
#' hm3 = process_hierarchical_model(
#'   noise_PP = NULL, noise_log_scale_bounds = NULL,
#'   scale_PP = NULL, scale_log_scale_bounds = NULL,
#'   range_PP = NULL, range_log_scale_bounds = NULL,
#'   observed_locs,
#'   matern_smoothness = 1.5,
#'   observed_field = observed_field,
#'   covariates = covariates,
#'   anisotropic = T) 
#' 
#' # there is an error, it's normal
#' hm4 = process_hierarchical_model(
#'   noise_PP = NULL, noise_log_scale_bounds = c(1, 3),
#'   scale_PP = PP, scale_log_scale_bounds = c(1, 3),
#'   range_PP = PP, range_log_scale_bounds = c(1, 3),
#'   observed_locs,
#'   matern_smoothness = 1.5,
#'   observed_field = observed_field,
#'   covariates = covariates,
#'   anisotropic = T) 
#' hm5 = process_hierarchical_model(
#'   noise_PP = PP, noise_log_scale_bounds = c(1, 3),
#'   scale_PP = NULL, scale_log_scale_bounds = c(1, 3),
#'   range_PP = PP, range_log_scale_bounds = c(1, 3),
#'   observed_locs,
#'   matern_smoothness = 1.5,
#'   observed_field = observed_field,
#'   covariates = covariates,
#'   anisotropic = T) 
#' hm6 = process_hierarchical_model(
#'   noise_PP = PP, noise_log_scale_bounds = c(1, 3),
#'   scale_PP = PP, scale_log_scale_bounds = c(1, 3),
#'   range_PP = NULL, range_log_scale_bounds = c(1, 3),
#'   observed_locs,
#'   matern_smoothness = 1.5,
#'   observed_field = observed_field,
#'   covariates = covariates,
#'   anisotropic = T) 

process_hierarchical_model <- function(noise_PP, noise_log_scale_bounds,
                                       scale_PP, scale_log_scale_bounds,
                                       range_PP, range_log_scale_bounds,
                                       observed_locs,
                                       matern_smoothness,
                                       observed_field,
                                       covariates,
                                       anisotropic) {
  
  # Processing PP priors
  noise_log_scale_bounds <- process_PP_prior(noise_PP, noise_log_scale_bounds, "noise")
  scale_log_scale_bounds <- process_PP_prior(scale_PP, scale_log_scale_bounds, "scale")
  range_log_scale_bounds <- process_PP_prior(range_PP, range_log_scale_bounds, "range")

  # Making a guess for maximum and minimum reasonable values for the range intercept
  # using as upper bound the geographic space size
  # and as lower bound the minimal distance between two space points
  alpha_max = -0.5 * log(8 * matern_smoothness) + log(max(dist(vecchia_approx$locs[seq(min(10000, nrow(vecchia_approx$locs))), ])) / 8)
  alpha_min = -0.5 * log(8 * matern_smoothness) + log(median(FNN::get.knn(vecchia_approx$locs, k = 1)$nn.dist) * 3)
  
  # OLS to get residual variance to make a guess for maximum and minimum reasonable values 
  # for NNGP and noise variance
  naive_ols =  lm(observed_field ~ covariates$X$X - 1)
  sigma_max = log(.99 * var(naive_ols$residuals))
  sigma_min = log(.01 * var(naive_ols$residuals))
  
  # Default mean prior computed from a reasonable case.
  res = list(
    noise_PP = noise_PP,
    noise_log_scale_bounds = noise_log_scale_bounds,
    scale_PP = scale_PP,
    scale_log_scale_bounds = scale_log_scale_bounds,
    range_PP = range_PP, 
    range_log_scale_bounds = range_log_scale_bounds, 
    anisotropic= anisotropic,
    matern_smoothness= matern_smoothness,
    beta_priors= list(),
    PP= PP,
    range_beta0_mean= (alpha_max + alpha_min) / 2,
    range_beta0_var= ((alpha_max - alpha_min) / 8)^2,
    noise_beta0_mean= (sigma_max + sigma_min) / 2,
    noise_beta0_var= ((sigma_max - sigma_min) / 8)^2,
    scale_beta0_mean= (sigma_max + sigma_min) / 2,
    scale_beta0_var=((sigma_max - sigma_min) / 8)^2,
    naive_ols = naive_ols
  )
  return(res)
}


#' Initialize transition kernels
#'
#' @param init a numeric value
#' @returns a list
#' @examples
#' tk <- process_transition_kernels()
#' nobs = 10000
#' observed_locs = cbind(runif(nobs), runif(nobs))
#' observed_field = rnorm(nobs)
#' vecchia_approx = createVecchia(observed_locs)
#' PP = createPP(vecchia_approx)
#' X = as.data.frame(cbind(rnorm(nobs), runif(nobs)))
#' covariates = list()
#' covariates$X = process_covariates(X, vecchia_approx = vecchia_approx)
#' covariates$X_range = process_covariates(X, vecchia_approx = vecchia_approx)
#' covariates$X_noise = process_covariates(X, vecchia_approx = vecchia_approx)
#' covariates$X_scale = process_covariates(X, vecchia_approx = vecchia_approx)
#' 
#' hm = process_hierarchical_model(
#'   noise_PP = PP, noise_log_scale_bounds = NULL,
#'   scale_PP = PP, scale_log_scale_bounds = NULL,
#'   range_PP = PP, range_log_scale_bounds = NULL,
#'   observed_locs,
#'   matern_smoothness = 1.5,
#'   observed_field = observed_field,
#'   covariates = covariates,
#'   anisotropic = T) 
#' process_transition_kernels(hm = hm)

process_transition_kernels <- function(init=-4, hm){
  # Transition kernels ###################################################################
  # Transition kernel state
  # transition kernel variance is given as the log
  # can be used in both stationary and nonstationary cases respectively as a random walk Metropolis or MALA step size
  # have an ancillary and a sufficient version when applicable
  # range
  res = 
  list(
    range_log_scale_sufficient = init,
    range_log_scale_ancillary =  init,
    range_beta_sufficient = rep(init, (1 + 2*hm$anisotropic)*(1+!is.null(hm$range_PP))),
    range_beta_ancillary  = rep(init, (1 + 2*hm$anisotropic)*(1+!is.null(hm$range_PP))),
    
    scale_beta_sufficient_mala = rep(init, 1+!is.null(hm$scale_PP)),
    scale_beta_ancillary_mala  = rep(init, 1+!is.null(hm$scale_PP)),
    scale_log_scale_sufficient = init,
    scale_log_scale_ancillary =  init,
    
    noise_beta_mala = rep(init, 1+!is.null(hm$noise_PP)),
    noise_log_scale = init
  )
  return(res)
}

#' Title
#'
#' @param hm hierarchical model (from `process_hierarchical_model()`)
#' @param covariates TODO
#' @param observed_field TODO
#' @param vecchia_approx TODO
#' @param init_tk TODO
#'
#' @returns a list


set.seed(100)
nobs = 10000
observed_locs = cbind(runif(nobs), runif(nobs))
observed_field = rnorm(nobs)

X = as.data.frame(cbind(runif(nobs), rnorm(nobs), rpois(nobs, 5)))
vecchia_approx = createVecchia(observed_locs, 12)
PP = createPP(vecchia_approx)

covariates = list(
  X = process_covariates(X = X, vecchia_approx), 
  scale_X = process_covariates(X = X, vecchia_approx, PP),
  noise_X = process_covariates(X = X, vecchia_approx, PP),
  range_X = process_covariates(X = X, vecchia_approx, PP)
)
hm = process_hierarchical_model(
  noise_PP = PP, noise_log_scale_bounds = NULL,
  scale_PP = PP, scale_log_scale_bounds = NULL,
  range_PP = PP, range_log_scale_bounds = NULL,
  observed_locs,
  matern_smoothness = 1.5,
  observed_field = observed_field,
  covariates = covariates,
  anisotropic = T) 

process_states(
    hm,
    covariates,
    observed_field,
    vecchia_approx,
    init_tk = -4
)

process_states <- function(
    hm,
    covariates,
    observed_field,
    vecchia_approx,
    init_tk = -4
) {
  
  # initializing sub-lisits in the state
  # Actual parameters
  params = list()
  # Metropolis proposal distribution width for MCMC
  ker_var <- process_transition_kernels(init_tk, hm)
  # Momenta for HMC
  momenta <- list()
  # Useful stuff pre-computed from the parameters
  stuff <- list()
  
  
  # Linear regression coefficients  ############################################
  #starting points for regression coeffs
  perturb = t(chol(vcov(hm$naive_ols))) %*% 
    rnorm(length(hm$naive_ols$coefficients))
  # parameter value 
  params[["beta"]] = hm$naive_ols$coefficients + perturb
  row.names(params[["beta"]]) = colnames(covariates$X$X)
  # Residuals of the OLS model that have to be explained by the latent field and the noise
  stuff$lm_fit = as.vector(covariates$X$X %*% matrix(params[["beta"]], ncol = 1))
  stuff$lm_fit_locs = as.vector(covariates$X$X_locs %*%
                                  matrix(params[["beta"]][covariates$X$which_locs], ncol = 1))
  stuff$lm_residuals = as.vector(observed_field - stuff$lm_fit)
  
  # Range of the NNGP and stuff depending on it ################################
  # parameter format and value
  if (is.null(hm$range_PP)){
    params$range_beta = matrix(0, ncol(covariates$range_X$X), 1 + 2 * hm$anisotropic)
    row.names(params$range_beta) = colnames(covariates$range_X$X)
  }
  if (!is.null(hm$range_PP)){
    params$range_beta = matrix(0, ncol(covariates$range_X$X) + hm$range_PP$n_knots, 1 + 2 * hm$anisotropic)
    row.names(params$range_beta) = c(colnames(covariates$range_X$X), 
                                     paste("PP", seq(hm$PP$n_knots), sep = "_"))
    if (!hm$anisotropic)params$range_log_scale =   hm$range_log_scale_bounds[1]
    if (hm$anisotropic)params$range_log_scale = c(rep(hm$range_log_scale_prior[1], 3),rep(0, 3))
  }
  params$range_beta[1, 1] = hm$range_beta0_mean
  # momenta
  momenta$range_beta_ancillary = matrix(rnorm(length(params$range_beta)), nrow(params$range_beta))
  momenta$range_beta_sufficient = matrix(rnorm(length(params$range_beta)), nrow(params$range_beta))
  if (!is.null(hm$range_PP)) {
    momenta$range_log_scale_ancillary  = rnorm(length(params$range_log_scale))
    momenta$range_log_scale_sufficient = rnorm(length(params$range_log_scale))
  }
  # cholesky factors of precision matrices
  stuff$compressed_chol = compute_sparse_chol(
    range_beta = params$range_beta, vecchia_approx = vecchia_approx, 
    range_X = covariates$range_X$X_locs, 
    PP = hm$range_PP, matern_smoothness = hm$matern_smoothness, 
    compute_derivative = T, num_threads = min(5, max(parallel::detectCores()-1, 1))
      )
  stuff$sparse_chol = Matrix::sparseMatrix(
    i = vecchia_approx$sparse_chol_i, 
    p = vecchia_approx$sparse_chol_p, 
    x = stuff$compressed_chol[,1,][vecchia_approx$sparse_chol_x_reorder], 
    triangular = T
  ) 
  #plot_pointillist_painting(vecchia_approx$locs, as.vector(Matrix::solve(stuff$sparse_chol, rnorm(nrow(vecchia_approx$locs)))))
  
  # Noise variance  and stuff depending on it ##################################
  # parameter format and value
  if(is.null(hm$noise_PP)) {
    params$noise_beta = matrix(rep(0, ncol(covariates$noise_X$X)), ncol = 1) #random starting values
    row.names(params$noise_beta) = colnames(covariates$noise_X$X)
  }
  if(!is.null(hm$noise_PP)) {
    params$noise_beta = matrix(rep(0, ncol(covariates$noise_X$X) + hm$noise_PP$n_knots), ncol = 1) #random starting values
    row.names(params$noise_beta) = c(colnames(covariates$noise_X$X), 
                                     paste("PP", seq(hm$PP$n_knots), sep = "_"))
    params$noise_log_scale = hm$noise_log_scale_prior[1]
  }
  params$noise_beta[1] = hm$noise_beta0_mean
  # momenta
  momenta$noise_beta = rnorm(length(params$noise_beta))
  # effective variance field, shall be used in density computations
  stuff$noise_var = exp(X_PP_mult_right(
    X = covariates$noise_X$X, PP = hm$noise_PP,
    vecchia_approx = vecchia_approx, Y = params$noise_beta, 
    permutate_PP_to_obs = T
  ))
  
  # Marginal variance of the NNGP and stuff depending on it ####################
  if (is.null(hm$scale_PP)) {
    params$scale_beta    = matrix(0, ncol(covariates$scale_X$X_locs), ncol = 1)
    row.names(params$scale_beta) = colnames(covariates$scale_X$X_locs)
    momenta$scale_beta_ancillary =  rnorm(ncol(covariates$scale_X$X_locs))
    momenta$scale_beta_sufficient = rnorm(ncol(covariates$scale_X$X_locs))
  } 
  if (!is.null(hm$scale_PP)) {
    params$scale_beta    = matrix(0, ncol(covariates$scale_X$X_locs) + hm$scale_PP$n_knots, ncol = 1)
    row.names(params$scale_beta) = c(colnames(covariates$scale_X$X_locs),
                                     paste("PP", seq(hm$PP$n_knots), sep = "_"))
    params$scale_log_scale = hm$scale_log_scale_prior[1]
    momenta$scale_beta_ancillary =  rnorm(ncol(covariates$scale_X$X_locs) + 
                                            hm$scale_PP$n_knots)
    momenta$scale_beta_sufficient = rnorm(ncol(covariates$scale_X$X_locs) + 
                                            hm$scale_PP$n_knots)
  }
  params$scale_beta[1] = hm$scale_beta0_mean
  # effective variance field, shall be used in density computations
  stuff$field_sd = exp(.5*X_PP_mult_right(
    X = covariates$scale_X$X_locs, PP = hm$scale_PP,
    vecchia_approx = vecchia_approx, Y = params$scale_beta, 
    permutate_PP_to_obs = F
  ))
  
  # Latent field ###############################################################
  params$field = stuff$field_sd * as.vector(Matrix::solve(
    stuff$sparse_chol,
    rnorm(vecchia_approx$n_locs)
  ))
  
  return(list(
    "params" = params, #parameters of interest to the model
    "stuff" = stuff, #stuff that is useful for computations
    "momenta" = momenta, # HMC momenta
    "ker_var" = ker_var # Starting points for transition kernels, will be adaptively tuned
  ))
}

# S3 class GeoNonStat
#' Create an object of class GeoNonStat
#' @param observed_locs a matrix of spatial coordinates where observations are done
#' @param observed_field a vector of observations of the interest variable
#' @param X a data.frame of covariates explaining the interest variable through fixed linear effects
#' @param m number of nearest neighbors to do Vecchia's approximation
#' @param matern_smoothness Matern smoothness, either 0.5 or 1.5
#' @param anisotropic anisotropic covariance
#' @param PP TODO
#' @param n_chains TODO
#' @param noise_PP TODO
#' @param noise_X a data.frame of covariates explaining the Gaussian noise variance through fixed linear effects
#' @param noise_log_scale_prior  1 times 2 matrix for the prior on the log-variance of the noise PP field.
#' @param scale_PP TODO
#' @param scale_X  a data.frame of covariates explaining the Gaussian process marginal variance through fixed linear effects
#' @param scale_log_scale_prior  1 times 2 matrix for the prior on the log-variance of the scale PP field.
#' @param range_PP TODO
#' @param range_X  a data.frame of covariates explaining the Gaussian process range through fixed linear effects
#' @param range_log_scale_prior 1 times 2 matrix for the prior on the log-variance of the range PP field. #' In the case of anisotropic range, input an 3 times 2 matrix, indicating bounds for the eigenvalues of the trivariate log-variance matrix.
#' @param seed TODO
#'
#' @returns a list
#' @export
#'
#' @examples
#' set.seed(100)
#' locs = cbind(runif(100), runif(100))
#' myPP = createPP(
#'   observed_locs = locs, # spatial sites
#'   matern_range = .1,
#'   knots = 50, # number of knots
#'   m = 15 # number of NNGP parents
#' )
#' myobj = GeoNonStat(
#'   observed_locs = locs, 
#'   observed_field = rnorm(100),
#'   matern_smoothness = 1.5, n_chains = 5,
#'   range_PP = T, PP = myPP, # use PP for range
#'   anisotropic = T # Covariance will be anisotropic
#' )
GeoNonStat <- 
  function(
    vecchia_approx,
    observed_field,  #spatial locations
    X = NULL, # Response variable 
    # Covariates per observation
    m = 10, #number of Nearest Neighbors
    matern_smoothness = 1.5, #Matern smoothness
    anisotropic = FALSE, 
    PP = NULL,
    n_chains = 2,
    # number of MCMC chains
    noise_PP = F,
    noise_X = NULL,
    noise_log_scale_prior = NULL,
    scale_PP = F,
    scale_X = NULL,
    scale_log_scale_prior = NULL,
    range_PP = F,
    range_X = NULL,
    range_log_scale_prior = NULL,
    seed = 1)
  {
    # time
    t_begin = Sys.time()
    # seed
    set.seed(seed)
    # cleansing RAM
    gc()
    
    # Sanity checks #####################################################################
    # format
    
    if (!is.vector(observed_field)) stop("observed_field should be a vector")
    
    if (!is.data.frame(X) & !is.null(X)) stop("X should be a data.frame or NULL")
    if (!is.data.frame(noise_X) & !is.null(noise_X)) stop("noise_X should be a data.frame or NULL")
    if (!is.data.frame(scale_X) & !is.null(scale_X)) stop("scale_X should be a data.frame or NULL")
    if (!is.data.frame(range_X) & !is.null(range_X)) stop("range_X should be a data.frame or NULL")
    
    if ((is.null(PP)) & (noise_PP | range_PP | scale_PP))
      stop("either noise_PP, range_PP, or scale_PP is TRUE, while nothing was provided for PP")
    #length of observations
    if (!all(
      unique(c(length(observed_field), nrow(observed_locs), nrow(X), nrow(scale_X), nrow(noise_X), nrow(range_X), length(PP$idx)))
        %in% c(0, length(observed_field)))) {
        stop(
          paste(
            "Lengths are not matching : observed_field has", length(observed_field), "observations,",
            "observed_locs has", nrow(observed_locs), "rows,",
            "X has", nrow(X), "rows,",
            "scale_X has", nrow(scale_X), "rows (can only be either 0 or the length of the observations),",
            "noise_X has", nrow(noise_X), "rows (can only be either 0 or the length of the observations),",
            "range_X has", nrow(range_X), "rows (can only be either 0 or the length of the observations),",
            "PP has", length(PP$idx), "locations (can only be either 0 or the length of the observations)"
          )
        )
     }
    
    # smoothness
    if (!matern_smoothness %in% c(1.5, .5))
      stop("only matern_smoothness = 1.5 or matern_smoothness = 0.5")
    
    # Re-ordering ############################################################################
    # remove duplicated locations
    # Vecchia approximation ##########################################################################
    # This object gathers the NNarray table used by GpGp package and related objects

    
    # covariates #########################################################
    
    covariates = list()
    # fixed effects for response
    covariates$X = process_covariates(X, observed_locs, vecchia_approx)
    # explicit PP basis
    explicit_PP_basis = NULL
    if (!is.null(PP)) {
      explicit_PP_basis = X_PP_mult_right(PP = PP,
                                          Y = diag(1, nrow(PP$knots), nrow(PP$knots)))
    }
    # fixed effects and PP for range
    covariates$range_X = process_covariates(range_X,
                                            observed_locs,
                                            vecchia_approx,
                                            range_PP)
    if (!identical(covariates$range_X$which_locs, seq(ncol(covariates$range_X$X_locs))))
      stop("The covariates range_X cannot vary within one spatial location of observed_locs")
    # fixed effects and PP for scale
    covariates$scale_X = process_covariates(scale_X,
                                            observed_locs,
                                            vecchia_approx,
                                            scale_PP)
    if (!identical(covariates$scale_X$which_locs, seq(ncol(covariates$scale_X$X))))
      stop("The covariates scale_X cannot vary within one spatial location of observed_locs")
    # fixed effects and PP for noise
    covariates$noise_X = process_covariates(noise_X,
                                            observed_locs,
                                            vecchia_approx,
                                            noise_PP)
    
    # Info about hierarchical model ##############################################################
    hierarchical_model <- process_hierarchical_model(
      PP,
      noise_PP, noise_log_scale_prior,
      scale_PP, scale_log_scale_prior,
      range_PP, range_log_scale_prior,
      locs=vecchia_approx$locs,
      matern_smoothness = matern_smoothness,
      observed_field = observed_field,
      covariates = covariates,
      anisotropic = anisotropic)
    
    naive_ols <- hierarchical_model$naive_ols
    
    # Chain states #################################################################
    # for each chain, creating sub-lists in order to stock all the stuff that is related to one chain, including :
    # transition_kernel_sd : a list that stocks the (current) automatically-tuned transition kernels standard deviations
    # params : a list that stocks the (current) parameters of the model, including covariance parameters, the value of the sampled field, etc
    # record : record of the MCMC iterations
    state = process_states(
      hm = hierarchical_model,
      covariates = covariates,
      observed_field = observed_field,
      anisotropic = anisotropic,
      vecchia_approx = vecchia_approx,
      init_tk = -4
    ) 
    # Remove unecessary naive OLS
    hierarchical_model$naive_ols <- NULL
    
    # Chain records setup #########################################################################
    
    # records is a list that stocks the recorded parameters of the model, including covariance parameters, the value of the sampled field, etc. In terms of RAM, those are the biggest bit !
    # iteration is a 2-colums matrix that records the iteration at the end of each chains join and the associated CPU time
    records = list()
    records =  sapply(state$params, function(x) NULL)
    iterations = list()
    iterations$checkpoints =  matrix(c(0, as.numeric(Sys.time() - t_begin, unit = "mins")), ncol = 2)
    colnames(iterations$checkpoints) = c("iteration", "time")
    iterations$thinning = c()
    
    # Result ####################################################################
    res <- structure(
      list(
      "data" =
        list(
          "locs" = locs,
          "observed_field" = observed_field,
          "observed_locs" = observed_locs,
          "covariates" = covariates
        ),
      "hierarchical_model" = hierarchical_model,
      "vecchia_approx" = vecchia_approx,
      "states" = state,
      "records" = records,
      "seed" = seed,
      "iterations" = iterations
    ), class = "GeoNonStat")
  message(paste(
    "Setup done,",
    as.numeric(Sys.time() - t_begin, units = "secs"),
    "s elapsed"
  ))
  return(res)
}

#' Print a 'GeoNonStat' object
#'
#' @param x an object of class \code{GeoNonStat}
#' @param ... Unused yet
#' @export
print.GeoNonStat <- function(x, ...) {
  cat("Object of class 'GeoNonStat'\n")
  print(paste(length(x$data$observed_field), "observed fields"))
  print(paste(nrow(x$data$observed_locs), "observed locations"))
  print(paste("Currently", length(x$states$params$beta), "states"))
}

#' Summary of a 'GeoNonStat' object 
#'
#' @param object an object of class \code{object}
#' @param ... additional arguments
#' @export
summary.GeoNonStat <- function(object, ...) {
  
  cat("### data ###")
  detailed_summary(object$data)
  cat("### hierarchical_model ###")
  detailed_summary(object$hierarchical_model)
  cat("### vecchia_approx ###",
    summary(object$vecchia_approx), 
    "### states ###",
    summary(object$states), 
    "### records ###", 
    summary(object$records), 
    "### seed ###", 
    object$seed, 
    "### iterations ###", 
    paste(
      nrow(object$iterations$checkpoints),
      "iterations for a total time of", 
      round(sum(object$iterations$checkpoints[,'time']),5), "seconds"
    ),
    sep="\n"
  )
}

detailed_summary <- function(partobject){
  sumdata <- summary(partobject)
  sumdata <- cbind(sumdata, 
                   "Value" = as.character(sapply(dimnames(sumdata)[[1]],
                                                 function(x) {
                                                   res <- ""
                                                   if(sumdata[x,"Length"] == 1 & sumdata[x,"Mode"] %in% c("numeric", "character"))
                                                     res <- partobject[[x]]
                                                   res
                                                 }
                   )))
  return(sumdata)
}