process_covariates = function(X, observed_locs, vecchia_approx, explicit_PP_basis = NULL, use_PP = F)
{
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
    res$X = matrix(model.matrix(~., as.data.frame(rep(1, nrow(observed_locs))))[,-1], nrow(observed_locs))
    colnames(res$X) = "(Intercept)"
  }
  X_ = res$X
  
  if(use_PP) X_ = cbind(res$X, explicit_PP_basis)
  # pre- computing XTX
  crossprod_X = crossprod(X_)
  res$chol_crossprod_X = chol(crossprod_X)
  res$n_regressors = ncol(res$X)
  # identifying  which X do not vary within location
  res$which_locs = c()
  for(i in seq(ncol(res$X))) 
  {
    if(all(duplicated(cbind(observed_locs, res$X[,i])) == vecchia_approx$duplicated_locs)) res$which_locs = c(res$which_locs, i)
  }
  res$X_locs = matrix(res$X[vecchia_approx$hctam_scol_1,res$which_locs], ncol = length(res$which_locs))
  colnames(res$X_locs) = colnames(res$X)[res$which_locs]
  X_locs_ = res$X_locs
  if(use_PP)X_locs_ = cbind(X_locs_, explicit_PP_basis[vecchia_approx$hctam_scol_1,])
  res$crossprod_X_locs = crossprod(X_locs_)
  res$chol_crossprod_X_locs = chol(res$crossprod_X_locs)
  #res$chol_crossprod_X_locs = (eigen(res$crossprod_X_locs)$val^.5) * t(eigen(res$crossprod_X_locs)$vec)
  res
}


#' Initializes a list containing the observations, the hierarchical model, and the MCMC chains.
#' 
#' @param observed_locs a matrix of spatial coordinates where observations are done
#' @param observed_field a vector of observations of the interest variable
#' @param X a data.frame of covariates explaining the interest variable through fixed linear effects
#' @param m number of nearest neighbors to do Vecchia's approximation
#' @param nu Matern smoothness, either 0.5 or 1.5
#' @param anisotropic anisotropic covariance
#' @param sphere Boolean, indicating lon-lat data
#' @param noise_X a data.frame of covariates explaining the Gaussian noise variance through fixed linear effects
#' @param range_X a data.frame of covariates explaining the Gaussian process range through fixed linear effects
#' @param noise_beta_mean vector indicating the prior mean for the regression coefficients explaining the Gaussian noise variance through fixed linear effects. Filled automatically if NULL. The number of rows is equal to the number of variables in X_noise after adding an intercept and expanding the factors, the number of columns is 1.
#' @param range_beta_mean vector (3-columns matrix in anisotropy case) indicating the prior mean for the regression coefficients explaining the Gaussian process range through fixed linear effects. Filled automatically if NULL.The number of rows is equal to the number of variables in X_range after adding an intercept, the number of columns is 1 if isotropic function or 3 if anisotropic function.
#' @param noise_beta_precision matrix for the prior precision for the regression coefficients explaining the Gaussian noise variance through fixed linear effects. Filled automatically if NULL.
#' @param range_beta_precision matrix for the prior precision for the regression coefficients explaining the Gaussian process range through fixed linear effects. Filled automatically if NULL. In anisotropic case, the matrix has size 3*n_var, each 3-block indicating the precision for a determinant-direction-direction 3-uplet: intecrept-det, intercept-dir, intercept-dir, V1-det, V1-dir, V1-dir, etc...
#' @param noise_log_var_prior 1 times 2 matrix for the prior on the log-variance of the noise PP field. 
#' @param range_log_var_prior 1 times 2 matrix for the prior on the log-variance of the range PP field. 
#' In the case of anisotropic range, input an 3 times 2 matrix, indicating bounds for the eigenvalues of the trivariate log-variance matrix. 
 

observed_locs = cbind(runif(10000), runif(10000))
observed_field = rnorm(10000)
X = NULL # Covariates per observation
m = 10 #number of Nearest Neighbors
nu =1.5 #Matern smoothness
anisotropic = F 
sphere = F
PP = NULL
n_chains = 2  # number of MCMC chains
noise_PP = F
noise_X = NULL
noise_beta_mean = NULL
noise_beta_precision = NULL
noise_log_var_prior = NULL
range_PP = F
range_X = NULL
range_beta_mean = NULL
range_beta_precision = NULL
range_log_var_prior = NULL 
seed = 1
log_gp_var_prior = NULL

mcmc_nngp_initialize_nonstationary = 
  function(observed_locs, #spatial locations
           observed_field, # Response variable
           X = NULL, # Covariates per observation
           m = 10, #number of Nearest Neighbors
           nu =1.5, #Matern smoothness
           anisotropic = F, 
           sphere = F, 
           PP = NULL, 
           n_chains = 2,  # number of MCMC chains
           noise_PP = F, noise_X = NULL, noise_beta_mean = NULL, noise_beta_precision = NULL, noise_log_var_prior = NULL, 
           range_PP = F, range_X = NULL, range_beta_mean = NULL, range_beta_precision = NULL, range_log_var_prior = NULL, 
           log_gp_var_prior = NULL, 
           seed = 1
  )
  {
    
    # time  
    t_begin = Sys.time()
    # seed
    set.seed(seed)
    # cleansing RAM
    gc()
    #################
    # Sanity checks #
    #################
    
    # format
    if(!is.matrix(observed_locs))stop("observed_locs should be a matrix")
    if(!is.vector(observed_field))stop("observed_field should be a vector")
    
    if(!is.data.frame(X) & !is.null(X))stop("X should be a data.frame or NULL")
    
    if(!is.data.frame(noise_X) & !is.null(noise_X))stop("noise_X should be a data.frame or NULL")
    if(!is.data.frame(range_X) & !is.null(range_X))stop("range_X should be a data.frame or NULL")
    
    if((is.null(PP)) & (noise_PP | range_PP))stop("either noise_PP or range_PP is TRUE, while nothing was provided for PP")
    
    if(!is.null(noise_beta_mean))if(!is.matrix(noise_beta_mean))stop("noise_beta_mean should be a matrix or NULL")
    if(!is.null(range_beta_mean))if(!is.matrix(range_beta_mean))stop("range_beta_mean should be a matrix or NULL")
    
    if(!is.null(noise_beta_precision))if(!is.matrix(noise_beta_precision))stop("noise_beta_precision should be a matrix or NULL")
    if(!is.null(range_beta_precision))if(!is.matrix(range_beta_precision))stop("range_beta_precision should be a matrix or NULL")
    
    if(!is.null(log_gp_var_prior)){
      if(length(log_gp_var_prior)!=2)  stop("log_gp_var_prior must be a numeric vector of length 2 or NULL")
      if(!is.numeric(log_gp_var_prior))stop("log_gp_var_prior must be a numeric vector of length 2 or NULL")
    }
    
    #length of observations
    if(
      !all(unique(c(
        length(observed_field),
        nrow(observed_locs),
        nrow(X),
        nrow(noise_X),
        nrow(range_X),
        length(PP$idx)
        )) %in% c(0, length(observed_field))
      )) stop(
        paste("Lengths are not matching : observed_field has", length(observed_field), "observations,",
              "observed_locs has", nrow(observed_locs), "rows,", 
              "X has", nrow(X), "rows,", 
              "noise_X has", nrow(noise_X), "rows (can only be either 0 or the length of the observations),", 
              "range_X has", nrow(range_X), "rows (can only be either 0 or the length of the observations),", 
              "PP has", length(PP$idx), "locations (can only be either 0 or the length of the observations)"
        )
      )
    
    # smoothness
    if (!nu %in% c(1.5, .5))stop("only nu = 1.5 or nu = 0.5")

  
    ###############
    # Re-ordering #
    ###############
    # remove duplicated locations
    duplicated_locs = duplicated (observed_locs)
    locs = observed_locs[duplicated_locs==F,]
    locs_reordering = order(runif(nrow(locs))); locs_reordering[seq(min(nrow(locs), 100000))] = locs_reordering[GpGp::order_maxmin(locs[locs_reordering[seq(min(nrow(locs), 100000))],])]
    locs = locs[locs_reordering,]
    # extracting number of locations as shortcut
    n = nrow(locs)
    
    #########################
    # Vecchia approximation #
    #########################
    
    # This object gathers the NNarray table used by GpGp package and related objects
    
    vecchia_approx = list()
    # storing numbers
    vecchia_approx$n_locs = n
    vecchia_approx$n_obs = length(observed_field)
    # matching observed locations with reordered, unrepeated locations
    locs_match = match(split(observed_locs, row(observed_locs)), split(locs, row(locs)))
    vecchia_approx$locs_match = locs_match
    vecchia_approx$locs_match_matrix = Matrix::sparseMatrix(i = vecchia_approx$locs_match, j = seq(vecchia_approx$n_obs), x = 1)
    # doing reversed operation : for a given unrepeated location, tell which observations correspond
    vecchia_approx$hctam_scol = split(seq(vecchia_approx$n_obs), locs_match)
    vecchia_approx$hctam_scol_1 = sapply(vecchia_approx$hctam_scol, function(x)x[1])
    # count how many observations correspond to one location
    vecchia_approx$obs_per_loc = unlist(sapply(vecchia_approx$hctam_scol, length))
    #extracting NNarray =  nearest neighbours for Vecchia approximation
    vecchia_approx$NNarray = GpGp::find_ordered_nn(locs, m)
    
    #computations from vecchia_approx$NNarray in order to create sparse Cholesky using Matrix::sparseMatrix
    #non_NA indices from vecchia_approx$NNarray
    vecchia_approx$NNarray_non_NA = !is.na(vecchia_approx$NNarray)
    #column idx of the uncompressed sparse Cholesky factor
    vecchia_approx$sparse_chol_column_idx = vecchia_approx$NNarray[vecchia_approx$NNarray_non_NA]
    #row idx of the uncompressed sparse Cholesky factor
    vecchia_approx$sparse_chol_row_idx = row(vecchia_approx$NNarray)[vecchia_approx$NNarray_non_NA]
    ### # adjacency matrix of MRF
    ### vecchia_approx$MRF_adjacency_mat =  Matrix::crossprod(Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = 1))
    ### # stupid trick to coerce adjacency matrix format...
    ### vecchia_approx$MRF_adjacency_mat[1, 2] = 0
    ### vecchia_approx$MRF_adjacency_mat[1, 2] = 1
    ### vecchia_approx$MRF_adjacency_mat@x = rep(1, length(vecchia_approx$MRF_adjacency_mat@x))
    #vecchia_approx$coloring = naive_greedy_coloring(vecchia_approx$MRF_adjacency_mat)
    # duplicated locs
    vecchia_approx$duplicated_locs = duplicated_locs
    # partition of locs for field update
    cl = parallel::makeCluster(max(1, parallel::detectCores()-2))
    parallel::clusterExport(cl = cl, varlist = c("locs", "n"), envir = environment())
    vecchia_approx$locs_partition = parallel::parSapply(
      cl = cl, 
      round(seq(n/10000+1, 2*(n/10000)+1, length.out = min(20, ceiling(n/10000)))), function(i)
      {
        kmeans(locs, centers = i, iter.max = 200,  algorithm = "Hartigan-Wong")$cluster
      })
    parallel::stopCluster(cl)
    
    ##############
    # covariates #
    ##############
    covariates = list()
    # fixed effects for response
    covariates$X = process_covariates(X, observed_locs, vecchia_approx)  
    # explicit PP basis
    explicit_PP_basis = NULL
    if(!is.null(PP))explicit_PP_basis = X_PP_mult_right(PP = PP, use_PP = T, Y = diag(1, nrow(PP$knots), nrow(PP$knots)))
    # fixed effects and PP for range
    covariates$range_X = process_covariates(range_X, observed_locs, vecchia_approx, explicit_PP_basis, range_PP)
    if(!identical(covariates$range_X$which_locs, seq(ncol(covariates$range_X$X_locs))))stop("The covariates range_X cannot vary within one spatial location of observed_locs")
    # fixed effects and PP for noise
    covariates$noise_X = process_covariates(noise_X, observed_locs, vecchia_approx, explicit_PP_basis, noise_PP)
    # explicit PP basis removal
    remove(explicit_PP_basis)
    
    #################################
    # Info about hierarchical model #
    #################################
    
    hierarchical_model = list()
    hierarchical_model$anisotropic = anisotropic
    hierarchical_model$sphere = sphere
    hierarchical_model$nu = nu
    hierarchical_model$beta_priors = list()
    hierarchical_model$PP = PP
    if(is.null(hierarchical_model$PP)) hierarchical_model$PP = list("n_PP" = 0)
    hierarchical_model$noise_PP = noise_PP
    hierarchical_model$range_PP = range_PP
    
    if(is.null(noise_log_var_prior)&noise_PP)
      {
      message("noise_log_var_prior was automatically set to an uniform on (-6, 2)")
      noise_log_var_prior = c(-6, 2)
      }
    if(!is.null(noise_log_var_prior))hierarchical_model$noise_log_var_prior = matrix(noise_log_var_prior)
    if(is.null(range_log_var_prior)&range_PP)
      {
      message("range_log_var_prior was automatically set to an uniform on (-6, 2)")
      hierarchical_model$range_log_var_prior = c(-6, 2)
      }
    if(!is.null(range_log_var_prior))hierarchical_model$range_log_var_prior = matrix(range_log_var_prior)
    
    
    
    # OLS to get residual variance to make a guess 
    naive_ols =  lm(observed_field~covariates$X$X-1)
    lm_fit = as.vector(covariates$X$X%*%matrix(naive_ols$coefficients, ncol = 1))
    lm_residuals = as.vector(observed_field- lm_fit)
    
    if(is.null(log_gp_var_prior))hierarchical_model$log_gp_var_prior = c(log(var(lm_residuals))-log(100), log(var(lm_residuals)))
    hierarchical_model$log_gp_var_prior = sort(hierarchical_model$log_gp_var_prior) 
      
    hierarchical_model$beta_priors$noise_beta_mean = noise_beta_mean
    hierarchical_model$beta_priors$range_beta_mean = range_beta_mean
    # Default mean prior computed from a reasonable case. 
    # The intercept is set to reasonable value and the rest is set to 0
    if(is.null(noise_beta_mean)) hierarchical_model$beta_priors$noise_beta_mean = matrix(c(log(var(lm_residuals)) - log(2),                     rep(0, covariates$noise_X$n_regressors-1)), ncol=1)
    if(is.null(range_beta_mean)) hierarchical_model$beta_priors$range_beta_mean = matrix(0, covariates$range_X$n_regressors, 1 + 2*anisotropic)
    hierarchical_model$beta_priors$range_beta_mean[1,1] = c(log(max(dist(locs[seq(1000), seq(2)])))-log(50))
    hierarchical_model$beta_priors$noise_beta_precision = noise_beta_precision
    hierarchical_model$beta_priors$range_beta_precision = range_beta_precision
    if(is.null(noise_beta_precision)) hierarchical_model$beta_priors$noise_beta_precision =diag(.01, covariates$noise_X$n_regressors, covariates$noise_X$n_regressors)
    if(is.null(range_beta_precision)) 
      {
      hierarchical_model$beta_priors$range_beta_precision = diag(.01, covariates$range_X$n_regressors, covariates$range_X$n_regressors)
      hierarchical_model$beta_priors$range_beta_precision = hierarchical_model$beta_priors$range_beta_precision %x% diag(1, (1+2*anisotropic))
    }
    
    
    ################
    # Chain states #
    ################
    # for each chain, creating sub-lists in order to stock all the stuff that is related to one chain, including : 
    # transition_kernel_sd : a list that stocks the (current) automatically-tuned transition kernels standard deviations
    # params : a list that stocks the (current) parameters of the model, including covariance parameters, the value of the sampled field, etc
    # records : records of the MCMC iterations
    
    states = list()
    # creating sub-list : each sub-list is a chain
    for(i in seq(n_chains))
    {
      states[[paste("chain", i, sep = "_")]] = list()
      #parameters of interest to the model
      states[[i]]$params = list()
      #stuff that is useful for computations
      states[[i]]$sparse_chol_and_stuff = list()
      # HMC momenta
      states[[i]]$momenta = list()
      
      ######################
      # Transition kernels #  
      ######################
      
      # Starting points for transition kernels, will be adaptively tuned
      states[[i]]$transition_kernels = list()
      # Transition kernel state
      # transition kernel variance is given as the log
      # can be used in both stationary and nonstationary cases respectively as a random walk Metropolis or MALA step size 
      # have an ancillary and a sufficient version when applicable
      # range
      states[[i]]$transition_kernels$range_log_var_sufficient = -4
      states[[i]]$transition_kernels$range_log_var_ancillary =  -4
      if(anisotropic){
        states[[i]]$transition_kernels$aniso_sufficient = c(-4)
        states[[i]]$transition_kernels$aniso_ancillary  = c(-4)
      }
      # range-var ratio
      states[[i]]$transition_kernels$kappa_sufficient = c(-4, 0)
      states[[i]]$transition_kernels$kappa_ancillary  = c(-4, 0)
      # noise variance
      states[[i]]$transition_kernels$noise_beta    = -4
      states[[i]]$transition_kernels$noise_log_var = -4
      
      #########
      # Range #
      #########
      
      # mean log-range (only parameter in the case of stationary model)
      if((!anisotropic)&(!sphere)) states[[i]]$params$range_beta = matrix(  sample(log(max(       dist(locs[sample(seq(vecchia_approx$n_locs)        , 1000),           ]))     )-log(seq(50, 500, 1)), 1), nrow = 1)
      if((!anisotropic)&( sphere)) states[[i]]$params$range_beta = matrix(  sample(log(max(fields::rdist.earth(locs[sample(seq(vecchia_approx$n_locs), 1000),           ]))/6000)-log(seq(50, 500, 1)), 1), nrow = 1)
      if(( anisotropic)&(!sphere)) states[[i]]$params$range_beta = matrix(c(sample(log(max(       dist(locs[sample(seq(vecchia_approx$n_locs)        , 1000),           ]))     )-log(seq(50, 500, 1)), 1), rep(0, 2)) , nrow = 1)
      if(( anisotropic)&( sphere)) states[[i]]$params$range_beta = matrix(c(sample(log(max(fields::rdist.earth(locs[sample(seq(vecchia_approx$n_locs), 1000), c(1, 2)   ]))/6000)-log(seq(50, 500, 1)), 1), rep(0, 2)) , nrow = 1)
      # initiate null regressors for other covariates                                                                                                                                          
      states[[i]]$params$range_beta = rbind(states[[i]]$params$range_beta, matrix(0, ncol(covariates$range_X$X_locs)+hierarchical_model$PP$n_PP*range_PP-1, length(states[[i]]$params$range_beta)))
      if(range_PP)row.names(states[[i]]$params$range_beta) = c(colnames(covariates$range_X$X_locs), paste("PP", seq(hierarchical_model$PP$n_PP), sep = "_"))
      if(!range_PP)row.names(states[[i]]$params$range_beta) = c(colnames(covariates$range_X$X_locs))
      if(range_PP)
      {
        if(!anisotropic) states[[i]]$params$range_log_var =   runif(1, hierarchical_model$range_log_var_prior[1], hierarchical_model$range_log_var_prior[2])
        if( anisotropic) states[[i]]$params$range_log_var = c(runif(3, hierarchical_model$range_log_var_prior[1], hierarchical_model$range_log_var_prior[2]), rep(0,3))
        states[[i]]$momenta$range_log_var_ancillary  = rnorm(length(states[[i]]$params$range_log_var))
        states[[i]]$momenta$range_log_var_sufficient = rnorm(length(states[[i]]$params$range_log_var))
      }
      states[[i]]$momenta$range_beta_ancillary =  matrix(rnorm(length(states[[i]]$params$range_beta)), nrow(states[[i]]$params$range_beta))
      states[[i]]$momenta$range_beta_sufficient = matrix(rnorm(length(states[[i]]$params$range_beta)), nrow(states[[i]]$params$range_beta))
        
      ###################################
      # Linear regression coefficients  #
      ###################################
      # Naive OLS to determine starting values for beta, the field, 
      naive_ols =  lm(observed_field~covariates$X$X-1)
      #starting points for regression coeffs
      perturb = t(chol(vcov(naive_ols)))%*%rnorm(length(naive_ols$coefficients))
      states[[i]]$params[["beta"]] = naive_ols$coefficients + perturb
      row.names(states[[i]]$params[["beta"]]) = colnames(covariates$X$X)
      # Residuals of the OLS model that have to be explained by the latent field and the noise
      states[[i]]$sparse_chol_and_stuff$lm_fit = as.vector(covariates$X$X%*%matrix(states[[i]]$params[["beta"]], ncol = 1))
      states[[i]]$sparse_chol_and_stuff$lm_fit_locs = as.vector(covariates$X$X_locs%*%matrix(states[[i]]$params[["beta"]][covariates$X$which_locs], ncol = 1))
      states[[i]]$sparse_chol_and_stuff$lm_residuals = as.vector(observed_field-  states[[i]]$sparse_chol_and_stuff$lm_fit)
      ##################
      # Noise variance #
      ##################
      # beta is just an intercept in stationary case
      states[[i]]$params$noise_beta    = matrix(rep(0, ncol(covariates$noise_X$X)+noise_PP*hierarchical_model$PP$n_PP), ncol = 1) #random starting values
      states[[i]]$params$noise_beta[1] = log(var(states[[i]]$sparse_chol_and_stuff$lm_residuals)) - log(2) +rnorm(1, 0, .5) # setting sensible value for the intercept
      if(!noise_PP) row.names(states[[i]]$params$noise_beta) = colnames(covariates$noise_X$X)
      if(noise_PP ) row.names(states[[i]]$params$noise_beta) = c(colnames(covariates$noise_X$X), paste("PP", seq(hierarchical_model$PP$n_PP), sep = "_"))
      states[[i]]$momenta$noise_beta = rnorm(ncol(covariates$noise_X$X)+noise_PP*hierarchical_model$PP$n_PP)
      # effective variance field, shall be used in density computations
      states[[i]]$sparse_chol_and_stuff$noise = variance_field(beta = states[[i]]$params$noise_beta, PP = PP, use_PP = noise_PP, X = covariates$noise_X$X)
      # log var if latent field
      if(noise_PP)states[[i]]$params$noise_log_var = runif(1, hierarchical_model$noise_log_var_prior[1], hierarchical_model$noise_log_var_prior[2])
      
      
      ####################
      # NNGP sparse chol #
      ####################
      states[[i]]$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = 
        compute_sparse_chol(anisotropic = anisotropic,
                            sphere = sphere, 
                            range_X = covariates$range_X$X_locs, 
                            range_beta = states[[i]]$params$range_beta, 
                            PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, 
                            NNarray = vecchia_approx$NNarray, locs_idx = vecchia_approx$hctam_scol_1, 
                            locs = locs, nu = nu, num_threads = max(1, parallel::detectCores()-2))
      states[[i]]$sparse_chol_and_stuff$sparse_chol = Matrix::sparseMatrix(x =  states[[i]]$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, triangular = T)
      states[[i]]$sparse_chol_and_stuff$precision_diag = as.vector((states[[i]]$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
      ################
      # Latent field #
      ################
      states[[i]]$params$field = exp(states[[i]]$params$log_gp_var) * as.vector(Matrix::solve(states[[i]]$sparse_chol_and_stuff$sparse_chol, rnorm(vecchia_approx$n_locs)))
    }
    
    #######################
    # Chain records setup #
    #######################
    
    # records is a list that stocks the recorded parameters of the model, including covariance parameters, the value of the sampled field, etc. In terms of RAM, those are the biggest bit !
    # iteration is a 2-colums matrix that records the iteration at the end of each chains join and the associated CPU time
    records = list()
    for(i in seq(n_chains))
    {
      records[[paste("chain", i, sep = "_")]] = list()
      records[[i]] =  sapply(states[[i]]$params, function(x)NULL)
      
    }
    iterations = list()
    iterations$checkpoints =  matrix(c(0, as.numeric(Sys.time()-t_begin, unit = "mins")), ncol = 2)
    colnames(iterations$checkpoints) = c("iteration", "time")
    iterations$thinning = c()
    
    ##########
    
    message(paste("Setup done,", as.numeric(Sys.time()- t_begin, units = "secs"), "s elapsed" ))
    return(list("data" = list("locs" = locs, "observed_field" = observed_field, "observed_locs" = observed_locs, "covariates" = covariates), "hierarchical_model" = hierarchical_model, 
                "vecchia_approx" = vecchia_approx, "states" = states, "records" = records, "t_begin" = t_begin, "seed" = seed, "iterations" = iterations))
  }  