#' initialize TODO Title

#' @param observed_locs a matrix of spatial coordinates where observations are done 
#' @param observed_field a vector of observations of the interest variable
#' @param X a data.frame of covariates explaining the interest variable through fixed linear effects
#' @param m number of nearest neighbors to do Vecchia's approximation
#' @param nu Matern smoothness, either 0.5 or 1.5
#' @param anisotropic anisotropic covariance
#' @param sphere Boolean, indicating lon-lat data
#' @param PP TODO
#' @param n_chains TODO
#' @param noise_PP TODO
#' @param noise_X a data.frame of covariates explaining the Gaussian noise variance through fixed linear effects
#' @param noise_log_scale_prior  1 times 2 matrix for the prior on the log-variance of the noise PP field. 
#' @param scale_PP 
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
#' locs = cbind(runif(10000), runif(10000))
#' PP = createPP(
#'   observed_locs = locs[seq(10000),], # spatial sites
#'   matern_range = .1,
#'   knots = 50, # number of knots
#'   m = 15 # number of NNGP parents
#' )
#' range_beta = rbind(c(-4, 0, 0), # intercept
#'                    matrix(.5*rnorm(150), 50)
#'                    ) %*% diag(c(1, 2,.5))
#' NNarray_aniso = GpGp::find_ordered_nn(locs[seq(10000),], 10)
#' 
#' # getting the coefficients
#' chol_precision = compute_sparse_chol(
#'     range_beta = range_beta, # range parameters
#'     NNarray = NNarray_aniso, # Vecchia approx DAG
#'     locs = locs[seq(10000),], # spatial coordinates
#'     range_X = matrix(1, 10000), # covariates for the range (just an intercept)
#'     PP = PP, use_PP = T, # predictive process
#'     nu = 1.5, # smoothness
#'     anisotropic = T # anisotropy
#'   )[[1]]
#' # putting coefficients in precision Cholesly
#' chol_precision = Matrix::sparseMatrix(
#'   x = chol_precision[!is.na(NNarray_aniso)], 
#'   i = row(NNarray_aniso)[!is.na(NNarray_aniso)], 
#'   j = (NNarray_aniso)[!is.na(NNarray_aniso)], 
#'   triangular = T
#' )
#' # sampling the anisotropic process
#' seed_vector = rnorm(10000)
#' aniso_latent_field = as.vector(Matrix::solve(chol_precision, seed_vector)) 
#' aniso_observed_field = aniso_latent_field + .8*rnorm(10000)
#' #plot_pointillist_painting(locs, aniso_latent_field)
#' MCMC_NNGP = initialize(
#'   observed_locs = locs[seq(10000),], observed_field = aniso_observed_field, 
#'   nu = 1.5, n_chains = 5,
#'   range_PP = T, PP = PP, # use PP for range
#'   anisotropic = T # Covariance will be anisotropic
#' )
initialize = 
  function(observed_locs = NULL, #spatial locations
           observed_field = NULL, # Response variable
           X = NULL, # Covariates per observation
           m = 10, #number of Nearest Neighbors
           nu =1.5, #Matern smoothness
           anisotropic = F, 
           sphere = F, 
           PP = NULL, 
           n_chains = 2,  # number of MCMC chains
           noise_PP = F, noise_X = NULL, noise_log_scale_prior = NULL, 
           scale_PP = F, scale_X = NULL, scale_log_scale_prior = NULL, 
           range_PP = F, range_X = NULL, range_log_scale_prior = NULL, 
           seed = 1
  )
  {
    
    # time  
    t_begin = Sys.time()
    # seed
    set.seed(seed)
    # cleansing RAM
    gc()

    # Sanity checks #####################################################################
    # format
    if(!is.matrix(observed_locs))stop("observed_locs should be a matrix")
    if(!is.vector(observed_field))stop("observed_field should be a vector")
    
    if(!is.data.frame(X) & !is.null(X))stop("X should be a data.frame or NULL")
    
    if(!is.data.frame(noise_X) & !is.null(noise_X))stop("noise_X should be a data.frame or NULL")
    if(!is.data.frame(scale_X) & !is.null(scale_X))stop("scale_X should be a data.frame or NULL")
    if(!is.data.frame(range_X) & !is.null(range_X))stop("range_X should be a data.frame or NULL")
    
    if((is.null(PP)) & (noise_PP | range_PP | scale_PP))stop("either noise_PP, range_PP, or scale_PP is TRUE, while nothing was provided for PP")
    #length of observations
    if(
      !all(unique(c(
        length(observed_field),
        nrow(observed_locs),
        nrow(X),
        nrow(scale_X),
        nrow(noise_X),
        nrow(range_X),
        length(PP$idx)
      )) %in% c(0, length(observed_field))
      )) stop(
        paste("Lengths are not matching : observed_field has", length(observed_field), "observations,",
              "observed_locs has", nrow(observed_locs), "rows,", 
              "X has", nrow(X), "rows,", 
              "scale_X has", nrow(scale_X), "rows (can only be either 0 or the length of the observations),", 
              "noise_X has", nrow(noise_X), "rows (can only be either 0 or the length of the observations),", 
              "range_X has", nrow(range_X), "rows (can only be either 0 or the length of the observations),", 
              "PP has", length(PP$idx), "locations (can only be either 0 or the length of the observations)"
        )
      )
    
    # smoothness
    if (!nu %in% c(1.5, .5))stop("only nu = 1.5 or nu = 0.5")
    
    
    # Re-ordering ############################################################################
    # remove duplicated locations
    duplicated_locs = duplicated (observed_locs)
    locs = observed_locs[duplicated_locs==F,]
    locs_reordering = order(runif(nrow(locs))); locs_reordering[seq(min(nrow(locs), 100000))] = locs_reordering[GpGp::order_maxmin(locs[locs_reordering[seq(min(nrow(locs), 100000))],])]
    locs = locs[locs_reordering,]
    # extracting number of locations as shortcut
    n = nrow(locs)
    

    # Vecchia approximation ##########################################################################
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
    
    # covariates #########################################################
    
    covariates = list()
    # fixed effects for response
    covariates$X = process_covariates(X, observed_locs, vecchia_approx)  
    # explicit PP basis
    explicit_PP_basis = NULL
    if(!is.null(PP))explicit_PP_basis = X_PP_mult_right(PP = PP, Y = diag(1, nrow(PP$knots), nrow(PP$knots)))
    # fixed effects and PP for range
    covariates$range_X = process_covariates(range_X, observed_locs, vecchia_approx, explicit_PP_basis, range_PP)
    if(!identical(covariates$range_X$which_locs, seq(ncol(covariates$range_X$X_locs))))stop("The covariates range_X cannot vary within one spatial location of observed_locs")
    # fixed effects and PP for scale
    covariates$scale_X = process_covariates(scale_X, observed_locs, vecchia_approx, explicit_PP_basis, scale_PP)
    if(!identical(covariates$scale_X$which_locs, seq(ncol(covariates$scale_X$X))))stop("The covariates scale_X cannot vary within one spatial location of observed_locs")
    # fixed effects and PP for noise
    covariates$noise_X = process_covariates(noise_X, observed_locs, vecchia_approx, explicit_PP_basis, noise_PP)
    # explicit PP basis removal
    remove(explicit_PP_basis)
    
    # Info about hierarchical model ##############################################################
    
    hierarchical_model = list()
    hierarchical_model$anisotropic = anisotropic
    hierarchical_model$sphere = sphere
    hierarchical_model$nu = nu
    hierarchical_model$beta_priors = list()
    hierarchical_model$PP = PP
    if(is.null(hierarchical_model$PP)) hierarchical_model$PP = list("n_PP" = 0)
    hierarchical_model$noise_PP = noise_PP
    hierarchical_model$scale_PP = scale_PP
    hierarchical_model$range_PP = range_PP
    
    if(is.null(noise_log_scale_prior)&noise_PP)
    {
      message("noise_log_scale_prior was automatically set to an uniform on (-6, 2)")
      noise_log_scale_prior = c(-6, 2)
    }
    if(!is.null(noise_log_scale_prior))hierarchical_model$noise_log_scale_prior = matrix(noise_log_scale_prior)
    if(is.null(scale_log_scale_prior)&scale_PP)
    {
      message("scale_log_scale_prior was automatically set to an uniform on (-6, 2)")
      scale_log_scale_prior = c(-6, 2)
    }
    if(!is.null(scale_log_scale_prior))hierarchical_model$scale_log_scale_prior = matrix(scale_log_scale_prior)
    if(is.null(range_log_scale_prior)&range_PP)
    {
      message("range_log_scale_prior was automatically set to an uniform on (-6, 2)")
      hierarchical_model$range_log_scale_prior = c(-6, 3)
    }
    
    # OLS to get residual variance to make a guess 
    naive_ols =  lm(observed_field~covariates$X$X-1)
    
    alpha_max = - 0.5 * log(8*nu) + log(max(dist(locs[seq(1000),]))/4)
    alpha_min = - 0.5 * log(8*nu) + log(median(FNN::get.knn(locs, k = 1)$nn.dist)*3)
    hierarchical_model$range_beta0_mean = (alpha_max + alpha_min)/2
    hierarchical_model$range_beta0_var = ((alpha_max - alpha_min)/8)^2
    sigma_max = log(var(naive_ols$residuals))
    sigma_min = log(var(naive_ols$residuals)/1000)
    hierarchical_model$noise_beta0_mean = (sigma_max + sigma_min)/2
    hierarchical_model$noise_beta0_var  = ((sigma_max - sigma_min)/8)^2
    hierarchical_model$scale_beta0_mean = (sigma_max + sigma_min)/2
    hierarchical_model$scale_beta0_var  = ((sigma_max - sigma_min)/8)^2
    # Default mean prior computed from a reasonable case. 
    
    
    # Chain states #################################################################
    # for each chain, creating sub-lists in order to stock all the stuff that is related to one chain, including : 
    # transition_kernel_sd : a list that stocks the (current) automatically-tuned transition kernels standard deviations
    # params : a list that stocks the (current) parameters of the model, including covariance parameters, the value of the sampled field, etc
    # record : record of the MCMC iterations
    
    state = list()
    #parameters of interest to the model
    state$params = list()
    #stuff that is useful for computations
    state$sparse_chol_and_stuff = list()
    # HMC momenta
    state$momenta = list()
    # Starting points for transition kernels, will be adaptively tuned
    state$transition_kernels = list()
    
    
    # Transition kernels ###################################################################  
    # Transition kernel state
    # transition kernel variance is given as the log
    # can be used in both stationary and nonstationary cases respectively as a random walk Metropolis or MALA step size 
    # have an ancillary and a sufficient version when applicable
    # range
    state$transition_kernels$range_log_scale_sufficient = -4
    state$transition_kernels$range_log_scale_ancillary =  -4
    state$transition_kernels$range_beta_sufficient = c(-4, -4)
    state$transition_kernels$range_beta_ancillary  = c(-4, -4)
    # scale
    state$transition_kernels$scale_beta_sufficient_mala = -4
    state$transition_kernels$scale_beta_ancillary_mala  = -4
    state$transition_kernels$scale_log_scale_sufficient = -4
    state$transition_kernels$scale_log_scale_ancillary =  -4
    # noise variance
    state$transition_kernels$noise_beta_mala = -4
    state$transition_kernels$noise_log_scale = -4
    
    # Linear regression coefficients  ###################################################################
    #starting points for regression coeffs
    perturb = t(chol(vcov(naive_ols)))%*%rnorm(length(naive_ols$coefficients))
    state$params[["beta"]] = naive_ols$coefficients + perturb
    row.names(state$params[["beta"]]) = colnames(covariates$X$X)
    # Residuals of the OLS model that have to be explained by the latent field and the noise
    state$sparse_chol_and_stuff$lm_fit = as.vector(covariates$X$X%*%matrix(state$params[["beta"]], ncol = 1))
    state$sparse_chol_and_stuff$lm_fit_locs = as.vector(covariates$X$X_locs%*%matrix(state$params[["beta"]][covariates$X$which_locs], ncol = 1))
    state$sparse_chol_and_stuff$lm_residuals = as.vector(observed_field-  state$sparse_chol_and_stuff$lm_fit)
    
    # Range #######################################################################
    
    # range beta
    state$params$range_beta = matrix(0, ncol(covariates$range_X$X) + range_PP * PP$n_PP, 1 + 2*anisotropic)
    if(!range_PP) row.names(state$params$range_beta) = colnames(covariates$range_X$X)
    if(range_PP ) row.names(state$params$range_beta) = c(colnames(covariates$range_X$X), paste("PP", seq(hierarchical_model$PP$n_PP), sep = "_"))
    state$params$range_beta[1,1] = hierarchical_model$range_beta0_mean
    state$momenta$range_beta_ancillary = matrix(rnorm(length(state$params$range_beta)), nrow(state$params$range_beta))
    state$momenta$range_beta_sufficient = matrix(rnorm(length(state$params$range_beta)), nrow(state$params$range_beta))
    # range log scale
    if(range_PP)
    {
      if(!anisotropic) state$params$range_log_scale =   hierarchical_model$range_log_scale_prior[1]
      if( anisotropic) state$params$range_log_scale = c(rep(hierarchical_model$range_log_scale_prior[1], 3), rep(0, 3))
      state$momenta$range_log_scale_ancillary  = rnorm(length(state$params$range_log_scale))
      state$momenta$range_log_scale_sufficient = rnorm(length(state$params$range_log_scale))
    }
    
    # Noise variance #######################################################
    
    # beta is just an intercept in stationary case
    state$params$noise_beta    = matrix(rep(0, ncol(covariates$noise_X$X)+noise_PP*hierarchical_model$PP$n_PP), ncol = 1) #random starting values
    if(!noise_PP) row.names(state$params$noise_beta) = colnames(covariates$noise_X$X)
    if(noise_PP ) row.names(state$params$noise_beta) = c(colnames(covariates$noise_X$X), paste("PP", seq(hierarchical_model$PP$n_PP), sep = "_"))
    state$params$noise_beta[1] = hierarchical_model$noise_beta0_mean
    state$momenta$noise_beta = rnorm(ncol(covariates$noise_X$X)+noise_PP*hierarchical_model$PP$n_PP)
    # noise log scale
    if(noise_PP)state$params$noise_log_scale = hierarchical_model$noise_log_scale_prior[1]
    # effective variance field, shall be used in density computations
    state$sparse_chol_and_stuff$noise = variance_field(beta = state$params$noise_beta, PP = PP, use_PP = noise_PP, X = covariates$noise_X$X)
    
    # Scale ################################################################
    
    state$params$scale_beta    = matrix(0, ncol(covariates$scale_X$X_locs) + scale_PP * hierarchical_model$PP$n_PP, ncol = 1) 
    if(!scale_PP)row.names(state$params$scale_beta) = colnames(covariates$scale_X$X_locs)
    if(scale_PP)row.names(state$params$scale_beta) = c(colnames(covariates$scale_X$X_locs),  paste("PP", seq(hierarchical_model$PP$n_PP), sep = "_"))
    state$params$scale_beta[1] = hierarchical_model$scale_beta0_mean
    state$momenta$scale_beta_ancillary =  rnorm(ncol(covariates$scale_X$X_locs) + scale_PP * hierarchical_model$PP$n_PP)
    state$momenta$scale_beta_sufficient = rnorm(ncol(covariates$scale_X$X_locs) + scale_PP * hierarchical_model$PP$n_PP)
    # variance
    if(scale_PP)state$params$scale_log_scale = hierarchical_model$scale_log_scale_prior[1]
    # effective variance field, shall be used in density computations
    state$sparse_chol_and_stuff$scale = variance_field(beta = state$params$scale_beta, PP = PP, use_PP = scale_PP, X = covariates$scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
    
    # NNGP sparse chol #############################################
    
    state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = 
      compute_sparse_chol(anisotropic = anisotropic,
                                  sphere = sphere, 
                                  range_X = covariates$range_X$X_locs, 
                                  range_beta = state$params$range_beta, 
                                  PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, 
                                  NNarray = vecchia_approx$NNarray, locs_idx = vecchia_approx$hctam_scol_1, 
                                  locs = locs, nu = nu, num_threads = max(1, parallel::detectCores()-2))
    state$sparse_chol_and_stuff$sparse_chol = Matrix::sparseMatrix(x =  state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, triangular = T)
    state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
    
    # Latent field ################################################################
    
    state$params$field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(state$sparse_chol_and_stuff$sparse_chol, rnorm(vecchia_approx$n_locs)))
    
    
    # Chain records setup #########################################################################
    
    # records is a list that stocks the recorded parameters of the model, including covariance parameters, the value of the sampled field, etc. In terms of RAM, those are the biggest bit !
    # iteration is a 2-colums matrix that records the iteration at the end of each chains join and the associated CPU time
    records = list()
    records =  sapply(state$params, function(x)NULL)
    iterations = list()
    iterations$checkpoints =  matrix(c(0, as.numeric(Sys.time()-t_begin, unit = "mins")), ncol = 2)
    colnames(iterations$checkpoints) = c("iteration", "time")
    iterations$thinning = c()
    
    # Result ####################################################################
    res = list("data" = 
                 list("locs" = locs, "observed_field" = observed_field, "observed_locs" = observed_locs, "covariates" = covariates), 
               "hierarchical_model" = hierarchical_model, 
               "vecchia_approx" = vecchia_approx, 
               "states" = state, "records" = records, 
               "t_begin" = t_begin, "seed" = seed, "iterations" = iterations)
    class(res) = "MCMC_NNGP"
    message(paste("Setup done,", as.numeric(Sys.time()- t_begin, units = "secs"), "s elapsed" ))
    return(res)
  }  