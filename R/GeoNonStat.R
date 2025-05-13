
#' Vecchia approximation setup
#'
#' @param observed_locs a matrix of spatial coordinates where observations are done
#' @param m number of nearest neighbors to do Vecchia's approximation
#'
#' @returns a list
#' @examples
#'   set.seed(100)
#'   size <- 2000
#'   observed_locs = cbind(runif(size), runif(size))
#'   res <- process_vecchia(observed_locs, m=10) 
process_vecchia <- function(observed_locs, m){
  message("building DAGs and indices for Vecchia approximation...")
  # Vecchia approximation ##########################################################################
  # This object gathers the NNarray table used by GpGp package and related objects
  
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
  #computations from process_vecchia$NNarray in order to create sparse Cholesky using Matrix::sparseMatrix
  #non_NA indices from process_vecchia$NNarray
  NNarray_non_NA = !is.na(NNarray)
  
  # partition of locs for field update
  cl = parallel::makeCluster(5)
  parallel::clusterExport(cl = cl, varlist = c("locs", "n"), envir = environment())
  locs_partition = parallel::parSapply(
    cl = cl, 
    round(seq(n / 10000 + 1, 2 * (n / 10000) + 1, length.out = min(10, ceiling(n / 10000)))), 
    function(i) { kmeans(locs, centers = i, iter.max = 200, algorithm = "Hartigan-Wong")$cluster}
  )
  colnames(locs_partition) = paste(round(seq(n / 10000 + 1, 2 * (n / 10000) + 1, length.out = min(20, ceiling(n / 10000)))), "clust")
  parallel::stopCluster(cl)
  
  return(list(
    n_locs = n_locs,
    n_obs = n_obs,
    locs = locs, 
    locs_match = locs_match,
    locs_match_matrix = locs_match_matrix,
    hctam_scol = hctam_scol,
    hctam_scol_1 = sapply(hctam_scol, function(x) x[1]),
    obs_per_loc = unlist(sapply(hctam_scol, length)), # count how many observations correspond to one location
    NNarray = NNarray,
    NNarray_non_NA = !is.na(NNarray),
    sparse_chol_column_idx = NNarray[NNarray_non_NA], #column idx of the uncompressed sparse Cholesky factor
    sparse_chol_row_idx = row(NNarray)[NNarray_non_NA], #row idx of the uncompressed sparse Cholesky factor
    locs_partition = locs_partition
  ))
}

#' Initialize hierarchical model
#'
#' @param PP TODO
#' @param noisePP TODO
#' @param scale_PP TODO
#' @param range_PP TODO
#' @param noise_log_scale_prior TODO
#' @param scale_log_scale_prior TODO
#' @param range_log_scale_prior TODO
#' @param locs TODO
#' @param nu TODO
#' @param observed_field TODO
#' @param covariates TODO
#' @param anisotropic TODO
#'
#' @returns a list
#' @examples
#' #TODO
process_hierarchical_model <- function(PP,
                                       noise_PP = NULL, noise_log_scale_prior,
                                       scale_PP = NULL, scale_log_scale_prior,
                                       range_PP = NULL, range_log_scale_prior,
                                       locs,
                                       nu,
                                       observed_field,
                                       covariates,
                                       anisotropic) {
  
  # Info about hierarchical model ##############################################################
  
  if (is.null(noise_log_scale_prior) & !is.null(noise_PP))
  {
    message("the prior for the marginal variance of the PP who describes the noise parameters was automatically set to an uniform on (-6, 3)")
    noise_log_scale_prior = c(-6, 3)
  }
  if (is.null(scale_log_scale_prior) & scale_PP)
  {
    message("the prior for the marginal variance of the PP who describes the scale parameters was automatically set to an uniform on (-6, 3)")
    scale_log_scale_prior = c(-6, 3)
  }
  if (is.null(range_log_scale_prior) & range_PP)
  {
    message("the prior for the marginal variance of the PP who describes the range parameters was automatically set to an uniform on (-6, 3)")
    range_log_scale_prior = c(-6, 3)
  }
  
  if (!is.null(noise_log_scale_prior))
    noise_log_scale_prior = matrix(noise_log_scale_prior)
  if (!is.null(scale_log_scale_prior))
    scale_log_scale_prior = matrix(scale_log_scale_prior)
  
  alpha_max = -0.5 * log(8 * nu) + log(max(dist(locs[seq(min(10000, nrow(locs))), ])) / 4)
  alpha_min = -0.5 * log(8 * nu) + log(median(FNN::get.knn(locs, k = 1)$nn.dist) * 3)
  
  # OLS to get residual variance to make a guess
  naive_ols =  lm(observed_field ~ covariates$X$X - 1)
  sigma_max = log(var(naive_ols$residuals))
  sigma_min = log(var(naive_ols$residuals) / 1000)
  # Default mean prior computed from a reasonable case.
  res <- list(
    anisotropic= anisotropic,
    nu= nu,
    beta_priors= list(),
    PP= PP,
    noise_PP= noise_PP,
    scale_PP= scale_PP,
    range_PP= range_PP,
    noise_log_scale_prior = noise_log_scale_prior,
    scale_log_scale_prior = scale_log_scale_prior,
    range_log_scale_prior = range_log_scale_prior,
    range_beta0_mean= (alpha_max + alpha_min) / 2,
    range_beta0_var= ((alpha_max - alpha_min) / 8)^2,
    noise_beta0_mean= (sigma_max + sigma_min) / 2,
    noise_beta0_var= ((sigma_max - sigma_min) / 8)^2,
    scale_beta0_mean= (sigma_max + sigma_min) / 2,
    scale_beta0_var=((sigma_max - sigma_min) / 8)^2,
    naive_ols = naive_ols
  )
  if(!noise_PP) res[["noise_log_scale_prior"]] <- NULL
  if(!scale_PP) res[["scale_log_scale_prior"]] <- NULL
  if(!range_PP) res[["range_log_scale_prior"]] <- NULL
  return(res)
}


#' Initialize transition kernels
#'
#' @param init a numeric value
#'
#' @returns a list
#'
#' @examples
#' tk <- process_transition_kernels()
process_transition_kernels <- function(init=-4){
  # Transition kernels ###################################################################
  # Transition kernel state
  # transition kernel variance is given as the log
  # can be used in both stationary and nonstationary cases respectively as a random walk Metropolis or MALA step size
  # have an ancillary and a sufficient version when applicable
  # range
  return(list(
    range_log_scale_sufficient = init,
    range_log_scale_ancillary =  init,
    range_beta_sufficient = c(init, init),
    range_beta_ancillary  = c(init, init),
    
    scale_beta_sufficient_mala = init,
    scale_beta_ancillary_mala  = init,
    scale_log_scale_sufficient = init,
    scale_log_scale_ancillary =  init,
    
    noise_beta_mala = init,
    noise_log_scale = init
  ))
}

#' Title
#'
#' @param hm hierarchical model (from `process_hierarchical_model()`)
#' @param covariates TODO
#' @param observed_field TODO
#' @param anisotropic TODO
#' @param vecchia_setup TODO
#' @param init_tk TODO
#'
#' @returns a list
process_states <- function(
    hm,
    covariates,
    observed_field,
    anisotropic,
    vecchia_setup,
    init_tk = -4
) {
  # Transition kernels ###################################################################
  # Starting points for transition kernels, will be adaptively tuned
  # Transition kernel state
  # transition kernel variance is given as the log
  # can be used in both stationary and nonstationary cases respectively as a random walk Metropolis or MALA step size
  # have an ancillary and a sufficient version when applicable
  # range
  transition_kernels <- process_transition_kernels(init_tk)
  
  # Linear regression coefficients  ###################################################################
  #starting points for regression coeffs
  perturb = t(chol(vcov(hm$naive_ols))) %*% 
    rnorm(length(hm$naive_ols$coefficients))
  
  params = list("beta"= NULL, 
                "range_beta" = NULL, 
                "range_log_scale" = NULL,
                "scale_beta" = NULL,
                "field" = NULL)
  
  params[["beta"]] = hm$naive_ols$coefficients + perturb
  row.names(params[["beta"]]) = colnames(covariates$X$X)
  
  sparse_chol_and_stuff <- list()
  # Residuals of the OLS model that have to be explained by the latent field and the noise
  sparse_chol_and_stuff$lm_fit = as.vector(covariates$X$X %*% matrix(params[["beta"]], ncol = 1))
  sparse_chol_and_stuff$lm_fit_locs = as.vector(covariates$X$X_locs %*%
                                                        matrix(params[["beta"]][covariates$X$which_locs], ncol = 1))
  sparse_chol_and_stuff$lm_residuals = as.vector(observed_field - sparse_chol_and_stuff$lm_fit)
  
  # Range #######################################################################
  # range beta
  params$range_beta = matrix(0,
                             ncol(covariates$range_X$X) + 
                               hm$range_PP * hm$PP$n_knots,
                             1 + 2 * anisotropic)
  if (!hm$range_PP)
    row.names(params$range_beta) = colnames(covariates$range_X$X)
  if (hm$range_PP)
    row.names(params$range_beta) = c(colnames(covariates$range_X$X), 
                                     paste("PP", seq(hm$PP$n_knots), sep = "_"))
  params$range_beta[1, 1] = hm$range_beta0_mean
  
  momenta <- list()
  momenta$range_beta_ancillary = matrix(rnorm(length(params$range_beta)), nrow(params$range_beta))
  momenta$range_beta_sufficient = matrix(rnorm(length(params$range_beta)), nrow(params$range_beta))
  # range log scale
  if (hm$range_PP) {
    if (!anisotropic)
      params$range_log_scale =   hm$range_log_scale_prior[1]
    if (anisotropic)
      params$range_log_scale = c(rep(hm$range_log_scale_prior[1], 3),
                                       rep(0, 3))
    momenta$range_log_scale_ancillary  = rnorm(length(params$range_log_scale))
    momenta$range_log_scale_sufficient = rnorm(length(params$range_log_scale))
  }
  
  # Noise variance #######################################################
  
  # beta is just an intercept in stationary case
  params$noise_beta = matrix(rep(0, ncol(covariates$noise_X$X) + 
                                   hm$noise_PP * hm$PP$n_knots), ncol = 1) #random starting values
  if (!hm$noise_PP)
    row.names(params$noise_beta) = colnames(covariates$noise_X$X)
  if (hm$noise_PP)
    row.names(params$noise_beta) = c(colnames(covariates$noise_X$X), 
                                     paste("PP", seq(hm$PP$n_knots), sep = "_"))
  params$noise_beta[1] = hm$noise_beta0_mean
  momenta$noise_beta = rnorm(ncol(covariates$noise_X$X) + hm$noise_PP *hm$PP$n_knots)
  # noise log scale
  if (hm$noise_PP)
    params$noise_log_scale = hm$noise_log_scale_prior[1]
  # effective variance field, shall be used in density computations
  sparse_chol_and_stuff$noise = variance_field(
    beta = params$noise_beta,
    PP = hm$PP,
    X = covariates$noise_X$X
  )
  
  # Scale ################################################################
  params$scale_beta    = matrix(0,
                                      ncol(covariates$scale_X$X_locs) + hm$scale_PP * hm$PP$n_PP,
                                      ncol = 1)
  if (!hm$scale_PP)
    row.names(params$scale_beta) = colnames(covariates$scale_X$X_locs)
  if (hm$scale_PP)
    row.names(params$scale_beta) = c(colnames(covariates$scale_X$X_locs),
                                           paste("PP", seq(hm$PP$n_knots), sep = "_"))
  params$scale_beta[1] = hm$scale_beta0_mean
  momenta$scale_beta_ancillary =  rnorm(ncol(covariates$scale_X$X_locs) + 
                                          hm$scale_PP * hm$PP$n_knots)
  momenta$scale_beta_sufficient = rnorm(ncol(covariates$scale_X$X_locs) + 
                                          hm$scale_PP * hm$PP$n_knots)
  
  # variance
  if (hm$scale_PP)
    params$scale_log_scale = hm$scale_log_scale_prior[1]
  # effective variance field, shall be used in density computations
  sparse_chol_and_stuff$scale = variance_field(
    beta = params$scale_beta,
    PP = hm$PP,
    X = covariates$scale_X$X_locs,
    locs_idx = vecchia_setup$hctam_scol_1
  )
  
  # NNGP sparse chol #############################################
  sparse_chol_and_stuff$compressed_sparse_chol_and_grad =
    compute_sparse_chol(
      anisotropic = anisotropic,
      range_X = covariates$range_X$X_locs,
      range_beta = params$range_beta,
      PP = hm$PP,
      use_PP = hm$range_PP,
      NNarray = vecchia_setup$NNarray,
      locs_idx = vecchia_setup$hctam_scol_1,
      locs = vecchia_setup$locs,
      smoothness = hm$nu,
      num_threads = max(1, parallel::detectCores() - 2)
    )
  
  sparse_chol_and_stuff$sparse_chol = Matrix::sparseMatrix(
    x =  sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_setup$NNarray_non_NA],
    i = vecchia_setup$sparse_chol_row_idx,
    j = vecchia_setup$sparse_chol_column_idx,
    triangular = T
  )
  sparse_chol_and_stuff$precision_diag = as.vector((
    sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_setup$NNarray_non_NA]^2
  ) %*% Matrix::sparseMatrix(
    i = seq(length(
      vecchia_setup$sparse_chol_column_idx
    )),
    j = vecchia_setup$sparse_chol_column_idx,
    x = rep(1, length(vecchia_setup$sparse_chol_row_idx))
  )
  )
  
  # Latent field ################################################################
  
  params$field = sqrt(sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(
    sparse_chol_and_stuff$sparse_chol,
    rnorm(vecchia_setup$n_locs)
  ))
  
  return(list(
    "params" = params, #parameters of interest to the model
    "sparse_chol_and_stuff" =sparse_chol_and_stuff, #stuff that is useful for computations
    "momenta" = momenta, # HMC momenta
    "transition_kernels" = transition_kernels # Starting points for transition kernels, will be adaptively tuned
  ))
}

# S3 class GeoNonStat
#' Create an object of class GeoNonStat
#' @param observed_locs a matrix of spatial coordinates where observations are done
#' @param observed_field a vector of observations of the interest variable
#' @param X a data.frame of covariates explaining the interest variable through fixed linear effects
#' @param m number of nearest neighbors to do Vecchia's approximation
#' @param nu Matern smoothness, either 0.5 or 1.5
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
#'   nu = 1.5, n_chains = 5,
#'   range_PP = T, PP = myPP, # use PP for range
#'   anisotropic = T # Covariance will be anisotropic
#' )
GeoNonStat <- 
  function(observed_locs,
           observed_field,  #spatial locations
           X = NULL, # Response variable 
           # Covariates per observation
           m = 10, #number of Nearest Neighbors
           nu = 1.5, #Matern smoothness
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
    if (!is.matrix(observed_locs)) stop("observed_locs should be a matrix")
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
    if (!nu %in% c(1.5, .5))
      stop("only nu = 1.5 or nu = 0.5")
    
    # Re-ordering ############################################################################
    # remove duplicated locations
    # Vecchia approximation ##########################################################################
    # This object gathers the NNarray table used by GpGp package and related objects

    vecchia_setup <- process_vecchia(observed_locs, m)
    
    # covariates #########################################################
    
    covariates = list()
    # fixed effects for response
    covariates$X = process_covariates(X, observed_locs, vecchia_setup)
    # explicit PP basis
    explicit_PP_basis = NULL
    if (!is.null(PP)) {
      explicit_PP_basis = X_PP_mult_right(PP = PP,
                                          Y = diag(1, nrow(PP$knots), nrow(PP$knots)))
    }
    # fixed effects and PP for range
    covariates$range_X = process_covariates(range_X,
                                            observed_locs,
                                            vecchia_setup,
                                            explicit_PP_basis,
                                            range_PP)
    if (!identical(covariates$range_X$which_locs, seq(ncol(covariates$range_X$X_locs))))
      stop("The covariates range_X cannot vary within one spatial location of observed_locs")
    # fixed effects and PP for scale
    covariates$scale_X = process_covariates(scale_X,
                                            observed_locs,
                                            vecchia_setup,
                                            explicit_PP_basis,
                                            scale_PP)
    if (!identical(covariates$scale_X$which_locs, seq(ncol(covariates$scale_X$X))))
      stop("The covariates scale_X cannot vary within one spatial location of observed_locs")
    # fixed effects and PP for noise
    covariates$noise_X = process_covariates(noise_X,
                                            observed_locs,
                                            vecchia_setup,
                                            explicit_PP_basis,
                                            noise_PP)
    # explicit PP basis removal
    remove(explicit_PP_basis)
    
    # Info about hierarchical model ##############################################################
    hierarchical_model <- process_hierarchical_model(
      PP,
      noise_PP, noise_log_scale_prior,
      scale_PP, scale_log_scale_prior,
      range_PP, range_log_scale_prior,
      locs=vecchia_setup$locs,
      nu = nu,
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
      vecchia_setup = vecchia_setup,
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
      "vecchia_setup" = vecchia_setup,
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
  cat("### vecchia_setup ###",
    summary(object$vecchia_setup), 
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