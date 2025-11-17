#' Vecchia approximation setup
#'
#' @param observed_locs a matrix of spatial coordinates where observations are done
#' @param m number of nearest neighbors to do Vecchia's approximation
#' @param ncores number of cores to run in parallel, default to 5.
#'
#' @returns a list of Vecchia approximation objects and metadata.
#' \describe{
#'   \item{\code{n_locs}}{Number of unique locations}
#'   \item{\code{n_obs}}{Number of observed locations}
#'   \item{\code{observed_locs}}{Observed locations}
#'   \item{\code{locs}}{Unique locations}
#'   \item{\code{t_locs}}{t(locs)}  # TODO : A mon avis ça ne sert à rien de stocker ça.
#'   \item{\code{locs_match}}{locs_match}
#'   \item{\code{locs_match_matrix}}{locs_match_matrix}
#'   \item{\code{hctam_scol}}{hctam_scol}
#'   \item{\code{hctam_scol_1}}{sapply(hctam_scol, function(x) x[1])}
#'   \item{\code{obs_per_loc}}{unlist(sapply(hctam_scol, length))} # count how many observations correspond to one location
#'   \item{\code{NNarray}}{NNarray}
#'   \item{\code{NNarray_non_NA}}{!is.na(NNarray)}
#'   \item{\code{sparse_chol_i}}{sparse_mat at i+1}
#'   \item{\code{sparse_chol_p}}{sparse_mat at p}
#'   \item{\code{sparse_chol_x_reorder}}{sparse_chol_x_reorder}
#'   \item{\code{locs_partition}}{locs_partition}
#' }
#' @export
#'
#' @examples
#' set.seed(100)
#' size <- 20000
#' observed_locs <- cbind(runif(size), runif(size))
#' res <- createVecchia(observed_locs, m = 10, ncores = 1)
createVecchia <- function(observed_locs,
                          m = 12,
                          ncores = 1) {
  # message("building DAGs and indices for Vecchia approximation...")
  # Vecchia approximation ##########################################################################
  # This object gathers the NNarray table used by GpGp package and related objects
  if (!is.matrix(observed_locs))
    stop("observed_locs should be a matrix")
  if (ncol(observed_locs) != 2)
    stop("observed_locs should have 2 columns")
  
  # remove duplicates
  locs <- unique(observed_locs)
  
  # storing numbers
  n_locs <- nrow(locs)
  n_obs <- nrow(observed_locs)
  
  # re-ordering spatial locations (random then maxmin on subset)
  max_locs <- min(n_locs, 10000)
  neworder <- order(runif(n_locs))
  neworder[seq(max_locs)] <- neworder[GpGp::order_maxmin(locs[neworder[seq(max_locs)], ])]
  locs <- locs[neworder, ]
  
  # matching observed locations with reordered, unrepeated locations
  locs_match <- match(split(observed_locs, row(observed_locs)), split(locs, row(locs)))
  
  locs_match_matrix <- Matrix::sparseMatrix(i = locs_match, j = seq(n_obs), x = 1)
  
  # doing reversed operation : for a given unrepeated location, tell which observations correspond
  hctam_scol <- split(seq(n_obs), locs_match)
  
  # extracting NNarray =  nearest neighbours for Vecchia approximation
  NNarray <- t(GpGp::find_ordered_nn(locs, m))
  NNNoNA <- !is.na(NNarray)
  
  sparse_mat <- Matrix::sparseMatrix(
    x = seq(sum(NNNoNA)),
    i = col(NNarray)[NNNoNA],
    j = NNarray[NNNoNA],
    triangular = TRUE
  )
  
  sparse_chol_x_reorder <- seq_along(NNarray)[NNNoNA][match(sparse_mat@x, seq_len(sum(NNNoNA)))]
  
  # Partitioning locations using parallel kmeans for field update
  locs_partition <- generate_location_partitions(locs, n_locs, ncores = ncores)
  
  markov_mat <- Matrix::crossprod(sparse_mat)
  markov_mat@x[] <- 1
  locs_partition_coloring <- lapply(split(locs_partition, col(locs_partition)), function(x) {
    M <- Matrix::sparseMatrix(i = seq_along(x),
                              j = x,
                              x = 1)
    return(naive_greedy_coloring(Matrix::t(M) %*% markov_mat %*% M))
  })
  
  for (i in seq_len(ncol(locs_partition))) {
    locs_partition[, i] <- locs_partition_coloring[[i]][locs_partition[, i]]
  }
  
  return(
    list(
      n_locs = n_locs,
      n_obs = n_obs,
      observed_locs = observed_locs,
      locs = locs,
      t_locs = t(locs),
      # TODO : A mon avis ça ne sert à rien de stocker ça.
      locs_match = locs_match,
      locs_match_matrix = locs_match_matrix,
      hctam_scol = hctam_scol,
      hctam_scol_1 = sapply(hctam_scol, function(x)
        x[1]),
      obs_per_loc = unlist(sapply(hctam_scol, length)),
      # count how many observations correspond to one location
      NNarray = NNarray,
      NNarray_non_NA = NNNoNA,
      # TODO : idem à voir si c'est utile de trainer ça ?
      sparse_chol_i = sparse_mat@i + 1,
      sparse_chol_p = sparse_mat@p,
      sparse_chol_x_reorder = sparse_chol_x_reorder,
      locs_partition = locs_partition,
      locs_partition_coloring = locs_partition_coloring
    )
  )
}


#' Partition spatial locations using parallel k-means
#' @noRd
generate_location_partitions <- function(locs, n, ncores = 1) {
  clust_size <- 10000
  n_clusters <- ceiling(n / clust_size)
  centers_seq <- round(seq(n_clusters, max(5, 2 * (n_clusters)), length.out = min(10, max(5, 2 * (
    n_clusters
  )) - (n_clusters) + 1)))
  
  if (ncores > 1) {
    cl <- parallel::makeCluster(min(ncores, parallel::detectCores(logical = FALSE)))
    on.exit(parallel::stopCluster(cl))
    
    parallel::clusterExport(cl, varlist = c("locs", "n"), envir = environment())
    
    locs_partition <- parallel::parSapply(cl, centers_seq, function(k) {
      kmeans(
        locs,
        centers = k,
        iter.max = 200,
        algorithm = "Hartigan-Wong"
      )$cluster
    })
  } else {
    locs_partition <- sapply(centers_seq, function(k) {
      kmeans(
        locs,
        centers = k,
        iter.max = 200,
        algorithm = "Hartigan-Wong"
      )$cluster
    })
  }
  
  colnames(locs_partition) <- paste0(centers_seq, "_clust")
  return(locs_partition)
}

#' process_covariates : pre-processes covariates by adding an intercept,
#' creating useful indices, and pre-computing useful matrices and vectors
#'
#' @param X a data.frame with as many rows as vecchia_approx$observed_locs
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param PP an object of class PP,
#' @param one_obs_per_locs (logical, default to FALSE). Should the covariate be
#' constrained to not vary within a spatial location
#' (should be TRUE to range_X and scale_X)
#'
#' @returns a list
#' @export
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#' nlocs <- 5000
#' nobs <- 10000
#' unique_locs <- cbind(runif(nlocs), runif(nlocs))
#' observed_locs <- unique_locs[as.numeric(cut(runif(nobs), seq(0, 1, length.out = nlocs))), ]
#' X <- as.data.frame(cbind(runif(nobs), rnorm(nobs), rpois(nobs, 5)))
#' vecchia_approx <- createVecchia(observed_locs, 12, ncores = 1)
#' PP <- createPP(vecchia_approx, plot = FALSE)
#'
#' # Good cases ######################
#' # no PP, an X
#' res <- process_covariates(X = X, vecchia_approx, PP = NULL, one_obs_per_locs = FALSE)
#' # no PP, no X
#' res <- process_covariates(X = NULL, vecchia_approx, PP = NULL, one_obs_per_locs = FALSE)
#' # PP and X
#' res <- process_covariates(X = X, vecchia_approx, PP = PP, one_obs_per_locs = FALSE)
#' # one obs of x per loc
#' X <- as.data.frame(cbind(observed_locs, observed_locs[, 1]^2 + observed_locs[, 2]^2))
#' res <- process_covariates(X = X, vecchia_approx, PP = NULL, one_obs_per_locs = TRUE)
#'
#' # Errors ##########################
#' # one obs of x per loc
#' # X = as.data.frame(cbind(runif(nobs), rnorm(nobs), rpois(nobs, 5)))
#' # res = process_covariates(X = X, vecchia_approx, PP = NULL, one_obs_per_locs = TRUE)
#' # bad X number of rows
#' # X = X[-1,]
#' # res = process_covariates(X = X, vecchia_approx, PP = NULL, one_obs_per_locs = TRUE)
#' # bad X format
#' # X = (cbind(runif(nobs), rnorm(nobs), rpois(nobs, 5)))
#' # res = process_covariates(X = X, vecchia_approx, PP = NULL, one_obs_per_locs = TRUE)
#' # not independent covariates
#' # X = as.data.frame(cbind(observed_locs, observed_locs, observed_locs[,1]^2+ observed_locs[,2]^2))
#' # res = process_covariates(X = X, vecchia_approx, PP = NULL, one_obs_per_locs = TRUE)
process_covariates <- function(X,
                               vecchia_approx,
                               PP = NULL,
                               one_obs_per_locs = FALSE) {
  covariate_name <- as.list(match.call())$X
  if (is.null(X))
    covariate_name <- X
  res <- list()
  
  # covariates in the observed field #
  if (!is.null(X)) {
    if (!is.data.frame(X)) {
      stop(paste(covariate_name, "should be a data.frame or NULL"))
    }
    if (nrow(X) != nrow(vecchia_approx$observed_locs)) {
      stop(
        paste(
          covariate_name,
          "should have the same number of rows as the vecchia_approx$observed_locs"
        )
      )
    }
    
    # creating model matrix
    # extracting a model matrix and storing the original argument
    res$arg <- X
    res$X <- model.matrix( ~ ., X)
  }
  
  # extracting a model matrix with only intercept and storing a message about the lack of original argument if no X is provided
  if (is.null(X)) {
    res$arg <- "No covariates were provided"
    res$X <- matrix(1, vecchia_approx$n_obs)
  }
  colnames(res$X)[1] <- "(Intercept)"
  res$X <- as(res$X, "sparseMatrix")
  
  if (det(as.matrix(crossprod(res$X))) < 1e-10) {
    stop(
      covariate_name,
      " does not induce an independent collection of covariates ",
      "(det (XTX) < 1e-10). In particular, remember that the Intercept is ",
      "automatically added."
    )
  }
  
  X_ <- res$X
  if (!is.null(PP)) {
    X_ <- cbind(
      res$X,
      X_PP_mult_right(
        PP = PP,
        vecchia_approx = vecchia_approx,
        permutate_PP_to_obs = TRUE,
        Y = diag(1, PP$n_knots)
      )
    )
  }
  # pre- computing XTX
  crossprod_X <- as.matrix(crossprod(X_)) + diag(1, ncol(X_), ncol(X_))
  res$crossprod_X <- crossprod_X
  res$chol_crossprod_X <- chol(crossprod_X)
  res$n_regressors <- ncol(res$X)

  # identifying  which X do not vary within location
  res$which_locs <- c()
  duplicated_locs <- duplicated(vecchia_approx$observed_locs)
  for (i in seq_len(ncol(res$X))) {
    if (all(duplicated(cbind(
      vecchia_approx$observed_locs, res$X[, i]
    )) == duplicated_locs)) {
      res$which_locs <- c(res$which_locs, i)
    }
  }
  if (one_obs_per_locs) {
    if (!identical(res$which_locs, seq_len(ncol(res$X)))) {
      stop(
        covariate_name,
        " cannot vary within one spatial location of ",
        "vecchia_approx$duplicated_locs"
      )
    }
  }

  res$X_locs <- matrix(
    res$X[vecchia_approx$hctam_scol_1, res$which_locs], 
    ncol = length(res$which_locs)
  )
  res$X_locs <- as(res$X_locs, "sparseMatrix")

  colnames(res$X_locs) <- colnames(res$X)[res$which_locs]
  X_locs_ <- res$X_locs
  if (!is.null(PP)) {
    X_locs_ <- cbind(
      X_locs_,
      X_PP_mult_right(
        PP = PP,
        vecchia_approx = vecchia_approx,
        permutate_PP_to_obs = FALSE,
        Y = diag(1, PP$n_knots)
      )
    )
  }
  res$crossprod_X_locs <- crossprod(X_locs_) + diag(1e-10, ncol(X_locs_), ncol(X_locs_))
  res$chol_crossprod_X_locs <- chol(res$crossprod_X_locs)

  res$X_locs_ <- X_locs_
  if (one_obs_per_locs) {
    res$X <- NULL
    res$chol_crossprod_X <- NULL
  }
  return(res)
}

#' Safety checks and automatic treatment of PP objects and their marginal variance bounds
#'
#' @param PP a PP object, as given by create_PP
#' @param log_scale_bounds a numeric vector of size 2 indicating the lower and upper PP log-variance bounds
#' @param parameter_name a character string indicating the number of the parameter, used for prints
#'
#' @export
#' @examples
#' vecchia_approx <- createVecchia(cbind(runif(100), runif(100)), 10, ncores = 1)
#' myPP <- createPP(vecchia_approx)
#' PPPP <- process_PP_prior(myPP, c(1, 2), "example")
process_PP_prior <- function(PP = NULL,
                             log_scale_bounds = NULL,
                             parameter_name) {
  if (is.null(PP) && is.null(log_scale_bounds)) {
    return(log_scale_bounds)
  }
  if (is.null(PP) && !is.null(log_scale_bounds)) {
    stop(
      paste(
        "No PP object was provided for the",
        parameter_name,
        "parameters, but log - marginal variance bounds were provided"
      )
    )
  }
  if (!inherits(PP, "PP")) {
    stop(paste(
      "PP who describes the",
      parameter_name,
      "parameters must be of class PP"
    ))
  }
  if (!is.null(PP) && is.null(log_scale_bounds)) {
    message(
      paste(
        "The log - marginal variance (log_scale) bounds for the PP who describes the",
        parameter_name,
        "parameters were automatically set to (-6, 0)"
      )
    )
    log_scale_bounds <- c(-6, 0)
  }
  if (!is.numeric(log_scale_bounds) || length(log_scale_bounds) != 2)
    stop(
      paste(
        "The log - marginal variance (log_scale) bounds for the PP who describes the",
        parameter_name,
        "parameters must be a numeric vector of length 2"
      )
    )
  log_scale_bounds <- sort(log_scale_bounds)
  return(log_scale_bounds)
}

#' Initialize hierarchical model
#'
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param noise_PP either an object of class PP used to model the field of log-noise parameters, or NULL in which case the model is stationary
#' @param noise_log_scale_bounds either a vector containing two numeric values bounding Uniform prior for the log-marginal variance of the noise's PP, or NULL in which case the bounds are set automatically
#' @param range_PP TODO
#' @param range_log_scale_bounds TODO
#' @param observed_locs TODO
#' @param matern_smoothness numerical value, a Matérn smoothness parameter, 
#' either 0.5 (aka ``exponential kernel'') or 1.5
#' @param observed_field a vector of observations.
#' @param covariates The list of covariates obtained with `process_covariates`
#' @param anisotropic logical, default to FALSE. Is the covariance anisotropic ?
#' @export
#'
#' @returns a list
#' @examples
#' nobs <- 10000
#' observed_locs <- cbind(runif(nobs), runif(nobs))
#' observed_field <- rnorm(nobs)
#' vecchia_approx <- createVecchia(observed_locs, ncores = 1)
#' PP <- createPP(vecchia_approx)
#' X <- as.data.frame(cbind(observed_locs, observed_locs^2))
#' covariates <- list()
#' covariates$X <- process_covariates(X, vecchia_approx = vecchia_approx)
#'
#' hm1 <- process_hierarchical_model(
#'   vecchia_approx = vecchia_approx,
#'   noise_PP = PP, noise_log_scale_bounds = NULL,
#'   range_PP = NULL,
#'   observed_locs = observed_locs,
#'   matern_smoothness = 1.5,
#'   observed_field = observed_field,
#'   covariates = covariates,
#'   anisotropic = TRUE
#' )
process_hierarchical_model <- function(vecchia_approx,
                                       noise_PP = NULL,
                                       noise_log_scale_bounds = NULL,
                                       range_PP = NULL,
                                       range_log_scale_bounds = NULL,
                                       observed_locs,
                                       matern_smoothness,
                                       observed_field,
                                       covariates,
                                       anisotropic) {
  if (!matern_smoothness %in% c(1.5, .5)) {
    stop("only matern_smoothness = 1.5 or matern_smoothness = 0.5")
  }
  # Processing PP priors
  noise_log_scale_bounds <- process_PP_prior(noise_PP, noise_log_scale_bounds, "noise")
  range_log_scale_bounds <- process_PP_prior(range_PP, range_log_scale_bounds, "range")
  
  # Making a guess for maximum and minimum reasonable values for the range intercept
  # using as upper bound the geographic space size
  # and as lower bound the minimal distance between two space points
  alpha_max <- -0.5 * log(8 * matern_smoothness) + log(max(dist(vecchia_approx$locs[seq(min(10000, nrow(vecchia_approx$locs))), ])) / 8)
  alpha_min <- -0.5 * log(8 * matern_smoothness) + log(median(FNN::get.knn(vecchia_approx$locs, k = 1)$nn.dist) * 3)
  # OLS to get residual variance to make a guess for maximum and minimum reasonable values
  # for NNGP and noise variance
  naive_ols <- lm(observed_field ~ as.matrix(covariates$X$X) - 1)
  sigma_max <- log(.99 * var(naive_ols$residuals))
  sigma_min <- log(.01 * var(naive_ols$residuals))
  # Default mean prior computed from a reasonable case.
  res <- list(
    noise = list(
      PP = noise_PP,
      log_scale_bounds = noise_log_scale_bounds,
      beta0_mean = (sigma_max + sigma_min) / 2,
      beta0_sd = (sigma_max - sigma_min) / 8
    ),
    scale = list(
      beta0_mean = (sigma_max + sigma_min) / 2,
      beta0_sd = (sigma_max - sigma_min) / 8
    ),
    range = list(
      PP = range_PP,
      log_scale_bounds = range_log_scale_bounds,
      beta0_mean = (alpha_max + alpha_min) / 2,
      beta0_sd = (alpha_max - alpha_min) / 8
    ),
    anisotropic = anisotropic,
    matern_smoothness = matern_smoothness,
    naive_ols = naive_ols
  )
  return(res)
}


#' Initialize transition kernels
#'
#' @param init a numeric value
#' @param hm a hierarchical model (obtained with `process_hierarchical_model()`)
#' @returns a list
#' @export
#' @description
#' TODO : revoir la description
#' Transition kernel state
#' transition kernel variance is given as the log
#' can be used in both stationary and nonstationary cases respectively
#' as a random walk Metropolis or MALA step size have an ancillary and
#' a sufficient version when applicable range description
#' @examples
#' nobs <- 10000
#' observed_locs <- cbind(runif(nobs), runif(nobs))
#' observed_field <- rnorm(nobs)
#' vecchia_approx <- createVecchia(observed_locs, ncores = 1)
#' PP <- createPP(vecchia_approx)
#' X <- as.data.frame(cbind(rnorm(nobs), runif(nobs)))
#' covariates <- list()
#' covariates$X <- process_covariates(X, vecchia_approx = vecchia_approx, one_obs_per_locs = FALSE)
#' covariates$X_range <- process_covariates(X, vecchia_approx = vecchia_approx, one_obs_per_locs = TRUE)
#' covariates$X_noise <- process_covariates(X, vecchia_approx = vecchia_approx, one_obs_per_locs = FALSE)
#'
#' hm <- process_hierarchical_model(
#'   vecchia_approx = vecchia_approx,
#'   noise_PP = PP,
#'   observed_locs = observed_locs,
#'   matern_smoothness = 1.5,
#'   observed_field = observed_field,
#'   covariates = covariates,
#'   anisotropic = TRUE
#' )
#' process_transition_kernels(hm = hm)
process_transition_kernels <- function(init = -4, hm) {
  return(
    list(
      range_log_scale_sufficient = rep(init, 1 + hm$anisotropic),
      range_log_scale_ancillary = rep(init, 1 + hm$anisotropic),
      range_beta_sufficient = rep(init, 1 + hm$anisotropic),
      range_beta_ancillary = rep(init, 1 + hm$anisotropic),
      range_log_scale_estimate = rep(hm$range$log_scale_bounds[1], 1 + hm$anisotropic),
      scale_beta_sufficient = init,
      scale_beta_ancillary = init,
      scale_log_scale_sufficient = init,
      scale_log_scale_ancillary = init,
      field_log_var_ancillary = init,
      noise_beta_mala = init,
      noise_log_scale = init
    )
  )
}

#' Title
#'
#' @param hm a hierarchical model (obtained with `process_hierarchical_model()`)
#' @param covariates The list of covariates obtained with `process_covariates`
#' @param observed_field TODO
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param init_tk numerical value to initialize transition kernels
#' @param seed integer value, seed used for reproducibility purposes. Default to 1
#' @export
#' @keywords internal
#'
#' @returns a list
#' @examples
#' set.seed(100)
#' nobs <- 10000
#' observed_locs <- cbind(runif(nobs), runif(nobs))
#' observed_field <- rnorm(nobs)
#'
#' X <- as.data.frame(cbind(runif(nobs), rnorm(nobs), rpois(nobs, 5)))
#' vecchia_approx <- createVecchia(observed_locs, 12, ncores = 1)
#' PP <- createPP(vecchia_approx)
#'
#' covariates <- list(
#'   X = process_covariates(X = X, vecchia_approx),
#'   noise_X = process_covariates(X = X, vecchia_approx, PP),
#'   range_X = process_covariates(X = X, vecchia_approx)
#' )
#' hm <- process_hierarchical_model(
#'   vecchia_approx = vecchia_approx,
#'   noise_PP = PP, noise_log_scale_bounds = NULL,
#'   range_PP = PP, range_log_scale_bounds = NULL,
#'   observed_locs,
#'   matern_smoothness = 1.5,
#'   observed_field = observed_field,
#'   covariates = covariates,
#'   anisotropic = TRUE
#' )
#'
#' process_states(
#'   hm,
#'   covariates,
#'   observed_field,
#'   vecchia_approx,
#'   init_tk = -4
#' )
process_states <- function(hm,
                           covariates,
                           observed_field,
                           vecchia_approx,
                           init_tk = -4,
                           seed = 1) {
  set.seed(seed)
  # initializing sub-lisits in the state
  # Actual parameters
  params <- list()
  # Metropolis proposal distribution width for MCMC
  ker_var <- process_transition_kernels(init_tk, hm)
  # Momenta for HMC
  momenta <- list()
  # Useful stuff pre-computed from the parameters
  stuff <- list()
  
  
  # Linear regression coefficients  ############################################
  # starting points for regression coeffs
  perturb <- t(chol(vcov(hm$naive_ols))) %*%
    rnorm(length(hm$naive_ols$coefficients))
  # parameter value
  params[["beta"]] <- hm$naive_ols$coefficients + perturb
  row.names(params[["beta"]]) <- colnames(covariates$X$X)
  # Residuals of the OLS model that have to be explained by the latent field and the noise
  stuff$lm_fit <- as.vector(covariates$X$X %*% matrix(params[["beta"]], ncol = 1))
  stuff$lm_residuals <- observed_field - stuff$lm_fit
  # Range of the NNGP and stuff depending on it ################################
  # parameter format and value
  if (is.null(hm$range$PP)) {
    params$range_beta <- matrix(
      data = 0,
      nrow = ncol(covariates$range_X$X_locs),
      ncol = 1 + 2 * hm$anisotropic
    ) # random starting values
    row.names(params$range_beta) <- colnames(covariates$range_X$X_locs)
  }
  
  if(!is.null(hm$range$PP)) {
    params$range_beta = matrix(data =  0, 
                               nrow = ncol(covariates$range_X$X_locs) + hm$range$PP$n_knots, 
                               ncol = 1 + 2 * hm$anisotropic) #random starting values
    row.names(params$range_beta) = c(colnames(covariates$range_X$X_locs), 
                                     paste("PP", seq(hm$range$PP$n_knots), sep = "_"))
    params$range_log_scale = matrix(runif(1 + hm$anisotropic, hm$range$log_scale_bounds[1], hm$range$log_scale_bounds[2]))
    row.names(params$range_log_scale) = c("range", "aniso")[seq(1 + hm$anisotropic)]
    momenta$range_log_scale_sufficient = rnorm(length(params$range_log_scale))
    momenta$range_log_scale_ancillary = rnorm(length(params$range_log_scale))
  }
  params$range_beta[-1] = rnorm(length(params$range_beta)-1, 0, .2)
  colnames(params$range_beta) = c("range", "aniso1", "aniso2")[seq(1 + 2 * hm$anisotropic)]
    #row.names(params$range_beta) = c(colnames(covariates$range_X$X))
  params$range_beta[1,1] = hm$range$beta0_mean + hm$range$beta0_sd * rnorm(1)
  
  # momenta
  momenta$range_beta_ancillary <- matrix(rnorm(length(params$range_beta)), nrow(params$range_beta))
  momenta$range_beta_sufficient <- matrix(rnorm(length(params$range_beta)), nrow(params$range_beta))
  # cholesky factors of precision matrices
  stuff$compressed_chol <- compute_sparse_chol(
    range_beta = params$range_beta,
    vecchia_approx = vecchia_approx,
    range_X = covariates$range_X,
    PP = hm$range$PP,
    matern_smoothness = hm$matern_smoothness,
    compute_derivative = TRUE,
    num_threads = min(5, max(parallel::detectCores() - 1, 1))
  )
  stuff$sparse_chol <- decompress_chol(vecchia_approx, stuff$compressed_chol)
  # plot_pointillist_painting(vecchia_approx$locs, as.vector(Matrix::solve(stuff$sparse_chol, rnorm(nrow(vecchia_approx$locs)))))
  
  # Noise variance  and stuff depending on it ##################################
  # parameter format and value
  if (is.null(hm$noise$PP)) {
    # random starting values
    params$noise_beta <- 
      matrix(rep(0, ncol(covariates$noise_X$X)), ncol = 1) 
    row.names(params$noise_beta) <- colnames(covariates$noise_X$X)
  } else {
    # random starting values
    params$noise_beta <- 
      matrix(rep(0, ncol(covariates$noise_X$X) + hm$noise$PP$n_knots), ncol = 1) 
    row.names(params$noise_beta) <- 
      c(colnames(covariates$noise_X$X),
        paste("PP", seq(hm$noise$PP$n_knots), sep = "_"))
    params$noise_log_scale <- hm$noise$log_scale_bounds[1]
  }
  
  if(!is.null(hm$noise$PP)) {
    params$noise_beta = matrix(rep(0, ncol(covariates$noise_X$X) + hm$noise$PP$n_knots), ncol = 1) #random starting values
    row.names(params$noise_beta) = c(colnames(covariates$noise_X$X), 
                                     paste("PP", seq(hm$noise$PP$n_knots), sep = "_"))
    params$noise_log_scale = matrix(runif(1, hm$noise$log_scale_bounds[1], hm$noise$log_scale_bounds[2]))
    row.names(params$noise_log_scale) = " "
  }
  params$noise_beta[1] = hm$noise$beta0_mean + hm$noise$beta0_sd * rnorm(1)
  params$noise_beta[-1] = rnorm(length(params$noise_beta[-1]), 0, .2)
  
  # momenta
  momenta$noise_beta <- rnorm(length(params$noise_beta))
  # effective variance field, shall be used in density computations
  stuff$noise_var = as.vector(exp(X_PP_mult_right(
    X = covariates$noise_X$X, PP = hm$noise$PP,
    vecchia_approx = vecchia_approx, Y = params$noise_beta, 
    permutate_PP_to_obs = T
  )))
  
  # Marginal variance of the NNGP and stuff depending on it ####################
  params$field_log_var = matrix(hm$scale$beta0_mean + hm$scale$beta0_sd * rnorm(1))
  row.names(params$field_log_var) = " "
  
  # Latent field ###############################################################
  params$field = c(exp(.5 * params$field_log_var)) * 
    as.vector(Matrix::solve(stuff$sparse_chol,rnorm(vecchia_approx$n_locs)))
  
  return(list(
    "params" = params,
    # parameters of interest to the model
    "stuff" = stuff,
    # stuff that is useful for computations
    "momenta" = momenta,
    # HMC momenta
    "ker_var" = ker_var # Starting points for tk, will be adaptively tuned
  ))
}

# S3 class GeoNonStat
#' Create an object of class GeoNonStat
#'
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param observed_field a vector of observations of the interest variable
#' @param X a data.frame of covariates explaining the interest variable 
#' through fixed linear effects
#' @param matern_smoothness numerical value, a Matérn smoothness parameter, 
#' either 0.5 (aka ``exponential kernel'') or 1.5
#' @param anisotropic logical, default to FALSE. Is the covariance anisotropic ?
#' @param n_chains number of MCMC chains
#' @param range_X a data.frame of covariates explaining the Gaussian process 
#' range through fixed linear effects
#' @param range_PP TODO
#' @param range_log_scale_bounds TODO
#' @param noise_X a data.frame of covariates explaining the Gaussian noise 
#' variance through fixed linear effects
#' @param noise_PP either an object of class PP used to model the field of 
#' log-noise parameters, or NULL in which case the model is stationary
#' @param noise_log_scale_bounds either a vector containing two numeric values 
#' bounding Uniform prior for the log-marginal variance of the noise's PP, 
#' or NULL in which case the bounds are set automatically
#' @param seed integer value, seed used for reproducibility purposes. Default to 1
#'
#' @returns an object of class `GeoNonStat`
#' @export
#'
#' @examples
#' set.seed(100)
#' nobs <- 10000
#' observed_locs <- cbind(runif(5000), runif(5000))
#' observed_locs <- observed_locs[sample(seq_len(5000), nobs, replace = TRUE), ]
#' observed_field <- rnorm(nobs)
#' vecchia_approx <- createVecchia(observed_locs, ncores = 1)
#' myPP <- createPP(
#'   vecchia_approx = vecchia_approx,
#'   matern_range = .1,
#'   knots = 600
#' )
#' noise_PP <- myPP
#' noise_log_scale_bounds <- c(-8, 3)
#'
#' X <- as.data.frame(cbind(runif(nobs), rnorm(nobs), rpois(nobs, 3)))
#' range_X <- as.data.frame(observed_locs)
#' noise_X <- X
#'
#'
#' myobj <- GeoNonStat(
#'   vecchia_approx = vecchia_approx,
#'   observed_field = observed_field, # spatial locations
#'   X = X, # Response variable
#'   # Covariates per observation
#'   matern_smoothness = 1.5, # Matern smoothness
#'   anisotropic = FALSE,
#'   n_chains = 5,
#'   noise_X = noise_X, noise_PP = noise_PP, noise_log_scale_bounds = NULL,
#'   range_X = range_X,
#'   seed = 1
#' )
#' summary(myobj)
#' myobj
GeoNonStat <- function(vecchia_approx,
                       observed_field,
                       # spatial locations
                       X = NULL,
                       # Response variable
                       # Covariates per observation
                       matern_smoothness = 1.5,
                       # Matern smoothness
                       anisotropic = FALSE,
                       n_chains = 2,
                       range_X = NULL,
                       range_PP = NULL,
                       range_log_scale_bounds = NULL,
                       noise_X = NULL,
                       noise_PP = NULL,
                       noise_log_scale_bounds = NULL,
                       seed = 1) {
  # time
  t_begin <- Sys.time()
  # seed
  set.seed(seed)
  # cleansing RAM
  gc()
  
  # Sanity checks #####################################################################
  
  # covariates #########################################################
  covariates <- list()
  # fixed effects for response
  covariates$X <- process_covariates(
    X = X,
    vecchia_approx = vecchia_approx,
    PP = NULL,
    one_obs_per_locs = FALSE
  )
  # fixed effects and PP for range
  covariates$range_X <- process_covariates(
    X = range_X,
    vecchia_approx = vecchia_approx,
    PP = range_PP,
    one_obs_per_locs = TRUE
  )
  # fixed effects and PP for noise
  covariates$noise_X <- process_covariates(
    X = noise_X,
    vecchia_approx = vecchia_approx,
    PP = noise_PP,
    one_obs_per_locs = FALSE
  )
  
  # Info about hierarchical model ##############################################################
  hierarchical_model <- process_hierarchical_model(
    vecchia_approx = vecchia_approx,
    noise_PP = noise_PP,
    noise_log_scale_bounds = noise_log_scale_bounds,
    range_PP = range_PP,
    range_log_scale_bounds = range_log_scale_bounds,
    observed_locs = vecchia_approx$observed_locs,
    matern_smoothness = matern_smoothness,
    observed_field = observed_field,
    covariates = covariates,
    anisotropic = anisotropic
  )
  
  
  # Chain states #################################################################
  # cl = parallel::makeCluster(min(parallel::detectCores()-1, n_chains))
  # parallel::clusterExport(cl, c("hierarchical_model", "covariates",
  #                              "vecchia_approx", "observed_field", "seed"))
  states <-
    # parallel::parLapply(cl =  cl,
    lapply(
      seq(n_chains), function(chain_number) {
        process_states(
          seed = chain_number + seed,
          hm = hierarchical_model,
          covariates = covariates,
          observed_field = observed_field,
          vecchia_approx = vecchia_approx,
          init_tk = -3
        )
      }
    )
  names(states) <- paste("chain", seq(n_chains), sep = "_")
  # parallel::stopCluster(cl)

  # Chain records setup #######################################################
  # records is a list that stocks the recorded parameters of the model,
  # including covariance parameters, the value of the sampled field, etc. 
  # In terms of RAM, those are the biggest bit !
  records <- lapply(
    states, 
    function(x) list(x$params)
  )
  # iteration is a 2-colums matrix that records the iteration at the end of 
  # each chains join and the associated CPU time
  checkpoints <- matrix(
    c(0, as.numeric(Sys.time() - t_begin, unit = "mins")), 
    ncol = 2
  )
  colnames(checkpoints) <- c("iteration", "time")
  
  # Result ####################################################################
  res <- structure(
    list(
      "covariates" = covariates,
      "observed_field" = observed_field,
      "hierarchical_model" = hierarchical_model,
      "vecchia_approx" = vecchia_approx,
      "states" = states,
      "records" = records,
      "seed" = seed,
      "checkpoints" = checkpoints
    ),
    class = "GeoNonStat"
  )
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
#' @param ... additional arguments (unused)
#' @rdname GeoNonStat
#' @export
#' @method print GeoNonStat
print.GeoNonStat <- function(x, ...) {
  cat("Object of class 'GeoNonStat'\n")
  print(paste(length(x$data$observed_field), "observed fields"))
  print(paste(nrow(x$data$observed_locs), "observed locations"))
  print(paste("Currently", length(x$states$params$beta), "states"))
}


#' summary of a part of a GeoNonStat object
#'
#' @param partobject object to summarize
detailed_summary <- function(partobject) {
  sumdata <- summary(partobject)
  return(
    cbind(
      sumdata, 
      "Value" = 
        as.character(
          sapply(
            dimnames(sumdata)[[1]], function(x) {
              res <- ""
              if (sumdata[x, "Length"] == 1 &&
                  sumdata[x, "Mode"] %in% c("numeric", "character")) {
                res <- partobject[[x]]
              }
              res
            }
          )
        )
    )
  )
}

#' Summary of a 'GeoNonStat' object
#'
#' @param object an object of class \code{GeoNonStat}
#' @param ... additional arguments (unused)
#' @rdname GeoNonStat
#' @export
#' @method summary GeoNonStat
summary.GeoNonStat <- function(object, ...) {
  cat("### data ###")
  detailed_summary(object$data)
  cat("### hierarchical_model ###")
  detailed_summary(object$hierarchical_model)
  cat(
    "### vecchia_approx ###",
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
      round(sum(object$iterations$checkpoints[, "time"]), 5),
      "seconds"
    ),
    sep = "\n"
  )
}
