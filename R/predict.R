extendVecchia <- function(vecchia_approx, new_locs) {
  observed_locs <- rbind(vecchia_approx$observed_locs, new_locs)
  locs <- unique(rbind(vecchia_approx$locs, new_locs))
  n_obs  <- nrow(observed_locs)
  new_n_obs  <- nrow(new_locs)
  n_locs <- nrow(locs)

  new_locs_unique <- locs[-seq_len(vecchia_approx$n_locs), , drop = FALSE]

  locs_match <- match(
    split(observed_locs, row(observed_locs)),
    split(locs, row(locs))
  )

  new_locs_match <- match(
    split(new_locs, row(new_locs)),
    split(locs, row(locs))
  )

  locs_match_matrix <- Matrix::sparseMatrix(
    i = locs_match,
    j = seq_len(n_obs),
    x = 1
  )

  hctam_scol <- split(seq_len(n_obs), locs_match)
  hctam_scol_1 <- sapply(hctam_scol, function(x) x[1])

  nn_new <- FNN::get.knnx(
    data  = vecchia_approx$locs,
    query = new_locs_unique,
    k     = nrow(vecchia_approx$NNarray) - 1
  )$nn.index
  
  NNarray <- vecchia_approx$NNarray
  if(length(nn_new)>0){
    NNarray <- cbind(
      vecchia_approx$NNarray,
      t(cbind(seq(vecchia_approx$n_locs + 1, n_locs), nn_new))
    )
  }
  NNNoNA <- !is.na(NNarray)
  sparse_mat <- Matrix::sparseMatrix(
    x = seq_len(sum(NNNoNA)),
    i = col(NNarray)[NNNoNA],
    j = NNarray[NNNoNA],
    triangular = TRUE
  )

  sparse_chol_x_reorder <- seq_along(NNarray)[NNNoNA][
    match(sparse_mat@x, seq_len(sum(NNNoNA)))
  ]

  list(
    new_observed_locs     = new_locs,
    observed_locs         = observed_locs,
    previous_n_obs        = vecchia_approx$n_obs,
    previous_n_locs       = vecchia_approx$n_locs,
    new_locs_match        = new_locs_match,
    n_obs                 = n_obs,
    new_n_obs             = new_n_obs,
    locs                  = locs,
    new_locs              = new_locs_unique,
    n_locs                = n_locs,
    t_locs                = t(locs),
    locs_match            = locs_match,
    locs_match_matrix     = locs_match_matrix,
    hctam_scol            = hctam_scol,
    hctam_scol_1          = hctam_scol_1,
    NNarray               = NNarray,
    NNNoNA                = NNNoNA,
    sparse_chol_x_reorder = sparse_chol_x_reorder
  )
}

# extended_PP = extendPP(PP, extended_vecchia_approx)
extendPP <- function(PP, extended_vecchia_approx) {
  if (is.null(PP)) {
    return(NULL)
  }
  createPP(
    vecchia_approx = extended_vecchia_approx,
    matern_range = PP$matern_range, knots = PP$knots, reorder_knots = F,
    plot = F
  )
}

predict1Field <- function(range_beta,
                          field,
                          extended_PP_range,
                          extended_vecchia_approx,
                          complete_range_X_locs,
                          hierarchical_model,
                          num_threads) {
  compressed_chol <-
    array(0, c(nrow(extended_vecchia_approx$NNarray), extended_vecchia_approx$n_locs, 1))
  log_range <- computeLogRange(
    range_beta = range_beta,
    PP = extended_PP_range,
    vecchia_approx = extended_vecchia_approx,
    range_X = complete_range_X_locs
  )
  res <- list()
  vecchia_(
    start_idx = extended_vecchia_approx$n_locs - 1,
    log_range = t(log_range), locs = extended_vecchia_approx$t_locs,
    NNarray = extended_vecchia_approx$NNarray,
    smoothness = hierarchical_model$matern_smoothness,
    compute_derivative = FALSE, num_threads = num_threads, result = compressed_chol
  )
  pred_field_sample <-
    c(
      field,
      (
        -GpGp::Linv_mult(
          Linv = t(compressed_chol[, , 1]),
          z = c(field, rep(0, nrow(extended_vecchia_approx$new_locs))),
          NNarray = t(extended_vecchia_approx$NNarray)
        )[-seq(extended_vecchia_approx$previous_n_locs)]
        + rnorm(extended_vecchia_approx$n_locs - extended_vecchia_approx$previous_n_locs)
      ) / compressed_chol[1, -seq(extended_vecchia_approx$previous_n_locs), 1]
    )
  res <- list(
    "log_range" = log_range[extended_vecchia_approx$new_locs_match, ],
    "field" = pred_field_sample[extended_vecchia_approx$new_locs_match]
  )

  return(res)
}


extendCovariate <- function(
    existing_covariate, new_X, new_X_name, one_obs_per_locs, 
    extended_vecchia_approx, extended_PP){
  if (!is.null(new_X) & !is.data.frame(new_X)) {
    stop(paste(new_X_name, "must be a data.frame or NULL"))
  }
  if (!identical(
    colnames(new_X),
    colnames(existing_covariate$arg)
  )) {
    stop(paste(c(
      paste("The variables names of", new_X_name, "must match those of noise_X used to train the model, who are"),
      colnames(existing_covariate$arg)
    ), collapse = ", "))
  }
  return(
    processCovariates(
    X = rbind(existing_covariate$arg, new_X), one_obs_per_locs = one_obs_per_locs,
    PP = extended_PP, vecchia_approx = extended_vecchia_approx)
  )
}

# TODO a paralleliser  (noise semble déjà marcher ??)
#' Predict a 'GeoNonStat' object on new locations
#'
#' @param object an object of class \code{GeoNonStat}
#' @param new_data_list a list of new data. The list should have 
#' a 'locs' entry, a data.frame (or 2 column matrix) of new locations,
#' a 'X' entry, for new covariates
#' a 'noise_X' entry (optional) for new X noise parameters
#' a 'range_X' entry (optional) for new X range parameters
#' @param burn_in numeric, value, between 0 and 1. Gives the proportion of
#' first records that should be removed. Default to 0.1
#' @param num_threads number of threads. Integer.
#' @param ... additional arguments (unused)
#' @rdname GeoNonStat
#' @export
#' @method predict GeoNonStat
#' @examples
#' # Recreating new data with less observations but same structure as example data
#' locs <- cbind(rnorm(200), rnorm(200))
#' X <- data.frame(matrix(rnorm(200*4), ncol=4))
#' colnames(X) <- paste0("V", 1:4)
#' testdata <- list("locs" = locs, "X" = X, noise_X=X)  
#' predict(processedGnsDemo, new_data_list=testdata)
predict.GeoNonStat <-  function(
    object, 
    new_data_list,
    burn_in=0.2, 
    num_threads = 5,
    ...
    ) {
  if(mcmcCount(object)<=0) {
    objname <- deparse(substitute(object))
    stop(paste0("MCMC doesn't seem to have run on '", 
                objname, 
                "',\ntry again after a call to ",
                "'GeoNonStatMcmc' or 'automaticMcmc' function"))
  }
  
  # setup ########
  samples <- list()
  # getting elements
  if(!("locs" %in% names(new_data_list))) {
    stop("new_data_list should have a 'locs' element.")
  }
  new_locs <- new_data_list$locs
  new_X <- new_noise_X <- new_range_X <- NULL
  if("X" %in% names(new_data_list))  new_X <- new_data_list$X
  if("range_X" %in% names(new_data_list))  new_range_X <- new_data_list$range_X
  if("noise_X" %in% names(new_data_list))  new_noise_X <- new_data_list$noise_X
  
  # extending Vecchia approx and new PP
  extended_vecchia_approx <- extendVecchia(object$vecchia_approx, new_locs)
  extended_PP_noise <- extendPP(object$hierarchical_model$noise$PP, extended_vecchia_approx)
  extended_PP_range <- extendPP(object$hierarchical_model$range$PP, extended_vecchia_approx)
  # Extending covariates
  extended_noise_X <- extendCovariate(
    new_X = new_noise_X, extended_PP = extended_PP_noise, 
    new_X_name = "new_noise_X", existing_covariate = object$covariates$noise_X,
    one_obs_per_locs = F,
    extended_vecchia_approx = extended_vecchia_approx)
  extended_range_X <- extendCovariate(
    new_X = new_range_X, extended_PP = extended_PP_range, 
    new_X_name = "new_range_X", existing_covariate = object$covariates$range_X,
    one_obs_per_locs = T,
    extended_vecchia_approx = extended_vecchia_approx)
  # relevant MCMC samples
  records <- aggregateRecords(
    object$records, 
    keep=c("all"),
    keep_separate_chains = FALSE, 
    burn_in = burn_in)
  # predicting new noise ########
  samples$noise_log_var <- matrix(
    NA,
    ncol = extended_vecchia_approx$new_n_obs,
    nrow = nrow(records$noise_beta)
  )
  for (i in seq_len(nrow(records$noise_beta))) {
    samples$noise_log_var[i,] <- log(noiseVar(
      X = extended_noise_X$X, noise_beta = records$noise_beta[i,], 
      vecchia_approx = extended_vecchia_approx, PP = extended_PP_noise
    )[-seq_len(extended_vecchia_approx$previous_n_obs)])
  }
  gc()
  # predicting new field, and new range ########
  samples$field <- samples$noise_log_var
  samples$`range_beta range` <- samples$noise_log_var
  if(object$hierarchical_model$anisotropic){
    samples$`range_beta aniso1` <- samples$noise_log_var
    samples$`range_beta aniso2` <- samples$noise_log_var
  }
  prediction_vecchia <- array(0, c(nrow(extended_vecchia_approx$NNarray), extended_vecchia_approx$n_locs, 1))
  prediction_chol <- Matrix::sparseMatrix(
    i = col(extended_vecchia_approx$NNarray)[!is.na(extended_vecchia_approx$NNarray)],
    j = extended_vecchia_approx$NNarray[!is.na(extended_vecchia_approx$NNarray)], 
    triangular = T, x = 1.0
  )
  log_range <- matrix(0, ncol(extended_vecchia_approx$NNarray), 1 + 2*object$hierarchical_model$anisotropic)
  range_beta <- object$records$chain_1[[1]]$range_beta
  # TODO tenter de paralleliser avec future. 
  for (i in seq_len(nrow(records$field))) {
    # getting range high level parameters
    range_beta[,1] <- (records$`range_beta range`[i,])
    if(object$hierarchical_model$anisotropic){
      range_beta[,2] <- records$`range_beta aniso1`[i,]
      range_beta[,3] <- records$`range_beta aniso2`[i,]
    }
    # getting spatial range 
    log_range[] <- computeLogRange(
      range_beta = range_beta,
      PP = extended_PP_range,
      vecchia_approx = extended_vecchia_approx,
      range_X = extended_range_X
    )
    # range at predicted locations
    samples$`range_beta range`[i,] <-  log_range[extended_vecchia_approx$new_locs_match,1]
    if(object$hierarchical_model$anisotropic){
      samples$`range_beta aniso1`[i,] <- log_range[extended_vecchia_approx$new_locs_match,2]
      samples$`range_beta aniso2`[i,] <- log_range[extended_vecchia_approx$new_locs_match,3]
    }
    # computing Vecchia at both observed and new observations
    vecchia_(
      start_idx = 1, 
      log_range = t(log_range), locs = extended_vecchia_approx$t_locs, 
      NNarray = extended_vecchia_approx$NNarray, 
      smoothness = object$hierarchical_model$matern_smoothness, 
      compute_derivative = F,
      num_threads = num_threads, 
      result = prediction_vecchia
    )
    prediction_chol@x <- prediction_vecchia[,,1][extended_vecchia_approx$sparse_chol_x_reorder]
    # predicting field
    field_pred <- c(records$field[i,], rep(0, extended_vecchia_approx$n_locs - extended_vecchia_approx$previous_n_locs))
    field_pred <- prediction_chol %*% field_pred
    field_pred[-seq(extended_vecchia_approx$previous_n_locs)] <- rnorm(extended_vecchia_approx$n_locs - extended_vecchia_approx$previous_n_locs)
    field_pred <- Matrix::solve(prediction_chol, field_pred)
    samples$field[i,] <- field_pred[extended_vecchia_approx$new_locs_match]
  }
  gc()
  # getting denoised signal at locations
  samples$fixed_effects <- as.matrix(records$beta %*% t(cbind(rep(1, nrow(new_locs)), new_X)))
  samples$denoised <- samples$field + samples$fixed_effects
  # getting summaries
  summaries <- list()
  for(name in names(samples)) summaries[[name]] <- summarizeRecords(samples[[name]])
  res <- list(samples = samples, summaries = summaries)
  return(res)
}
