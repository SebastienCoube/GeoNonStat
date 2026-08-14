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
    PP = extended_PP_range, vecchia_approx = extended_vecchia_approx)
  )
}

# TODO a paralleliser  (noise semble déjà marcher ??)
#' Predict a 'GeoNonStat' object on new locations
#'
#' @param object an object of class \code{GeoNonStat}
#' @param new_locs new locations
#' @param new_X new X
#' @param new_noise_X new X for noise parameters
#' @param new_range_X new X for range parameters
#' @param burn_in numeric, value, between 0 and 1. Gives the proportion of
#' first records that should be removed. Default to 0.1
#' @param num_threads number of threads. Integer.
#' @param ... additional arguments (unused)
#' @rdname GeoNonStat
#' @export
#' @method predict GeoNonStat
predict.GeoNonStat <-  function(
    object, 
    new_locs = NULL, 
    new_X = NULL, 
    new_noise_X = NULL, 
    new_range_X = NULL, 
    burn_in=0.2, 
    num_threads = 5,
    return_samples = T,
    ...
    ) {
  # extending Vecchia approx and new PP
  extended_vecchia_approx <- extendVecchia(object$vecchia_approx, new_locs)
  extended_PP_noise <- extendPP(object$hierarchical_model$noise$PP, extended_vecchia_approx)
  extended_PP_range <- extendPP(object$hierarchical_model$range$PP, extended_vecchia_approx)
  # Extending covariates
  extended_X <- extendCovariate(
    new_X = new_X, extended_PP = NULL, 
    new_X_name = "new_X", existing_covariate = object$covariates$X,
    one_obs_per_locs = F,
    extended_vecchia_approx = extended_vecchia_approx)
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
  # predicting new noise
  pred_noise_var <- matrix(
    NA,
    nrow = extended_vecchia_approx$new_n_obs,
    ncol = nrow(records$noise_beta)
  )
  for (i in seq_len(nrow(records$noise_beta))) {
    pred_noise_var[, i] <- noiseVar(
      X = extended_noise_X$X, noise_beta = records$noise_beta[i,], 
      vecchia_approx = extended_vecchia_approx, PP = extended_PP_noise
    )[-seq_len(extended_vecchia_approx$previous_n_obs)]
  }
  # predicting new field, and new range
  predicted_field <- matrix(
    NA,
    nrow = extended_vecchia_approx$new_n_obs,
    ncol = nrow(records$noise_beta)
  )
  predicted_range <- predicted_field
  if(object$hierarchical_model$anisotropic){
    predicted_aniso1 <- predicted_field
    predicted_aniso2 <- predicted_field
  }
  
  for (i in seq_len(nrow(records$noise_beta))) {
    pred_noise_var[, i] <- noiseVar(
      X = extended_noise_X$X, noise_beta = records$noise_beta[i,], 
      vecchia_approx = extended_vecchia_approx, PP = extended_PP_noise
    )[-seq_len(extended_vecchia_approx$previous_n_obs)]
  }
  
   pred_log_range <- lapply(predicted_latent, function(x)x$log_range)
   
   res = list(
     noise_log_var = t(summarizeRecords(t(predicted_log_variance_noise))), 
     fixed_effects = t(summarizeRecords(t(new_X %*% t(estimates$beta)))), 
     field         = t(summarizeRecords(t(sapply(predicted_latent, function(x)x$field))))
   )
   

   res$range         = summarizeRecords()
   if(object$hierarchical_model$anisotropic){
     res$pred_log_range_aniso1 <- t(summarizeRecords(t(sapply(pred_log_range, function(x)x[,2]))))
     res$pred_log_range_aniso2 <- t(summarizeRecords(t(sapply(pred_log_range, function(x)x[,3]))))
   }
   
  rm(predicted_latent)
  rm(latent_samples)
  rm(predicted_log_variance_noise)
  rm(complete_range_X_locs)
  rm(pred_log_range)
  invisible(gc())
   return(res)
}
