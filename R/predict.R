extendVecchia <- function(vecchia_approx, new_locs) {
  observed_locs <- rbind(vecchia_approx$observed_locs, new_locs)
  locs <- unique(rbind(vecchia_approx$locs, new_locs))
  n_obs <- nrow(observed_locs)
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

  NNarray <- cbind(
    vecchia_approx$NNarray,
    t(cbind(seq(vecchia_approx$n_locs + 1, n_locs), nn_new))
  )
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


predict1LogVarNoise <- function(noise_beta, extended_PP_noise, extended_vecchia_approx, new_noise_X) {
  PP_effect <- xPPMultRight(
    X = NULL, PP = extended_PP_noise,
    vecchia_approx = extended_vecchia_approx,
    permutate_PP_to_obs = T,
    Y = noise_beta[-seq_len(ncol(new_noise_X))]
  )

  PP_effect <- PP_effect[-seq_len(extended_vecchia_approx$previous_n_obs)]
  X_effect <- new_noise_X %*% noise_beta[seq_len(ncol(new_noise_X))]
  res <- PP_effect + X_effect
  return(res)
}

# predictNoise <- function(geo_non_stat, extended_vecchia_approx, new_noise_X, burn_in){
#   extended_PP_noise <- extendPP(geo_non_stat$hierarchical_model$noise$PP)
#   new_noise_X <- cbind(rep(1, extended_vecchia_approx$n_obs -
#                              extended_vecchia_approx$previous_n_obs),
#                        new_noise_X)
#   # apply to MCMC samples of noise_beta from geo_non_stat
# }

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
predict.GeoNonStat <- function(
  object,
  new_locs,
  new_X = NULL,
  new_noise_X = NULL,
  new_range_X = NULL,
  burn_in = 0.2,
  num_threads = 5,
  ...
) {
  # extending Vecchia approx and new PP
  extended_vecchia_approx <- extendVecchia(object$vecchia_approx, new_locs)
  extended_PP_noise <- extendPP(object$hierarchical_model$noise$PP, extended_vecchia_approx)
  extended_PP_range <- extendPP(object$hierarchical_model$range$PP, extended_vecchia_approx)

  # setup for new_noise_X
  if (!is.null(new_noise_X) & !is.data.frame(new_noise_X)) {
    stop("new_noise_X must be a data.frame or NULL")
  }
  if (!identical(
    colnames(new_noise_X),
    colnames(object$covariates$noise_X$arg)
  )) {
    stop(paste(c(
      "The variables names of new_noise_X must match those of noise_X used to train the model, who are",
      colnames(object$covariates$noise_X$arg)
    ), collapse = ", "))
  }
  new_noise_X <- model.matrix(rep(1, nrow(new_locs)) ~ ., new_noise_X)

  # setup for new_X
  if (!is.null(new_X) & !is.data.frame(new_X)) {
    stop("new_X must be a data.frame or NULL")
  }
  if (!identical(
    colnames(new_X),
    colnames(object$covariates$X$arg)
  )) {
    stop(paste(c(
      "The variables names of new_X must match those of X used to train the model, who are",
      colnames(object$covariates$X$arg)
    ), collapse = ", "))
  }
  new_X <- model.matrix(rep(1, nrow(new_locs)) ~ ., new_X)

  # setup for new_range_X
  if (!identical(
    colnames(new_range_X),
    colnames(object$covariates$range_X$arg)
  )) {
    stop(paste(c(
      "The variables names of new_range_X must match those of range_X, who are",
      colnames(object$covariates$range_X$arg)
    ), collapse = ", "))
  }
  complete_range_X_locs <- processCovariates(
    X = rbind(object$covariates$range_X$arg, new_range_X), one_obs_per_locs = T,
    PP = extended_PP_range, vecchia_approx = extended_vecchia_approx
  )

  # getting MCMC samples
  estimates <- aggregateRecords(object,
    keep = c("beta", "noise_beta"),
    keep_separate_chains = FALSE,
    burn_in = burn_in
  )
  predicted_log_variance_noise <- apply(estimates$noise_beta, 1, function(x) {
    predict1LogVarNoise(
      x,
      extended_PP_noise,
      extended_vecchia_approx,
      new_noise_X
    )
  })
  gc()

  # field + range
  latent_samples <- filterRecords(object$records,
    keep = c("range_beta", "field"),
    burn_in = burn_in
  )
  latent_samples <- unlist(latent_samples, recursive = FALSE)
  predicted_latent <- lapply(
    latent_samples, function(x) {
      predict1Field(
        range_beta = x$range_beta, field = x$field,
        extended_PP_range = extended_PP_range,
        extended_vecchia_approx = extended_vecchia_approx,
        complete_range_X_locs = complete_range_X_locs,
        object$hierarchical_model,
        num_threads = num_threads
      )
    }
  )
  gc()
  pred_log_range <- lapply(predicted_latent, function(x) x$log_range)

  samples <- list(
    noise_log_var = predicted_log_variance_noise,
    fixed_effects = new_X %*% t(estimates$beta),
    field = sapply(predicted_latent, function(x) x$field),
    range = sapply(pred_log_range, function(x) x[, 1])
  )
  if (object$hierarchical_model$anisotropic) {
    samples$pred_log_range_aniso1 <- sapply(pred_log_range, function(x) x[, 2])
    samples$pred_log_range_aniso2 <- sapply(pred_log_range, function(x) x[, 3])
  }
  gc()
  return(list(
    samples = samples,
    summaries = lapply(samples, function(x) t(summarizeRecords(t(x))))
  ))
}
