rbind_by_name <- function(x) {
  nms <- names(x[[1]])
  res <- lapply(nms, function(nm) {
    elems <- lapply(x, `[[`, nm)
    is_null <- vapply(elems, is.null, logical(1))
    if (all(is_null)) return(NULL)
    if (any(is_null)) {
      stop("Element '", nm, "' is missing from some of the data lists")
    }
    if (all(vapply(elems, function(e) is.matrix(e) || is.data.frame(e), logical(1)))) {
      do.call(rbind, unname(elems))
    } else {
      do.call(c, unname(elems))
    }
  })
  names(res) <- nms
  res
}


PPPredict <- function(PP, vecchia_approx_chunk) {
  if (is.null(PP)) return(NULL)
  return(createPP(
    vecchia_approx = vecchia_approx_chunk, matern_range = PP$matern_range,
    knots = PP$knots, reorder_knots = F, plot = F, verbose = F))
}


#' Predict a 'GeoNonStat' object on new locations
#'
#' @param object an object of class \code{GeoNonStat}
#' @param new_data_lists a list of new data. The list should have 
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
#' testpredictions <- predict(processedGnsDemo, list(testdata))
predict.GeoNonStat <-  function(
    object, 
    new_data_lists,
    burn_in = 0.2,
    num_threads = 5,
    ...
) {
  # setup ########
  if (mcmcCount(object) <= 0) {
    objname <- deparse(substitute(object))
    stop(paste0("MCMC doesn't seem to have run on '",
                objname,
                "',\ntry again after a call to ",
                "'GeoNonStatMcmc' or 'automaticMcmc' function"))
  }
  aniso <- isTRUE(object$hierarchical_model$anisotropic)
  
  # putting NAs if no observed field
  for (i in seq_along(new_data_lists)) {
    if (is.null(new_data_lists[[i]]$observed_field)) {
      new_data_lists[[i]]$observed_field <- rep(NA_real_, nrow(new_data_lists[[i]]$locs))
    }
  }
  
  # merging data lists ####
  data_list_idx <- rep(seq_along(new_data_lists), sapply(new_data_lists, function(x) nrow(x$locs)))
  data_list_pred <- rbind_by_name(new_data_lists)
  n_all <- length(data_list_idx)
  
  # instantiating result
  empty_summary <- matrix(NA_real_, n_all, 5)
  colnames(empty_summary) <- c("mean", "sd", "q 2.5%", "median", "q 97.5%")
  predicted_stuff_names <-
    c("field", "noise_log_var", "fixed", "field_and_fixed", "range")
  if (aniso) predicted_stuff_names <- c(predicted_stuff_names, "aniso1", "aniso2")
  summaries <- lapply(predicted_stuff_names, function(x) empty_summary)
  names(summaries) <- predicted_stuff_names
  scores <- matrix(NA_real_, n_all, 4)
  colnames(scores) <- c("95% coverage", "CRPS", "predictive density", "squared error")
  
  
  # relevant MCMC samples
  records <- aggregateRecords(
    object$records,
    keep = c("all"),
    keep_separate_chains = FALSE,
    burn_in = burn_in)
  gc()
  
  # predicting
  # exact matching of duplicated locations (much faster than matching lists of rows)
  pred_locs <- as.matrix(data_list_pred$locs)
  locs_key <- do.call(paste, c(
    lapply(seq_len(ncol(pred_locs)), function(j) sprintf("%.17g", pred_locs[, j])),
    sep = "|"))
  first_occurrence <- !duplicated(locs_key)
  unique_pred_locs <- pred_locs[first_occurrence, , drop = FALSE]
  pred_locs_match <- match(locs_key, locs_key[first_occurrence])
  
  # chunk of non duplicated locations
  chunk_size <- 1000
  n_unique <- nrow(unique_pred_locs)
  chunks <- split(seq_len(n_unique), ceiling(seq_len(n_unique) / chunk_size))
  
  n_mcmc_iters <- nrow(records$field)
  
  # range coefficients at one MCMC iteration (isotropic models have no aniso coefficients)
  getRangeBeta <- function(mcmc_iter) {
    if (aniso) {
      cbind(
        records$`range_beta range`[mcmc_iter, ],
        records$`range_beta aniso1`[mcmc_iter, ],
        records$`range_beta aniso2`[mcmc_iter, ])
    } else {
      cbind(records$`range_beta range`[mcmc_iter, ])
    }
  }
  
  # log-range at training locations, computed once for all chunks
  # (array of size n_train x n_range_params x n_mcmc_iters)
  n_range_params <- 1 + 2 * aniso
  n_train <- object$vecchia_approx$n_locs
  log_range_train_all <- array(0, dim = c(n_train, n_range_params, n_mcmc_iters))
  for (mcmc_iter in seq_len(n_mcmc_iters)) {
    log_range_train_all[, , mcmc_iter] <- computeLogRange(
      range_beta = getRangeBeta(mcmc_iter),
      range_X = object$covariates$range_X,
      PP = object$hierarchical_model$range$PP,
      vecchia_approx = object$vecchia_approx
    )[]
  }

  # subsetting potentially duplicated locations using non-duplicated locations
  for (chunk in chunks) {
    # setup
    # matching potentially duplicated locations to unique location
    which_obs_in_chunk_idx <- which(pred_locs_match %in% chunk)
    # extending Vecchia approx and PP at chunk of observed locations
    vecchia_approx_chunk <- createVecchia(pred_locs[which_obs_in_chunk_idx, , drop = FALSE])
    range_PP_chunk <- PPPredict(vecchia_approx_chunk = vecchia_approx_chunk, PP = object$hierarchical_model$range$PP)
    noise_PP_chunk <- PPPredict(vecchia_approx_chunk = vecchia_approx_chunk, PP = object$hierarchical_model$noise$PP)
    # Processing covariates at chunk of observed locations
    range_X_chunk <- processCovariates(
      X = data_list_pred$range_X[which_obs_in_chunk_idx, , drop = FALSE], vecchia_approx = vecchia_approx_chunk,
      PP = range_PP_chunk, one_obs_per_locs = T)
    noise_X_chunk <- processCovariates(
      X = data_list_pred$noise_X[which_obs_in_chunk_idx, , drop = FALSE], vecchia_approx = vecchia_approx_chunk,
      PP = noise_PP_chunk, one_obs_per_locs = F)
    X_chunk <- processCovariates(
      X = data_list_pred$X[which_obs_in_chunk_idx, , drop = FALSE], vecchia_approx = vecchia_approx_chunk,
      one_obs_per_locs = F)
    # allocating MCMC samples at chunk of observed locations
    chunk_samples <- lapply(predicted_stuff_names, function(x) matrix(0, n_mcmc_iters, length(which_obs_in_chunk_idx)))
    names(chunk_samples) <- predicted_stuff_names

    # predicting fixed effects
    chunk_samples$fixed[] <- (records$beta) %*% t(X_chunk$X)

    # Initializing log-range arrays for prediction
    log_range_chunk <- matrix(0, length(chunk), n_range_params)
    # Initializing Vecchia approximation for prediction
    m <- nrow(object$vecchia_approx$NNarray) - 1
    NNarray <- t(FNN::knnx.index(
      data  = object$vecchia_approx$locs,
      query = vecchia_approx_chunk$locs, k = m))
    vecchia_chunk <- matrix(0, length(chunk), m + 1)
    vecchia_chunk_expanded <- Matrix::sparseMatrix(
      i = c(t(col(NNarray))), j = c(t((NNarray))),
      x = c(vecchia_chunk[, -1]), repr = "T",
      dims = c(length(chunk), object$vecchia_approx$n_locs)
    )

    for (mcmc_iter in seq_len(n_mcmc_iters)) {
      # predicting noise log variance
      chunk_samples$noise_log_var[mcmc_iter, ] <- log(noiseVar(
        X = noise_X_chunk$X, PP = noise_PP_chunk,
        noise_beta = records$noise_beta[mcmc_iter, ], vecchia_approx = vecchia_approx_chunk))

      # log-range at training locations, precomputed before the chunk loop
      log_range_train <- matrix(log_range_train_all[, , mcmc_iter], n_train, n_range_params)
      log_range_chunk[] <- computeLogRange(
        range_beta = getRangeBeta(mcmc_iter),
        range_X = range_X_chunk, PP = range_PP_chunk,
        vecchia_approx = vecchia_approx_chunk
      )[]

      # Vecchia approximation
      vecchia_predict(
        log_range = t(log_range_train),
        predicted_log_range = t(log_range_chunk),
        predicted_locs = t(vecchia_approx_chunk$locs),
        locs = object$vecchia_approx$t_locs,
        NNarray = NNarray, num_threads = num_threads,
        result = vecchia_chunk,
        smoothness = object$hierarchical_model$matern_smoothness)
      vecchia_chunk_expanded@x <- c(vecchia_chunk[, -1])
      chunk_samples$field[mcmc_iter, ] <- (
        as.vector(
          - (vecchia_chunk_expanded) %*% records$field[mcmc_iter, ] +
            rnorm(vecchia_approx_chunk$n_locs)
        ) / vecchia_chunk[, 1]
      )[vecchia_approx_chunk$locs_match]
      # storing samples of log-range
      if (!aniso) {
        chunk_samples$range[mcmc_iter, ] <-  log_range_chunk[vecchia_approx_chunk$locs_match, 1]/2
      }else{
        # reparametrizing to range / aniso format
        log_range_chunk <- log_range_chunk %*% solve(matrix(c(1 / sqrt(2), 1 / sqrt(2), 0, 1 / sqrt(2), -1 / sqrt(2), 0, 0, 0, 1), 3) * sqrt(2) * 2)
        chunk_samples$range[mcmc_iter, ] <-  log_range_chunk[vecchia_approx_chunk$locs_match, 1]
        chunk_samples$aniso1[mcmc_iter, ] <- log_range_chunk[vecchia_approx_chunk$locs_match, 2]
        chunk_samples$aniso2[mcmc_iter, ] <- log_range_chunk[vecchia_approx_chunk$locs_match, 3]
      }
    }
    chunk_samples$field_and_fixed[] <-
      chunk_samples$fixed[] + chunk_samples$field[]
    for (name in predicted_stuff_names) {
      summaries[[name]][which_obs_in_chunk_idx, ] <- summarizeRecords(chunk_samples[[name]])
    }
    scores[which_obs_in_chunk_idx, ] <- allScores(
      true_y = data_list_pred$observed_field[which_obs_in_chunk_idx],
      sd_samples   = t(exp(.5 * chunk_samples$noise_log_var)),
      mean_samples = t(chunk_samples$field_and_fixed)
    )$scores_per_obs
  }

  res <- list()
  for (i in seq_along(new_data_lists)) {
    idx_i <- which(data_list_idx == i)
    res[[i]] <- list()
    res[[i]]$summaries <- lapply(summaries, function(summa) summa[idx_i, , drop = FALSE])
    res[[i]]$scores <- list()
    res[[i]]$scores$per_obs <- scores[idx_i, , drop = FALSE]
    res[[i]]$scores$aggregated <- signif(apply(res[[i]]$scores$per_obs, 2, mean, na.rm = TRUE), digits = 2)
    res[[i]]$locs <- new_data_lists[[i]]$locs
  }
  names(res) <- names(new_data_lists)
  return(res)
}