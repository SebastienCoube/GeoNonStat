#' Filter records of a GeoNonStat object
#'
#' @param records a list containing records from a `GeoNonStat` object.
#' @param burn_in numeric, value, between 0 and 1. Gives the percentage of
#' first records that should be removed If burn_in = 0.1 : the first 10% of
#' records will be removed from the result.
#' @param keep character vector. What params should be kept. Can also be
#' `"all"`, to keep all parameters, or `"all_even_the_field"` to keep all
#' parameters, event the `field.`
#'
#' @returns a filtered list or records
#' @export
#' @keywords internal
#'
#' @examples
#' redCords <- filterRecords(processedGnsDemo$records, burn_in = 0.2)
filterRecords <- function(records, burn_in = 0, keep = "all") {
  # Reduce paramaters
  namesparam <- names(records[[1]][[1]])
  if (length(keep) == 1 && keep == "all_even_the_field") keep <- namesparam
  if (length(keep) == 1 && keep == "all") keep <- setdiff(namesparam, "field")
  if (any(sapply(keep, function(x) !(x %in% namesparam)))) {
    stop(
      "Try to keep a parameter that doesn't exists. Acceptable values for",
      "`keep` are 'all', 'all_even_the_field' or a character vector with ",
      "parameters to keep"
    )
  }
  records <- lapply(records, function(x) {
    lapply(x, function(y) y[keep])
  })

  # Reduce iterations
  n_iter <- length(records[[1]])
  start_iter <- max(0, floor(n_iter * burn_in))
  if (n_iter - start_iter <= 0) {
    stop("Not enough iterations, reduce burn_in if necessary")
  }
  if (!all(lengths(records) == n_iter)) {
    stop(
      "Something is wrong with your object, ",
      "all the MCMC chains should have the same number of iterations"
    )
  }
  keep_iters <- (start_iter + 1):n_iter

  return(
    lapply(records, function(x) {
      return(x[keep_iters])
    })
  )
}

#' Transpose a list
#'
#' @param list a list containing sublists (all having the same structure).
#'
#' @returns a list, transposed
#' @export
#' @keywords internal
#'
#' @examples
#' A <- list("R1" = "a1", "R2" = "a2", "R3" = "a3")
#' B <- list("R1" = "b1", "R2" = "b2", "R3" = "b3")
#' C <- list("R1" = "c1", "R2" = "c2", "R3" = "c3")
#' transposeList(list("A" = A, "B" = B, "C" = C))
transposeList <- function(list) {
  nameslist <- names(list[[1]])
  tl <- lapply(
    seq_along(list[[1]]),
    function(i) lapply(list, `[[`, i)
  )
  names(tl) <- nameslist
  return(tl)
}

#' Aggregate records of a GeoNonStat object, for 1 chain
#'
#' @param chain a list of records containing for each iteration,
#' a list of parameters
#' @returns a list containing matrices of parameters,
#' aggregated over the iterations.
#' @export
#' @keywords internal
#'
#' @examples
#' tmp <- transposeList(processedGnsDemo$records[[1]])
#' agg <- aggregateRecordsChain(tmp)
aggregateRecordsChain <- function(chain) {
  res <- lapply(chain, function(param) {
    # Numeric case
    if (!is.matrix(param[[1]])) {
      param <- lapply(param, function(x) {
        dim(x) <- c(length(x), 1)
        return(x)
      })
    }

    if (ncol(param[[1]]) == 1) {
      bindres <- do.call(rbind, lapply(param, as.numeric))
      colnames(bindres) <- rownames(param[[1]])
    } else {
      bindres <- do.call(
        rbind,
        lapply(param, function(m) {
          v <- as.vector(m)
          names(v) <- paste(
            rep(rownames(m), times = ncol(m)),
            rep(colnames(m), each = nrow(m)),
            sep = "_"
          )
          v
        })
      )
    }
    return(bindres)
  })
  return(res)
}


#' Aggregate records of a GeoNonStat object, for 1 chain
#'
#' @param object a GeoNonStat object containing records
#' @param burn_in numeric, value, between 0 and 1. Gives the percentage of
#' first records that should be removed If burn_in = 0.1 : the first 10% of
#' records will be removed from the result.
#' @param keep character vector. What params should be kept. Can also be
#' `"all"`, to keep all parameters, or `"all_even_the_field"` to keep all
#' parameters, event the `field.`
#'
#' @returns a list containing matrices of parameters,
#' aggregated over the iterations.
#' @export
#' @keywords internal
#'
#' @examples
#' GeoNonStat:::aggregateRecords(processedGnsDemo)
aggregateRecords <- function(object, burn_in = 0.1, keep = "all", keep_separate_chains = FALSE) {
  res <- filterRecords(object$records, burn_in = burn_in, keep = keep)
  namesparams <- names(res[[1]][[1]])
  res <- lapply(res, transposeList)
  res <- lapply(res, aggregateRecordsChain)
  if (keep_separate_chains) {
    res <- lapply(res, function(x) {
      names(x) <- namesparams
      return(x)
    })
  } else {
    res <- do.call(Map, c(list(rbind), res))
    names(res) <- namesparams
  }
  return(res)
}

#' Summarize a matrix by column
#'
#' @param mat a numeric matrix
#' @param quant named numerical vector. What quantiles to give
#'
#' @returns a matrix
#' @export
#'
#' @examples
#' tt <- matrix(rnorm(200), ncol = 20)
#' summarizeRecords(tt)
summarizeRecords <- function(mat, quant = c("q 2.5%" = .025, "median" = .5, "q 97.5%" = .975)) {
  matstat <- apply(mat, 2, FUN = function(v) {
    vstat <- c(
      mean(v),
      sd(v),
      quantile(v, quant, names = FALSE)
    )
    return(signif(vstat, 3))
  })
  rownames(matstat) <- c("mean", "sd", names(quant))
  return(matstat)
}


#' Estimate the records of a GeoNonStat object.
#'
#' @param object a GeoNonStat object containing records
#' @param burn_in numeric, value, between 0 and 1. Gives the percentage of
#' first records that should be removed If burn_in = 0.1 : the first 10% of
#' records will be removed from the result.
#' @param keep character vector. What params should be kept. Can also be
#' `"all"`, to keep all parameters, or `"all_even_the_field"` to keep all
#' parameters, event the `field.`
#'
#' @returns a list of numerical summary by parameters.
#' @export
#'
#' @examples
#' estimate(processedGnsDemo)
estimate <- function(object, burn_in = .1, keep = "all") {
  res <- aggregateRecords(object, burn_in = burn_in, keep = "all")
  res <- sapply(res, summarizeRecords, simplify = FALSE, USE.NAMES = TRUE)
  return(res)
}


#' Title
#'
#' @param X a matrix of samples
#' @param beta_samples a matrix of beta
#' @param field_samples a matrix of field
#' @param noise_var_samples a matrix of variance noise
#' @param observed_field a matrix, the observed field
#'
#' @returns a list
elppd <- function(X, beta_samples, field_samples, noise_var_samples, observed_field) {
  elppd_per_obs <- (
    -log((sqrt(2) * pi)) # 1/sqrt(2 pi) passed to the log
    - .5 * log(noise_var_samples) # variance penalty of Gaussian density
      - .5 * (
        (X %*% t(beta_samples) + t(field_samples)) # mean
        - observed_field # observation
      )^2 / noise_var_samples # variance
  )
  elppd_per_obs <- apply(elppd_per_obs, 1, mean)
  return(list("elppd_per_obs" = elppd_per_obs, "mean_elppd" = mean(elppd_per_obs)))
}

#' Compute the log-dens of the observations
#'
#' @param object a GeoNonStat object
#' @param burn_in numeric, value, between 0 and 1. Gives the percentage of
#' first records that should be removed If burn_in = 0.1 : the first 10% of
#' records will be removed from the result.
#'
#' @returns a list
#' @export
#'
#' @examples
#' lppd <- elppdTrain(processedGnsDemo)
elppdTrain <- function(object, burn_in = .1) {
  aggregated_records <- aggregateRecords(object, burn_in, keep = c("beta", "noise_beta", "field"))
  noise_var <-
    apply(aggregated_records$noise_beta, 1, function(x) {
      as.vector(exp(xPPMultRight(
        X = object$covariates$noise_X$X, PP = object$hierarchical_model$noise$PP,
        vecchia_approx = object$vecchia_approx, Y = x,
        permutate_PP_to_obs = T
      )))
    })
  # Gaussian log-dens of the observations
  return(elppd(
    object$covariates$X$X,
    beta_samples = aggregated_records$beta,
    field_samples = aggregated_records$field[, object$vecchia_approx$locs_match],
    noise_var_samples = noise_var,
    observed_field = object$observed_field
  ))
}
