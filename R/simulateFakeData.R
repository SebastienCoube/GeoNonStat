#' Ingredients needed to simulate observations from a GeoNonStat model: 
#' explanatory variables, spatial process, and noise process. 
#' 
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param X a data.frame of covariates explaining the interest variable
#' through fixed linear effects
#' @param matern_smoothness numerical value, a Matérn smoothness parameter,
#' either 0.5 (aka ``exponential kernel'') or 1.5
#' @param noise_X a data.frame of covariates explaining the Gaussian noise
#' variance through fixed linear effects
#' @param noise_PP either an object of class PP used to model the field of
#' log-noise parameters, or NULL in which case the model is stationary
#' @param range_X a data.frame of covariates explaining the Gaussian process
#' range through fixed linear effects
#' @param range_PP TODO
#' @param anisotropic logical, default to FALSE. Is the covariance anisotropic ?
#'
#' @returns a list
#' @export
#'
#' @examples
#' set.seed(2)
#' nlocs <- 15000
#' observed_locs <- cbind(runif(nlocs), runif(nlocs))
#' vecchia_approx <- createVecchia(observed_locs, m = 6)
#'
#' # fixed effects
#' X <- data.frame(
#'   V1 = observed_locs[, 1],
#'   V2 = rnorm(nrow(observed_locs)),
#'   V3 = rnorm(nrow(observed_locs)),
#'   V4 = rnorm(nrow(observed_locs))
#' )
#'
#' PP_range <- createPP(vecchia_approx, knots = 16, matern_range = .3)
#' PP_noise <- createPP(vecchia_approx, knots = 500, matern_range = .06)
#'
#' gns_simulator <- createGnsSimulator(
#'   vecchia_approx = vecchia_approx,
#'   X = X,
#'   noise_X = X, noise_PP = PP_noise,
#'   range_X = NULL, range_PP = PP_range,
#'   anisotropic = TRUE
#' )
createGnsSimulator <- function(vecchia_approx,
                               X = NULL, # Response variable
                               matern_smoothness = 1.5,
                               # Matern smoothness
                               noise_X = NULL,
                               noise_PP = NULL,
                               range_X = NULL,
                               range_PP = NULL,
                               anisotropic = FALSE) {
  covariates <- list()
  covariates$X <- processCovariates(
    X = X,
    vecchia_approx = vecchia_approx,
    PP = NULL,
    one_obs_per_locs = FALSE
  )
  # fixed effects and PP for range
  covariates$range_X <- processCovariates(
    X = range_X,
    vecchia_approx = vecchia_approx,
    PP = range_PP,
    one_obs_per_locs = TRUE
  )
  # fixed effects and PP for noise
  covariates$noise_X <- processCovariates(
    X = noise_X,
    vecchia_approx = vecchia_approx,
    PP = noise_PP,
    one_obs_per_locs = FALSE
  )
  return(
    list(
      "covariates" = covariates,
      "vecchia_approx" = vecchia_approx,
      "noise_PP" = noise_PP,
      "range_PP" = range_PP,
      "matern_smoothness" = matern_smoothness,
      "anisotropic" = anisotropic
    )
  )
}


#' @title From a GeoNonStatSimulator object, creates parameters at the right format
#'
#' @param gns_simulator an object created by `createGnsSimulator()`
#' @param range_intercept TODO
#' @param range_PP_log_var TODO
#' @param noise_intercept TODO
#' @param noise_PP_log_var TODO
#' @param field_log_var TODO
#' @param verbose logical, default to TRUE. Active verbose mode ?
#' @param seed integer value, seed used for reproducibility purposes. Default to 1
#'
#' @returns a list
#' @export
#'
#' @examples
#' set.seed(2)
#' nlocs <- 15000
#' observed_locs <- cbind(runif(nlocs), runif(nlocs))
#' vecchia_approx <- createVecchia(observed_locs, m = 6)
#'
#' # fixed effects
#' X <- data.frame(
#'   V1 = observed_locs[, 1],
#'   V2 = rnorm(nrow(observed_locs)),
#'   V3 = rnorm(nrow(observed_locs)),
#'   V4 = rnorm(nrow(observed_locs))
#' )
#'
#' # X_range = as.data.frame(observed_locs)
#'
#' PP_range <- createPP(vecchia_approx, knots = 16, matern_range = .3)
#' PP_noise <- createPP(vecchia_approx, knots = 500, matern_range = .06)
#'
#'
#' # full monty
#' gns_simulator <- createGnsSimulator(
#'   vecchia_approx = vecchia_approx,
#'   X = X,
#'   noise_X = X, noise_PP = PP_noise,
#'   range_X = NULL, range_PP = PP_range,
#'   anisotropic = TRUE
#' )
#' gns_simulator_parameters <- createGnsSimulatorParameters(gns_simulator)
createGnsSimulatorParameters <- function(gns_simulator,
                                         range_intercept = NULL,
                                         range_PP_log_var = NULL,
                                         noise_intercept = NULL,
                                         noise_PP_log_var = NULL,
                                         field_log_var = NULL,
                                         verbose = TRUE,
                                         seed = 1) {
  set.seed(seed)
  covs <- gns_simulator$covariates
  res <- list(
    "X_coeff" = matrix(
      0,
      nrow = ncol(covs$X$X),
      dimnames = list(colnames(covs$X$X))
    ),
    "noise_X_coeff" = matrix(
      0,
      nrow = ncol(covs$noise_X$X),
      dimnames = list(colnames(covs$noise_X$X))
    ),
    "range_X_coeff" = matrix(
      0,
      nrow = ncol(covs$range_X$X_locs),
      dimnames = list(colnames(covs$range_X$X_locs), c("det", "aniso", "aniso")[seq(1 + 2 * gns_simulator$anisotropic)]),
      ncol = 1 + 2 * gns_simulator$anisotropic
    ),
    "noise_PP_coeff" = NULL,
    "field_log_var" = NULL
  )

  if (verbose) {
    message("Creating parameters to simulate from GeoNonStat model")
  }

  # fixed effects for the response
  res$X_coeff[] <- rnorm(length(res$X_coeff))
  if (verbose) {
    message("------------------")
    message("LINEAR REGRESSION PARAMETERS")
    message("The coefficients in X_coeff have been set to random normal coeffs")
  }
  # field log variance
  if (!is.null(field_log_var)) {
    res$field_log_var <- field_log_var
  } else {
    res$field_log_var <- 0
    if (verbose) {
      message("------------------")
      message("FIELD LOG VARIANCE PARAMETER")
      message("The field log-variance has been set to 0")
    }
  }
  if (verbose) {
    message("------------------")
    message("RANGE PARAMETERS")
    message(
      paste(
        "Note : the model is",
        c(
          "isotropic so range_X_coeff and range_PP_coeff (if applicable) have 1 column",
          "anisotropic so range_X_coeff and range_PP_coeff (if applicable) have 3 columns"
        )[1 + gns_simulator$anisotropic]
      )
    )
  }
  # fixed effects for the range
  if (!is.null(range_intercept)) {
    res$range_X_coeff[1, 1] <- range_intercept
    if (verbose) {
      message(
        paste(
          "The range intercept (range_X_coeff[1,1]) has been set to the chosen value of",
          range_intercept
        )
      )
    }
  } else {
    res$range_X_coeff[1, 1] <-
      -0.5 * log(8 * gns_simulator$matern_smoothness) +
      log(max(dist(gns_simulator$vecchia_approx$locs[seq(min(nrow(gns_simulator$vecchia_approx$locs), 10000)), ])) / 20)
    if (verbose) {
      message(
        paste(
          "The log range intercept (range_X_coeff[1,1]) has been set automatically to",
          round(res$range_X_coeff[1, 1], 3)
        )
      )
    }
  }
  if (verbose) {
    message(
      "Except for the range intercept (range_X_coeff[1,1]), ",
      "all the coefficients in range_X_coeff are still set to 0"
    )
  }
  # PP effects for the range
  if (!is.null(gns_simulator$range_PP)) {
    # log_var
    if (is.null(range_PP_log_var)) {
      range_PP_log_var <- rep(-6, 1 + gns_simulator$anisotropic)
      if (verbose) {
        message(
          "range_PP_log_var has been set to c(",
          do.call(paste, as.list(range_PP_log_var)),
          ")"
        )
      }
    }
    if (length(range_PP_log_var) != 1 + gns_simulator$anisotropic) {
      stop("range_PP_log_var must have length 1 if isotropic, 2 if anisotropic")
    }
    if (!is.numeric(range_PP_log_var)) {
      stop("range_PP_log_var must be numeric")
    }
    # PP coeffs
    if (!is.null(gns_simulator$range_PP)) {
      res$range_PP_coeff <- matrix(
        0,
        nrow = gns_simulator$range_PP$n_knots,
        ncol = 1 + 2 * gns_simulator$anisotropic,
        dimnames = list(paste0("PP", seq(gns_simulator$range_PP$n_knots)))
      )
      res$range_PP_coeff[, 1] <- exp(.5 * range_PP_log_var[1]) * rnorm(nrow(res$range_PP_coef))
      if (verbose) {
        message(
          "The first column of range_PP_coeff (who parametrizes brute range) ",
          "has been filled with independent Normal coefficients with ",
          "log-variance range_PP_log_var[1]=", range_PP_log_var[1]
        )
      }
      if (gns_simulator$anisotropic) {
        res$range_PP_coeff[, c(2, 3)] <- exp(.5 * range_PP_log_var[2]) * rnorm(2 * nrow(res$range_PP_coeff))
        if (verbose) {
          message(
            "The second and third columns of range_PP_coeff (who parametrize anisotropy) have been filled with independent Normal coefficients with log-variance range_PP_log_var[2]=",
            range_PP_log_var[2]
          )
        }
      }
    }
  }

  if (verbose) {
    message("------------------")
  }
  if (verbose) {
    message("NOISE PARAMETERS")
  }
  # fixed effects for the noise
  if (!is.null(noise_intercept)) {
    res$noise_X_coeff[1, 1] <- noise_intercept
  } else {
    res$noise_X_coeff[1, 1] <- 0
    if (verbose) {
      message("The noise intercept (noise_X_coeff[1]) has been set to 0")
    }
  }
  if (verbose) {
    message("All the rest of the coefficients in noise_X_coeff are still set to 0")
  }

  # PP effects for the noise
  if (!is.null(gns_simulator$noise_PP)) {
    # log_var
    if (is.null(noise_PP_log_var)) {
      noise_PP_log_var <- -6
      if (verbose) {
        message("noise_PP_log_var has been set to ", noise_PP_log_var)
      }
    }
    if (length(noise_PP_log_var) != 1) {
      stop("noise_PP_log_var must have length 1")
    }
    if (!is.numeric(noise_PP_log_var)) {
      stop("noise_PP_log_var must be numeric")
    }
    # PP coeffs
    res$noise_PP_coeff <- matrix(
      0,
      nrow = gns_simulator$noise_PP$n_knots,
      ncol = 1,
      dimnames = list(paste(
        "PP", seq(gns_simulator$noise_PP$n_knots),
        sep = ""
      ))
    )
    res$noise_PP_coeff[, 1] <- exp(.5 * noise_PP_log_var[1]) * rnorm(nrow(res$noise_PP_coeff))
    if (verbose) {
      message(
        "The first column of noise_PP_coeff has been filled with independent Normal coefficients with log-variance noise_PP_log_var[1] = ",
        noise_PP_log_var[1]
      )
    }
  }
  if (verbose) {
    message("------------------")
    message("You can disable this message by setting 'verbose' to FALSE")
  }
  res
}


#' @title Simulates fake data
#' @description From a GeoNonStatSimulator object and
#' parameters at the right format,
#' simulates fake data
#'
#' @param gns_simulator  an object created by `createGnsSimulator()`
#' @param gns_params an object created by `createGnsSimulatorParameters()`
#' @param num_threads integer, number of treads to use. Default to 1.
#'
#' @returns A list containing the simulated data
#' @export
#'
#' @examples
#' set.seed(2)
#' nlocs <- 15000
#' observed_locs <- cbind(runif(nlocs), runif(nlocs))
#' vecchia_approx <- createVecchia(observed_locs, m = 6)
#'
#' # fixed effects
#' X <- data.frame(
#'   V1 = observed_locs[, 1],
#'   V2 = rnorm(nrow(observed_locs)),
#'   V3 = rnorm(nrow(observed_locs)),
#'   V4 = rnorm(nrow(observed_locs))
#' )
#'
#' PP_range <- createPP(vecchia_approx, knots = 16, matern_range = .3)
#' PP_noise <- createPP(vecchia_approx, knots = 500, matern_range = .06)
#'
#' # full monty
#' gns_simulator <- createGnsSimulator(
#'   vecchia_approx = vecchia_approx,
#'   X = X,
#'   noise_X = X, noise_PP = PP_noise,
#'   range_X = NULL, range_PP = PP_range,
#'   anisotropic = TRUE
#' )
#' gns_params <- createGnsSimulatorParameters(gns_simulator)
#' gns_data <- simulateGnsData(gns_simulator = gns_simulator, gns_params = gns_params)
simulateGnsData <- function(gns_simulator,
                            gns_params,
                            num_threads = 5) {
  log_range_field <- computeLogRange(
    range_beta = rbind(
      gns_params$range_X_coeff,
      gns_params$range_PP_coeff
    ),
    vecchia_approx = gns_simulator$vecchia_approx,
    range_X = gns_simulator$covariates$range_X,
    PP = gns_simulator$range_PP
  )


  log_noise_field <- xPPMultRight(
    X = gns_simulator$covariates$noise_X$X,
    PP = gns_simulator$noise_PP,
    vecchia_approx = gns_simulator$vecchia_approx,
    Y = rbind(
      gns_params$noise_X_coeff,
      gns_params$noise_PP_coeff
    ),
    permutate_PP_to_obs = TRUE
  )
  compressed_chol <- array(0,
    dim = c(
      nrow(gns_simulator$vecchia_approx$NNarray),
      gns_simulator$vecchia_approx$n_locs, 1
    )
  )
  vecchia_(
    log_range = t(log_range_field),
    locs = gns_simulator$vecchia_approx$t_locs,
    NNarray = gns_simulator$vecchia_approx$NNarray,
    smoothness = gns_simulator$matern_smoothness,
    compute_derivative = FALSE,
    num_threads = max(7, parallel::detectCores() - 1),
    result = compressed_chol
  )
  sparse_chol <- decompressChol(
    vecchia_approx = gns_simulator$vecchia_approx,
    compressed_sparse_chol = compressed_chol
  )

  latent_field <- as.vector(
    exp(.5 * gns_params$field_log_var) * Matrix::solve(sparse_chol, rnorm(nrow(sparse_chol)))
  )
  noise <- as.vector(exp(.5 * log_noise_field) * rnorm(length(log_noise_field)))
  fixed_effects <- as.vector(gns_simulator$covariates$X$X %*% gns_params$X_coeff)

  observed_field <-
    as.vector(t(gns_simulator$vecchia_approx$locs_match_matrix) %*% latent_field) +
    noise +
    fixed_effects

  res <- list(
    "log_range_field" = log_range_field,
    "log_noise_var_field" = log_noise_field,
    "sparse_chol" = sparse_chol,
    "noise" = noise,
    "fixed_effects" = fixed_effects,
    "latent_field" = latent_field,
    "observed_field" = observed_field
  )
}
