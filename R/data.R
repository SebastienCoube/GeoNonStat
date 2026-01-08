#' @title
#' GeoNonStat demo data
#'
#' @description
#' A S3 GeoNonStat object TODO
#' Created with the following code :
#' 
#' set.seed(2)
#' nlocs = 3000
#' observed_locs = cbind(runif(nlocs), runif(nlocs))
#' vecchia_approx = createVecchia(observed_locs, m = 6)
#' 
#' # fixed effects
#' X = as.data.frame(cbind(observed_locs[,1], rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs))))
#' 
#' PP_range = createPP(vecchia_approx, knots = 16, matern_range = .3)
#' PP_noise = createPP(vecchia_approx, knots = 30, matern_range = .25)
#' 
#' gns_simulator = createGnsSimulator(
#'   vecchia_approx = vecchia_approx,
#'   X = X, 
#'   noise_X = X, noise_PP = PP_noise, 
#'   range_X = NULL,range_PP = PP_range, 
#'   anisotropic = T
#' )
#' 
#' field_log_var=  2
#' # case with nice non-stationarity and little noise. 
#' range_log_scale = c(-1,-1)
#' noise_intercept = 1
#' noise_PP_log_var = -1
#' coeff_list = createGnsSimulatorParameters(
#'   gns_simulator, 
#'   range_PP_log_var =range_log_scale, 
#'   noise_intercept = noise_intercept, 
#'   noise_PP_log_var = noise_PP_log_var, 
#'   field_log_var = field_log_var)
#' coeff_list$range_X_coeff[1] = -5
#' 
#' fakeData = simulateGnsData(
#' gns_simulator = gns_simulator, 
#' gns_params = coeff_list)
#' 
#' gnsDemo = GeoNonStat(
#'   vecchia_approx = vecchia_approx, 
#'   observed_field = fakeData$observed_field, X = X, 
#'   matern_smoothness = 1.5, anisotropic = T, 
#'   n_chains = 3, 
#'   noise_X = NULL, range_X = NULL,
#'   noise_PP = PP_noise,
#'   range_PP = PP_range
#' )
#' @format
#' #TODO
#'
#' @usage gnsDemo
#'
#' @return
#' A \code{\link{GeoNonStat}}.
#'
#' @examples
#' gnsDemo
"gnsDemo"


#' @title
#' GeoNonStat demo data
#'
#' @description
#' A S3 GeoNonStat object TODO
#' set.seed(123)
#' processedGnsDemo = run_parallel_version_goret(
#'   object = gnsDemo, n_chains_in_parallel = 1,
#'   n_threads_per_chain = 6, n_iterations = 15, seed = 1)
#' @format
#' #TODO
#'
#' @usage gnsDemo
#'
#' @return
#' A \code{\link{GeoNonStat}} with processed MCMC.
#'
#' @examples
#' processedGnsDemo
"processedGnsDemo"