#' @title
#' GeoNonStat demo data
#'
#' @description
#' A S3 GeoNonStat object with nice non-stationarity and little noise.
#' generated for demonstration and examples purpose.
#' @format a GeoNonStat object with 3000 observed locations and field.
#' No MCMC has run on this object.
#'
#' @details
#' The data have been generated with the following code :
#'
#' \preformatted{
#' set.seed(2)
#' nlocs = 3000
#' observed_locs = cbind(runif(nlocs), runif(nlocs))
#' vecchia_approx = createVecchia(observed_locs, m = 6)
#'
#' # fixed effects
#' X = data.frame(
#'  V1 = observed_locs[,1],
#'  V2 = rnorm(nrow(observed_locs)),
#'  V3 = rnorm(nrow(observed_locs)),
#'  V4 = rnorm(nrow(observed_locs))
#' )
#'
#' PP_range = createPP(vecchia_approx,
#'                     knots = 16,
#'                     matern_range = .3,
#'                     plot=FALSE)
#' PP_noise = createPP(vecchia_approx,
#'                     knots = 30,
#'                     matern_range = .25,
#'                     plot=FALSE)
#'
#' gns_simulator = createGnsSimulator(
#'   vecchia_approx = vecchia_approx,
#'   X = X,
#'   noise_X = X, noise_PP = PP_noise,
#'   range_X = NULL,range_PP = PP_range,
#'   anisotropic = TRUE
#' )
#'
#' coeff_list = createGnsSimulatorParameters(
#'   gns_simulator,
#'   range_PP_log_var =c(-1,-1),
#'   noise_intercept = 1,
#'   noise_PP_log_var = -1,
#'   field_log_var = 2,
#'   verbose = FALSE)
#' coeff_list$range_X_coeff[1] = -5
#'
#' fake_data = simulateGnsData(
#'   gns_simulator = gns_simulator,
#'   gns_params = coeff_list)
#'
#' gnsDemo = GeoNonStat(
#'   vecchia_approx = vecchia_approx,
#'   observed_field = fake_data$observed_field, X = X,
#'   matern_smoothness = 1.5, anisotropic = T,
#'   n_chains = 3,
#'   noise_X = NULL, range_X = NULL,
#'   noise_PP = PP_noise,
#'   range_PP = PP_range
#' )
#' }
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
#' A S3 GeoNonStat object with nice non-stationarity and little noise.
#' generated for demonstration and examples purpose. MCMC has already run on
#' this object.
#' @format a GeoNonStat object with 3000 observed locations and field.
#' 15 iterations of MCMC have run on this object.
#'
#' @details
#' The data have been generated with the following code :
#'
#' \preformatted{
#'   set.seed(123)
#'   processedGnsDemo = runParallelGeoNonStatMcmc(
#'     object = gnsDemo,
#'     n_chains_in_parallel = 1,
#'     n_threads_per_chain = 1,
#'     n_iterations = 15,
#'     seed = 1
#'   )
#' }
#'
#' @return
#' A \code{\link{GeoNonStat}} with processed MCMC.
#'
#' @examples
#' processedGnsDemo
"processedGnsDemo"
