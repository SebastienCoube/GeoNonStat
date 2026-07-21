#' Run the MCMC on one chain of a `GeoNonStat` object
#'
#' @param covariates  The list of covariates obtained with `processCovariates()`
#' @param observed_field a vector of observations of the interest variable
#' @param hierarchical_model a hierarchical model (obtained with `processHierarchicalModel()`)
#' @param vecchia_approx an object created by `vecchiaApprox()`
#' @param state state (from a chain of GeoNonstat object) to update
#' @param n_iterations numeric value, number of iterations of MCMC. Default to 100
#' @param num_threads number of threads. Integer.
#' @param iter_start integer, iteration where to start (if some iterations are 
#' already computed, begin after the last one)
#'
#' @returns a list containaing last state and records of parameters for each
#' MCMC iteration
#' @export
#' @keywords internal
#'
#' @examples
#' chain_id <- 1
#' iter_start <- length(gnsDemo$records$chain_1)
#' processedChain <-
#'   runOneChainMcmc(
#'     covariates = gnsDemo$covariates,
#'     observed_field = gnsDemo$observed_field,
#'     hierarchical_model = gnsDemo$hierarchical_model,
#'     vecchia_approx = gnsDemo$vecchia_approx,
#'     state = gnsDemo$states[[chain_id]],
#'     n_iterations = 20,
#'     num_threads = 1,
#'     iter_start = iter_start
#'   )
#' names(processedChain)
runOneChainMcmc <- function(
  covariates,
  observed_field,
  hierarchical_model,
  vecchia_approx,
  state,
  n_iterations = 100,
  num_threads = 1,
  iter_start
) {
  # Initialization of parameters (empty whith right structure)
  params_records <- lapply(
    seq_len(n_iterations),
    function(i) {
      initStateParams(state$params)
    }
  )

  for (iter in seq_len(n_iterations)) {
    # Regression coefficients #################################
    res <- updateBeta(state$params, state$stuff, covariates$X, vecchia_approx, observed_field)
    state$params <- res[["params"]]
    state$stuff$lm_fit <- res[["stuff"]]$lm_fit
    state$stuff$lm_residuals <- res[["stuff"]]$lm_residuals
    # Latent field ###############################
    res <- updateLatentField(
      state = state, hierarchical_model = hierarchical_model,
      vecchia_approx = vecchia_approx, observed_field = observed_field,
      iter = iter, num_threads = num_threads
    )
    state$params$field <- res$field

    if (iter + iter_start > 30) {
       # Field variance ###############################
      res <- updateFieldLogVar(state, hierarchical_model$scale, vecchia_approx, iter, iter_start)
      state$params$field <- res$params$field
      state$params$field_log_var <- res$params$field_log_var
      state$ker_var <- res$ker_var
    }
    if (iter + iter_start > 60) {
       # Range beta ###############################
       res <- updateRangeBetaAndLogVarMALA(
         state, hierarchical_model, vecchia_approx, 
         range_X = covariates$range_X, iter, 
         iter_start, num_threads)
       state <- res[["state"]]
       # Variance of the  range PP ###############################
       if(!is.null(hierarchical_model$range$PP)){
         state$params$range_log_scale = updateVarPPSuff(
           hm4params = hierarchical_model$range, 
           beta4params = state$params$range_beta, current_range_log_scale = state$params$range_log_scale)
         
         #state <- updateVarPPAncillaryXSufficient(
         #  state, 
         #  hierarchical_model, 
         #  vecchia_approx, range_X = covariates$range_X, iter, iter_start, num_threads)
         
         state$params$range_log_scale = updateVarPPSuff(
           hm4params = hierarchical_model$range, 
           beta4params = state$params$range_beta, current_range_log_scale = state$params$range_log_scale)
       }
    }

    # Noise ###############################
    # Noise beta ###############################

    res <- updateNoiseBetaFisher(state, hierarchical_model$noise, covariates$noise_X, vecchia_approx, iter, iter_start)
    state <- res[["state"]]

    # Noise log scale ###############################
    if (!is.null(hierarchical_model$noise$PP)) {
      state$params$noise_log_scale <- updateVarPPSuff(
        hm4params = hierarchical_model$noise,
        beta4params = state$params$noise_beta,
        current_range_log_scale = state$params$noise_log_scale
      )
    }
    # Storing the samples ###############################
    params_records[[iter]] <- state$params
  }
  return(list("state" = state, "params_records" = params_records))
}


#' Run the MCMC on a `GeoNonStat` object, with parallel execution on chains
#'
#' @param object an object of class `GeoNonStat`
#' @param n_chains_in_parallel numeric, number of chains in parallel, default to NULL
#' @param n_threads_per_chain numeric, number of threads by markov chain, default to 5
#' @param n_iterations numeric value, number of iterations of MCMC. Default to 100
#' @details The MCMC is run using [future::plan()].
#'
#' \preformatted{
#'  if (n_chains_in_parallel == 1) {
#'    # For no parallelisation
#'    future::plan(future::sequential)
#'  } else {
#'    # For multisession or multicore (prefered):
#'    if (parallelly::supportsMulticore()) {
#'      future::plan(future::multicore, workers = n_chains_in_parallel)
#'    } else {
#'      future::plan(future::multisession, workers = n_chains_in_parallel)
#'    }
#'  }
#' }
#' You can see `?future::plan` to define your plan more precisely.
#'
#' @returns a GeoNonStat object, with updated state (last state by chain)
#' and records (params records for all iterations, by chain)
#' @export
#'
#' @examples
#' processedGns <- GeoNonStatMcmc(
#'   gnsDemo,
#'   n_chains_in_parallel = 1,
#'   n_threads_per_chain = 1, # example in sequential
#'   n_iterations = 15
#' )
GeoNonStatMcmc <- function(
  object,
  n_chains_in_parallel = NULL,
  n_threads_per_chain = 5,
  n_iterations = 100
) {
  if (is.null(n_chains_in_parallel)) n_chains_in_parallel <- length(object$states)
  iter_start <- length(object$records$chain_1)
  if (n_iterations + iter_start < 60) {
    warning(
      "The algorithm is implemented to update the spatial range parameters ",
      "after at least 60 iterations. You should increase n_iterations."
    )
  }
  usedPlan <- utils::capture.output({
    print(future::plan())
  })[1]
  cat("-------- The parallelism on chains is", usedPlan, "--------\n")
  cat("         you can change it by using                        \n")
  cat("         future::plan(multicore, workers=x)    # if supported \n")
  cat("         future::plan(multisession, workers=x) # on Windows or RStudio\n")
  cat("         future::plan(sequential)              # no parallelism\n")
  cat("-------- MCMC Running..... --------\n")
  chains <- seq_along(object$states)

  covariates_local <- object$covariates
  observed_field_local <- object$observed_field
  hierarchical_model_local <- object$hierarchical_model
  vecchia_approx_local <- object$vecchia_approx
  states_local <- object$states

  res <- future.apply::future_mapply(
    FUN = function(states_local) {
      runOneChainMcmc(
        state = states_local,
        covariates = covariates_local,
        observed_field = observed_field_local,
        hierarchical_model = hierarchical_model_local,
        vecchia_approx = vecchia_approx_local,
        iter_start = iter_start,
        n_iterations = n_iterations,
        num_threads = n_threads_per_chain
      )
    },
    object$states,
    SIMPLIFY = FALSE,
    future.seed = TRUE,
    future.globals = list(
      covariates_local = covariates_local,
      observed_field_local = observed_field_local,
      hierarchical_model_local = hierarchical_model_local,
      vecchia_approx_local = vecchia_approx_local,
      iter_start = iter_start,
      n_iterations = n_iterations,
      n_threads_per_chain = n_threads_per_chain,
      runOneChainMcmc = runOneChainMcmc
    )
  )

  # Returning object processed, with for each chain : last state  and records for all iterations
  new_object <- object
  new_object$states <- lapply(res, function(x) x[["state"]])
  newrecords <- lapply(res, function(x) x[["params_records"]])
  new_object$records <- mapply(function(x1, x2) c(x1, x2), object$records, newrecords, SIMPLIFY = FALSE)
  return(new_object)
}
