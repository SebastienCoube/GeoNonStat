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
  # settings for MCMC 
  MALA_range_beta <- 3
  begin_range_update <- 50
  fisher = list(begin_learn = 100, begin_introduce = 150, end_introduce = 180)
  begin_var_update <- 30
  end_var_update <- 250
  for (iter in seq_len(n_iterations)) {
    # Latent field ###############################
    state$params$field <- updateLatentField(
      state = state, hierarchical_model = hierarchical_model,
      vecchia_approx = vecchia_approx, observed_field = observed_field,
      iter = iter, num_threads = num_threads
    )$field
    # Regression coefficients #################################
    res <- updateBeta(state$params, state$stuff, covariates$X, vecchia_approx, observed_field)
    state$params <- res$params
    state$stuff$lm_fit <- res$stuff$lm_fit; state$stuff$lm_residuals <- res$stuff$lm_residuals
    
    # Latent field ###############################
    state$params$field <- updateLatentField(
      state = state, hierarchical_model = hierarchical_model,
      vecchia_approx = vecchia_approx, observed_field = observed_field,
      iter = iter, num_threads = num_threads
    )$field
    if ((iter + iter_start > begin_var_update) & (iter + iter_start > end_var_update)) {
       # Field variance ###############################
      res <- updateFieldLogVar(state, hierarchical_model$scale, vecchia_approx, iter, iter_start)
      state$params$field[] <- res$params$field[]
      state$params$field_log_var <- res$params$field_log_var
      state$ker_var <- res$ker_var
    }
    if (iter + iter_start > begin_range_update) {
       # Range beta ###############################
      state <- rangeBetaS(
         state, hierarchical_model, vecchia_approx, 
         range_X = covariates$range_X, iter, iter_start, num_threads,
         fisher = fisher, n_MALA = MALA_range_beta
         )
      state <- rangeBetaA(
         state, hierarchical_model, vecchia_approx, 
         range_X = covariates$range_X, iter, iter_start, num_threads,
         fisher = fisher, n_MALA = MALA_range_beta
         )
       # Variance of the  range PP ###############################
       if(!is.null(hierarchical_model$range$PP)){
         state <- rangeLogScaleS(
           state, hierarchical_model, vecchia_approx, 
           range_X = covariates$range_X, iter, iter_start, num_threads)
         state$params$range_log_scale <- updateVarPPSuff(
           hm4params = hierarchical_model$range, 
           beta4params = state$params$range_beta, current_range_log_scale = state$params$range_log_scale)
         state <- rangeLogScaleA(
           state, hierarchical_model, vecchia_approx, 
           range_X = covariates$range_X, iter, iter_start, num_threads)
         state$params$range_log_scale <- updateVarPPSuff(
           hm4params = hierarchical_model$range, 
           beta4params = state$params$range_beta, current_range_log_scale = state$params$range_log_scale)
       }
    }

    # Latent field ###############################
    state$params$field <- updateLatentField(
      state = state, hierarchical_model = hierarchical_model,
      vecchia_approx = vecchia_approx, observed_field = observed_field,
      iter = iter, num_threads = num_threads
    )$field
    
    # Noise ###############################
    # Noise beta ###############################
    state <- noiseBeta(state, hm_noise = hierarchical_model$noise, 
                       noise_X = covariates$noise_X, 
                       vecchia_approx, iter, iter_start, 
                       fisher = fisher)
    
    # Noise log scale ###############################
    if (!is.null(hierarchical_model$noise$PP)) {
      state$params$noise_log_scale <- updateVarPPSuff(
        hm4params = hierarchical_model$noise,
        beta4params = state$params$noise_beta,
        current_range_log_scale = state$params$noise_log_scale
      )
      state <- noiseLogScaleAncillary(
        state, hierarchical_model, vecchia_approx, 
        noise_X = covariates$noise_X, iter, iter_start)
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
  n_iterations = 100, 
  verbose = T
) {
  if (is.null(n_chains_in_parallel)) n_chains_in_parallel <- length(object$states)
  iter_start <- length(object$records$chain_1)
  usedPlan <- utils::capture.output({
    print(future::plan())
  })[1]
  if(verbose){
  cat("-------- The parallelism on chains is", usedPlan, "--------\n")
  cat("         you can change it by using                        \n")
  cat("         future::plan(multicore, workers=x)    # if supported \n")
  cat("         future::plan(multisession, workers=x) # on Windows or RStudio\n")
  cat("         future::plan(sequential)              # no parallelism\n")
  cat("-------- MCMC Running..... --------\n")
  }
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



#' @export
automaticMcmc <- function(
    object,
    n_chains_in_parallel = NULL,
    n_threads_per_chain = 5, 
    satisfying_ESS = 100,
    satisfying_Gelman_Rubin = 1.1, 
    iter_per_step = 500,
    burn_in = .3, 
    verbose = T
){
  satisfying_ESS <- max(satisfying_ESS, 0)
  satisfying_Gelman_Rubin <- max(satisfying_Gelman_Rubin, 1.01)
  if(mcmcCount(object)>0){
    mcmc_diags <- mcmcDiags(object, burn_in, verbose = F)
    worst_gr <- mcmc_diags$worst$value[match("Point est. Gelman", mcmc_diags$worst$criterium)]
    worst_ess <- mcmc_diags$worst$value[match("ess", mcmc_diags$worst$criterium)]
  }
  if(mcmcCount(object)==0){
    worst_ess <- 0 
    worst_gr <- Inf 
  }
  if(verbose)cat("Starting auto MCMC... ")
  while(worst_gr > satisfying_Gelman_Rubin | worst_ess < satisfying_ESS){
    cat(paste(mcmcCount(object), "iterations already done, worst Gelman Rubin R-hat =", 
                  worst_gr, ", worst ESS =", worst_ess, ", going on for", iter_per_step, "more MCMC iterations."))
    object <- GeoNonStatMcmc(
      object,
      n_chains_in_parallel = NULL,
      n_threads_per_chain = 5,
      n_iterations = iter_per_step, 
    )
    mcmc_diags <- mcmcDiags(object, burn_in, verbose = F)
    worst_gr <- mcmc_diags$worst$value[match("Point est. Gelman", mcmc_diags$worst$criterium)]
    worst_ess <- mcmc_diags$worst$value[match("ess", mcmc_diags$worst$criterium)]
  }
  cat(paste("Stopping at", mcmcCount(object), "iterations, worst Gelman Rubin R-hat =", 
                worst_gr, ", worst ESS =", worst_ess))
  return(object)
}
