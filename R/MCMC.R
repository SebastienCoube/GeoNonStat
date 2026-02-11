#' Run the MCMC on one chain of a `GeoNonStat` object
#'
#' @param covariates  The list of covariates obtained with `processCovariates`
#' @param observed_field TODO
#' @param hierarchical_model a hierarchical model (obtained with `process_hierarchical_model()`)
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param state TODO
#' @param n_iterations_update TODO
#' @param num_threads TODO
#' @param iter_start TODO
#' @param seed integer value, seed used for reproducibility purposes. Default to 1
#'
#' @returns a list containaing last state and records of parameters for each 
#' MCMC iteration
#' @export
#' 
#' @examples
#' chain_id <- 1
#' seed <- 1
#' iter_start <- length(gnsDemo$records$chain_1)
#' processedChain <- 
#'     runGeoNonStatMcmc(
#'         covariates = gnsDemo$covariates, 
#'         observed_field = gnsDemo$observed_field, 
#'         hierarchical_model = gnsDemo$hierarchical_model, 
#'         vecchia_approx = gnsDemo$vecchia_approx, 
#'         state = gnsDemo$states[[chain_id]], 
#'         n_iterations_update = 20, 
#'         num_threads = 1, 
#'         iter_start = iter_start, 
#'         seed = iter_start + seed + chain_id
#'         )
#' names(processedChain)
runGeoNonStatMcmc = function(
    covariates,
    observed_field,
    hierarchical_model,
    vecchia_approx,
    state,
    n_iterations_update = 100, 
    num_threads = 1, 
    iter_start,
    seed=123 
) {
  set.seed(seed)
  # Initialization of parameters (empty whith right structure)
  params_records = lapply(seq_len(n_iterations_update), 
                          function(i) {initStateParams(state$params)}
  )
  
  for(iter in seq_len(n_iterations_update)){
    
    if(iter/10 == iter %/% 10)cat("iter = ", iter, "\n")
    # Regression coefficients ###############################
    # Regression coefficients #################################
    res <- updateBeta(state$params, state$stuff, covariates$X, vecchia_approx, observed_field)
    state$params <- res[["params"]]
    state$stuff <- res[["stuff"]]
    
    # Latent field ###############################
    state$params$field <- updateLatentField(state$params, state$stuff, vecchia_approx, observed_field, iter, num_threads)
    
    # Field log var ###############################
    state <- updateFieldLogVar(state, hierarchical_model$scale, vecchia_approx, iter, iter_start)
    
    if(iter + iter_start > 100){
      # Range beta ###############################
      res <- updateRangeBeta(state, hierarchical_model, vecchia_approx, range_X =  covariates$range_X, iter, iter_start, num_threads)
      state <- res[["state"]]
      
      # Variance of the  range PP ###############################
      state <- updateVarianceRangepp(
        state = state, hierarchical_model = hierarchical_model, range_X = covariates$range_X, 
        vecchia_approx = vecchia_approx, iter = iter, 
        iter_start = iter_start, num_threads = num_threads)
    }
    
    # Noise ###############################
    # Noise beta ###############################
    
    res <- updateNoiseBeta(state, hierarchical_model$noise, covariates$noise_X, vecchia_approx, iter, iter_start)
    state <- res[["state"]]
    
    # Noise log scale ###############################
    state <- updateNoiseLogScale(state, hierarchical_model$noise, covariates$noise_X, vecchia_approx, res[["squared_residuals"]], iter, iter_start)
    
    # Storing the samples ###############################
    params_records[[iter]] = state$params
  }
  return(list("state" = state, "params_records" = params_records))
}


#' Run the MCMC on a `GeoNonStat` object, with parallel execution on chains
#'
#' @param object an object of class `GeoNonStat`
#' @param n_chains_in_parallel numeric, number of chains in parallel, default to NULL
#' @param n_threads_per_chain numeric, number of threads by markov chain, default to 5
#' @param n_iterations numeric value, number of iterations. Default to 100
#' @param seed integer value, seed used for reproducibility purposes. Default to 1
#'
#' @returns a GeoNonStat object, with updated state (last state by chain) 
#' and records (params records for all iterations, by chain)
#' @export
#'
#' @examples
#' processedGns <- runParallelGeoNonStatMcmc(
#'      gnsDemo, 
#'      n_chains_in_parallel = 1,
#'      n_threads_per_chain = 1, # example in sequential
#'      n_iterations = 15)
runParallelGeoNonStatMcmc = function(
    object, 
    n_chains_in_parallel = NULL, 
    n_threads_per_chain = 5, 
    n_iterations = 100, 
    seed = 1){
  if(is.null(n_chains_in_parallel)) n_chains_in_parallel = length(object$states)
  iter_start = length(object$records$chain_1)
  
  # TODO : ça n'est pas recommandé par la communauté ici. 
  # Il faut laisser l'utilisateur choisir. 
  # On peut mettre le code ci-dessous dans les exemples
  # Et afficher le plan utiliser pour info. 
  if(n_chains_in_parallel == 1){
    future::plan(future::sequential)
  } else {
    if(parallelly::supportsMulticore()) {
      future::plan(future::multicore, workers=n_chains_in_parallel)
    } else {
      future::plan(future::multisession, workers=n_chains_in_parallel)
    }
  }
  usedPlan <- utils::capture.output({print(future::plan())})[1]
  cat("-------- The parallelism on chains is", usedPlan,"--------\n")
  cat("You can change it by using:\n")
  cat("       # For no parallelisation\n")
  cat("       future::plan(future::sequential) # No parallelisation\n")
  cat("       # For multisession or multicore (prefered):\n")
  cat("       if(parallelly::supportsMulticore()) {\n")
  cat("         future::plan(future::multicore, workers=n_chains_in_parallel)\n")
  cat("       } else {\n")
  cat("         future::plan(future::multisession, workers=n_chains_in_parallel)\n")
  cat("       }\n")
  cat("See 'help(future::plan)' for more info\n\n")
  cat("-------- MCMC Running..... --------\n")
  chains <- seq_along(object$states)
  
  covariates_local <- object$covariates
  observed_field_local <- object$observed_field
  hierarchical_model_local <- object$hierarchical_model
  vecchia_approx_local <- object$vecchia_approx
  states_local <- object$states
  
  tmpfun <- function(chain, ...) {
    runGeoNonStatMcmc(
      state = states_local[[chain]],
      seed = iter_start + seed + chain, 
      ...
    )
  }
    
  res <- future.apply::future_lapply(
    chains,
    FUN = tmpfun,
    covariates = covariates_local,
    observed_field = observed_field_local,
    hierarchical_model = hierarchical_model_local,
    vecchia_approx = vecchia_approx_local,
    iter_start = iter_start,
    n_iterations_update = n_iterations,
    num_threads = n_threads_per_chain,
    future.seed = TRUE
  )
  
  # Returning object processed, with for each chain : last state  and records for all iterations
  new_object <- object
  new_object$states <- lapply(res, function(x) x[["state"]])
  new_object$records <- lapply(res, function(x) x[["params_records"]])
  return(new_object)
}

