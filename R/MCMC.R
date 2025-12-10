#' Title TODO
#'
#' @param covariates  The list of covariates obtained with `process_covariates`
#' @param observed_field TODO
#' @param hierarchical_model a hierarchical model (obtained with `process_hierarchical_model()`)
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param state TODO
#' @param n_iterations_update TODO
#' @param num_threads TODO
#' @param iter_start TODO
#' @param seed integer value, seed used for reproducibility purposes. Default to 1
#'
#' @returns TODO
#' @export
#'
#' @examples
#' #TODO
#' #TODO
MessyMessyMcmc = function(
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
  
  # Regression coefficients ###############################
  # Regression coefficients #################################
  res <- update_beta(state$params, state$stuff, covariates$X, vecchia_approx, observed_field)
  state$params <- res[["params"]]
  state$stuff <- res[["stuff"]]
  
  # Latent field ###############################
  state$params$field <- update_latent_field(state$params, state$stuff, vecchia_approx, observed_field, iter, num_threads)
  
  # Field log var ###############################
  state <- update_field_log_var(state, hierarchical_model$scale, vecchia_approx, iter, iter_start)
  
  # Range beta ###############################
  
  res <- update_range_beta(state, hierarchical_model, covariates$range_X, iter, iter_start, num_threads)
  state <- res[["state"]]
  
  # Variance of the  range PP ###############################
  state <- update_variance_rangepp(state, hierarchical_model, covariates$range_X, vecchia_approx, res[["p"]], iter, iter_start, num_threads)
  
  
  # Noise ###############################
  # Noise beta ###############################
  
  res <- update_noise_beta (state, hierarchical_model$noise, covariates$noise_X, vecchia_approx, iter, iter_start)
  state <- res[["state"]]
  
  # Noise log scale ###############################
  state <- update_noise_log_scale(state, hierarchical_model$noise, covariates$noise_X, vecchia_approx, res[["squared_residuals"]], iter, iter_start)
  
  # Storing the samples ###############################
  params_records[[iter]] = state$params
}
  return(list("state" = state, "params_records" = params_records))
}


#' Title
#'
#' @param object an object of class `GeoNonStat`
#' @param n_chains_in_parallel numeric, number of chains in parallel, default to NULL
#' @param n_threads_per_chain numeric, number of threads by markov chain, default to 10
#' @param n_iterations numeric value, number of iterations. Default to 100
#' @param seed integer value, seed used for reproducibility purposes. Default to 1
#'
#' @returns a list
#' @export
#'
#' @examples
#' # TODO
run_parallel_version_goret = function(
    object, 
    n_chains_in_parallel = NULL, 
    n_threads_per_chain = 10, 
    n_iterations = 100, 
    seed = 1){
  if(is.null(n_chains_in_parallel)) n_chains_in_parallel = length(object$states)
  iter_start = length(object$records$chain_1)
  if(n_chains_in_parallel == 1){
    res <- lapply(seq_along(object$states), function(chain_idx) {
      list("state" = object$states[[chain_idx]],
           "params" = lapply(
             seq_len(n_iterations), 
             function(i) initStateParams(object$states[[chain_idx]][["params"]])
           )
      )
    })
    
    for(chain_idx in seq_along(object$states)){
      res[[chain_idx]] <- MessyMessyMcmc(
        covariates = object$covariates, 
        observed_field = object$observed_field, 
        hierarchical_model = object$hierarchical_model, 
        vecchia_approx = object$vecchia_approx, 
        state = object$states[[chain_idx]], 
        n_iterations_update = n_iterations, 
        num_threads = n_threads_per_chain, 
        iter_start = iter_start, 
        seed = iter_start + seed + chain_idx
      )
    }
  } else {
    cl = parallel::makeCluster(n_chains_in_parallel)
    parallel::clusterExport(cl = cl, varlist = c("iter_start", "seed", "object", "n_iterations"), envir = environment())
    res = parallel::parLapply(
      cl = cl, 
      X = seq(length(object$states)), function(chain_idx){
        MessyMessyMcmc(
          covariates = object$covariates, observed_field = object$observed_field, 
          hierarchical_model = object$hierarchical_model, vecchia_approx = object$vecchia_approx, 
          state = object$states[[chain_idx]], 
          n_iterations_update = n_iterations, 
          num_threads = n_threads_per_chain, 
          iter_start = iter_start, 
          seed = iter_start + seed + chain_idx
        )
      }
    )
  }
  return(res)
}

