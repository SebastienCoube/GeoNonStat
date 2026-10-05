

#' Aggregate records of a GeoNonStat object, for 1 chain
#'
#' @param object a GeoNonStat object containing records
#' @param burn_in numeric, value, between 0 and 1. Gives the proportion of
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
#' aggregateRecords(processedGnsDemo$records)
aggregateRecords <- function(records, burn_in = .3, keep = "all", keep_separate_chains = FALSE){
  n_iter <- length(records[[1]])
  if(burn_in<0)stop("burn_in must either be a proportion between 0 and 1 or a number of iterations greater than 1")
  if(burn_in>1)burn_in <- burn_in/n_iter
  # separate chains
   if(keep_separate_chains){
     res <- list()
     for(i in seq(length(records))){
       res[[i]] <- aggregateRecords(records[i], burn_in, keep)
     }
     names(res) <- names(records)
     return(res)
   }
  # pooling chains
  n_chains <- length(records)
  n_iter_kept <- ceiling(n_iter*(1-burn_in))
  iter_start <- n_iter - n_iter_kept 
  res = list()
  
  params_to_aggregate <- names(records[[1]][[1]])
  if(keep[1] == "nofield") keep <- setdiff(params_to_aggregate, "field")
  if(keep[1] == "all") keep <- names(records[[1]][[1]])
  params_to_aggregate <- intersect(keep, params_to_aggregate)
  # putting log scales at end
  params_to_aggregate <- c(params_to_aggregate[!grepl("log_scale", params_to_aggregate)], params_to_aggregate[grepl("log_scale", params_to_aggregate)]) 
  
  for(param_name in  params_to_aggregate){
    for(aniso_dim in seq(dim(records[[1]][[1]][[param_name]])[2])){
      param_dim <- dim(records[[1]][[1]][[param_name]])[1]
      param_dimnames <- row.names((records[[1]][[1]][[param_name]]))
      resname <- trimws(paste(param_name, dimnames(records[[1]][[1]][[param_name]])[[2]][aniso_dim]))
      res[[resname]] <- matrix(
        0, nrow = n_iter_kept * n_chains, 
        ncol = param_dim, 
        dimnames = list(NULL, param_dimnames)
      )
      for(chain_idx in seq_len(n_chains)){
       for(iter in seq_len(n_iter_kept)){
         res[[resname]][(chain_idx-1)*n_iter_kept + iter,] <-
           records[[chain_idx]][[iter + iter_start]][[param_name]][,aniso_dim]
       }
      }
    }
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
#' @keywords internal
#'
#' @examples
#' tt <- matrix(rnorm(200), ncol = 20)
#' summarizeRecords(tt)
summarizeRecords <- function(mat, quant = c("q 2.5%" = .025, "median" = .5, "q 97.5%" = .975)) {
  means <- matrixStats::colMeans2(mat)
  sds   <- matrixStats::colSds(mat)
  qs    <- matrixStats::colQuantiles(mat, probs = quant)  # ncol(mat) x length(quant)
  if(ncol(mat)>1) {
    matstat <- rbind(means, sds, t(qs))
  } else {
    matstat <- matrix(c(means, sds, qs), nrow = 2 + length(quant), ncol = ncol(mat))
    
  }
  matstat <- signif(matstat, 3)
  colnames(matstat) <- colnames(mat)
  rownames(matstat) <- c("mean", "sd", names(quant))
  return(matstat)
}


#' Estimation of the model parameters using the MCMC records of a GeoNonStat object.
#'
#' @param object a GeoNonStat object containing records
#' @param burn_in numeric, value, between 0 and 1. Gives the proportion of
#' first records that should be removed. Default to 0.1
#' @param keep character vector. What params should be kept. Can also be
#' `"all"`, to keep all parameters, or `"all_even_the_field"` to keep all
#' parameters, event the `field.`
#'
#' @returns a list of numerical summary by parameters.
#' @export
#'
#' @examples
#' estimate(processedGnsDemo)
estimate <- function(object, burn_in = 0.3) {
  summaries <- list()
  scores <- list()
  scores$per_obs <- matrix(0, 4, object$vecchia_approx$n_obs)
  # getting MCMC samples
  samples <- aggregateRecords(object$records, burn_in = burn_in, keep ="all")
  # getting samples of noise PP at unique locations
  if(!is.null(object$hierarchical_model$noise$PP)){
    samples$noise_PP <- samples$field
    for(i in seq_len(nrow(samples$noise_beta))){
      samples$noise_PP[i,] <- xPPMultRight(
        X = NULL, vecchia_approx = object$vecchia_approx, 
        PP = object$hierarchical_model$noise$PP, 
        Y = t(samples$noise_beta[i, -seq(object$covariates$noise_X$n_regressors),drop=F]),
        permutate_PP_to_obs = F
      )
    }
    gc()
  }
  # getting summaries of MCMC samples
  summaries$hyperparameters <- lapply(samples, summarizeRecords)
  summaries$at_unique_locs <- list()
  summaries$at_unique_locs$field <- summaries$hyperparameters$field; summaries$hyperparameters$field <- NULL
  if(!is.null(object$hierarchical_model$noise$PP)){
    summaries$at_unique_locs$noise_PP <- summaries$hyperparameters$noise_PP; summaries$hyperparameters$noise_PP <- NULL
  }
  
  # working on observed locations
  summaries$at_observed_locs <- list()
  summaries$at_observed_locs$field <- summaries$at_unique_locs$field[,object$vecchia_approx$locs_match]
  # initializing summaries at observed locs
  summaries$at_observed_locs$field_and_fixed <- summaries$at_observed_locs$field
  summaries$at_observed_locs$noise_log_var <- summaries$at_observed_locs$field
  # computing over chunks to avoid RAM filling
  chunk_size <- 10000
  chunks <- split(
    seq_len(object$vecchia_approx$n_obs), 
    seq_len(object$vecchia_approx$n_obs)%/%chunk_size)
  for(chunk in chunks){
    fixed_effect_sample <- samples$beta %*% t(object$covariates$X$X[chunk,])
    latent_field_sample <- samples$field[,object$vecchia_approx$locs_match[chunk]]
    field_and_fixed_sample <- fixed_effect_sample + latent_field_sample
    log_noise_variance_sample <-
      samples$noise_beta[,seq(object$covariates$noise_X$n_regressors)] %*%
      t(object$covariates$noise_X$X[chunk,])
    if(!is.null(object$hierarchical_model$noise$PP)){
      log_noise_variance_sample <-
        log_noise_variance_sample +
        samples$noise_PP[,object$vecchia_approx$locs_match[chunk]]
    }
    summaries$at_observed_locs$field_and_fixed[,chunk] <-
      summarizeRecords(field_and_fixed_sample)
    summaries$at_observed_locs$noise_log_var[,chunk] <-
      summarizeRecords(log_noise_variance_sample)
    chunk_scores <- t(allScores(
      true_y = object$observed_field[chunk], 
      mean_samples = t(field_and_fixed_sample), 
      sd_samples = t(exp(.5 * log_noise_variance_sample))
    )$scores_per_obs)
    scores$per_obs[,chunk] <- chunk_scores
  }
  row.names(scores$per_obs) <- row.names(chunk_scores)
  scores$aggregated <- apply(scores$per_obs, 1, mean)
  gc()

  return(list(summaries = summaries, scores = scores))
}


