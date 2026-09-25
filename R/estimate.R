

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
#' GeoNonStat::aggregateRecords(processedGnsDemo$records)
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
estimate <- function(object, burn_in = 0.2) {
  samples <- aggregateRecords(object$records, burn_in = burn_in, keep ="all")
  samples$field <- samples$field[, object$vecchia_approx$locs_match]
  samples$fixed_effects <- as.matrix(samples$beta %*% t(object$covariates$X$X))
  samples$denoised <- samples$field + samples$fixed_effects
  samples$noise_log_var <- samples$denoised
  # TODO tenter de paralleliser avec future. 
  for(i in seq_len(nrow(samples$noise_beta))){
    samples$noise_log_var[i,] <- log(noiseVar(
      X = object$covariates$noise_X$X, 
      PP = object$hierarchical_model$noise$PP, 
      vecchia_approx = object$vecchia_approx, 
      noise_beta = samples$noise_beta[i,]))
  }
  
  summaries <- list()
  for(name in names(samples)) summaries[[name]] <- GeoNonStat::summarizeRecords(samples[[name]])
  return(list(samples = samples, summaries = summaries))
}


