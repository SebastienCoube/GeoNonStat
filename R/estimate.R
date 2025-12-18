reduceRecords <- function(records, burn_in = 0, keep="all") {
  # Reduce paramaters
  namesparam <- names(records[[1]][[1]])
  keep <- switch(keep,
                "all" = setdiff(namesparam, "field"),
                "all_even_the_field" = namesparam,
                keep)
  records <- lapply(records, function(x) 
    lapply(x, function(y) y[keep]))
  
  # Reduce iterations
  n_iter <- length(records[[1]])
  start_iter <- max(0, floor(n_iter * burn_in))
  if(n_iter - start_iter <= 0)
    stop("Not enough iterations, reduce burn_in if necessary")
  if(!all(lengths(records) == n_iter)) 
    stop("Something is wrong with your object, ",
         "all the MCMC chains should have the same number of iterations")
  keep_iters <- (start_iter+1):n_iter
  
  return(
    lapply(records, function(x) {
      return(x[keep_iters])
      }
    )
  )
}

transposeList <- function(list) {
  lapply(seq_along(list[[1]]),
         function(i) lapply(list, `[[`, i))
}

AggregateRecordsPerChain = function(chain){
  res = lapply(chain, function(param){
    # Numeric case
    if(!is.matrix(param[[1]])) {
      param <- lapply(param, function(x) {
        dim(x) <- c(length(x), 1)
        return(x)
      })
    }
    
    if(ncol(param[[1]]) == 1) {
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

summarize = function(v){
  if(is.matrix(v)) return(apply(v,2,summarize))
  return(signif(c("mean" = mean(v), "sd" = sd(v), 
                  "q 2.5%" = unname(quantile(v, .025)), 
                  "median" = unname(quantile(v, .5)), 
                  "q 97.5%" = unname(quantile(v, .975))
  ), 3))
}

summarizevect = function(v, quant){
  vstat <- c(mean(v),
             sd(v),
            quantile(v, quant, names=FALSE))
  return(signif(vstat, 3))
}

summarizemat = function(mat, quant=c("q 2.5%" = .025, "median" = .5, "q 97.5%" = .975)){
  matstat <- apply(mat,2,summarizevect, quant)
  rownames(matstat) <- c("mean", "sd", names(quant))
  return(matstat)
}


Estimate <- function(object, burn_in = .1, keep = "all_even_the_field"){
  res <- reduceRecords(object$records, burn_in = burn_in, keep=keep)
  res <- lapply(res, transposeList)
  namesparams <- names(res)
  res <- lapply(res, AggregateRecordsPerChain)
  res <- do.call(Map, c(list(rbind), res))
  names(res) <- namesparams
  res <-  lapply(res, summarizemat)
  return(res)
}



Elppd = function(X, beta_samples, field_samples, noise_var_samples, observed_field){
  elppd_per_obs = (
    - log((sqrt(2)*pi)) # 1/sqrt(2 pi) passed to the log
    -.5*  log(noise_var_samples) # variance penalty of Gaussian density
    - .5 * (
      (X %*% t(beta_samples) + t(field_samples)) #mean
      - observed_field # observation
    )^2 / noise_var_samples # variance 
  )
  elppd_per_obs = apply(elppd_per_obs, 1, mean)
  return(list("elppd_per_obs" = elppd_per_obs, "mean_elppd" = mean(elppd_per_obs)))
}

ElppdTrain = function(geo_non_stat, burn_in = .1){
  aggregated_records = AggregateRecords(geo_non_stat, burn_in, who = c("beta", "noise_beta", "field"))
  noise_var = 
    apply(aggregated_records$noise_beta, 1, function(x)
        as.vector(exp(X_PP_mult_right(
        X = geo_non_stat$covariates$noise_X$X, PP = geo_non_stat$hierarchical_model$noise$PP,
        vecchia_approx = geo_non_stat$vecchia_approx, Y = x, 
        permutate_PP_to_obs = T
      )))
      )
  # Gaussian log-dens of the observations
  return(Elppd(
    geo_non_stat$covariates$X$X, 
    beta_samples = aggregated_records$beta, 
    field_samples = aggregated_records$field[,geo_non_stat$vecchia_approx$locs_match], 
    noise_var_samples = noise_var, 
    observed_field = geo_non_stat$observed_field))
}
