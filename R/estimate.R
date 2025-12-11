summarize = function(v){
  if(is.matrix(v))return(apply(v,2,summarize))
  return(signif(c("mean" = mean(v), "sd" = sd(v), 
           "q 2.5%" = unname(quantile(v, .025)), "median" = unname(quantile(v, .5)), "q 97.5%" = unname(quantile(v, .975))
           ), 3))
}

AggregateRecordsPerChain = function(records, burn_in = .1, who = "all"){
  if(identical(who, "all")) who = setdiff(names(records[[1]][[1]]), "field")
  if(identical(who, "all_even_the_field")) who = names(records[[1]][[1]])
  res = lapply(records, function(record){
    res = list()
    for(name in who){
      res[[name]] = matrix(0, ceiling(length(record)*(1-burn_in)), length(record[[1]][[name]]))
      if(!is.null( row.names(record[[1]][[name]]))){
        if(ncol(record[[1]][[name]])==1)colnames(res[[name]]) = row.names(record[[1]][[name]])
        if(ncol(record[[1]][[name]])>1)colnames(res[[name]]) = c(outer(row.names(record[[1]][[name]]), colnames(record[[1]][[name]]), function(x, y)paste(x, y, sep = "_")))
      }
      for(iter in seq(max(1, floor(length(record) * burn_in)), length(record))) {
        res[[name]][iter - floor(length(record) * burn_in),] = record[[iter]][[name]]
      }
    }
    res
  })
}

AggregateRecords = function(records, burn_in = .1, who = "all"){
  aggregated_records = AggregateRecordsPerChain(records = records, burn_in = burn_in, who = who)
  res=  list()
  for(name in names(aggregated_records[[1]]))res[[name]] = do.call(rbind, lapply(aggregated_records, function(x)x[[name]]))
  return(res)
}


Estimate = function(geo_non_stat, burn_in = .1, who = "all_even_the_field"){
  res = AggregateRecords(geo_non_stat, burn_in, who)
  res=  lapply(res, summarize)
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
