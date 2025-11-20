ExtendVecchia = function(vecchia_approx, new_locs){
  
}

Predict = function(geo_non_stat, burn_in = .1, 
                   new_locs, 
                   new_X = NULL,
                   new_noise_X = NULL,
                   new_range_X = NULL
                   ){
  aggregated_records = AggregateRecords(geo_non_stat, burn_in)
  
  res = lapply(res, summarize)
  return(res)
}


geo_non_stat$covariates$X$