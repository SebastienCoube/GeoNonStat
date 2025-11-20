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
      for(iter in seq( max(1, floor(length(record) * burn_in)), length(record))){
        res[[name]][iter - floor(length(record) * burn_in),] = record[[iter]][[name]]
      }
    }
    res
  })
}

AggregateRecords = function(records, burn_in = .1, who = "all"){
  aggregated_records = AggregateRecordsPerChain(records = geo_non_stat$records, burn_in = burn_in, who = who)
  res=  list()
  for(name in names(aggregated_records[[1]]))res[[name]] = do.call(rbind, lapply(aggregated_records, function(x)x[[name]]))
  return(res)
}

PrintMcmcDiags = function(geo_non_stat, burn_in = .1, min_ESS = 50, max_Gelman_upper = 1.05){
  aggregated_records = AggregateRecordsPerChain(geo_non_stat$records, burn_in)
  ess = lapply(aggregated_records, function(x)lapply(x, function(x)coda::effectiveSize(x)))
  ess = lapply(ess, unlist)
  ess = Reduce("+", ess)
  ess = signif(ess, 2)
  
  aggregated_records_ = lapply(aggregated_records, function(x)do.call("cbind", x))
  aggregated_records_ = lapply(aggregated_records_, coda::as.mcmc)
  aggregated_records_ = coda::as.mcmc.list(aggregated_records_)
  gelman_diags = t(coda::gelman.diag(aggregated_records_, multivariate = F)[[1]])
  gelman_diags = signif(gelman_diags, 2)
  res = t(rbind(ess, gelman_diags)) 
  worst = c(min(res[,1]), max(res[,2]), max(res[,3])) 
  res = rbind("which.worst" = c(row.names(res)[which.min(res[,1])], row.names(res)[which.max(res[,2])], row.names(res)[which.max(res[,3])]), res)
  res = rbind("worst" = worst, res)
  colnames(res) = c("ess", "Point est. Gelman", "Upper C.I. Gelman")
  print(res)
}



TracePlots = function(geo_non_stat, burn_in = .1, who = "all", mfrow = c(3,3)){
  aggregated_records = AggregateRecordsPerChain(geo_non_stat$records, burn_in, who = who)
  plotted_iters = seq(
    ceiling(length(geo_non_stat$records[[1]])*(burn_in)), 
    length(geo_non_stat$records[[1]])
    )
  for(name in names(aggregated_records[[1]])){
    par(mfrow = mfrow)
    for(i in seq(ncol(aggregated_records[[1]][[name]]))){
      to_be_plotted = sapply(aggregated_records, function(x)x[[name]][,i])
      plot(0, 0, ylim = c(min(to_be_plotted), max(to_be_plotted)), 
           xlim = c(min(plotted_iters), max(plotted_iters)), 
           type = "n", main  = paste("trace plot of", name, colnames(aggregated_records[[1]][[name]])[i]))
      for(j in seq(ncol(to_be_plotted))){
        lines(plotted_iters, to_be_plotted[,j])
      }
    }
  }
}

