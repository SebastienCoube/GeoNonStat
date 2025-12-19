PrintMcmcDiags = function(object, burn_in = .1, min_ESS = 50, max_Gelman_upper = 1.05){
  records <- aggregateRecords(object, burn_in = burn_in, keep="all", keep_separate_chains = TRUE)
  ess = lapply(records, function(x)lapply(x, function(x)coda::effectiveSize(x)))
  ess = lapply(ess, unlist)
  ess = Reduce("+", ess)
  ess = signif(ess, 2)
  
  records = lapply(records, function(x) do.call("cbind", x))
  records = lapply(records, coda::as.mcmc)
  records = coda::as.mcmc.list(records)
  gelman_diags = t(coda::gelman.diag(records, multivariate = F)[[1]])
  gelman_diags = signif(gelman_diags, 2)
  res = t(rbind(ess, gelman_diags)) 
  worst = c(min(res[,1]), max(res[,2]), max(res[,3])) 
  which.worst = matrix(
    c(row.names(res)[c(which.min(res[,1]), which.max(res[,2]), which.max(res[,3]))]),
    nrow=1)
  res = rbind("worst" = worst, res)
  colnames(res)  <- colnames(which.worst) <- c("ess", "Point est. Gelman", "Upper C.I. Gelman")
  print(which.worst)
  print(res)
  return(list("diags"=res, "which.worst"=which.worst))
}

TracePlots = function(object, burn_in = .1, keep = "all", mfrow = c(3,2)){
  par(mfrow = mfrow)
  n_iters <- length(object$records[[1]])
  records = aggregateRecords(object, burn_in, keep = keep, keep_separate_chains = TRUE)
  records <- transposeList(records)
  plotted_iters = seq(max(0, floor(n_iters * burn_in))+1,n_iters)
  xlim <- c(min(plotted_iters), max(plotted_iters))
  nbchains <- length()
  for(c in seq_len(ncol(records))){
     for(name in names(records)[[1]]){
       to_be_plotted
      to_be_plotted = records[[c]][[name]]
      plot(0, 0, 
           ylim = c(min(to_be_plotted), max(to_be_plotted)), 
           xlim = xlim, 
           type = "n", 
           main  = paste("trace plot of", name, colnames(records[[name]])))
      for(j in seq_len(ncol(to_be_plotted))){
        lines(plotted_iters, to_be_plotted[,j])
      }
    }
  }
}
