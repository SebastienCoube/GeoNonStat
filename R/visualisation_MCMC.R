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

TraceOneParam = function(name, param, jstart, jend, xlim, color_lines){
  dimparam <- jend - jstart + 1
  ncol <- min(dimparam, ifelse(dimparam>9, 4, 3))
  nrow <- ceiling(dimparam / ncol)
  par(mfrow = c(nrow, ncol))
  
  # Iterate over columns
  for(col in seq(jstart, jend)) { 
    namecol <- colnames(param[[1]])[col]
    to_be_plotted <- do.call("cbind", lapply(param, function(x) x[,col]))
    matplot(to_be_plotted, 
            type="l", 
            lty = 1,
            lwd=2,
            main  = paste0(name, "\n", namecol),
            col = color_lines)
  }
}

split_intervals <- function(total_size, starts, nmax =  16) {
  intervals <- list()
  ends <- c(starts[-1] - 1, total_size)
  for (i in seq_along(starts)) {
    s <- starts[i]
    e <- ends[i]
    
    while ((e - s + 1) > nmax) {
      intervals <- append(intervals, list(c(s, s + nmax - 1)))
      s <- s + nmax
    }
    intervals <- append(intervals, list(c(s, e)))
  }
  intervals
}

TracePlots = function(object, burn_in = .1, keep = "all", mfrow = c(3,2)){
  ormfrow <- par("mfrow")
  omar <- par("mar")
  par(mar=c(2.1, 3.1, 3.1, 1.1))
  n_iters <- length(object$records[[1]])
  records = aggregateRecords(object, burn_in, keep = keep, keep_separate_chains = TRUE)
  records <- transposeList(records)
  names_params <- names(records)
  n_chain <- length(records[[1]])
  plotted_iters = seq(max(0, floor(n_iters * burn_in))+1,n_iters)
  
  xlim <- c(min(plotted_iters), max(plotted_iters))
  color_lines <- getColorsCat(3)
  
  for(name in names_params){
    print(name)
    if(name=="noise_beta"){
      print("toto")
    }
    refparams <- records[[name]][[1]]
    n <- dim(refparams)[2]
    # For parms with several intercepts, split representation !
    intercepts <- grep("Intercept", colnames(refparams))
    columns <- list()
    
    if(length(intercepts)>0) {
      starts <- intercepts
    } else {
      starts <- 1
    }
    columns <- split_intervals(n, starts, nmax = 16)
   
    for(i in length(columns)) {
      TraceOneParam(name, records[[name]], columns[[i]][1], columns[[i]][2], xlim = xlim, color_lines = color_lines)
    }
  }
    
  # Resetting original parameters mfrow
  par("mfrow"= ormfrow)
  par("mar" = omar)
  return(invisible(NULL))
}
