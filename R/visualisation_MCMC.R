

#  #' MCMC diagnostics 
#  #'
#  #' @param object 
#  #' @param burn_in 
#  #' @param min_ESS 
#  #' @param max_Gelman_upper 
#  #' @description
#  #' The MCMC chain(s) must explore properly the posterior distribution in order to be able to compute sensible estimates. 
#  #' The catch is that the samples are correlated with each other, because we have a Markov chain. How to know when the exploration is done properly ? We use two classical methods from the literature. 
#  #' 
#  #' The first is the so-called Effective Sample Size (ESS). 
#  #' The ESS adjusts the number of iterations by the auto-correlation function: the higher the auto-correlation, the lower the ESS. 
#  #' Think of this number as the number of samples from the posterior distribution you would like to compute estimates and histograms from. 
#  #' The higher the ESS, the better the estimates will be. 
#  #' https://www.rdocumentation.org/packages/sns/versions/1.2.2/topics/ess
#  #' https://mc-stan.org/docs/2_18/reference-manual/effective-sample-size-section.html
#  #' 
#  #' The second is the Gelman-Rubin diagnostic. It is some kind of ANOVA between two parallel chains. 
#  #' If all the chains are exploring the same posterior distribution, their samples are difficult to tell apart.
#  #' On the other hand, if each chain is stuck in its own place and has not reached the posterior distribution, then their samples are very easy to tell apart. 
#  #' Note that Gelman-Rubin diagnostics can be "tricked" by chains whose starting points are not well separated, which is the case in our implementation, and give values close to 1 even if the chains have not converged yet. 
#  #' As a consequence, wait for a couple hundred iterations before paying attention to the Gelman-Rubin diagnostics. 
#  #' The closer to 1, the better. 
#  #' 
printMcmcDiags = function(object, burn_in = .1, min_ESS = 50, max_Gelman_upper = 1.05){
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
  res <- data.frame(res)
  colnames(res)  <- c("ess", "Point est. Gelman", "Upper C.I. Gelman")
  
  regle <- c("min", "max", "max")
  worst <- data.frame(
    criterium = colnames(res), 
    worst=regle,
    value = mapply(function(col, r) {
      if (r == "min") min(col) else max(col)
    }, res, regle, USE.NAMES=FALSE),
    which_worst = mapply(function(col, r) {
      i <- if (r == "min") which.min(col) else which.max(col)
      rownames(res)[i]
    }, res, regle, USE.NAMES=FALSE)
  )
  
  cat("Summary of diagnostics\n")
  print(summary(res))
  cat("----------------------\n")
  print(worst)
  return(list("diags"=res, "worst" = worst))
}

traceOneParam = function(name, param, start, end, nrow, ncol, xlim, color_lines){
  par(mfrow = c(nrow, ncol))
  
  # Iterate over columns
  for(col in seq(start, end)) { 
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

splitIntervals <- function(total_size, starts, nmax =  16) {
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
  
  dims <- max(sapply(intervals, function(x) x[2]-x[1] +1))
  ncol <- min(dims, ifelse(dims>9, 4, 3))
  nrow <- ceiling(dims / ncol)
  intervals <- lapply(intervals, function(x) c(x, nrow, ncol))
  intervals
}

#' trace the records of a processed GeoNonStat object
#'
#' @param object a processed GeoNonStat object
#' @param burn_in percentage of first records to remove
#' @param keep names of pramaters to plot
#'
#' @returns NULL
#' @export
#'
#' @examples
#' tracePlots(processedGnsDemo)
tracePlots = function(object, burn_in = .1, keep = "all"){
  # Checks
  if(length(object$records[[1]])<=1) {
    warning( length(object$records[[1]]), 
             " record found, has the MCMC been executed on this object ?",
             " No plot produced.")
    return(NULL)
  }
  # Prepare data
  records = aggregateRecords(object, 
                             burn_in, 
                             keep = keep, 
                             keep_separate_chains = TRUE)
  records <- transposeList(records)
  
  # How to plot
  ormfrow <- par("mfrow")
  omar <- par("mar")
  par(mar=c(2.1, 3.1, 3.1, 1.1))
  n_iters <- length(object$records[[1]])
  names_params <- names(records)
  plotted_iters = seq(max(0, floor(n_iters * burn_in))+1,n_iters)
  
  xlim <- c(min(plotted_iters), max(plotted_iters))
  color_lines <- getColorsCat(3)
  
  for(name in names_params){
    refparams <- records[[name]][[1]]
    n <- dim(refparams)[2]
    
    # For params with several intercepts, split representation !
    columns <- list()
    starts <- 1
    intercepts <- grep("Intercept", colnames(refparams))
    if(length(intercepts)>0) 
      starts <- intercepts
    
    columns <- splitIntervals(n, starts, nmax = 16)
    cat("Plot", name)
    if(length(columns)>1) cat(" (on", length(columns), "plots)")
    cat("\n")
    for(i in 1:length(columns)) {
      traceOneParam(name, 
                    records[[name]], 
                    columns[[i]][1], 
                    columns[[i]][2], 
                    columns[[i]][3], 
                    columns[[i]][4], 
                    xlim = xlim, 
                    color_lines = color_lines)
    }
  }
    
  # Resetting original parameters mfrow
  par("mfrow"= ormfrow)
  par("mar" = omar)
  return(invisible(NULL))
}
