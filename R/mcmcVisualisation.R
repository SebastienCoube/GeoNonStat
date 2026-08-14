getRecordVarNames <- function(records)unname(sapply(unlist(mapply(paste, names(records[[1]]), sapply(records[[1]], colnames))), trimws))

#' summary of MCMC diagnostics
#'
#' @param object a GeoNonStat object
#' @param burn_in numeric, value, between 0 and 1. Gives the proportion of
#' first records that should be removed. Default to 0.1
#' @return a list of two data.frames
#'
#' @export
#' @keywords internal
#' @examples
#' MCMC_diags <- mcmcDiags(processedGnsDemo)
#' # TODO - quid des NA ?
mcmcDiags <- function(object, burn_in = 0.3, verbose = TRUE, keep = "nofield") {
  # only two significant digits are kept
  records <- aggregateRecords(object$records, burn_in = burn_in, keep = keep, keep_separate_chains = TRUE)
  res <- matrix(0, sum(sapply(records[[1]], ncol)), 3)
  colnames(res) <- c("ess", "Point est. Gelman", "Upper C.I. Gelman")
  row.names(res) <- getRecordVarNames(records)
  i <- 1 
  for(var_name in names(records[[1]])){
    for (var_num in seq(ncol(records[[1]][[var_name]]))){
      samples <- list() 
      for(j in seq_len(length(records))){
        samples[[j]] <- coda::mcmc(records[[j]][[var_name]][,var_num,drop=F])
        res[i,1] <- res[i,1] + coda::effectiveSize(samples[[j]])
      }
      res[i,c(2,3)] <- coda::gelman.diag(coda::mcmc.list(samples))$psrf
      i <- i+1
    }
  }
  
  res <- signif(res, 3)
  res <- data.frame(res)
  if (anyNA(res)) {
    warning(
      "Some parameters haven't been updated yet by the MCMC.",
      "You should increase n_iterations and re-run MCMC before requesting diags."
    )
    return(res)
  }

  regle <- c("min", "max", "max")
  worst <- data.frame(
    criterium = colnames(res),
    worst = regle,
    value = mapply(function(col, r) {
      if (r == "min") min(col) else max(col)
    }, res, regle, USE.NAMES = FALSE),
    which_worst = mapply(function(col, r) {
      i <- if (r == "min") which.min(col) else which.max(col)
      rownames(res)[i]
    }, res, regle, USE.NAMES = FALSE)
  )

  if (verbose) {
    cat("Summary of diagnostics\n")
    print(summary(res))
    cat("----------------------\n")
    print(worst)
  }
  return(list("diags" = res, "worst" = worst))
}

traceOneParam <- function(name, records, start, end, nrow, ncol, xlim, color_lines) {
  par(mfrow = c(nrow, ncol))

  to_be_plotted <- matrix(0, nrow(records[[1]][[1]]), length(records))
  # Iterate over columns
  for (col in seq(start, end)) {
    namecol <- colnames(records[[1]][[name]])[col]
    for(i in seq(length(records))){
      to_be_plotted[,i] <- records[[i]][[name]][,col]
    }
    graphics::matplot(
      to_be_plotted,
      type = "l",
      lty = 1,
      lwd = 2,
      main = paste0(name, "\n", namecol),
      col = color_lines
    )
  }
}

splitIntervals <- function(total_size, starts, nmax = 16) {
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

  dims <- max(sapply(intervals, function(x) x[2] - x[1] + 1))
  ncol <- min(dims, ifelse(dims > 9, 4, 3))
  nrow <- ceiling(dims / ncol)
  intervals <- lapply(intervals, function(x) c(x, nrow, ncol))
  intervals
}

#' Trace the records of a processed GeoNonStat object
#'
#' @param object a processed GeoNonStat object
#' @param burn_in numeric, value, between 0 and 1. Gives the proportion of
#' first records that should be removed. Default to 0.1
#' @param keep names of paramaters to plot
#'
#' @returns NULL
#' @export
#'
#' @examples
#' tracePlots(processedGnsDemo)
tracePlots <- function(object, burn_in = .3, keep = "nofield") {
  # Checks
  if (length(object$records[[1]]) <= 1) {
    warning(
      length(object$records[[1]]),
      " record found, has the MCMC been executed on this object ?",
      " No plot produced."
    )
    return(NULL)
  }
  # Prepare data
  records <- aggregateRecords(object$records,
    burn_in,
    keep = keep,
    keep_separate_chains = TRUE
  )

  # How to plot
  ormfrow <- par("mfrow")
  omar <- par("mar")
  par(mar = c(2.1, 3.1, 3.1, 1.1))
  n_iters <- length(object$records[[1]])
  names_params <- names(records[[1]])
  plotted_iters <- seq(max(0, floor(n_iters * burn_in)) + 1, n_iters)

  xlim <- c(min(plotted_iters), max(plotted_iters))
  color_lines <- getColorsCat(3)
  nplots <- 0
  for (name in names_params) {
    refparams <- records[[1]][[name]]
    n <- dim(refparams)[2]

    # For params with several intercepts, split representation !
    columns <- list()
    starts <- 1
    intercepts <- grep("Intercept", colnames(refparams))
    if (length(intercepts) > 0) {
      starts <- intercepts
    }

    columns <- splitIntervals(n, starts, nmax = 9)
    cat("Plot", name)
    if (length(columns) > 1) cat(" (on", length(columns), "plots)")
    cat("\n")
    for (i in 1:length(columns)) {
      traceOneParam(
        name,
        records, 
        start = columns[[i]][1],
        end =   columns[[i]][2],
        nrow =  columns[[i]][3],
        ncol =  columns[[i]][4],
        xlim = xlim,
        color_lines = color_lines
      )
      nplots <- nplots + 1
    }
  }
  cat("...", nplots, "plots produced. ")
  # Resetting original parameters mfrow
  par("mfrow" = ormfrow)
  par("mar" = omar)
  return(invisible(NULL))
}
