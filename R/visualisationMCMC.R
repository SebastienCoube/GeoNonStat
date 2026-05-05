#' @export
printMcmcDiags <- function(object, burn_in = .1, min_ESS = 50, max_Gelman_upper = 1.05) {
  # only two significant digits are kept
  records <- aggregateRecords(object, burn_in = burn_in, keep = "all", keep_separate_chains = TRUE)
  ess <- lapply(records, function(x) lapply(x, function(x) coda::effectiveSize(x)))
  ess <- lapply(ess, unlist)
  ess <- Reduce("+", ess)
  ess <- signif(ess, 2)

  records <- lapply(records, function(x) do.call("cbind", x))
  records <- lapply(records, coda::as.mcmc)
  records <- coda::as.mcmc.list(records)
  gelman_diags <- t(coda::gelman.diag(records, multivariate = F)[[1]])
  gelman_diags <- signif(gelman_diags, 2)
  res <- t(rbind(ess, gelman_diags))
  res <- data.frame(res)
  colnames(res) <- c("ess", "Point est. Gelman", "Upper C.I. Gelman")

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

  cat("Summary of diagnostics\n")
  print(summary(res))
  cat("----------------------\n")
  print(worst)
  return(list("diags" = res, "worst" = worst))
}

traceOneParam <- function(name, param, start, end, nrow, ncol, xlim, color_lines) {
  par(mfrow = c(nrow, ncol))

  # Iterate over columns
  for (col in seq(start, end)) {
    namecol <- colnames(param[[1]])[col]
    to_be_plotted <- do.call("cbind", lapply(param, function(x) x[, col]))
    graphics::matplot(to_be_plotted,
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
tracePlots <- function(object, burn_in = .1, keep = "all") {
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
  records <- aggregateRecords(object,
    burn_in,
    keep = keep,
    keep_separate_chains = TRUE
  )
  records <- transposeList(records)

  # How to plot
  ormfrow <- par("mfrow")
  omar <- par("mar")
  par(mar = c(2.1, 3.1, 3.1, 1.1))
  n_iters <- length(object$records[[1]])
  names_params <- names(records)
  plotted_iters <- seq(max(0, floor(n_iters * burn_in)) + 1, n_iters)

  xlim <- c(min(plotted_iters), max(plotted_iters))
  color_lines <- getColorsCat(3)

  for (name in names_params) {
    refparams <- records[[name]][[1]]
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
      traceOneParam(name,
        records[[name]],
        columns[[i]][1],
        columns[[i]][2],
        columns[[i]][3],
        columns[[i]][4],
        xlim = xlim,
        color_lines = color_lines
      )
    }
  }

  # Resetting original parameters mfrow
  par("mfrow" = ormfrow)
  par("mar" = omar)
  return(invisible(NULL))
}
