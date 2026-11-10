#' Coverage intervals of the response variable
#' @param true_y a vector of true observations of the response variable
#' @param y_samples samples of y observations
#'
#' @returns a vector of Booleans indicating if each true observation is within the predicted coverage interval
posterior_coverage <- function(true_y, y_samples) {
  predicted_quantiles <- rowQuantiles(y_samples, probs = c(.025, .975))
  inBounds(predicted_quantiles, true_y)
}

#' Predictive Density
#' @param true_y a vector of true observations of the response variable
#' @param mean_samples a matrix of samples of mean from the model
#' @param sd_samples a matrix of samples of standard deviation from the model
#' @returns A predictive density for each true observation
#' @export
#' @keywords internal
#' @examples
#' true_y <- rnorm(100)
#' mean_samples <- matrix(rep(0, 1000 * length(true_y)), length(true_y))
#' sd_samples <- matrix(rep(1, 1000 * length(true_y)), length(true_y))
#' predictiveDensity(true_y, mean_samples, sd_samples)
predictiveDensity <- function(true_y, mean_samples, sd_samples) {
  # Gaussian log-dens of the observations
  per_obs <- log(apply(dnorm(true_y, mean_samples, sd_samples), 1, mean))
  return(per_obs)
}

#' Squared prediction error
#' @param true_y a vector of true observations of the response variable
#' @param mean_samples a matrix of samples of mean from the model
#' @returns A squared error for each true observation
#' @export
#' @keywords internal
#' @examples
#' true_y <- rnorm(100)
#' mean_samples <- matrix(rep(0, 1000 * length(true_y)), length(true_y))
#' squaredError(true_y, mean_samples)
squaredError <- function(true_y, mean_samples) {
  per_obs <- (true_y - apply(mean_samples, 1, mean))^2
  return(per_obs)
}

#' Return averaged and per-observation model scoring
#' @param true_y a vector of true observations of the response variable
#' @param mean_samples a matrix of samples of mean from the model
#' @param sd_samples a matrix of samples of standard deviation from the model
#' @export
#' @keywords internal
#' @examples
#' true_y <- rnorm(100)
#' mean_samples <- matrix(rep(0, 1000 * length(true_y)), length(true_y))
#' sd_samples <- matrix(rep(1, 1000 * length(true_y)), length(true_y))
#' allScores(true_y, mean_samples, sd_samples)
allScores <- function(true_y, mean_samples, sd_samples) {
  y_samples <- mean_samples + sd_samples * rnorm(length(sd_samples))
  per_obs <- matrix(NA_real_, length(true_y), 4)
  colnames(per_obs) <- c("95% coverage", "CRPS", "predictive density", "squared error")
  
  ok <- !is.na(true_y)
  if (any(ok)) {
    per_obs[ok, "95% coverage"] <- posterior_coverage(true_y = true_y[ok], y_samples = y_samples[ok,])
    per_obs[ok, "CRPS"] <- scoringRules::crps_sample(true_y[ok], y_samples[ok,])
    per_obs[ok, "predictive density"] <- predictiveDensity(true_y[ok], mean_samples[ok,], sd_samples[ok,])
    per_obs[ok, "squared error"] <- squaredError(true_y[ok], mean_samples[ok,])
  }
  
  res <- list()
  res$scores <- signif(apply(per_obs, 2, mean, na.rm = TRUE), 2)
  res$scores_per_obs <- per_obs
  return(res)
}

#' @title Find scattered spatial locations for train-test partition
#' @description A max-min heuristic is used to scatter the test locations across space,
#'  so that test locations are independent-ish from each other, and surrounded
#'  by training locations.
#' @param locs a matrix of spatial locations.
#' @param prop_test the proportion of test locations. Some locations can
#' be observed several times, so the proportion of test indices can be
#' greater than prop_test. If NULL, set to 0.1.
#' @returns A list of train and test indices
#' @export
#' @keywords internal
#' 
#' @examples
#' # duplicated observed locations
#' observed_locs <- cbind(runif(2000), runif(2000))
#' observed_locs <- rbind(observed_locs, observed_locs[1:1000, ])
#' idx <- getIndexesClose(observed_locs, prop_test = 0.1)
getIndexesClose <- function(locs, prop_test = 0.1) {
  if (!inBounds(c(0, 1), prop_test, includeBounds = FALSE)) {
    stop("prop_clust and prop_test must be >0 and <1")
  }
  ulocs <- unique(locs)
  ulocs <- ulocs[order(runif(nrow(ulocs))),]
  ulocs[seq(min(10000, nrow(ulocs))),] <- ulocs[GpGp::order_maxmin(ulocs[seq(min(10000, nrow(ulocs))),]),]
  n_test <- floor(nrow(ulocs) * prop_test)
  test <- seq(n_test)
  locs_match <- match(
    split(locs, row(locs)),
    split(ulocs, row(ulocs))
  )
  idx_test <- which(locs_match %in% test)
  idx_train <- which(!(locs_match %in% test))
  return(list(
    "train" = idx_train,
    "test" = idx_test, 
    "n_locs_test" = n_test,  
    "n_obs_test"  = length(idx_test),  
    "n_locs_train"= nrow(ulocs) - n_test,  
    "n_obs_train" = length(idx_train)   
  ))
}

#' Find scattered spatial locations for train-test partition and remove
#' train locations who are too close to test locations in order to assess long
#' range prediction.
#' @description
#'  A K-means algorithm is run in order to parition the spatial locations into
#'  clusters. A max-min heuristic is used on the clusters to scatter the
#'  test clusters across space, so that test clusters are independent-ish from
#'  each other, and surrounded by training locations. The observation that is
#'  closest to the center of each cluster is kept for validation and the rest of
#'  the cluster is discarded, in order to forbid immediate interpolation of
#'  the test data.
#' @param locs a matrix of spatial locations.
#' @param n_clust the number of clusters. If NULL, will be set to 1\% of the data.
#' @param prop_test the proportion of clusters used as test locations.
#' Default to 0.01.
#' @returns A list of train and test indices
#' @export
#' @keywords internal
#'
#' @examples
#' observed_locs <- round(cbind(runif(10000), runif(10000)), 3)
#' idx <- getIndexesFar(observed_locs)
getIndexesFar <- function(locs, n_clust = NULL, prop_test = 0.1) {
  if (is.null(n_clust)) n_clust <- floor(nrow(locs) / 100)
  if (!inBounds(c(0, nrow(locs)), n_clust, includeBounds = FALSE) ||
    !inBounds(c(0, 1), prop_test, includeBounds = FALSE)) {
    stop("prop_clust and prop_test must be >0 and <1")
  }
  ulocs <- unique(locs)
  K <- kmeans(ulocs, centers = n_clust)
  n_test <- floor(n_clust * prop_test)
  if(n_test < 3 ) {
    warning("Not enough observations to produce far test sample")
    return(NULL)
  }
  if(n_test<=0) stop("n_clust * prop_test should be > 0")
  test_clust <- GpGp::order_maxmin(K$centers)[seq(n_test)]
  test <- FNN::knnx.index(data = ulocs, k = 1, query = K$centers[test_clust, ])
  discarded <- setdiff(which(K$cluster %in% test_clust), test)
  train <- which(!K$cluster %in% test_clust)
  locs_match <- match(
    split(locs, row(locs)),
    split(ulocs, row(ulocs))
  )
  idx_test <- which(locs_match %in% test)
  idx_train <- which(!(locs_match %in% discarded) & !(locs_match %in% test))
  idx_discarded <- which(locs_match %in% discarded)
  return(list(
    "train" = idx_train,
    "test" = idx_test, 
    "discarded" = idx_discarded, 
    "n_obs_test" = length(idx_test), 
    "n_locs_test" = nrow(test), 
    "n_obs_discarded" = length(which(locs_match %in% discarded)),
    "n_locs_discarded" =  length(discarded)
  ))
}

#' Plot the summary of a data split
#'
#' @param list_locs a list of length n, each entries having a 'locs' component
#' that should be a matrix (or data.frame or array) with 2 columns containing
#' spatial locations
#' @return a list
#' @export
#' 
#' @examples
#' locs <- cbind(runif(2000), runif(2000))
#' locslist <- list("train" = list("locs"=locs[1:1900,]), 
#'                  "test" = list(locs=locs[1901:2000,]))
#' splitSummary(locslist)
splitSummary <- function(list_locs) {
  info <- attr(list_locs, "info")
  print(paste("Train test split of a spatial data set with", info["total", "observations"], 
              "observations on", info["total", "spatial locations"], "distinct spatial locations"))
  print(info)
}

#' Plot the locations of a data split
#'
#' @param list_locs a list of length n, each entries having a 'locs' component
#' that should be a matrix (or data.frame or array) with 2 columns containing
#' spatial locations
#' @param pch a vector of pch for points shape of length 1 (will be repeated) or n. 
#' Passed to `plot`
#' @param cex a numeric vector for points size of length 1 (will be repeated) or n. 
#' Passed to `plot`
#' @param col a vector for points color of length n. Passed to `plot`
#' @returns a plot
#' @export
#' 
#' @examples
#' locs <- cbind(runif(2000), runif(2000))
#' locslist <- list("train" = list("locs"=locs[1:1900,]), 
#'                  "test" = list(locs=locs[1901:2000,]))
#' plotSplit(locslist)
plotSplit <- function(list_locs,
                      col = c("gray", "red", "orange", "blue"),
                      cex = c(1, .5, .5, 1),
                      pch = list(".", 15, 15, ".")) {
  if(!is.list(list_locs)) stop("list_locs should be a list")
  if(!all(sapply(list_locs, function(x) "locs" %in% names(x))))
    stop("Except info, all list_locs entries should have a 'locs' element")
  list_locs <- lapply(list_locs, function(x) x[["locs"]])
  nblist <- length(list_locs)
  nameslocs <- names(list_locs)
  
  replaceparam <- function(arg, user, is_missing) {
    def <- eval(formals(plotSplit)[[arg]])
    out <- rep_len(def, nblist)
    names(out) <- nameslocs
      if (is.null(names(user))) { # unnammed -> recycled.
        out <- rep_len(user, nblist)
        names(out) <- nameslocs
      } else { # named -> replaced
        for (k in intersect(names(user), nameslocs)) 
          out[[k]] <- user[[k]]
      }
    out
  }
  col <- replaceparam("col", col, missing(col))
  cex <- replaceparam("cex", cex, missing(cex))
  pch <- replaceparam("pch", pch, missing(pch))
  
  lims <- sapply(list_locs, function(m) apply(m, 2, range))
  xlim <- range(lims[1:2, ])
  ylim <- range(lims[3:4, ])
  
  xlab <- "locs[,1]"
  ylab <- "locs[,2]"
  if (!is.null(colnames(list_locs[[1]]))) {
    xlab <- colnames(list_locs[[1]])[1]
    ylab <- colnames(list_locs[[1]])[2]
  }
  
  oldpar <- par(c("mar", "xpd"))
  on.exit(par(mar=oldpar$mar, xpd=oldpar$xpd))
  if(def.par$mar[3]<=3) {
    par(mar = def.par$mar + c(0, 0, 3-def.par$mar[3], 0), xpd = NA)
  }
  
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
  
  for (k in nameslocs)
    points(list_locs[[k]], col = col[[k]], pch = pch[[k]], cex = cex[[k]])
  
  if (!is.null(names(list_locs))) {
    legend_labels <- nameslocs
    u <- par("usr")
    legend(
      x = u[1],
      y = u[4],
      legend = legend_labels,
      y.intersp = 1,
      yjust=0,
      xpd = NA,
      ncol = 2,
      fill=col,
      horiz = FALSE,
      bty = "n"
    )
  }
}
  

#' Split data between train and test
#' @description Split a GeoNonStat object between test and train samples
#'
#' @param locs a matrix of spatial locations.
#' @param data a list containing vectors or data.frames to split into
#' train/test
#' @param n_clust the number of clusters. If NULL, will be set to 1\% of the data.
#' @param prop_test_close the proportion of observations used for "close" 
#' validation, among "not far" data. Default to 0.1
#' @param prop_test_far the proportion of observations used for "far" 
#' validation, as a proportion of spatial clusters. Default to 0.1
#' @param n_clust_far number of spatial clusters for "far" validation. 
#' Default to NULL, set to 1% of the data
#' @param round_locs numeric, default to 0. rounding factor for locations coordinates. 
#' @export
#'
#' @examples
#' n <- 10000
#' observed_locs <- round(cbind(runif(n), runif(n)), 3)
#' X <- as.data.frame(matrix(rnorm(n * 4), ncol = 4))
#' range_X <- as.data.frame(matrix(rnorm(n * 4), ncol = 4))
#' noise_X <- as.data.frame(matrix(rnorm(n * 4), ncol = 4))
#' y <- rnorm(n)
#' datalist <- list("X" = X, "range_X" = range_X, "noise_X" = noise_X, "y" = y)
#' traintest <- splitData(observed_locs, datalist, prop_test_close = 0.1)
#' traintest <- splitData(observed_locs, datalist, prop_test_close = 0.1, prop_test_far=0.1)
splitData <- function(locs, 
                      datalist,
                      prop_test_close = 0.1,
                      prop_test_far = 0.1, 
                      n_clust_far = NULL,
                      round_locs = 0){
  if (round_locs > 0) locs <- round_locs * round(locs / round_locs)
  
  # Split "far" : désactivé si n_clust_far = 0 ou prop_test_far = 0
  idx_far <- NULL
  if ((n_clust_far > 0 || is.null(n_clust_far)) && prop_test_far > 0) {
    idx_far <- getIndexesFar(locs, n_clust = n_clust_far, prop_test = prop_test_far)
  }
  has_far <- !is.null(idx_far)
  
  # Create indexes for train, test far and test close
  idx_avail <- if (has_far) idx_far$train else seq_len(nrow(locs))
  idx_close <- getIndexesClose(locs = locs[idx_avail, ], prop_test = prop_test_close)
  
  final_split <- list(
    "train"      = idx_avail[idx_close$train],
    "test_close" = idx_avail[idx_close$test]
  )
  info <- rbind(
    "train"      = c(idx_close$n_locs_train, idx_close$n_obs_train),
    "test_close" = c(idx_close$n_locs_test,  idx_close$n_obs_test)
  )
  
  if (has_far) {
    final_split$test_far  <- idx_far$test
    final_split$discarded <- idx_far$discarded
    info <- rbind(
      "train"      = info["train", , drop = FALSE],
      "test_close" = info["test_close", , drop = FALSE],
      "test_far"   = c(idx_far$n_locs_test,      idx_far$n_obs_test),
      "discarded"  = c(idx_far$n_locs_discarded, idx_far$n_obs_discarded)
    )
  }
  info <- rbind(info, "total" = colSums(info))
  colnames(info) <- c("spatial locations", "observations")


  makeBlock <- function(idx) c(list(locs = locs[idx, , drop = FALSE]), subsetData(datalist, idx))
  output <- list("train" = makeBlock(final_split$train),
                 "test_close" = makeBlock(final_split$test_close))
  if (has_far) {
    output$test_far <- makeBlock(final_split$test_far)
    output$discarded <- makeBlock(final_split$discarded)
  }
  attr(output, "info") <- info

  cat("Summary of split data\n")
  splitSummary(output)
  plotSplit(output)
  return(output)
}
