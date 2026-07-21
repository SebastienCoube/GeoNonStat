#' Coverage intervals of the response variable
#' @param true_y a vector of true observations of the response variable
#' @param y_samples a matrix of pseudo-observations sampled following the model
#'
#' @returns a vector of Booleans indicating if each true observation is within the predicted coverage interval
#' @export
#' @keywords internal
#'
#' @examples
#' true_y <- rnorm(100)
#' y_samples <- matrix(rnorm(1000 * length(true_y)), length(true_y))
#' coverage(true_y, y_samples)
coverage <- function(true_y, y_samples) {
  quantiles <- apply(y_samples, 1, function(x) quantile(x, c(.025, .975)))
  per_obs <- mapply(
    inBounds,
    x = true_y,
    bounds = split(quantiles, col(quantiles))
  )
  return(per_obs)
}

#' \href{https://en.wikipedia.org/wiki/Scoring_rule#Examples_of_proper_scoring_rules}{Continuous ranked probability score} (CRPS)
#' @param true_y a vector of true observations of the response variable
#' @param y_samples a matrix of pseudo-observations sampled following the model
#' @returns A CRPS for each true observation
#' @export
#' @keywords internal
#'
#' @examples
#' true_y <- rnorm(100)
#' y_samples <- matrix(rnorm(1000 * length(true_y)), length(true_y))
#' CRPS(true_y, y_samples)
CRPS <- function(true_y, y_samples) {
  ecdfs <- apply(y_samples, 1, stats::ecdf)
  per_obs <- mapply(
    function(f, x) f(x),
    x = true_y,
    f = ecdfs
  )
  return(per_obs)
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
#' score(true_y, mean_samples, sd_samples)
score <- function(true_y, mean_samples, sd_samples) {
  y_samples <- mean_samples + sd_samples * rnorm(length(sd_samples))
  per_obs <- matrix(0, length(true_y), 4)
  colnames(per_obs) <- c("95% coverage", "CRPS", "predictive density", "squared error")
  per_obs[, "95% coverage"] <- coverage(true_y, y_samples)
  per_obs[, "CRPS"] <- CRPS(true_y, y_samples)
  per_obs[, "predictive density"] <- predictiveDensity(true_y, mean_samples, sd_samples)
  per_obs[, "squared error"] <- squaredError(true_y, mean_samples)
  res <- list()
  res$criteria <- signif(apply(per_obs, 2, mean), 2)
  names(res$criteria) <- c("95% coverage", "CRPS", "ELPPD", "MSE")
  res$criteria_per_obs <- per_obs
  return(res)
}


#' Model scoring for train and test data
#' @description
#'  Several scores are implemented: the \href{https://en.wikipedia.org/wiki/Coverage_probability}{Coverage by 95\% intervals},
#'  the \href{https://en.wikipedia.org/wiki/Scoring_rule#Examples_of_proper_scoring_rules}{Continuous ranked probability score} (CRPS),
#'  the \href{https://arxiv.org/abs/1507.04544}{Empirical log Predictive Pointwise Density} (ELPPD),
#'  and the \href{https://en.wikipedia.org/wiki/Mean_squared_error}{Mean Squared Error} (MSE).
#'  Train scores are computed from a GeoNonStat object, while test scores are computed from a prediction and a vector of test observations.
#' @param object, an object of class GeoNonStat in order to evaluate the scores on the training data set, or NULL
#' @param prediction, a prediction from a GeoNonStat object in order to evaluate the scores on the test data set, or NULL
#' @param new_y a vector on new observations from the interest variable in order to evaluate the scores on the test data set, or NULL
#' @export
scores <- function(object = NULL, prediction = NULL, new_y = NULL) {
  train <- !is.null(object) & is.null(prediction) & is.null(new_y) & (class(object) == "GeoNonStat")
  test <- is.null(object) & !is.null(prediction) & !is.null(new_y)
  if (!(train | test)) {
    stop(
      "Either: object must be from class GeoNonStat, prediction be NULL,",
      "and new_y be NULL; or geo_non_stat must be NULL, ",
      "prediction be a prediction from a GeoNonStat model, ",
      "and y be test observations"
    )
  }
  # TODO check new_y and prediction same length
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
#' @examples
#' # duplicated observed locations
#' observed_locs <- cbind(runif(2000), runif(2000))
#' observed_locs <- rbind(observed_locs, observed_locs[1:1000, ])
#' idx <- getIndexesClose(observed_locs, prop_test = 0.1)
#' plotSplit(observed_locs, idx)
getIndexesClose <- function(locs, prop_test = 0.1) {
  if (!inBounds(c(0, 1), prop_test, includeBounds = FALSE)) {
    stop("prop_clust and prop_test must be >0 and <1")
  }
  ulocs <- unique(locs)
  ulocs <- ulocs[order(runif(nrow(ulocs))),]
  ulocs[seq(10000),] <- ulocs[GpGp::order_maxmin(ulocs[seq(10000),]),]
  n_test <- floor(nrow(ulocs) * prop_test)
  test <- seq(n_test)
  locs_match <- match(
    split(locs, row(locs)),
    split(ulocs, row(ulocs))
  )
  idx_test <- which(locs_match %in% test)
  if (length(idx_test) > n_test) {
    message(
      "Due to duplicates, ", round(100 * length(idx_test) / nrow(ulocs), 2),
      "% indices returned for 'close' test sample"
    )
  }
  idx_train <- which(!(locs_match %in% test))
  return(list(
    "train" = idx_train,
    "test" = idx_test
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
#'
#' @examples
#' observed_locs <- round(cbind(runif(10000), runif(10000)), 3)
#' idx <- getIndexesFar(observed_locs)
#' plotSplit(observed_locs, idx)
getIndexesFar <- function(locs, n_clust = NULL, prop_test = 0.1) {
  if (is.null(n_clust)) n_clust <- floor(nrow(locs) / 100)
  if (!inBounds(c(0, nrow(locs)), n_clust, includeBounds = FALSE) ||
    !inBounds(c(0, 1), prop_test, includeBounds = FALSE)) {
    stop("prop_clust and prop_test must be >0 and <1")
  }
  ulocs <- unique(locs)
  K <- kmeans(ulocs, centers = n_clust)
  n_test <- floor(n_clust * prop_test)
  test_clust <- GpGp::order_maxmin(K$centers)[seq(n_test)]
  test <- FNN::knnx.index(data = ulocs, k = 1, query = K$centers[test_clust, ])
  discarded <- setdiff(which(K$cluster %in% test_clust), test)
  train <- which(!K$cluster %in% test_clust)
  locs_match <- match(
    split(locs, row(locs)),
    split(ulocs, row(ulocs))
  )
  idx_test <- which(locs_match %in% test)
  if (length(idx_test) > n_test) {
    message(
      "Due to duplicates, ",
      round(100 * length(idx_test) / n_clust, 2),
      "% indices returned for 'far' test sample"
    )
  }
  idx_train <- which(!(locs_match %in% discarded) & !(locs_match %in% test))
  return(list(
    "train" = idx_train,
    "test" = idx_test
  ))
}


#' Plot the locations of a data split
#'
#' @param locs a matrix (or data.frame or array) with 2 columns containing
#' spatial locations
#' @param indexes a list of length n, containing usually train and test indices,
#' among row indexes of `locs`
#' @param pch a vector of pch for points shape of length 1 or 2. Passed to `plot`
#' @param cex a numeric vector for points size of length 1 or 2. Passed to `plot`
#' @param col a vector for points color of length 1 or 2. Passed to `plot`
#' @returns a plot
#' @export
#'
#' @examples
#' locs <- cbind(runif(220), runif(220))
#' idx <- list("train" = 1:200, "test" = 201:220)
#' plotSplit(locs, idx)
plotSplit <- function(locs, indexes,
                      pch = c(16, 16),
                      cex = c(0.3, 1),
                      col = c("blue", "red")) {
  locs_name <- deparse(substitute(locs))
  n <- length(indexes)
  if (length(pch) != n) pch <- rep(pch[1], n)
  if (length(cex) != n) cex <- rep(cex[1], n)
  if (length(col) != n) col <- grDevices::rainbow(n)
  # plot locations
  plot(locs[indexes[[1]], ],
    pch = pch[1], cex = cex[1], col = col[1],
    xlab = paste0(locs_name, "[,1]"),
    ylab = paste0(locs_name, "[,2]")
  )
  for (i in 2:length(indexes)) {
    points(locs[indexes[[i]], ],
      pch = pch[i], cex = cex[i], col = col[i],
      xlab = "", ylab = ""
    )
  }
  if (!is.null(names(indexes))) {
    legend(
      x = "topleft",
      legend = names(indexes),
      pt.cex = 0,
      inset = c(0, -0.2),
      xpd = TRUE,
      text.col = col,
      horiz = TRUE,
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
#' @param prop_test the proportion of clusters used as test locations.
#' Default to 0.01.
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
#' traintest <- createSplitData(observed_locs, datalist, prop_test = c(0.1, 0.01))
createSplitData <- function(locs, data,
                            n_clust = NULL,
                            prop_test = 0.1, 
                            round_locs = 0){
  if (length(prop_test) == 1) prop_test <- rep(prop_test, 2)
  if (round_locs > 0) locs <- round_locs * round(locs / round_locs)
  # Create indexes for train, test far and test close
  idx_far <- getIndexesFar(locs, n_clust = n_clust, prop_test = prop_test[1])
  idx_close <- getIndexesClose(locs[idx_far$train, ], prop_test = prop_test[2])
  final_split <- list(
    "train" = idx_far$train[idx_close$train],
    "test_far" = idx_far$test,
    "test_close" = idx_far$train[idx_close$test]
  )
  plotSplit(locs, final_split,
    col = c("black", "red", "orange"),
    cex = c(0.1, 0.5, 0.5)
  )
  return(list(
    "train" = c(
      list(locs = locs[final_split$train, ]),
      subsetData(data, final_split$train)
    ),
    "test_far" = c(
      list(locs = locs[final_split$test_far, ]),
      subsetData(data, final_split$test_far)
    ),
    "test_close" = c(
      list(locs = locs[final_split$test_close, ]),
      subsetData(data, final_split$test_close)
    )
  ))
}
