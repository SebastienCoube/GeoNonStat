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
  per_obs <- matrix(0, length(true_y), 4)
  colnames(per_obs) <- c("95% coverage", "CRPS", "predictive density", "squared error")
  per_obs[, "95% coverage"] <- coverage(true_y, y_samples)
  per_obs[, "CRPS"] <- scoringRules::crps_sample(true_y, y_samples)
  per_obs[, "predictive density"] <- predictiveDensity(true_y, mean_samples, sd_samples)
  per_obs[, "squared error"] <- squaredError(true_y, mean_samples)
  res <- list()
  res$scores <- signif(apply(per_obs, 2, mean), 2)
  res$scores_per_obs <- per_obs
  return(res)
}




#' Model scoring for train data
#'
#' @param object a GeoNonStat object
#' @param burn_in numeric, value, between 0 and 1. Gives the proportion of
#' first records that should be removed. Default to 0.1
#'
#' @description
#' Several scores are implemented: the \href{https://en.wikipedia.org/wiki/Coverage_probability}{Coverage by 95\% intervals},
#' the \href{https://en.wikipedia.org/wiki/Scoring_rule#Examples_of_proper_scoring_rules}{Continuous ranked probability score} (CRPS),
#' the \href{https://arxiv.org/abs/1507.04544}{Empirical log Predictive Pointwise Density} (ELPPD),
#' and the \href{https://en.wikipedia.org/wiki/Mean_squared_error}{Mean Squared Error} (MSE).
#' Train scores are computed from a GeoNonStat object, while test scores are computed from a prediction and a vector of test observations.
#' @export
#' 
#' @examples
#' trainScores(processedGnsDemo)
trainScores <- function(object, burn_in){
  estimates <- estimate(object)
  allScores(
    true_y = object$observed_field, 
    mean_samples = t(estimates$samples$denoised), 
    sd_samples = t(exp(.5*estimates$samples$noise_log_var)))
}

#' Model scoring for test data
#'
#' @param object a GeoNonStat object
#' @param burn_in numeric, value, between 0 and 1. Gives the proportion of
#' first records that should be removed. Default to 0.1
#' @param test_datalist TODO 
#' @param num_threads numeric TODO
#'
#' @description
#' Several scores are implemented: the \href{https://en.wikipedia.org/wiki/Coverage_probability}{Coverage by 95\% intervals},
#' the \href{https://en.wikipedia.org/wiki/Scoring_rule#Examples_of_proper_scoring_rules}{Continuous ranked probability score} (CRPS),
#' the \href{https://arxiv.org/abs/1507.04544}{Empirical log Predictive Pointwise Density} (ELPPD),
#' and the \href{https://en.wikipedia.org/wiki/Mean_squared_error}{Mean Squared Error} (MSE).
#' Train scores are computed from a GeoNonStat object, while test scores are computed from a prediction and a vector of test observations.
#' @export
#'  
#' @examples
#' # Recreating new data with less observations but same structure as example data
#' locs <- cbind(rnorm(200), rnorm(200))
#' X <- data.frame(matrix(rnorm(200*4), ncol=4))
#' colnames(X) <- paste0("V", 1:4)
#' Yobs <- rnorm(200)
#' testdata <- list("locs" = locs, 
#'                  "X" = X, 
#'                  noise_X=X,
#'                  observed_field=Yobs)  
#' testScores(processedGnsDemo, test_datalist=testdata)
testScores <- function(object, 
                       burn_in = 0.1, 
                       test_datalist,
                       num_threads=5){
  
  # getting elements
  if(!("observed_field" %in% names(test_datalist))) {
    stop("test_datalist should have a 'observed_field' element.")
  }
  test_observed_field <- test_datalist$observed_field
  if(!("locs" %in% names(test_datalist))) {
    stop("test_datalist should have a 'locs' element.")
  }
  test_locs <- test_datalist$locs
  
  test_X <- test_noise_X <- test_range_X <- NULL
  if("X" %in% names(test_datalist))  test_X <- test_datalist$X
  if("range_X" %in% names(test_datalist))  test_range_X <- test_datalist$range_X
  if("noise_X" %in% names(test_datalist))  test_noise_X <- test_datalist$noise_X
  
  prediction <- predict.GeoNonStat(
    object = object, new_data_list=test_datalist, burn_in = burn_in, 
    num_threads = num_threads)
  allScores(
    true_y = test_observed_field, 
    mean_samples = t(prediction$samples$denoised), 
    sd_samples = t(exp(.5*prediction$samples$noise_log_var)))
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
  print(paste("Train test split of a spatial data set with", list_locs$info["total", "observations"], 
              "observations on", list_locs$info["total", "spatial locations"], "distinct spatial locations"))
  print(list_locs$info)
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
                      pch = 16,
                      cex_train = 1,
                      cex_discarded = 1,
                      cex_close = .5,
                      cex_far = .5,
                      col_train = "gray", 
                      col_discarded = "blue",
                      col_close = "orange",
                      col_far = "red", 
                      legend_position = "topleft") {
  if(!is.list(list_locs)) stop("list_locs should be a list")
  trainlocs = unique(list_locs$train$locs)
  discardedlocs = unique(list_locs$discarded$locs)
  closelocs = unique(list_locs$test_close$locs)
  farlocs = unique(list_locs$test_far$locs)
  plot(rbind(
    trainlocs,
    discardedlocs,
    closelocs,
    farlocs
  ), type = "n", xlab ="", ylab ="")
  points(trainlocs, col = col_train, pch = ".", cex = cex_train)
  points(discardedlocs, col = col_discarded, pch = ".", cex = cex_discarded)
  points(farlocs, col = col_far, pch = 15, cex = cex_far)
  points(closelocs, col = col_close, pch = 15, cex = cex_close)
  legend(legend_position, legend = c("train", "discarded", "test far", "test close"),
         fill = c(col_train, col_discarded, col_close, col_far))
}
  

#' Split data between train and test
#' @description Split a GeoNonStat object between test and train samples
#'
#' @param locs a matrix of spatial locations.
#' @param data a list containing vectors or data.frames to split into
#' train/test
#' @param n_clust the number of clusters. If NULL, will be set to 1\% of the data.
#' @param prop_test_locs the proportion observsations used for validation. 
#' For "close" validation, the proportion is the proportion of spatial locations from the train data set
#' For "far" validation, the proportion is the proportion of spatial clusters. 
#' Default to 0.01. Can be of length 2 to change proportions between far 
#' and close test locations (far first, then close)
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
#' traintest <- splitData(observed_locs, datalist, prop_test = c(0.1, 0.01))
splitData <- function(locs, data,
                            n_clust = NULL,
                            prop_test = 0.1, 
                            round_locs = 0){
  if (length(prop_test) == 1) prop_test <- rep(prop_test, 2)
  if (round_locs > 0) locs <- round_locs * round(locs / round_locs)
  # Create indexes for train, test far and test close
  idx_far <- getIndexesFar(locs, n_clust = n_clust, prop_test = prop_test[1])
  idx_close <- getIndexesClose(locs = locs[idx_far$train, ], prop_test = prop_test[2])
  final_split <- list(
    "train" = idx_far$train[idx_close$train],
    "test_far" = idx_far$test,
    "test_close" = idx_far$train[idx_close$test],
    "discarded" = idx_far$discarded
  )
  info = rbind(
    c(idx_close$n_locs_train,      idx_close$n_obs_train),
    c(idx_far$n_locs_test,         idx_far$n_obs_test),
    c(idx_close$n_locs_test,       idx_close$n_obs_test),
    c(idx_far$n_locs_discarded,    idx_far$n_obs_discarded))
  info = rbind(info, apply(info, 2, sum))
  row.names(info) = c("train", "far", "close", "discarded", "total")
  colnames(info) = c("spatial locations", "observations")
  output <- list(
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
    ), 
    "discarded" = c(
      list(locs = locs[final_split$discarded, ]),
      subsetData(data, final_split$discarded)
    ), 
    "info" = info
  )
  cat("Summary of split data\n")
  splitSummary(output)
  plotSplit(output)
  return(output)
}
