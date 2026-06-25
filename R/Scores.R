#' Coverage intervals of the response variable
#' @param true_y a vector of true observations of the response variable
#' @param y_samples a matrix of pseudo-observations sampled following the model
#'
#' @returns a vector of Booleans indicating if each true observation is within the predicted coverage interval
#' @export
#' @keywords internal
#'
#' @examples
#' true_y = rnorm(100)
#' y_samples = matrix(rnorm(1000*length(true_y)), length(true_y))
#' coverage(true_y, y_samples)
coverage <- function(true_y, y_samples){
  quantiles <- apply(y_samples, 1, function(x)quantile(x, c(.025, .975)))
  per_obs <- mapply(
    inBounds, 
    x = true_y, 
    bounds = split(quantiles, col(quantiles)))
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
#' y_samples <- matrix(rnorm(1000*length(true_y)), length(true_y))
#' CRPS(true_y, y_samples)
CRPS <- function(true_y, y_samples){
  ecdfs <- apply(y_samples, 1, ecdf)
  per_obs <- mapply(
    function(f, x)f(x), 
    x = true_y, 
    f = ecdfs)
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
#' mean_samples <- matrix(rep(0, 1000*length(true_y)), length(true_y))
#' sd_samples <- matrix(rep(1, 1000*length(true_y)), length(true_y))
#' predictiveDensity(true_y, mean_samples, sd_samples)
predictiveDensity <- function(true_y, mean_samples, sd_samples) {
  # Gaussian log-dens of the observations
  per_obs =  log(apply(dnorm(true_y, mean_samples, sd_samples), 1, mean))
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
#' mean_samples <- matrix(rep(0, 1000*length(true_y)), length(true_y))
#' MSE(true_y, mean_samples) 
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
#' mean_samples <- matrix(rep(0, 1000*length(true_y)), length(true_y))
#' sd_samples <- matrix(rep(1, 1000*length(true_y)), length(true_y))
#' score(true_y, mean_samples, sd_samples)

score <- function(true_y, mean_samples, sd_samples){
  y_samples = mean_samples + sd_samples * rnorm(length(sd_samples))
  per_obs = matrix(0, length(y_true), 4)
  colnames(per_obs) = c("95% coverage", "CRPS", "predictive density", "squared error")
  per_obs[,"95% coverage"] = coverage(true_y, y_samples)
  per_obs[,"CRPS"] = CRPS(true_y, y_samples)
  per_obs[,"predictive density"] = predictiveDensity(true_y, mean_samples, sd_samples)
  per_obs[,"squared error"] = squaredError(true_y, mean_samples)
  res = list()
  res$criteria <- signif(apply(per_obs, 2, mean), 2)
  names(res$criteria) = c("95% coverage", "CRPS", "ELPPD", "MSE")
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
#' @param geo_non_stat, an object of class GeoNonStat in order to evaluate the scores on the training data set, or NULL
#' @param prediction, a prediction from a GeoNonStat object in order to evaluate the scores on the test data set, or NULL
#' @param new_y a vector on new observations from the interest variable in order to evaluate the scores on the test data set, or NULL
#' @export

scores <- function(geo_non_stat = NULL, prediction = NULL, test_y = NULL){
  train <- !is.null(geo_non_stat)&is.null(prediction)&is.null(test_y)&(class(geo_non_stat) == "GeoNonStat")
  test <- is.null(geo_non_stat)&!is.null(prediction)&!is.null(test_y)
  if(!(train|test))stop("Either: geo_non_stat must be from class GeoNonStat, prediction be NULL, and test_y be NULL; or geo_non_stat must be NULL, prediction be a prediction from a GeoNonStat model, and y be test observations"
  # TODO check test_y and prediction same length
  )
}

#' Find scattered spatial locations for train-test partition
#' @description 
#'  A max-min heuristic is used to scatter the test locations across space,
#'  so that test locations are independent-ish from each other, and surrounded
#'  by training locations. 
#' @param observed_locs a matrix of spatial locations.
#' @param n_test the number of test locations. Some locations can be observed 
#' several times, so the number of test indices can be greater than n_test.  
#' @returns A list of train and test indices
#' @export
#' @example
#' 
#' # duplicated observed locations  
#' observed_locs <- cbind(runif(20000), runif(20000))
#'  observed_locs <- rbind(observed_locs, observed_locs)
#' # number of non-duplicated test locations
#'  n_test = 100
#' # train-test split
#' idx <- testIndicesClose(observed_locs, n_test = n_test)
#' plot(observed_locs, pch = ".")
#' # checking number of non-duplicated test locations
#' points(observed_locs[idx$test,], pch = 16, col=2)
#' length(idx$test)
#' sum(!duplicated(observed_locs[idx$test,]))
  
testIndicesClose <- function(observed_locs, n_test = NULL){
  locs  <- unique(observed_locs)
  if(is.null(n_test))n_test <- ceiling(nrow(locs)/10)
  test <- GpGp::order_maxmin(locs)[seq(n_test)]
  locs_match = match(split(observed_locs, row(observed_locs)), split(locs, row(locs)))
  return(list("train" = which(!locs_match%in%test) , "test_close" = which(locs_match%in%test)))
}

#' Find scattered spatial locations for train-test partition and remove 
#' train locations who are too close to test locations in order to assess long 
#' range prediction
#' @description 
#'  A K-means algorithm is run in order to parition the spatial locations into
#'  clusters. A max-min heuristic is used on the clusters to scatter the 
#'  test clusters across space, so that test clusters are independent-ish from 
#'  each other, and surrounded by training locations. The observation that is 
#'  closest to the center of each cluster is kept for validation and the rest of
#'  the cluster is discarded, in order to forbid immediate interpolation of 
#'  the test data.
#' @param observed_locs a matrix of spatial locations.
#' @param n_clust the number of test locations.
#' @param n_test the number of test clusters.
#' @export
#' 
#' @example
#'  observed_locs <- cbind(runif(20000), runif(20000))
#'  observed_locs <- rbind(observed_locs, observed_locs)
#'  plot(observed_locs, pch = ".")
#'  idx <- testIndicesFar(observed_locs)
#'  plot(observed_locs[idx$train,], pch = ".", xlab = "", ylab = "")
#'  points(observed_locs[idx$discarded,], pch = ".", xlab = "", ylab = "", col=  2)
#'  points(observed_locs[idx$test,], pch = 16, xlab = "", ylab = "", col=  4, cex = .5)
testIndicesFar <- function(observed_locs, n_clust = 200, n_test = 40){
  locs  <- unique(observed_locs)
  K = kmeans(locs, centers = n_clust)
  test_clust = GpGp::order_maxmin(K$centers)[seq(n_test)]
  test = FNN::knnx.index(data = locs, k = 1, query = K$centers[test_clust,])
  discarded = setdiff(which(K$cluster %in% test_clust), test)
  train = which(! K$cluster %in% test_clust)
  # intersect(train, discarded)
  # intersect(train, test)
  # intersect(test, discarded)
  locs_match = match(split(observed_locs, row(observed_locs)), split(locs, row(locs)))
  return(list("train" = which(locs_match%in%train), 
              "test_far" = which(locs_match%in%test), 
              "discarded" = which(locs_match%in%discarded)))
}

# TODO testIndicesCloseAndFar
# TODO plotTestIndices(test_indices, locs)
# TODO split(y, observed_locs, X, X_range, 
#            X_noise, train_test_indices)



#' Split data between train and test
#' @description
#'  
#' @param geo_non_stat, an object of class GeoNonStat in order to evaluate the scores on the training data set, or NULL
#' @param prediction, a prediction from a GeoNonStat object in order to evaluate the scores on the test data set, or NULL
#' @param new_y a vector on new observations from the interest variable in order to evaluate the scores on the test data set, or NULL
#' @export
