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


#' Return averaged and per-observation model comparisons 
#' @param true_y a vector of true observations of the response variable
#' @param mean_samples a matrix of samples of mean from the model
#' @param sd_samples a matrix of samples of standard deviation from the model
#' @export
#' @keywords internal
#' @examples
#' true_y <- rnorm(100)
#' mean_samples <- matrix(rep(0, 1000*length(true_y)), length(true_y))
#' sd_samples <- matrix(rep(1, 1000*length(true_y)), length(true_y))
#' doTheScores(true_y, mean_samples, sd_samples)

doTheScores <- function(true_y, mean_samples, sd_samples){
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
#'  and the \href{https://en.wikipedia.org/wiki/Mean_squared_error}{Mean Squared Error} (MSE).\n 
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