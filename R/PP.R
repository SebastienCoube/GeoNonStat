#' @title Generate spatial knots using k-means clustering
#' @description Selects a set of spatial knots using k-means clustering from the observed locations.
#' If the number of locations is large, a subsample of up to 10,000 points is used,
#' and a small random perturbation is added to break ties and improve cluster separation.
#'
#' @param knots_number Integer. The number of knots (i.e., clusters) to generate.
#' @param locs A matrix of spatial coordinates (typically with two columns for 2D locations).
#'
#' @return A matrix of size \code{knots_number} x ncol(\code{locs}) containing the spatial knot coordinates (cluster centers).
#' @export
#' @keywords internal
#' @examples
#' locs <- cbind(runif(5000), runif(5000))
#' knots <- knotsFromKmeans(100, locs)
#' plot(locs, col = "grey", pch = 16, cex = 0.5)
#' points(knots, col = "red", pch = 19)
knotsFromKmeans <- function(knots_number, locs) {
  n_sample <- min(nrow(locs), 50000)
  if (knots_number >= n_sample) {
    stop("49999 knots maximum allowed")
  }
  sampled_locs <- locs[sample(seq_len(nrow(locs)), n_sample, replace = FALSE), ]

  # TODO ici j'ai un warning Les etapes de transfer (quick-TRANSfer stage) ont depasse le maximum (= 2500000)
  # Si on met l'algorithme "Lloyd" ca resoud le pb.
  centers <- kmeans(sampled_locs,
    knots_number,
    algorithm = "Hartigan-Wong",
    iter.max = 50
  )$centers
  return(centers)
}

#' @title Create a Predictive Process (PP) object
#' @description
#' Creates an object of class `PP`, a low-rank predictive process used as a prior
#' to describe spatial variations of covariance parameters.
#'
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param matern_range numeric (optional). Matern range parameter for the PP.
#' If `NULL`, a default value is inferred from spatial locations.
#' @param knots Either:
#'   * a matrix of spatial knots,
#'   * a positive integer (number of knots, determined via k-means), or
#'   * `NULL` (default) - in which case a default number and placement of knots is used.
#' @param reorder_knots logical, default to FALSE. Should knots be reordered ?
#' If so, the knots will be reordered using `GpGp::order_maxmin`
#' @param seed integer value, seed used for reproducibility purposes. Default to 1
#' @param plot Logical, whether to produce diagnostic plots (default `TRUE`).
#'
#' @rdname PP
#' @aliases PP
#' @return An object of class `PP`.
#' @export
#'
#' @examples
#' vecchia_approx <- createVecchia(cbind(runif(10000), runif(10000)), 10)
#' # automatic
#' pepito <- createPP(vecchia_approx)
#' \dontrun{
#' # choosing manually Matern range, too small wrt number of knots
#' pepito <- createPP(vecchia_approx, matern_range = .1)
#' # choosing manually Matern range, way too small wrt number of knots
#' pepito <- createPP(vecchia_approx, matern_range = .01)
#' # choosing manually number of knots in order to adjust to Matern range
#' pepito <- createPP(vecchia_approx, knots = 1000, matern_range = .1)
#' # choosing manually number of knots, but picking too few for default Matern range
#' pepito <- createPP(vecchia_approx, knots = 20)
#' # choosing manually Matern range in order to adjust to the number of knots
#' pepito <- createPP(vecchia_approx, knots = 20, matern_range = .5)
#' # inputing an user-specified grid of knots
#' pepito <- createPP(vecchia_approx,
#'   knots = as.matrix(
#'     expand.grid(
#'       seq(-.05, 1.05, .05),
#'       seq(-.05, 1.05, .05)
#'     )
#'   )
#' )
#' # inputing an user-specified grid of knots in order to adjust to small Matern range
#' pepito <- createPP(vecchia_approx,
#'   knots = as.matrix(
#'     expand.grid(
#'       seq(-.05, 1.05, .05),
#'       seq(-.05, 1.05, .05)
#'     )
#'   ),
#'   matern_range = .1
#' )
#' pepito <- createPP(vecchia_approx,
#'   knots = as.matrix(
#'     expand.grid(
#'       seq(-.05, 1.05, .05),
#'       seq(-.05, 1.05, .05)
#'     )
#'   ),
#'   matern_range = .05
#' )
#' pepito <- createPP(vecchia_approx,
#'   knots = as.matrix(
#'     expand.grid(
#'       seq(-.05, 1.05, .025),
#'       seq(-.05, 1.05, .025)
#'     )
#'   ),
#'   matern_range = .05
#' )
#' }


createPP <- function(vecchia_approx,
                     matern_range = NULL,
                     knots = NULL,
                     reorder_knots = TRUE,
                     seed = 1234,
                     plot = TRUE) {
  if (!is.list(vecchia_approx)) {
    stop("Argument 'vecchia_approx' must be a list.")
  }
  if (!is.null(matern_range) &&
    matern_range <= 0) {
    stop("'matern_range' must be positive.")
  }
  if (!is.null(knots) &&
    is.numeric(knots) && is.vector(knots) == 1 && any(knots <= 0)) {
    stop("'knots' must be positive.")
  }
  if (is.matrix(knots) &&
    ncol(knots) != ncol(vecchia_approx$locs)) {
    stop("Matrix 'knots' must have the same number of columns as spatial locations.")
  }
  if (is.vector(knots) && knots > nrow(vecchia_approx$locs)) {
    warning("You requested more knots than spatial locations.")
  }

  # Generate knots
  if (is.null(knots)) {
    knots <- 25
    message(paste("number of knots set to", knots))
  }
  if (!is.matrix(knots)) {
    knots <- knotsFromKmeans(knots, vecchia_approx$locs)
    message("knot placement done by default using k-means")
  }

  # matern range
  if (is.null(matern_range)) {
    matern_range <- max(dist(knots)) * .25
    message(
      paste("Matern range set to ", signif(matern_range, 3)),
      " 25 % of the space pseudo-diameter"
    )
  }

  # knots order
  if (is.null(rownames(knots))) {
    rownames(knots) <- seq_len(nrow(knots))
  }
  if(reorder_knots) knots <- knots[GpGp::order_maxmin(knots), ]

  additional <- NULL
  if (nrow(vecchia_approx$NNarray) - 1 - nrow(knots) > 0) {
    additional <-
      GpGp::find_ordered_nn(vecchia_approx$locs,
        m = nrow(vecchia_approx$NNarray) - 1 - nrow(knots)
      )[, -1]
  }
  # NNarray
  NNarray <- rbind(
    cbind(
      GpGp::find_ordered_nn(knots, nrow(vecchia_approx$NNarray) - 1),
      matrix(NA, nrow(knots), max(
        0, nrow(vecchia_approx$NNarray) - nrow(knots)
      ))
    ),
    cbind(
      nrow(knots) + seq_len(vecchia_approx$n_locs),
      FNN::get.knnx(
        query = vecchia_approx$locs,
        data = knots,
        k = min(nrow(knots), nrow(vecchia_approx$NNarray) - 1)
      )$nn.index,
      additional + nrow(knots)
    )
  )

  # Cholesky matrix
  combined_locs <- rbind(knots, vecchia_approx$locs)
  Linv_vals <- GpGp::vecchia_Linv(
    covparms = c(1, matern_range, 1e-6),
    covfun_name = "matern15_isotropic",
    locs = combined_locs,
    NNarray = NNarray
  )
  notnaNNarray <- !is.na(NNarray)
  sparse_chol <- Matrix::sparseMatrix(
    i = row(NNarray)[notnaNNarray],
    j = NNarray[notnaNNarray],
    x = Linv_vals[notnaNNarray],
    triangular = TRUE
  )

  res <- structure(
    list(
      "knots" = knots,
      "matern_range" = matern_range,
      "sparse_chol" = sparse_chol,
      "t_sparse_chol" = Matrix::t(sparse_chol),
      "n_knots" = nrow(knots),
      "vecchia_locs" = vecchia_approx$locs
    ),
    class = "PP"
  )

  if (plot) {
    plot(res, mar_var_loss = TRUE)
  } else {
    # Just to print the diagnostic of var loss
    varloss <- varLossPP(res)
  }

  return(res)
}

#' Summary of a 'PP' object
#'
#' @param object an object of class \code{PP}
#' @param ... additional arguments, unused
#' @export
#' @rdname PP
#' @aliases PP
#' @method summary PP
#' @examples
#' vecchia_approx <- createVecchia(cbind(runif(100), runif(100)), 10)
#' pepito <- createPP(vecchia_approx, plot = FALSE)
#' summary(pepito)
summary.PP <- function(object, ...) {
  cat(
    "Object of class \'PP\' with",
    object$n_knots,
    "knots,",
    "based on",
    dim(object$vecchia_locs)[1],
    "locations,",
    "matern range =",
    object$matern_range
  )
  varLossPP(object)
  return(invisible(NULL))
}

#' @title Compute the percentage of marginal variance who is lost because of the use of a PP
#' @param x an object of class PP, create with `createPP`
#' @return a numeric vector
#' @export
#' @keywords internal
#' @examples
#' vecchia <- createVecchia(cbind(runif(1000), runif(1000)))
#' pepito <- createPP(vecchia, plot = FALSE)
#' varLossPP(pepito)
varLossPP <- function(x) {
  if (is.null(x$knots | is.null(x$sparse_chol))) {
    stop("x must contains 'knots' and 'sparse_chol'")
  }
  PP_mar_var <- apply(Matrix::solve(x$sparse_chol, Matrix::diag(
    nrow = nrow(x$sparse_chol), ncol = nrow(x$knots)
  )), 1, function(x) {
    sum(x^2)
  })
  # max(0) because of tiny numerical errors
  PP_mar_var <- (pmax(0, 1.000001 - PP_mar_var) / 1.000001) * 100
  mean_mar_var <- mean(PP_mar_var)
  msg <- if (mean_mar_var > 10) {
    "quite a bit of loss, and may be fixed by adding more knots or increasing the Matern range."
  } else if (mean_mar_var > 3) {
    "fairly good, but it might be improved by adding more knots or increasing the Matern range."
  } else {
    "great !"
  }
  message(
    round(mean_mar_var, 1),
    "% of marginal variance on average is lost with the use of a PP.\nThis is ",
    msg
  )

  return(PP_mar_var)
}


#' TODO : est-ce que cette fonction est utilisee ?
#'
#' @param PP an object of class `PP`
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param df logical, default to TRUE. Return result as a data.frame ?
getBasis <- function(PP, vecchia_approx, df = TRUE) {
  res <- xPPMultRight(
    X = NULL,
    PP = PP,
    vecchia_approx = vecchia_approx,
    Y = diag(1, PP$n_knots),
    permutate_PP_to_obs = TRUE
  )
  colnames(res) <- paste("Basis_", seq_len(ncol(res)), sep = "")
  if (any(res < .001)) {
    res[res < .001] <- 0
  }
  if (df) {
    res <- as.data.frame(res)
  } else {
    res <- as(res, "sparseMatrix")
  }
  return(res)
}
