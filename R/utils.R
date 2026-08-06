
#' Subset a list containing matrices or a data.frames objects.
#' If NULL, return NULL
#'
#' @param data a list containing matrices or data.frames objects.
#' @param idx a numeric vector of indexes to filter rows
#'
#' @returns a list
#' @keywords internal
subsetData <- function(data, idx) {
  lapply(data, function(x) {
    if (is.null(x)) return(NULL) 
    if(is.data.frame(x)) return(x[idx, , drop=FALSE])
    if(is.vector(x)) return(x[idx])
    return(stop("subsetData only works on data.frame or vector"))
  })
}

#' Exponential of a square matrix, adding a small numeric value on the diagonal
#'
#' @param coords a numeric vector, of length such as
#' `(sqrt(8 * length(coords) + 1) - 1) / 2` is an integer
#' @param eps a numeric value, default to .0001
#'
#' @returns a square matrix
#' @export
#' @keywords internal
#'
#' @examples
#' expmat(c(1, 2, 3, 4, 5, 6))
expmat <- function(coords, eps = 0.0001) {
  res <- expm::expm(symmat(coords))
  res + diag(eps, nrow(res), ncol(res))
}


#' Create symetric matrix from coordinates
#'
#' @param coords a numeric vector, of length such as
#' `(sqrt(8 * length(coords) + 1) - 1) / 2` is an integer
#'
#' @returns a matrix
#' @export
#' @keywords internal
#'
#' @examples
#' symmat(c(1, 2, 3, 4, 5, 6))
symmat <- function(coords) {
  # check if length of vector is compatible
  n <- (sqrt(8 * length(coords) + 1) - 1) / 2
  if (as.integer(n) != n || n <= 0) {
    stop("length of coords incompatible with a symetric matrix")
  }
  symmat <- matrix(0, nrow = n, ncol = n)
  diag(symmat) <- coords[1:n]
  symmat[lower.tri(symmat)] <- coords[-seq_len(n)]
  symmat[upper.tri(symmat)] <- symmat[lower.tri(symmat)]
  return(symmat)
}

#' Title
#'
#' @param M a SparseMatrix
#'
#' @returns a vector of colors, of length the number of rows of M
#' @export
#' @keywords internal
#' @examples
#' n <- 5
#' i <- c(rep(1, n), 2:n) # 1re ligne + 1re colonne sauf la 1re case
#' j <- c(1:n, rep(1, n - 1))
#' M <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n, n))
#' naiveGreedyColoring(M)
#' M[2, 3] <- 1
#' M[3, 2] <- 1
#' naiveGreedyColoring(M)
naiveGreedyColoring <- function(M) {
  # number of nodes
  if (!Matrix::isSymmetric(M, checkDN = FALSE)) {
    stop("M must be symmetric")
  }
  n_obs <- nrow(M)
  # deducting degrees
  degrees <- as.vector(rep(1, n_obs) %*% M)
  # getting adjacent nodes of a given node
  idx <- split(M@i + 1, rep(seq_along(diff(M@p)), diff(M@p)))
  # creating a color * node matrix of incompatibilities
  incompatibilities <- matrix(0, n_obs + 1, max(degrees))
  cols <- rep(0, n_obs)

  for (i in seq_len(n_obs))
  {
    cols[i] <- match(0, incompatibilities[i, ])
    incompatibilities[idx[[i]], cols[i]] <- 1
  }
  return(cols)
}

#' Title TODO
#'
#' @param vecchia_approx  an object created by `vecchia_approx()`
#' @param compressed_sparse_chol an object created with `compute_sparse_chol()`
#'
#' @returns a sparse triangular matrix
#' @export
#' @keywords internal
decompressChol <- function(vecchia_approx, compressed_sparse_chol) {
  n <- length(vecchia_approx$sparse_chol_p) - 1L
  return(
    new("dtCMatrix",
      i = vecchia_approx$sparse_chol_i - 1L,
      p = vecchia_approx$sparse_chol_p,
      x = compressed_sparse_chol[, , 1][vecchia_approx$sparse_chol_x_reorder],
      Dim = c(n, n),
      uplo = "L",
      diag = "N"
    )
  )
}

#' Computes a Vecchia sparse Cholesky factor and its derivatives
#'
#' @param range_beta parameter for the range.
#' If the covariance is anisotropic, it must have 3 columns. It the covariance is isotropic, it must have 1 column.
#' The first coefficients are multiplied with range_X
#' The last coefficients are multiplied with the spatial basis functions of PP
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param range_X covariates for range, treated using process_covariates
#' @param PP predictive process obtained through `createPP()`
#'
#' @returns a list
#' @export
#' @keywords internal
#'
#' @examples
#' locs <- cbind(runif(1000), runif(1000))
#' vecchia_approx <- createVecchia(locs)
#' X <- data.frame(rnorm(nrow(locs)))
#' range_X <- processCovariates(X, vecchia_approx)
#' range_beta <- matrix(c(1, -2))
#' PP <- suppressMessages(createPP(vecchia_approx, plot = FALSE))
#' lr <- computeLogRange(
#'   range_beta = range_beta,
#'   vecchia_approx = vecchia_approx,
#'   range_X = range_X,
#'   PP = NULL
#' )
#' # test with a PP
#' range_beta <- matrix(rnorm(ncol(range_X$X_locs) + PP$n_knots), ncol = 1)
#' lr_with_PP <- computeLogRange(
#'   range_beta = range_beta,
#'   vecchia_approx = vecchia_approx,
#'   range_X = range_X,
#'   PP = PP
#' )
#' range_beta <- matrix(rnorm(3 * (2 + PP$n_knots)), ncol = 3)
#' lr_range3col <- computeLogRange(
#'   range_beta = range_beta,
#'   vecchia_approx = vecchia_approx,
#'   range_X = range_X,
#'   PP = PP
#' )
computeLogRange <- function(range_beta,
                            vecchia_approx,
                            range_X,
                            PP = NULL) {
  if (ncol(range_beta) == 3) {
    Y <- range_beta %*% matrix(c(1 / sqrt(2), 1 / sqrt(2), 0, 1 / sqrt(2), -1 / sqrt(2), 0, 0, 0, 1), 3) * sqrt(2) * 2
  } else if (ncol(range_beta) == 1) {
    Y <- range_beta * 2
  } else {
    stop("range_beta is expected to have 1 (isotropic case) or 3 (anisotropic case) columns")
  }

  log_range <- as.matrix(
    xPPMultRight(
      vecchia_approx = vecchia_approx,
      X = range_X$X_locs,
      PP = PP,
      Y = Y,
      permutate_PP_to_obs = FALSE
    )
  )
  return(log_range)
}


#' @title Multiply the concatenation of a matrix of covariates and a PP by a matrix
#' @description
#' Multiply the concatenation of a matrix of covariates and a PP by a matrix
#' `(X|PP) %*% Y`
#'
#' @param X a matrix of covariates who will be multiplied by the first rows of Y
#' @param PP either a PP whose basis will be multiplied by the last columns of Y
#' @param Y the matrix who multiplies the covariates and the PP
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param permutate_PP_to_obs logical, default to FALSE. TODO
#'
#' @returns a matrix
#' @export
#' @keywords internal
#'
#' @examples
#' locs <- cbind(runif(1000), runif(1000))
#' locs <- rbind(locs, locs)
#' vecchia_approx <- createVecchia(locs, 12)
#' PP <- suppressMessages(createPP(vecchia_approx, plot = FALSE))
#' covariate_coefficients <- c(4, 1, 1, .5)
#' knots_coeffs <- rnorm(PP$n_knots)
#' X <- cbind(1, vecchia_approx$locs, rnorm(nrow(vecchia_approx$locs)))
#'
#' # multiplying X alone
#' res1 <- xPPMultRight(X = X, Y = covariate_coefficients, vecchia_approx = vecchia_approx)
#' # multiplying PP alone
#' res2 <- xPPMultRight(PP = PP, Y = knots_coeffs, vecchia_approx = vecchia_approx)
#' # multiplying PP and matrix of covariate, one obs for each location
#' X_by_loc <- cbind(1, vecchia_approx$locs, rnorm(vecchia_approx$n_locs))
#' res3 <- xPPMultRight(
#'   PP = PP, X = X_by_loc,
#'   Y = c(covariate_coefficients, knots_coeffs),
#'   vecchia_approx = vecchia_approx
#' )
#' plotPointillistPainting(vecchia_approx$locs, res3,
#'   main = "PP + covariates,\n one covariate for each location", pch = 15
#' )
#'
#' # TODO Nun functional example
#' # # multiplying PP and matrix of covariates with an index
#' # X_by_obs = cbind(1, vecchia_approx$observed_locs, rnorm(vecchia_approx$n_obs))
#' # res4 <- xPPMultRight(X = X_by_obs, PP = PP,
#' #                        Y = c(covariate_coefficients, knots_coeffs),
#' #                        vecchia_approx = vecchia_approx)
#' # plotPointillistPainting(vecchia_approx$observed_locs,
#' #                           res4,
#' #                           main = "PP + covariates,\n  one covariate for each observation")
#'
#' # multiplying
#' res5 <- xPPMultRight(
#'   X = NULL, PP = PP,
#'   Y = diag(1, PP$n_knots),
#'   vecchia_approx = vecchia_approx,
#'   permutate_PP_to_obs = TRUE
#' )
xPPMultRight <- function(X = NULL,
                         PP = NULL,
                         vecchia_approx,
                         Y,
                         permutate_PP_to_obs = FALSE) {
  if (is.null(X) & is.null(PP)) stop("X and PP can't be both NULL")
  # Sanity checks
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  expected_rows <- 0
  if (!is.null(X)) expected_rows <- expected_rows + ncol(X)
  if (!is.null(PP)) expected_rows <- expected_rows + nrow(PP$knots)
  if (nrow(Y) != expected_rows) {
    stop("Y should have ", expected_rows, " rows it has ", nrow(Y))
  }

  if (permutate_PP_to_obs) {
    locs_idx <- vecchia_approx$locs_match
  } else {
    locs_idx <- seq_len(vecchia_approx$n_locs)
  }

  # Multiply X and Y
  xrow_offset <- 0
  res <- NULL

  if (!is.null(X)) {
    xrow_offset <- ncol(X)
    # using res[] is important for performance
    res <- X %*% Y[seq_len(xrow_offset), , drop = FALSE]
  }

  if (!is.null(PP)) {
    # remove X rows from Y if needed
    idxY_PP <- seq.int(xrow_offset + 1, nrow(Y))
    m <- length(idxY_PP)
    k <- nrow(PP$sparse_chol)
    V <- rbind(Y[idxY_PP, , drop = FALSE], matrix(0, k - m, ncol(Y)))
    # V <- matrix(0, nrow(PP$sparse_chol), ncol(Y))
    # V[seq_along(idxY_PP), ] <- Y[idxY_PP, , drop=FALSE]
    solved <- Matrix::solve(PP$sparse_chol, V, triangular = TRUE)
    PP_result <- solved[-seq_len(nrow(PP$knots)), , drop = FALSE]
    # using res[] is important for performance
    if (is.null(res)) {
      res <- PP_result[locs_idx, , drop = FALSE]
    } else {
      res[] <- res[] + PP_result[locs_idx, , drop = FALSE]
    }
  }

  if (!is.matrix(res)) {
    res <- as.matrix(res)
  }
  if (ncol(res) == 3) {
    colnames(res) <- c("det", "aniso1", "aniso2")
  } else if (ncol(res) == 1) {
    colnames(res) <- "det"
  }
  return(res)
}

#' @title Do the cross-product of the concatenation of a matrix of covariates and a PP and a matrix
#' @description Do the cross-product of the concatenation of a matrix of covariates and a PP and a matrix
#' `t(X|PP) %*% Y`
#' @param X a matrix of covariates who will be multiplied by the first rows of Y
#' @param PP either a PP whose basis will be multiplied by the last columns of Y
#' @param Y the matrix who multiplies the covariates and the PP
#' @param vecchia_approx an object created by `vecchia_approx()`
#' @param permutate_PP_to_obs logical, default to FALSE. TODO
#'
#' @returns a matrix
#' @export
#' @keywords internal
#'
#' @examples
#' set.seed(123)
#' locs <- cbind(runif(50), runif(50))
#' vecchia_approx <- createVecchia(locs)
#' PP <- createPP(vecchia_approx, plot = FALSE)
#' X <- matrix(rnorm(100), 50)
#' Y <- matrix(rnorm(30 * nrow(X)), nrow(X))
#'
#' # just surrogate of crossprod
#' res1 <- xPPCrossprod(X = X, Y = Y)
#' identical(crossprod(X, Y), res1)
#'
#' # crossprod + PP with observations of X on the locs
#' res2 <- xPPCrossprod(
#'   X = X, PP = PP, Y = Y,
#'   vecchia_approx = vecchia_approx,
#'   permutate_PP_to_obs = FALSE
#' )
#'
#' # crossprod + PP with observations of X on the obs
#' res3 <- xPPCrossprod(
#'   X = X, PP = PP, Y = Y,
#'   vecchia_approx = vecchia_approx,
#'   permutate_PP_to_obs = TRUE
#' )
xPPCrossprod <- function(X,
                         PP = NULL,
                         Y,
                         vecchia_approx = NULL,
                         permutate_PP_to_obs = FALSE) {
  if(!is.null(X)){if (nrow(X) != nrow(Y)) {
    stop("X and Y should have the same number of rows")
  }}
  if (permutate_PP_to_obs && is.null(vecchia_approx)) {
    stop("To permutate PP to observed values vecchia_approx needs to be provided.")
  }
  if (!is.null(PP)) {
    if(!is.null(X)){if ((nrow(X) != nrow(PP$vecchia_locs)) & (!permutate_PP_to_obs)) {
      stop("X should have the same number of rows as locations in vecchia (permutate_PP_to_obs = FALSE)")
    }}
    if(!is.null(X)){if ((nrow(X) != vecchia_approx$n_obs) & (permutate_PP_to_obs)) {
      stop("X should have the same number of rows as observations in vecchia (permutate_PP_to_obs = TRUE)")
    }}
    if (is.null(vecchia_approx)) {
      stop("PP provided, vecchia_approx should be provided too")
    }
    if (vecchia_approx$n_locs != nrow(PP$vecchia_locs)) {
      stop("vecchia_approx should have the same number of locations as the locations of PP")
    }
  }
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  res <- NULL
  if(!is.null(X)){res <- crossprod(x = X, y = Y)}
  if (!is.null(PP)) {
    if (permutate_PP_to_obs) {
      Y <- vecchia_approx$locs_match_matrix %*% Y
    }
    res <-
      rbind(res, Matrix::solve(PP$t_sparse_chol, rbind(matrix(
        0, nrow(PP$knots), ncol(Y)
      ), Y))[1:PP$n_knots, , drop = FALSE])
  }
  if (!is.matrix(res)) {
    res <- as.matrix(res)
  }
  return(res)
}


#' Transpose a list
#'
#' @param list a list containing sublists (all having the same structure).
#'
#' @returns a list, transposed
#' @export
#' @keywords internal
#'
#' @examples
#' A <- list("R1" = "a1", "R2" = "a2", "R3" = "a3")
#' B <- list("R1" = "b1", "R2" = "b2", "R3" = "b3")
#' C <- list("R1" = "c1", "R2" = "c2", "R3" = "c3")
#' transposeList(list("A" = A, "B" = B, "C" = C))
transposeList <- function(list) {
  nameslist <- names(list[[1]])
  tl <- lapply(
    seq_along(list[[1]]),
    function(i) lapply(list, `[[`, i)
  )
  names(tl) <- nameslist
  return(tl)
}
