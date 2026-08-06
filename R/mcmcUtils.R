renewMomentum <- function(momentum, kept_momentum = .9) {
  if (!inBounds(c(0, 1), kept_momentum, includeBounds = TRUE)) {
    stop("kept_momentum must be between 0 and 1")
  }
  momentum[] <- sqrt(kept_momentum) * momentum[] + sqrt(1 - kept_momentum) * rnorm(length(momentum[]))
  momentum
}

updateKernel <- function(iter,
                         iter_start,
                         kernel_value,
                         mult) {
  kernel_value <-
    kernel_value +
    length(kernel_value) * mult / sqrt(10 + iter + iter_start)
  kernel_value <- max(kernel_value, -30)
  kernel_value <- min(kernel_value, 2)
  kernel_value
}

#' Initiate a NA vector of matrix of the same size as the param given
#'
#' @param stateParams a list of vector or a matrix
#'
#' @returns a list of vector or a matrix
#' @export
#' @keywords internal
#'
#' @examples
#' initStateParams(list(
#'   rnorm(20),
#'   matrix(rnorm(20), ncol = 4)
#' ))
initStateParams <- function(stateParams) {
  lapply(stateParams, function(x) {
    tmp <- rep(NA_real_, length(x))
    if (is.matrix(x)) {
      dim(tmp) <- dim(x)
    }
    return(tmp)
  })
}

#' Creates a conditioning matrix from a default prior Fisher and an empirical Fisher matrix 
condMat <- function(empirical_fisher, prior_fisher, iter, iter_start){
  mix <- min(1, max((500 - (iter + iter_start)) / 200, .05))
  renorm <- function(M)M/(sqrt(sum(M^2))+1e-6) # Frobenius norm
  t(chol(solve(
    as(
      (1 - mix) * renorm(empirical_fisher) +
      (mix) *     renorm(prior_fisher),
      "sparseMatrix"
    )
  )))
}

# update a Fisher information matrix using gradients
updateFisher <- function(iter, iter_start, dens_grad, empirical_fisher){
  empirical_fisher <-
    ((iter+iter_start-1)/(iter+iter_start)) * empirical_fisher +
    tcrossprod((dens_grad)) / (iter+iter_start)
  empirical_fisher
}

# density of a momentum
momentumDens <- function(momentum)-.5*sum(momentum^2)

# MALA forward step
moveForward <- function(stepsize, cond_mat, dens_grad,
                        position, momentum){
  position +
    stepsize/2 *     as.vector((cond_mat %*% (crossprod(cond_mat, c(dens_grad))))) +
    sqrt(stepsize) * as.vector( cond_mat %*% c(momentum))
}

#' position <- rnorm(10)
#' momentum <- rnorm(10)
#' dens_grad <- rnorm(10)
#' dens_grad_back <- rnorm(10)
#' cond_mat <- t(chol(crossprod(as(matrix(rnorm(1000), 100), "sparseMatrix"))))
#' stepsize <- .0001
#' new_position <- moveForward(
#'   momentum = momentum, stepsize = stepsize, cond_mat = cond_mat,
#'   position = position, dens_grad = dens_grad)
#' plot(
#'   solve(
#'    sqrt(stepsize) * cond_mat,
#'   position - moveForward(
#'     momentum = 0*momentum, stepsize = stepsize, cond_mat = cond_mat,
#'     position = new_position, dens_grad = dens_grad_back)
#'   ),
#'   momentumBack(
#'     stepsize = stepsize, cond_mat = cond_mat, current_params = position,
#'     proposed_params = new_position, dens_grad_back = dens_grad_back
#'   ))
#' abline(a = 0, b = 1)

# Finding reverse momentum from a MALA step
momentumBack <- function(
    current_params, proposed_params,
    dens_grad_back,
    cond_mat, stepsize
){
  solve(
    cond_mat,
    c(
    c(current_params)
    - c(proposed_params)
    - as.vector(cond_mat %*% (crossprod(cond_mat, c(dens_grad_back)))) * stepsize / 2)
  )/sqrt(stepsize)
}



#' Computes density gradient of a log scale parameter using the Chain Rule
logScaleDensGrad <- function(PP_coeff, grad_PP_coeff){
  res <- .5 * PP_coeff * grad_PP_coeff
  if(ncol(PP_coeff)==3)return(c(sum(res[,1]), sum(res[,-1])))
  return(sum(res))
}

#' Extracts PP coefficients from regression coefficients
getPPCoeff <- function(range_beta, n_knots){
  return(range_beta[-seq_len(nrow(range_beta) - n_knots),])
}


#' Moves PP coefficients if a new log scale parameter is proposed
#' # Testing coherence of logScaleDensGrad formula and movePPCoeff
#' hierarchical_model <- list(anisotropic = T, range = list(PP = list(n_knots = 9)))
#' range_beta <- matrix(rnorm(30), 10)
#' range_log_scale <- c(-1, -3)
#' range_log_scale_ <- range_log_scale + 1e-6
#' range_beta_ <- movePPCoeff(range_beta, range_log_scale, range_log_scale_, hierarchical_model)
#' plot(
#' (range_beta_ -range_beta)[-1,]*1e6,
#' .5*range_beta[-1,]
#' )
movePPCoeff <- function(beta, log_scale, new_log_scale, n_knots){
  beta[-seq_len(nrow(beta) - n_knots),] <-
    beta[-seq_len(nrow(beta) - n_knots),] %*% 
    diag(exp(.5*(new_log_scale - log_scale)) [c(1, rep(2, 2*(ncol(beta)==3)))])
  return(beta)
}


#' Is a numeric include in an interval ?
#'
#' @param bounds Numeric of length 2, bounds of the interval.
#' @param x a numeric value to test (or a vector)
#' @param includeBounds logical (of length 1 or 2). Are the bounds included in the interval ?
#'
#' @returns a logical value, of same length as x
#' @export
#' @keywords internal
#'
#' @examples
#' inBounds(c(0, 12), 0, includeBounds = FALSE)
#' inBounds(c(0, 12), 12, includeBounds = TRUE)
inBounds <- function(bounds, x, includeBounds = FALSE) {
  if (length(bounds) != 2 || !is.numeric(bounds)) stop("bounds should be a numeric of length 2")
  if (length(includeBounds) == 1) includeBounds <- rep(includeBounds, 2)
  if (length(includeBounds) != 2 || !is.logical(includeBounds)) stop("includeBounds should be a logical of length 1 or 2")
  includeBounds <- 1 * includeBounds + 1
  res1 <- switch(includeBounds[1],
                 (x > bounds[1]),
                 (x >= bounds[1])
  )
  res2 <- switch(includeBounds[2],
                 (x < bounds[2]),
                 (x <= bounds[2])
  )
  res <- (res1 & res2)
  return(res)
}





#' Compute prior logarithmic density of a matrix
#'
#' @param beta a numeric matrix, with 1 (isotropic case) or 3 (anisotropic case) columns
#' @param n_PP number of PP
#' @param beta0_mean numeric value, mean of beta_0
#' @param beta0_var numeric value, mean of beta_0
#' @param log_scale numeric vector of length 1 for a 1 column beta matrix
#' and length 2 for a 3 column beta matrix
#' @description
#' beta follows a Normal distribution, with independent components.
#'
#' In the isotropic case, beta has one column.
#' It has 3 categories of rows :
#' 
#' - the first row, corresponding to the Intercept.
#'  
#' - the next rows, corresponding to the rest of the explanatory variables if there are any.
#'  
#' - the last rows, corresponding to the PP coefficients if there are any.
#' 
#' For example with \code{nrow(X)=5} and \code{n_PP=6} :
#' 
#' \preformatted{
#'                       The mean          The diagonal of the
#'                      is a vector        variance matrix is
#'                      with shape :      a vector with shape :
#'
#'    (Intercept)      beta_0_mean             beta_0_var
#'            X    |-       0                     0.01
#'            X    |        0                     0.01
#'            X   -|        0                     0.01
#'            X    |        0                     0.01
#'            X    |        0                     0.01
#'           PP    |        0            exp(range_log_scale)
#'           PP    |        0            exp(range_log_scale)
#'           PP   -|        0            exp(range_log_scale)
#'           PP    |        0            exp(range_log_scale)
#'           PP    |        0            exp(range_log_scale)
#'           PP    |        0            exp(range_log_scale)
#' }
#'
#' In the anisotropic case, beta has 3 columns, one for the Range, two for the Anisotropy
#' The categories of rows are the same as in the Isotropic case.
#'
#' The mean and the diagonal of the variance matrix are stored under
#' the shape of matrices with 3 columns, like beta.
#'
#' \preformatted{
#'                                     mean                                                 variance
#'
#'                      (range)      (aniso)    (aniso)
#'
#'    (Intercept)      beta_0_mean      0           0              beta_0_var                  0.01                      0.01
#'                 |-       0           0           0                 0.01                     0.01                      0.01
#'                 |        0           0           0                 0.01                     0.01                      0.01
#'            X   -|        0           0           0                 0.01                     0.01                      0.01
#'                 |        0           0           0                 0.01                     0.01                      0.01
#'                 |_       0           0           0                 0.01                     0.01                      0.01
#'                 |-       0           0           0        exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])
#'                 |        0           0           0        exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])
#'            PP  -|        0           0           0        exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])
#'                 |        0           0           0        exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])
#'                 |        0           0           0        exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])
#'                 |_       0           0           0        exp(range_log_scale[1])  exp(range_log_scale[2])   exp(range_log_scale[2])
#' }
#' @returns a numeric value
#' @export
#' @keywords internal
#'
#' @examples
#' betaPriorLogDens(
#'   beta = matrix(rnorm(300), 100, ncol = 3),
#'   n_PP = 90,
#'   beta0_mean = -5,
#'   beta0_var = 2,
#'   log_scale = c(-3, -2)
#' )
betaPriorLogDens <- function(beta,
                             n_PP,
                             beta0_mean,
                             beta0_var,
                             log_scale) {
  nrb <- nrow(beta)
  ncb <- ncol(beta)
  if (is.null(n_PP)) n_PP <- 0
  if (!ncb %in% c(1, 3)) stop("beta is expected to have 1 or 3 columns")
  if (n_PP > nrb + 2) stop("n_PP can't be greater than nrow(beta) + 2")
  
  mean_mat <- matrix(0, nrb, ncb)
  var_mat <- matrix(0.01, nrb, ncb)
  mean_mat[1, 1] <- beta0_mean # Intercept for range
  var_mat[1, 1] <- beta0_var # Intercept for range
  determinant_part = 0
  if(n_PP>0){
    var_mat[-seq_len(nrb - n_PP), 1] <- exp(log_scale[1])
    if (ncb == 3) {
      if (!length(log_scale) == 2) {
        stop("log_scale is supposed to be of length 2")
      }
      var_mat[-seq_len(nrb - n_PP), c(2, 3)] <- exp(log_scale[2])
    }
    determinant_part <- -0.5 * log_scale[1] * n_PP
    if (length(log_scale) == 2) {
      determinant_part <- determinant_part - 2 * 0.5 * log_scale[2] * n_PP
    }
  }
  quadratic_part <- -0.5 * sum((beta - mean_mat)^2 / var_mat)
  return(quadratic_part + determinant_part)
}

#' Compute gradient of logarithmic density prior
#'
#' @param beta a numeric matrix, with 1 (isotropic case) or 3 (anisotropic case) columns
#' @param n_PP number of prediction points.
#' @param beta0_mean numeric value, mean of beta_0
#' @param beta0_var numeric value, mean of beta_0
#' @param log_scale numeric vector of length 1 for a 1 column beta matrix
#' and length 2 for a 3 column beta matrix
#'
#' @returns an array
#' @export
#' @keywords internal
#' @examples
#' betaPriorLogDensDerivative(
#'   beta = matrix(rnorm(300), 100, 3),
#'   n_PP = 90,
#'   beta0_mean = -5,
#'   beta0_var = 2,
#'   log_scale = c(-3, -2)
#' )
betaPriorLogDensDerivative <- function(beta,
                                       n_PP,
                                       beta0_mean,
                                       beta0_var,
                                       log_scale) {
  nrb <- nrow(beta)
  ncb <- ncol(beta)
  if (!ncb %in% c(1, 3)) {
    stop("beta is expected to have 1 or 3 columns")
  }
  if(is.null(n_PP))n_PP <- 0
  if (n_PP > nrb + 2) {
    stop("n_PP can't be greater than nrow(beta) + 2")
  }
  mean_mat <- matrix(0, nrb, ncb)
  var_mat <- matrix(0.01, nrb, ncb)
  mean_mat[1, 1] <- beta0_mean
  var_mat[1, 1] <- beta0_var
  if(n_PP > 0){
    var_mat[-seq_len(nrb - n_PP), 1] <- exp(log_scale[1])
    if (ncb == 3) {
      if (!length(log_scale) == 2) {
        stop("log_scale is supposed to be of length 2")
      }
      var_mat[-seq_len(nrb - n_PP), c(2, 3)] <- exp(log_scale[2])
    }
  }
  return(-(beta - mean_mat) / var_mat)
}

