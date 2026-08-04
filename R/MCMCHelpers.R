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

updateBeta <- function(params, stuff, X, vecchia_approx, observed_field) {
  # centered parametrization of latent field
  beta_covmat <- solve(crossprod(X$X / stuff$noise_var, X$X))
  if (all(!is.infinite(beta_covmat) & !is.nan(beta_covmat))) {
    if (all(eigen(beta_covmat)$val > 0)) {
      beta_mean <- beta_covmat %*% crossprod(X$X, ((observed_field - params$field[vecchia_approx$locs_match]) / stuff$noise_var))
      params$beta[] <- (beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))[]
    }
  }
  # un-centered parametrization of latent field
  centered_field <- as.vector(params$field + X$X_locs %*% matrix(params$beta[X$which_locs], ncol = 1))
  sparse_chol_X <- as.matrix(stuff$sparse_chol %*% (X$X_locs)) / exp(.5 * params$field_log_var[1, 1])
  beta_precision <- crossprod(x = sparse_chol_X, y = sparse_chol_X)
  beta_covmat <- solve(beta_precision, tol = min(rcond(beta_precision), .Machine$double.eps))
  if (all(!is.infinite(beta_covmat) & !is.nan(beta_covmat))) {
    if (all(eigen(beta_covmat)$d > 0)) {
      beta_mean <- c(as.vector(stuff$sparse_chol %*% (centered_field / exp(.5 * params$field_log_var[1, 1]))) %*% sparse_chol_X %*% beta_covmat)
      params$beta[X$which_locs] <- as.vector(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
      params$field <- centered_field - as.vector(X$X_locs %*% matrix(params$beta[X$which_locs], ncol = 1))
    }
  }
  # updating state$stuff
  stuff$lm_fit[] <- as.vector(X$X %*% params$beta)
  stuff$lm_residuals <- observed_field - stuff$lm_fit
  return(list("params" = params, "stuff" = stuff))
}


updateLatentField <- function(state, vecchia_approx, hierarchical_model, observed_field, iter, num_threads) {
  # Field log var ###############################
  cluster_idx <- 1
  chosen_locs_partition_idx <- 1 + iter %% ncol(vecchia_approx$locs_partition)
  locs_partition <- vecchia_approx$locs_partition[, chosen_locs_partition_idx]
  precision_from_obs <- vecchia_approx$locs_match_matrix %*% (1 / state$stuff$noise_var)
  mean_from_obs <- as.vector(vecchia_approx$locs_match_matrix %*% ((observed_field - state$stuff$lm_fit) / state$stuff$noise_var))

  d <- 1 / exp(.5 * state$params$field_log_var[1, 1])
  unique_clusters <- unique(locs_partition)

  chol_list <- lapply(
    unique(locs_partition),
    function(cluster_idx) { # posterior_precision_subset
      idx <- which(locs_partition == cluster_idx)
      pbs <- Matrix::crossprod(state$stuff$sparse_chol[, idx])
      pbs <- Matrix::Diagonal(length(idx), d) %*% pbs %*% Matrix::Diagonal(length(idx), d)
      Matrix::diag(pbs) <- Matrix::diag(pbs) + as.vector(precision_from_obs[idx])
      pbs <- Matrix::expand(Matrix::Cholesky(pbs))
      return(pbs)
    }
  )

  for (asdf in seq_len(1)) {
    # field itself
    for (cluster_idx in unique(locs_partition)) {
      selected_idx <- which(locs_partition == cluster_idx)
      additional_mean_from_field <-
        as.vector(
          (1 / exp(.5 * state$params$field_log_var[1, 1])) *
            Matrix::crossprod(
              state$stuff$sparse_chol[, selected_idx],
              (state$stuff$sparse_chol %*%
                ((state$params$field * (locs_partition != cluster_idx)) / exp(.5 * state$params$field_log_var[1, 1]))
              )
            )
        )
      chololo <- chol_list[[match(cluster_idx, unique(locs_partition))]]
      state$params$field[selected_idx] <- as.vector(
        Matrix::t(chololo$P) %*%
          Matrix::solve(
            Matrix::t(chololo$L),
            rnorm(nrow(chololo$L)) +
              Matrix::solve(
                chololo$L, # inverse of precision matrix...
                chololo$P %*%
                  (-additional_mean_from_field + mean_from_obs[selected_idx])
              )
          )
      )
    }
  }
  return(list(field = state$params$field))
}


updateFieldLogVar <- function(state, scale, vecchia_approx, iter, iter_start) {
    # ancillary
    for (field_log_var_idx in seq_len(2)) {
      new_field_log_var <- state$params$field_log_var[1, 1] + exp(.5 * state$ker_var$field_log_var_ancillary) * rnorm(1)
      new_field <- state$params$field * exp(.5 * (new_field_log_var - state$params$field_log_var[1, 1]))
      current_U <-
        (
          -betaPriorLogDens(
            beta = as.matrix(state$params$field_log_var[1, 1]), n_PP = 0, log_scale = 0,
            beta0_mean = scale$beta0_mean,
            beta0_var = scale$beta0_sd^2
          ) # normal prior
          + .5 * sum((state$stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2 / state$stuff$noise_var) # observation ll
        )

      proposed_U <-
        (
          -betaPriorLogDens(
            beta = as.matrix(new_field_log_var), n_PP = 0, log_scale = 0,
            beta0_mean = scale$beta0_mean,
            beta0_var = scale$beta0_sd^2
          ) # normal prior
          + .5 * sum((state$stuff$lm_residuals - new_field[vecchia_approx$locs_match])^2 / state$stuff$noise_var) # observation ll
        )
      state$ker_var$field_log_var_ancillary <- updateKernel(
        iter_start = iter_start,
        kernel_value = state$ker_var$field_log_var_ancillary, iter = iter, mult = -.25
      )
      if (current_U - proposed_U > log(runif(1))) {
        state$ker_var$field_log_var_ancillary <- updateKernel(
          iter_start = iter_start,
          kernel_value = state$ker_var$field_log_var_ancillary, iter = iter, mult = 1
        )
        state$params$field_log_var[1, 1] <- new_field_log_var
        state$params$field <- new_field
      }
    }

    # Sufficient
    fieldT_cholT_chol_field <- sum((state$stuff$sparse_chol %*% state$params$field)^2)
    for(repp in seq(5)){
    for (field_log_var_idx in seq_len(5)) {
      new_field_log_var <- state$params$field_log_var[1, 1] + rnorm(1) * .5^(field_log_var_idx ) # exp(state$ker_var$field_log_var_sufficient)
      current_U <-
        (
          -betaPriorLogDens(
            beta = as.matrix(state$params$field_log_var[1, 1]), n_PP = 0, log_scale = 0,
            beta0_mean = scale$beta0_mean,
            beta0_var = scale$beta0_sd^2
          ) # normal prior
          + .5 * fieldT_cholT_chol_field / exp(state$params$field_log_var[1, 1]) # observation ll
            + vecchia_approx$n_locs * (.5 * state$params$field_log_var[1, 1]) # observation ll
        )

      proposed_U <-
        (
          -betaPriorLogDens(
            beta = as.matrix(new_field_log_var), n_PP = 0, log_scale = 0,
            beta0_mean = scale$beta0_mean,
            beta0_var = scale$beta0_sd^2
          ) # normal prior
          + .5 * fieldT_cholT_chol_field / exp(new_field_log_var) # observation ll
            + vecchia_approx$n_locs * (.5 * new_field_log_var) # observation ll
        )
      state$ker_var$field_log_var_sufficient <- updateKernel(
        iter_start = iter_start,
        kernel_value = state$ker_var$field_log_var_sufficient, iter = iter, mult = -.25
      )
      if (current_U - proposed_U > log(runif(1))) {
        state$ker_var$field_log_var_sufficient <- updateKernel(
          iter_start = iter_start,
          kernel_value = state$ker_var$field_log_var_sufficient, iter = iter, mult = 1
        )
        state$params$field_log_var[1, 1] <- new_field_log_var
      }
    }
    }
  
  return(state)
}


updateVarPPSuff <- function(hm4params, beta4params, current_range_log_scale) {
  for (i in seq_len(10))
  {
    q <- current_range_log_scale + rnorm(length(current_range_log_scale), 0, .05)
    if (all(inBounds(hm4params$log_scale_bounds, q)) &
      (
        +betaPriorLogDens(
          beta = beta4params, n_PP = hm4params$PP$n_knots,
          beta0_mean = hm4params$beta0_mean,
          beta0_var = hm4params$beta0_sd^2,
          log_scale = q
        )
        - betaPriorLogDens(
            beta = beta4params, n_PP = hm4params$PP$n_knots,
            beta0_mean = hm4params$beta0_mean,
            beta0_var = hm4params$beta0_sd^2,
            log_scale = current_range_log_scale
          )
        > log(runif(1))
      )
    ) {
      current_range_log_scale[] <- q
    }
  }
  return(current_range_log_scale)
}

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

