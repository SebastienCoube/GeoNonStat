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
  kernel_value <- max(kernel_value, -12)
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
  for (i in seq(2)) {
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
    for (field_log_var_idx in seq_len(5)) {
      new_field_log_var <- state$params$field_log_var[1, 1] + rnorm(1) * .5^(2 + field_log_var_idx %% 5) # exp(state$ker_var$field_log_var_sufficient)
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

updateRangeBetaAndLogVarMALA <- function(state, hierarchical_model, vecchia_approx,
                                         range_X, iter, iter_start, num_threads) {
  range_reparam_mat <- matrix(2)
  if (hierarchical_model$anisotropic) {
    range_reparam_mat <- matrix(c(
      2, 2, 0,
      2, -2, 0,
      0, 0, 2 * sqrt(2)
    ), 3)
  }

  ##########################
  # Range beta (ancillary) #
  ##########################
  
  # computing gradient
  dens_grad <- c(
    # log var
    (
      betaPriorLogDensDerivative(
        beta = as.matrix(state$params$field_log_var), n_PP = 0, log_scale = 0,
        beta0_mean = hierarchical_model$scale$beta0_mean,
        beta0_var = hierarchical_model$scale$beta0_sd^2
      )
      + .5 * sum(
        state$params$field[vecchia_approx$locs_match] *
          (state$stuff$lm_residuals - state$params$field[vecchia_approx$locs_match]) /
          state$stuff$noise_var
      ) # observation ll
    ),
    # field range
    -(
      -betaPriorLogDensDerivative(
        beta = state$params$range_beta, n_PP = hierarchical_model$range$PP$n_knots,
        beta0_mean = hierarchical_model$range$beta0_mean,
        beta0_var = hierarchical_model$range$beta0_sd^2,
        state$params$range_log_scale
      ) # normal prior
      + xPPCrossprod(
        X = range_X$X_locs, vecchia_approx = vecchia_approx,
        permutate_PP_to_obs = FALSE,
        PP = hierarchical_model$range$PP,
        Y = # Jacobian of range field wrt range_beta
          t(
            # natural gradient of obs likelihood wrt range field
            derivativeSandwiches_(
              vecchia = state$stuff$compressed_chol,
              left_vector = as.vector(
                Matrix::solve(
                  Matrix::t(state$stuff$sparse_chol),
                  -as.vector(vecchia_approx$locs_match_matrix %*% # gradient of  Gaussian observations ll wrt latent field
                               ((state$params$field[vecchia_approx$locs_match] - state$stuff$lm_residuals) / state$stuff$noise_var))
                  * exp(.5 * state$params$field_log_var[1, 1]) # part of sparse chol
                )
              ),
              right_vector = state$params$field / exp(.5 * state$params$field_log_var[1, 1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
              NNarray = vecchia_approx$NNarray,
              sauce_determinant_chef = FALSE,
              num_threads = num_threads
            )
          ) %*% range_reparam_mat
      )
    )
  )
  
  # conditioning matrix
  cond_mat <- solve(as(
    (1 - min(1, max((500 - (iter + iter_start)) / 200, .05))) *
      state$stuff$range_beta_conditioning_a / sqrt(sum(state$stuff$range_beta_conditioning_a^2))
      + Matrix::bdiag(
        1,
        min(1, max((500 - (iter + iter_start)) / 200, .05)) *
          (diag(rep(1, 1 + 2 * hierarchical_model$anisotropic)) %x% as.matrix(range_X$crossprod_X)) /
          sqrt(3 * sum(as.matrix(range_X$crossprod_X)^2))
      ),
    "sparseMatrix"
  ))
  cond_mat <- cond_mat / sqrt(sum(cond_mat^2))
  cond_mat <- t(chol(cond_mat))

  for (asdf in seq(1)) {
    # updating conditioning matrix
    if (iter + iter_start > 200) {
      state$stuff$range_beta_conditioning_a <-
        ((iter + iter_start) / (iter + iter_start - 1)) * state$stuff$range_beta_conditioning_a +
        tcrossprod(c(dens_grad)) / (iter + iter_start)
    }
    stepsize <- exp(state$ker_var$range_beta_ancillary[1])
    # updating MALA innovation with auto-correlation
    state$momenta$field_log_var_ancillary_grouped <- renewMomentum(state$momenta$field_log_var_ancillary_grouped)
    state$momenta$range_beta_ancillary <- renewMomentum(state$momenta$range_beta_ancillary)
    # making MALA step
    move_forward <-
      as.vector(stepsize * (cond_mat %*% (crossprod(cond_mat, c(dens_grad)))) / 2) +
      sqrt(stepsize) * as.vector(cond_mat %*% c(state$momenta$field_log_var_ancillary_grouped, state$momenta$range_beta_ancillary))
    new_range_beta <- state$params$range_beta + move_forward[-1]
    new_field_log_var <- state$params$field_log_var + move_forward[1]

    # re-computing Vecchia approx
    log_range <- computeLogRange(
      range_beta = new_range_beta,
      PP = hierarchical_model$range$PP,
      vecchia_approx = vecchia_approx,
      range_X = range_X
    )
    vecchia_(
      start_idx = 1,
      log_range = t(log_range), 
      locs = vecchia_approx$t_locs,
      NNarray = vecchia_approx$NNarray,
      smoothness = hierarchical_model$matern_smoothness,
      compute_derivative = TRUE, 
      num_threads = num_threads, 
      result = state$stuff$proposed_compressed_chol
    )
    state$stuff$proposed_sparse_chol@x <- state$stuff$proposed_compressed_chol[, , 1][vecchia_approx$sparse_chol_x_reorder]
    new_field <-
      as.vector(Matrix::solve(state$stuff$proposed_sparse_chol, state$stuff$sparse_chol %*% (state$params$field))) *
        exp(.5 * (new_field_log_var[1, 1] - state$params$field_log_var[1, 1]))

    # computing gradient
    dens_grad_back <- c(
      # log var
      (
        betaPriorLogDensDerivative(
          beta = as.matrix(new_field_log_var), n_PP = 0, log_scale = 0,
          beta0_mean = hierarchical_model$scale$beta0_mean,
          beta0_var = hierarchical_model$scale$beta0_sd^2
        )
        + .5 * sum(
            new_field[vecchia_approx$locs_match] *
              (state$stuff$lm_residuals - new_field[vecchia_approx$locs_match]) /
              state$stuff$noise_var
          ) # observation ll
      ),
      # range
      -(
        -betaPriorLogDensDerivative(
          beta = new_range_beta, n_PP = hierarchical_model$range$PP$n_knots,
          beta0_mean = hierarchical_model$range$beta0_mean,
          beta0_var = hierarchical_model$range$beta0_sd^2,
          state$params$range_log_scale
        ) # normal prior
        + xPPCrossprod(
            X = range_X$X_locs, vecchia_approx = vecchia_approx,
            permutate_PP_to_obs = FALSE,
            PP = hierarchical_model$range$PP,
            Y = # Jacobian of range field wrt range_beta
              t(
                # natural gradient of obs likelihood wrt range field
                derivativeSandwiches_(
                  vecchia = state$stuff$proposed_compressed_chol,
                  left_vector = as.vector(
                    Matrix::solve(
                      Matrix::t(state$stuff$proposed_sparse_chol),
                      -as.vector(vecchia_approx$locs_match_matrix %*% # gradient of  Gaussian observations ll wrt latent field
                        ((new_field[vecchia_approx$locs_match] - state$stuff$lm_residuals) / state$stuff$noise_var))
                      * exp(.5 * new_field_log_var[1, 1]) # part of sparse chol
                    )
                  ),
                  right_vector = new_field / exp(.5 * new_field_log_var[1, 1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                  NNarray = vecchia_approx$NNarray,
                  sauce_determinant_chef = FALSE,
                  num_threads = num_threads
                )
              )
          ) %*% range_reparam_mat
      )
    )

    # forward proposal density
    q_proposal <- -.5 * sum(c(
      state$momenta$field_log_var_ancillary_grouped,
      state$momenta$range_beta_ancillary
    )^2)
    # computing backward proposal
    innov_back <- solve(
      sqrt(stepsize) * cond_mat,
      c(state$params$field_log_var, state$params$range_beta)
      - c(new_field_log_var, new_range_beta)
        - as.vector(cond_mat %*% (crossprod(cond_mat, c(dens_grad_back)))) * stepsize / 2
    )
    # backward proposal density
    q_back <- -.5 * sum(innov_back^2)

    # current negated log density
    current_U <-
      (
        # field log var prior
        -betaPriorLogDens(
          beta = as.matrix(state$params$field_log_var), n_PP = 0, log_scale = 0,
          beta0_mean = hierarchical_model$scale$beta0_mean,
          beta0_var = hierarchical_model$scale$beta0_sd^2
        )
        # field range prior
        - betaPriorLogDens(
            beta = state$params$range_beta,
            n_PP = hierarchical_model$range$PP$n_knots,
            beta0_mean = hierarchical_model$range$beta0_mean,
            beta0_var = hierarchical_model$range$beta0_sd^2,
            state$params$range_log_scale
          )
          # data ll
          + .5 * sum((state$stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2 / state$stuff$noise_var) # observation ll
      )
    # proposed negated log density
    proposed_U <-
      (
        # field log var prior
        -betaPriorLogDens(
          beta = as.matrix(new_field_log_var), n_PP = 0, log_scale = 0,
          beta0_mean = hierarchical_model$scale$beta0_mean,
          beta0_var = hierarchical_model$scale$beta0_sd^2
        )
        # field range prior
        - betaPriorLogDens(
            beta = new_range_beta,
            n_PP = hierarchical_model$range$PP$n_knots,
            beta0_mean = hierarchical_model$range$beta0_mean,
            beta0_var = hierarchical_model$range$beta0_sd^2,
            state$params$range_log_scale
          )
          # data ll
          + .5 * sum((state$stuff$lm_residuals - new_field[vecchia_approx$locs_match])^2 / state$stuff$noise_var) # observation ll
      )
    
    #par(mfrow = c(2,2))
    #if(iter %/% 1 == iter / 1){
    #  empirical_grad = c()
    #  sparse_chol_grad_test = state$stuff$sparse_chol
    #  compressed_chol_grad_test = state$stuff$compressed_chol
    #  for(i in seq(length(state$params$range_beta)+1)){
    #    range_beta_grad_test <- state$params$range_beta 
    #    field_log_var_grad_test <- state$params$field_log_var 
    #    if(i==1) field_log_var_grad_test = field_log_var_grad_test + 1e-6
    #    if(i>1) range_beta_grad_test[i-1] = range_beta_grad_test[i-1] + 1e-6
    #    
    #    # re-computing Vecchia approx
    #    log_range <- computeLogRange(
    #      range_beta = range_beta_grad_test,
    #      PP = hierarchical_model$range$PP,
    #      vecchia_approx = vecchia_approx,
    #      range_X = range_X
    #    )
    #    vecchia_(start_idx = 1,
    #             log_range = t(log_range), locs = vecchia_approx$t_locs,
    #             NNarray = vecchia_approx$NNarray,
    #             smoothness = hierarchical_model$matern_smoothness,
    #             compute_derivative = TRUE, num_threads = num_threads, result = compressed_chol_grad_test
    #    )
    #    sparse_chol_grad_test@x <- compressed_chol_grad_test[,,1][vecchia_approx$sparse_chol_x_reorder]
    #    field_grad_test <- 
    #      as.vector(Matrix::solve(sparse_chol_grad_test, state$stuff$sparse_chol %*% (state$params$field))) * 
    #      exp(.5 * (field_log_var_grad_test[1,1] - state$params$field_log_var[1, 1]))
    #    
    #    # proposed negated log density
    #    test_grad_U <-
    #      (
    #        # field log var prior
    #        - betaPriorLogDens(
    #          beta = as.matrix(field_log_var_grad_test), n_PP = 0, log_scale = 0,
    #          beta0_mean = hierarchical_model$scale$beta0_mean,
    #          beta0_var = hierarchical_model$scale$beta0_sd^2)
    #        # field range prior
    #        -betaPriorLogDens(
    #          beta = range_beta_grad_test,
    #          n_PP = hierarchical_model$range$PP$n_knots,
    #          beta0_mean = hierarchical_model$range$beta0_mean,
    #          beta0_var = hierarchical_model$range$beta0_sd^2,
    #          state$params$range_log_scale
    #        ) 
    #        # data ll
    #        + .5 * sum((state$stuff$lm_residuals - field_grad_test[vecchia_approx$locs_match])^2 / state$stuff$noise_var) # observation ll
    #      )
    #    empirical_grad = c(empirical_grad, 
    #                       -(test_grad_U - current_U)*1e6
    #    )
    #  }
    #  plot(dens_grad, empirical_grad,  main  = paste("ancillary iter", iter))
    #  abline(a=0, b=1)
    #  abline(h=0)
    #  abline(v=0)
    #}
    #if((iter %/% 1) == (iter/ 1)
    #){
    #  empirical_grad = c()
    #  sparse_chol_grad_test = state$stuff$proposed_sparse_chol
    #  compressed_chol_grad_test = state$stuff$proposed_compressed_chol
    #  for(i in seq(length(state$params$range_beta)+1)){
    #    range_beta_grad_test <- new_range_beta 
    #    field_log_var_grad_test <- new_field_log_var 
    #    if(i==1) field_log_var_grad_test = field_log_var_grad_test + 1e-6
    #    if(i>1) range_beta_grad_test[i-1] = range_beta_grad_test[i-1] + 1e-6
    #    
    #    # re-computing Vecchia approx
    #    log_range <- computeLogRange(
    #      range_beta = range_beta_grad_test,
    #      PP = hierarchical_model$range$PP,
    #      vecchia_approx = vecchia_approx,
    #      range_X = range_X
    #    )
    #    vecchia_(start_idx = 1,
    #             log_range = t(log_range), locs = vecchia_approx$t_locs,
    #             NNarray = vecchia_approx$NNarray,
    #             smoothness = hierarchical_model$matern_smoothness,
    #             compute_derivative = TRUE, num_threads = num_threads, result = compressed_chol_grad_test
    #    )
    #    sparse_chol_grad_test@x <- compressed_chol_grad_test[,,1][vecchia_approx$sparse_chol_x_reorder]
    #    field_grad_test <- 
    #      as.vector(Matrix::solve(sparse_chol_grad_test, state$stuff$proposed_sparse_chol %*% (new_field))) * 
    #      exp(.5 * (field_log_var_grad_test[1,1] - new_field_log_var[1, 1]))
    #    
    #    # proposed negated log density
    #    test_grad_U <-
    #      (
    #        # field log var prior
    #        - betaPriorLogDens(
    #          beta = as.matrix(field_log_var_grad_test), n_PP = 0, log_scale = 0,
    #          beta0_mean = hierarchical_model$scale$beta0_mean,
    #          beta0_var = hierarchical_model$scale$beta0_sd^2)
    #        # field range prior
    #        -betaPriorLogDens(
    #          beta = range_beta_grad_test,
    #          n_PP = hierarchical_model$range$PP$n_knots,
    #          beta0_mean = hierarchical_model$range$beta0_mean,
    #          beta0_var = hierarchical_model$range$beta0_sd^2,
    #          state$params$range_log_scale
    #        ) 
    #        # data ll
    #        + .5 * sum((state$stuff$lm_residuals - field_grad_test[vecchia_approx$locs_match])^2 / state$stuff$noise_var) # observation ll
    #      )
    #    empirical_grad = c(empirical_grad, 
    #                       -(test_grad_U - proposed_U)*1e6
    #    )
    #  }
    #  plot(dens_grad_back, empirical_grad, main  = paste("ancillary iter", iter))
    #  abline(a=0, b=1)
    #  abline(h=0)
    #  abline(v=0)
    #}
    
    # updating kernel down
    state$ker_var$range_beta_ancillary[1] <- updateKernel(
      iter = iter, iter_start = iter_start,
      kernel_value = state$ker_var$range_beta_ancillary[1], mult = -2
    )

    if (!is.nan(current_U - proposed_U - q_proposal + q_back)) {
      if (log(runif(1)) < current_U - proposed_U - q_proposal + q_back) {
        # updating kernel up
        state$ker_var$range_beta_ancillary[1] <- updateKernel(
          iter = iter, iter_start = iter_start,
          kernel_value = state$ker_var$range_beta_ancillary[1], mult = 3
        )
        # replacing current momentum
        state$momenta$range_beta_ancillary <- innov_back[-1]
        state$momenta$field_log_var_ancillary_grouped <- innov_back[1]
        # replacing current gradient
        dens_grad <- dens_grad_back
        # replacing state$stuff
        state$params$field <- new_field
        # invert positions of proposed and current matrix to avoid recreating objects.
        stuff_compressed <- match(c("compressed_chol", "proposed_compressed_chol"), names(state$stuff))
        stuff_compressed <- stuff_compressed[!is.na(stuff_compressed)]
        names(state$stuff)[stuff_compressed] <- names(state$stuff)[rev(stuff_compressed)]
        stuff_sparse <- match(c("sparse_chol", "proposed_sparse_chol"), names(state$stuff))
        stuff_sparse <- stuff_sparse[!is.na(stuff_sparse)]
        names(state$stuff)[stuff_sparse] <- names(state$stuff)[rev(stuff_sparse)]
        state$params$range_beta[] <- new_range_beta
        state$params$field_log_var[] <- new_field_log_var
        # print("ANCILLAAAAARRRYYYYYY")
      }
    }
    # negating momentum
    state$momenta$range_beta_ancillary <- -state$momenta$range_beta_ancillary
    state$momenta$field_log_var_ancillary_grouped <- -state$momenta$field_log_var_ancillary_grouped
  }
  ###########################
  # Range beta (sufficient) #
  ###########################
  dens_grad <- c(
    # log var
    betaPriorLogDensDerivative(
      beta = as.matrix(state$params$field_log_var), n_PP = 0, log_scale = 0,
      beta0_mean = hierarchical_model$scale$beta0_mean,
      beta0_var = hierarchical_model$scale$beta0_sd^2
    )
    + .5 * sum((state$stuff$sparse_chol %*% state$params$field)^2) / exp(state$params$field_log_var[1, 1]) # observation ll
      - .5 * vecchia_approx$n_locs,
    # range
    +betaPriorLogDensDerivative(
      beta = state$params$range_beta,
      n_PP = hierarchical_model$range$PP$n_knots,
      beta0_mean = hierarchical_model$range$beta0_mean,
      beta0_var = hierarchical_model$range$beta0_sd^2,
      log_scale = state$params$range_log_scale
    ) # normal prior
    # normal prior derivative
    - xPPCrossprod(
      X = range_X$X_locs, PP = hierarchical_model$range$PP,
      permutate_PP_to_obs = FALSE,
      vecchia_approx = vecchia_approx,
      Y = # Jacobian of range field wrt range_beta
        t( # natural gradient of obs likelihood wrt range field
          derivativeSandwiches_(
            vecchia = state$stuff$compressed_chol, # derivative of the (unscaled) NNGP factor
            left_vector = as.vector(state$stuff$sparse_chol %*% (state$params$field / exp(.5 * state$params$field_log_var[1, 1]))), # left vector = whitened latent field
            right_vector = state$params$field / exp(.5 * state$params$field_log_var[1, 1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
            NNarray = vecchia_approx$NNarray,
            sauce_determinant_chef = TRUE,
            num_threads = num_threads
          )
        )
    ) %*% range_reparam_mat
  )
  
  #if(iter %/% 1 == iter /1
  #   ){
  #empirical_grad = c()
  #sparse_chol_grad_test = state$stuff$proposed_sparse_chol
  #compressed_chol_grad_test = state$stuff$proposed_compressed_chol
  #for(i in seq(length(state$params$range_beta)+1)){
  #  range_beta_grad_test <- state$params$range_beta 
  #  field_log_var_grad_test <- state$params$field_log_var 
  #  if(i==1) field_log_var_grad_test = field_log_var_grad_test + 1e-6
  #  if(i>1) range_beta_grad_test[i-1] = range_beta_grad_test[i-1] + 1e-6
  #  
  #  # re-computing Vecchia approx
  #  log_range <- computeLogRange(
  #    range_beta = range_beta_grad_test,
  #    PP = hierarchical_model$range$PP,
  #    vecchia_approx = vecchia_approx,
  #    range_X = range_X
  #  )
  #  vecchia_(start_idx = 1,
  #    log_range = t(log_range), locs = vecchia_approx$t_locs,
  #    NNarray = vecchia_approx$NNarray,
  #    smoothness = hierarchical_model$matern_smoothness,
  #    compute_derivative = TRUE, num_threads = num_threads, result = compressed_chol_grad_test
  #  )
  #  sparse_chol_grad_test@x <- compressed_chol_grad_test[,,1][vecchia_approx$sparse_chol_x_reorder]
  #  
  #  
  #  current_U <-       (
  #    # field log var prior
  #    -betaPriorLogDens(
  #      beta = as.matrix(state$params$field_log_var), n_PP = 0, log_scale = 0,
  #      beta0_mean = hierarchical_model$scale$beta0_mean,
  #      beta0_var = hierarchical_model$scale$beta0_sd^2)
  #    # range
  #    -betaPriorLogDens(
  #      beta = state$params$range_beta,
  #      n_PP = hierarchical_model$range$PP$n_knots,
  #      beta0_mean = hierarchical_model$range$beta0_mean,
  #      beta0_var = hierarchical_model$range$beta0_sd^2,
  #      log_scale = state$params$range_log_scale
  #    )
  #    # data ll
  #    + .5 * sum((state$stuff$sparse_chol %*% (state$params$field / exp(.5 * state$params$field_log_var[1, 1])))^2)
  #    - sum(log(Matrix::diag(state$stuff$sparse_chol)))
  #    + vecchia_approx$n_locs * (.5 * state$params$field_log_var[1, 1])
  #  )
  #  proposed_U <-       (
  #    # field log var prior
  #    -betaPriorLogDens(
  #      beta = as.matrix(field_log_var_grad_test), n_PP = 0, log_scale = 0,
  #      beta0_mean = hierarchical_model$scale$beta0_mean,
  #      beta0_var = hierarchical_model$scale$beta0_sd^2)
  #    # range
  #    -betaPriorLogDens(
  #      beta = range_beta_grad_test,
  #      n_PP = hierarchical_model$range$PP$n_knots,
  #      beta0_mean = hierarchical_model$range$beta0_mean,
  #      beta0_var = hierarchical_model$range$beta0_sd^2,
  #      state$params$range_log_scale
  #    )
  #    # data ll
  #    + .5 * sum((sparse_chol_grad_test %*% (state$params$field / exp(.5 * field_log_var_grad_test[1,1])))^2)
  #    - sum(log(Matrix::diag(sparse_chol_grad_test)))
  #    + vecchia_approx$n_locs * (.5 * field_log_var_grad_test[1,1])
  #  )
  #  empirical_grad = c(empirical_grad, 
  #                     -(proposed_U - current_U)*1e6
  #  )
  #}
  #plot(dens_grad, empirical_grad,  main  = paste("sufficient iter", iter))
  #abline(a=0, b=1)
  #abline(h=0)
  #abline(v=0)
  #}
  #cond_mat <- solve(state$stuff$range_beta_conditioning_s)
  cond_mat <- solve(
    as(
      (1 - min(1, max((500 - (iter + iter_start))/200, .05))) *
        state$stuff$range_beta_conditioning_s / sqrt(sum(state$stuff$range_beta_conditioning_s^2)) 
      + Matrix::bdiag(
        1, 
        min(1, max((500 - (iter + iter_start))/200, .05)) *
          (diag(rep(1, 1+ 2*hierarchical_model$anisotropic)) %x% as.matrix(range_X$crossprod_X)) / 
          sqrt(3*sum(as.matrix(range_X$crossprod_X)^2)))
      , 
      "sparseMatrix"
    )
  )
  cond_mat <- cond_mat / sqrt(sum(cond_mat^2))
  cond_mat <- t(chol(cond_mat))
  
  for(asdf in seq(1)){
    if(iter + iter_start >200){
      state$stuff$range_beta_conditioning_s <- 
        ((iter+iter_start)/(iter+iter_start-1)) * state$stuff$range_beta_conditioning_s + 
        tcrossprod(c(dens_grad)) / (iter+iter_start)
    }    
    stepsize <- exp(state$ker_var$range_beta_sufficient[1])
    state$momenta$range_beta_sufficient = renewMomentum(state$momenta$range_beta_sufficient)
    state$momenta$field_log_var_sufficient_grouped = renewMomentum(state$momenta$field_log_var_sufficient_grouped)
    move_forward <-
      as.vector(stepsize * (cond_mat %*% (crossprod(cond_mat, c(dens_grad))))/2) + 
      sqrt(stepsize) * as.vector(cond_mat %*% c( state$momenta$field_log_var_sufficient_grouped, state$momenta$range_beta_sufficient))
    new_range_beta <- state$params$range_beta + move_forward[-1]
    new_field_log_var <- state$params$field_log_var + move_forward[1]
    
    log_range <- computeLogRange(
      range_beta = new_range_beta,
      PP = hierarchical_model$range$PP,
      vecchia_approx = vecchia_approx,
      range_X = range_X
    )
    vecchia_(start_idx = 1,
      log_range = t(log_range), locs = vecchia_approx$t_locs,
      NNarray = vecchia_approx$NNarray,
      smoothness = hierarchical_model$matern_smoothness,
      compute_derivative = TRUE, num_threads = num_threads, result = state$stuff$proposed_compressed_chol
    )
    state$stuff$proposed_sparse_chol@x <- state$stuff$proposed_compressed_chol[,,1][vecchia_approx$sparse_chol_x_reorder]
    
    dens_grad_back <- c(
      # log var
       betaPriorLogDensDerivative(
         beta = as.matrix(new_field_log_var), n_PP = 0, log_scale = 0,
         beta0_mean = hierarchical_model$scale$beta0_mean,
         beta0_var = hierarchical_model$scale$beta0_sd^2
       )
      + .5 * sum((state$stuff$proposed_sparse_chol %*% state$params$field)^2) / exp(new_field_log_var[1, 1]) # observation ll
      - .5 * vecchia_approx$n_locs, 
      # range
      betaPriorLogDensDerivative(
        beta = new_range_beta,
        n_PP = hierarchical_model$range$PP$n_knots,
        beta0_mean = hierarchical_model$range$beta0_mean,
        beta0_var = hierarchical_model$range$beta0_sd^2,
        log_scale = state$params$range_log_scale) 
      - xPPCrossprod(
        X = range_X$X_locs, PP = hierarchical_model$range$PP,
        permutate_PP_to_obs = FALSE,
        vecchia_approx = vecchia_approx,
        Y = # Jacobian of range field wrt range_beta
          t( # natural gradient of obs likelihood wrt range field
            derivativeSandwiches_(
              vecchia = state$stuff$compressed_chol, # derivative of the (unscaled) NNGP factor
              left_vector = as.vector(state$stuff$sparse_chol %*% (state$params$field / exp(.5 * state$params$field_log_var[1, 1]))), # left vector = whitened latent field
              right_vector = state$params$field / exp(.5 * state$params$field_log_var[1, 1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
              NNarray = vecchia_approx$NNarray,
              sauce_determinant_chef = TRUE,
              num_threads = num_threads
            )
          )
      ) %*% range_reparam_mat
  )
  
    #if(iter %/% 1 == iter / 1
    #   ){
    # empirical_grad = c()
    # sparse_chol_grad_test = state$stuff$proposed_sparse_chol
    # compressed_chol_grad_test = state$stuff$proposed_compressed_chol
    # for(i in seq(length(state$params$range_beta)+1)){
    #   range_beta_grad_test <- new_range_beta 
    #   field_log_var_grad_test <- new_field_log_var 
    #   if(i==1) field_log_var_grad_test = field_log_var_grad_test + 1e-6
    #   if(i>1) range_beta_grad_test[i-1] = range_beta_grad_test[i-1] + 1e-6
    #   
    #   # re-computing Vecchia approx
    #   log_range <- computeLogRange(
    #     range_beta = range_beta_grad_test,
    #     PP = hierarchical_model$range$PP,
    #     vecchia_approx = vecchia_approx,
    #     range_X = range_X
    #   )
    #   vecchia_(start_idx = 1,
    #     log_range = t(log_range), locs = vecchia_approx$t_locs,
    #     NNarray = vecchia_approx$NNarray,
    #     smoothness = hierarchical_model$matern_smoothness,
    #     compute_derivative = TRUE, num_threads = num_threads, result = compressed_chol_grad_test
    #   )
    #   sparse_chol_grad_test@x <- compressed_chol_grad_test[,,1][vecchia_approx$sparse_chol_x_reorder]
    #   
    #   
    #   current_U <-       (
    #     # field log var prior
    #     -betaPriorLogDens(
    #       beta = as.matrix(new_field_log_var), n_PP = 0, log_scale = 0,
    #       beta0_mean = hierarchical_model$scale$beta0_mean,
    #       beta0_var = hierarchical_model$scale$beta0_sd^2)
    #     # range
    #     -betaPriorLogDens(
    #       beta = new_range_beta,
    #       n_PP = hierarchical_model$range$PP$n_knots,
    #       beta0_mean = hierarchical_model$range$beta0_mean,
    #       beta0_var = hierarchical_model$range$beta0_sd^2,
    #       log_scale = state$params$range_log_scale
    #     )
    #     # data ll
    #     + .5 * sum((state$stuff$proposed_sparse_chol %*% (state$params$field / exp(.5 * new_field_log_var[1, 1])))^2)
    #     - sum(log(Matrix::diag(state$stuff$proposed_sparse_chol)))
    #     + vecchia_approx$n_locs * (.5 * new_field_log_var[1, 1])
    #   )
    #   proposed_U <-       (
    #     # field log var prior
    #     -betaPriorLogDens(
    #       beta = as.matrix(field_log_var_grad_test), n_PP = 0, log_scale = 0,
    #       beta0_mean = hierarchical_model$scale$beta0_mean,
    #       beta0_var = hierarchical_model$scale$beta0_sd^2)
    #     # range
    #     -betaPriorLogDens(
    #       beta = range_beta_grad_test,
    #       n_PP = hierarchical_model$range$PP$n_knots,
    #       beta0_mean = hierarchical_model$range$beta0_mean,
    #       beta0_var = hierarchical_model$range$beta0_sd^2,
    #       state$params$range_log_scale
    #     )
    #     # data ll
    #     + .5 * sum((sparse_chol_grad_test %*% (state$params$field / exp(.5 * field_log_var_grad_test[1,1])))^2)
    #     - sum(log(Matrix::diag(sparse_chol_grad_test)))
    #     + vecchia_approx$n_locs * (.5 * field_log_var_grad_test[1,1])
    #   )
    #   empirical_grad = c(empirical_grad, 
    #                      -(proposed_U - current_U)*1e6
    #   )
    # }
    # plot(dens_grad_back, empirical_grad,  main  = paste("sufficient iter", iter))
    # abline(a=0, b=1)
    # abline(h=0)
    # abline(v=0)
    #}
     
    q_proposal <- -.5 * sum(c(
      state$momenta$field_log_var_sufficient_grouped, 
      state$momenta$range_beta_sufficient)^2)
    # computing backward proposal
    innov_back <- solve(
      sqrt(stepsize) * cond_mat,
      c(state$params$field_log_var, state$params$range_beta)
      - c(new_field_log_var, new_range_beta)
      - as.vector(cond_mat %*% (crossprod(cond_mat, c(dens_grad_back)))) * stepsize / 2
    )
    # backward proposal density
    q_back <- -.5 * sum(innov_back^2)
    vecchia_(
      start_idx = 1,
      log_range = t(log_range), locs = vecchia_approx$t_locs,
      NNarray = vecchia_approx$NNarray,
      smoothness = hierarchical_model$matern_smoothness,
      compute_derivative = TRUE, num_threads = num_threads, result = state$stuff$proposed_compressed_chol
    )
    state$stuff$proposed_sparse_chol@x <- state$stuff$proposed_compressed_chol[, , 1][vecchia_approx$sparse_chol_x_reorder]

    dens_grad_back <- c(
      # log var
      betaPriorLogDensDerivative(
        beta = as.matrix(new_field_log_var), n_PP = 0, log_scale = 0,
        beta0_mean = hierarchical_model$scale$beta0_mean,
        beta0_var = hierarchical_model$scale$beta0_sd^2
      )
      + .5 * sum((state$stuff$proposed_sparse_chol %*% state$params$field)^2) / exp(new_field_log_var[1, 1]) # observation ll
        - .5 * vecchia_approx$n_locs,
      # range
      betaPriorLogDensDerivative(
        beta = new_range_beta,
        n_PP = hierarchical_model$range$PP$n_knots,
        beta0_mean = hierarchical_model$range$beta0_mean,
        beta0_var = hierarchical_model$range$beta0_sd^2,
        log_scale = state$params$range_log_scale
      )
      - xPPCrossprod(
          X = range_X$X_locs, PP = hierarchical_model$range$PP,
          permutate_PP_to_obs = FALSE,
          vecchia_approx = vecchia_approx,
          Y = # Jacobian of range field wrt range_beta
            t( # natural gradient of obs likelihood wrt range field
              derivativeSandwiches_(
                vecchia = state$stuff$proposed_compressed_chol, # derivative of the (unscaled) NNGP factor
                left_vector = as.vector(state$stuff$proposed_sparse_chol %*% (state$params$field / exp(.5 * new_field_log_var[1, 1]))), # left vector = whitened latent field
                right_vector = state$params$field / exp(.5 * new_field_log_var[1, 1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                NNarray = vecchia_approx$NNarray,
                sauce_determinant_chef = TRUE,
                num_threads = num_threads
              )
            )
        ) %*% range_reparam_mat
    )
  
    q_proposal <- -.5 * sum(c(
      state$momenta$field_log_var_sufficient_grouped,
      state$momenta$range_beta_sufficient
    )^2)
    innov_back <- solve(
      sqrt(stepsize) * cond_mat,
      c(state$params$field_log_var, state$params$range_beta)
      - c(new_field_log_var, new_range_beta)
        - as.vector(cond_mat %*% (crossprod(cond_mat, c(dens_grad_back)))) * stepsize / 2
    )
    q_back <- -.5 * sum(innov_back^2)
  
    current_U <- (
      # field log var prior
      -betaPriorLogDens(
        beta = as.matrix(state$params$field_log_var), n_PP = 0, log_scale = 0,
        beta0_mean = hierarchical_model$scale$beta0_mean,
        beta0_var = hierarchical_model$scale$beta0_sd^2
      )
      # range
      - betaPriorLogDens(
          beta = state$params$range_beta,
          n_PP = hierarchical_model$range$PP$n_knots,
          beta0_mean = hierarchical_model$range$beta0_mean,
          beta0_var = hierarchical_model$range$beta0_sd^2,
          log_scale = state$params$range_log_scale
        )
        # data ll
        + .5 * sum((state$stuff$sparse_chol %*% (state$params$field / exp(.5 * state$params$field_log_var[1, 1])))^2)
        - sum(log(Matrix::diag(state$stuff$sparse_chol)))
        + vecchia_approx$n_locs * (.5 * state$params$field_log_var[1, 1])
    )
    proposed_U <- (
      # field log var prior
      -betaPriorLogDens(
        beta = as.matrix(new_field_log_var), n_PP = 0, log_scale = 0,
        beta0_mean = hierarchical_model$scale$beta0_mean,
        beta0_var = hierarchical_model$scale$beta0_sd^2
      )
      # range
      - betaPriorLogDens(
          beta = new_range_beta,
          n_PP = hierarchical_model$range$PP$n_knots,
          beta0_mean = hierarchical_model$range$beta0_mean,
          beta0_var = hierarchical_model$range$beta0_sd^2,
          state$params$range_log_scale
        )
        # data ll
        + .5 * sum((state$stuff$proposed_sparse_chol %*% (state$params$field / exp(.5 * new_field_log_var[1, 1])))^2)
        - sum(log(Matrix::diag(state$stuff$proposed_sparse_chol)))
        + vecchia_approx$n_locs * (.5 * new_field_log_var[1, 1])
    )
  
    state$ker_var$range_beta_sufficient[1] <- updateKernel(
      iter = iter, iter_start = iter_start,
      kernel_value = state$ker_var$range_beta_sufficient[1], mult = -2
    )
    if (!is.nan(current_U - proposed_U - q_proposal + q_back)) {
      if (log(runif(1)) < current_U - proposed_U - q_proposal + q_back) {
        state$ker_var$range_beta_sufficient[1] <- updateKernel(
          iter = iter, iter_start = iter_start,
          kernel_value = state$ker_var$range_beta_sufficient[1], mult = 3
        )
        state$momenta$range_beta_sufficient <- innov_back[-1]
        state$momenta$field_log_var_sufficient_grouped <- innov_back[1]
        dens_grad <- dens_grad_back
        # invert positions of proposed and current matrix to avoid recreating objects.
        stuff_compressed <- match(c("compressed_chol", "proposed_compressed_chol"), names(state$stuff))
        stuff_compressed <- stuff_compressed[!is.na(stuff_compressed)]
        names(state$stuff)[stuff_compressed] <- names(state$stuff)[rev(stuff_compressed)]
        stuff_sparse <- match(c("sparse_chol", "proposed_sparse_chol"), names(state$stuff))
        stuff_sparse <- stuff_sparse[!is.na(stuff_sparse)]
        names(state$stuff)[stuff_sparse] <- names(state$stuff)[rev(stuff_sparse)]
        state$params$range_beta[] <- new_range_beta
        state$params$field_log_var[] <- new_field_log_var
        #print("SUFFICIEEEENT")
      }
    }
    # negating momentum
    state$momenta$range_beta_sufficient <- -state$momenta$range_beta_sufficient
    state$momenta$field_log_var_sufficient_grouped <- -state$momenta$field_log_var_sufficient_grouped
  }


  return(list(state = state))
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

updateNoiseBeta <- function(state, noise, noise_X, vecchia_approx, iter, iter_start) {
  # VEWY IMPOWTANT don't remove or comment
  squared_residuals <- as.matrix(state$stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2
  dens_grad <- (
    -betaPriorLogDensDerivative(
      beta = state$params$noise_beta, n_PP = noise$PP$n_knots,
      beta0_mean = noise$beta0_mean,
      beta0_var = noise$beta0_sd^2,
      log_scale = state$params$noise_log_scale
    ) # normal prior
    + xPPCrossprod(
        X = noise_X$X, noise$PP, vecchia_approx = vecchia_approx, permutate_PP_to_obs = TRUE,
        Y =
          (
            +.5 # determinant part of normal likelihood
            - (squared_residuals / state$stuff$noise_var) / 2 # exponential part of normal likelihood
          )
      )
  )
  for (asdf in seq(25)) {
    # HMC update
    q <- noise_X$L_minus_one %*% state$params$noise_beta
    state$momenta$noise_beta <- renewMomentum(state$momenta$noise_beta)
    p <- state$momenta$noise_beta
    # Make a half step for momentum at the beginning
    exp_noise_mala <- exp(state$ker_var$noise_beta_mala)
    p <- p - exp_noise_mala * crossprod(noise_X$L, dens_grad) / 2

    n_hmc_steps <- 1
    for (hmc_step in seq_len(n_hmc_steps)) {
      # Make a full step for the position
      q <- q + exp_noise_mala * p
      new_noise_beta <- noise_X$L %*% q
      new_noise_var <- as.vector(exp(xPPMultRight(
        X = noise_X$X, PP = noise$PP,
        vecchia_approx = vecchia_approx, Y = new_noise_beta,
        permutate_PP_to_obs = TRUE
      )))
      # Make a half step for momentum at the end
      new_dens_grad <- (
        -betaPriorLogDensDerivative(
          beta = new_noise_beta, n_PP = noise$PP$n_knots,
          beta0_mean = noise$beta0_mean,
          beta0_var = noise$beta0_sd^2,
          log_scale = state$params$noise_log_scale
        ) # normal prior
        + xPPCrossprod(
            X = noise_X$X,
            noise$PP,
            vecchia_approx = vecchia_approx,
            permutate_PP_to_obs = TRUE,
            Y = (+.5 # determinant part of normal likelihood
            - (squared_residuals / new_noise_var) / 2 # exponential part of normal likelihood
            )
          )
      )
      p <- p - exp_noise_mala * crossprod(noise_X$L, new_dens_grad) / (1 + (hmc_step == n_hmc_steps))
    }

    # Evaluate potential and kinetic energies at start and end of trajectory
    current_U <- (
      -betaPriorLogDens(
        beta = state$params$noise_beta,
        n_PP = noise$PP$n_knots,
        beta0_mean = noise$beta0_mean,
        beta0_var = noise$beta0_sd^2,
        log_scale = state$params$noise_log_scale
      ) # normal prior
      + .5 * sum(log(state$stuff$noise_var)) # det
        + .5 * sum(squared_residuals / state$stuff$noise_var) # observations
    )
    current_K <- sum(state$momenta$noise_beta^2) / 2
    proposed_U <- (
      -betaPriorLogDens(
        beta = new_noise_beta,
        n_PP = noise$PP$n_knots,
        beta0_mean = noise$beta0_mean,
        beta0_var = noise$beta0_sd^2,
        log_scale = state$params$noise_log_scale
      ) # normal prior
      + .5 * sum(log(new_noise_var)) # det
        + .5 * sum(squared_residuals / new_noise_var) # observations
    )
    proposed_K <- sum(p^2) / 2

    state$ker_var$noise_beta_mala <- updateKernel(
      iter = iter,
      iter_start = iter_start,
      kernel_value = state$ker_var$noise_beta_mala,
      mult = -.8
    )
    if (!is.nan(current_U - proposed_U + current_K - proposed_K)) {
      if (log(runif(1)) < current_U - proposed_U + current_K - proposed_K) {
        state$ker_var$noise_beta_mala <- updateKernel(
          iter = iter,
          iter_start = iter_start,
          kernel_value = state$ker_var$noise_beta_mala,
          mult = 1
        )
        dens_grad <- new_dens_grad
        state$momenta$noise_beta <- -p
        state$params$noise_beta[] <- new_noise_beta
        state$stuff$noise_var <- new_noise_var
      }
    }
    # negating momentum
    state$momenta$noise_beta <- -state$momenta$noise_beta
  }

  return(list("state" = state, "squared_residuals" = squared_residuals))
}


updateNoiseBetaFisher <- function(state, noise, noise_X, vecchia_approx, iter, iter_start) {
  # VEWY IMPOWTANT don't remove or comment
  squared_residuals <- as.matrix(state$stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2
  dens_grad <- (
    -betaPriorLogDensDerivative(
      beta = state$params$noise_beta, n_PP = noise$PP$n_knots,
      beta0_mean = noise$beta0_mean,
      beta0_var = noise$beta0_sd^2,
      log_scale = state$params$noise_log_scale
    ) # normal prior
    + xPPCrossprod(
        X = noise_X$X, noise$PP, vecchia_approx = vecchia_approx, permutate_PP_to_obs = TRUE,
        Y =
          (
            +.5 # determinant part of normal likelihood
            - (squared_residuals / state$stuff$noise_var) / 2 # exponential part of normal likelihood
          )
      )
  )

  cond_mat <- solve(
    as(
      (1 - min(1, max((500 - (iter + iter_start)) / 200, .05))) *
        state$stuff$noise_beta_conditioning / sqrt(sum(state$stuff$noise_beta_conditioning^2))
        +
        min(1, max((500 - (iter + iter_start)) / 200, .05)) *
          as.matrix(noise_X$crossprod_X) /
          sqrt(3 * sum(as.matrix(noise_X$crossprod_X)^2)),
      "sparseMatrix"
    )
  )
  cond_mat <- cond_mat / sqrt(sum(cond_mat^2))
  cond_mat <- t(chol(cond_mat))

  # for conditioning matrix update
  grad_record_for_Fisher <- list()

  n_mala <- 30
  for (asdf in seq(n_mala)) {
    # for conditioning matrix update
    grad_record_for_Fisher[[length(grad_record_for_Fisher) + 1]] <- as.vector(dens_grad)

    # HMC update
    q <- solve(cond_mat, state$params$noise_beta)
    state$momenta$noise_beta <- renewMomentum(state$momenta$noise_beta)
    p <- state$momenta$noise_beta
    # Make a half step for momentum at the beginning
    exp_noise_mala <- exp(state$ker_var$noise_beta_mala)
    p <- p - exp_noise_mala * crossprod(cond_mat, dens_grad) / 2

    n_hmc_steps <- 1
    for (hmc_step in seq_len(n_hmc_steps)) {
      # Make a full step for the position
      q <- q + exp_noise_mala * p
      new_noise_beta <- cond_mat %*% q
      new_noise_var <- as.vector(exp(xPPMultRight(
        X = noise_X$X, PP = noise$PP,
        vecchia_approx = vecchia_approx, Y = new_noise_beta,
        permutate_PP_to_obs = TRUE
      )))
      # Make a half step for momentum at the end
      new_dens_grad <- (
        -betaPriorLogDensDerivative(
          beta = new_noise_beta, n_PP = noise$PP$n_knots,
          beta0_mean = noise$beta0_mean,
          beta0_var = noise$beta0_sd^2,
          log_scale = state$params$noise_log_scale
        ) # normal prior
        + xPPCrossprod(
            X = noise_X$X,
            noise$PP,
            vecchia_approx = vecchia_approx,
            permutate_PP_to_obs = TRUE,
            Y = (+.5 # determinant part of normal likelihood
            - (squared_residuals / new_noise_var) / 2 # exponential part of normal likelihood
            )
          )
      )
      p <- p - exp_noise_mala * crossprod(cond_mat, new_dens_grad) / (1 + (hmc_step == n_hmc_steps))
    }

    # Evaluate potential and kinetic energies at start and end of trajectory
    current_U <- (
      -betaPriorLogDens(
        beta = state$params$noise_beta,
        n_PP = noise$PP$n_knots,
        beta0_mean = noise$beta0_mean,
        beta0_var = noise$beta0_sd^2,
        log_scale = state$params$noise_log_scale
      ) # normal prior
      + .5 * sum(log(state$stuff$noise_var)) # det
        + .5 * sum(squared_residuals / state$stuff$noise_var) # observations
    )
    current_K <- sum(state$momenta$noise_beta^2) / 2
    proposed_U <- (
      -betaPriorLogDens(
        beta = new_noise_beta,
        n_PP = noise$PP$n_knots,
        beta0_mean = noise$beta0_mean,
        beta0_var = noise$beta0_sd^2,
        log_scale = state$params$noise_log_scale
      ) # normal prior
      + .5 * sum(log(new_noise_var)) # det
        + .5 * sum(squared_residuals / new_noise_var) # observations
    )
    proposed_K <- sum(p^2) / 2

    state$ker_var$noise_beta_mala <- updateKernel(
      iter = iter,
      iter_start = iter_start,
      kernel_value = state$ker_var$noise_beta_mala,
      mult = -.8
    )
    if (!is.nan(current_U - proposed_U + current_K - proposed_K)) {
      if (log(runif(1)) < current_U - proposed_U + current_K - proposed_K) {
        state$ker_var$noise_beta_mala <- updateKernel(
          iter = iter,
          iter_start = iter_start,
          kernel_value = state$ker_var$noise_beta_mala,
          mult = 1
        )
        dens_grad <- new_dens_grad
        state$momenta$noise_beta <- -p
        state$params$noise_beta[] <- new_noise_beta
        state$stuff$noise_var <- new_noise_var
      }
    }
    # negating momentum
    state$momenta$noise_beta <- -state$momenta$noise_beta
  }

  # updating conditioning matrix
  if (iter + iter_start > 150) {
    state$stuff$noise_beta_conditioning[] <-
      ((iter + iter_start) / (iter + iter_start - 1)) * state$stuff$noise_beta_conditioning[] +
      tcrossprod(do.call(cbind, grad_record_for_Fisher))[] / (iter + iter_start)
  }

  return(list("state" = state, "squared_residuals" = squared_residuals))
}


convertBetaAncillary <- function(PP, old_log_scale, new_log_scale, beta){
  nscale <- length(old_log_scale)
  if(nscale != length(new_log_scale) || !(nscale %in% c(1,2)))
    stop("wrong format")
  if(nscale==1){
    mult_mat <- matrix(exp(new_log_scale-old_log_scale))
  }
  if(nscale==2){
    mult_mat <- diag(exp(new_log_scale-old_log_scale)[c(1,2,2)])
  }
  rows <- seq(nrow(beta)-PP$n_knots)
  beta[-rows,] =  beta[-rows,] %*% mult_mat
  return(beta)
}


updateVarPPAncillaryXSufficient <- function(state, hierarchical_model, vecchia_approx, range_X, iter, iter_start, num_threads){
  new_range_log_scale <- state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, exp(state$ker_var$range_log_scale_sufficient))
      state$ker_var$range_log_scale_sufficient <- updateKernel(
        iter = iter, iter_start = iter_start,
        kernel_value = state$ker_var$range_log_scale_sufficient, mult = -1
      )
  if (all(inBounds(hierarchical_model$range$log_scale_bounds, new_range_log_scale))) {
    new_range_beta <- convertBetaAncillary(
      beta = state$params$range_beta, new_log_scale = new_range_log_scale, 
      old_log_scale = state$params$range_log_scale, 
      PP = hierarchical_model$range$PP
    )
    log_range <- computeLogRange(
      range_beta = new_range_beta,
      PP = hierarchical_model$range$PP,
      vecchia_approx = vecchia_approx,
      range_X = range_X
    )
    vecchia_(start_idx = 1,
             log_range = t(log_range), locs = vecchia_approx$t_locs,
             NNarray = vecchia_approx$NNarray,
             smoothness = hierarchical_model$matern_smoothness,
             compute_derivative = TRUE, num_threads = num_threads, result = state$stuff$proposed_compressed_chol
    )
    state$stuff$proposed_sparse_chol@x <- state$stuff$proposed_compressed_chol[,,1][vecchia_approx$sparse_chol_x_reorder]
    
    current_U <-       (
      + .5 * sum((state$stuff$sparse_chol %*% (state$params$field / exp(.5 * state$params$field_log_var[1, 1])))^2)
      - sum(log(Matrix::diag(state$stuff$sparse_chol)))
    )
    proposed_U <-       (
      + .5 * sum((state$stuff$proposed_sparse_chol %*% (state$params$field / exp(.5 * state$params$field_log_var[1,1])))^2)
      - sum(log(Matrix::diag(state$stuff$proposed_sparse_chol)))
    )
    
    if (log(runif(1)) < current_U - proposed_U) {
      state$ker_var$range_log_scale_sufficient <- updateKernel(
        iter = iter, iter_start = iter_start,
        kernel_value = state$ker_var$range_log_scale_sufficient, mult = 3
      )
      # invert positions of proposed and current matrix to avoid recreating objects.
      stuff_compressed <- match(c("compressed_chol", "proposed_compressed_chol"), names(state$stuff))
      stuff_compressed <- stuff_compressed[!is.na(stuff_compressed)]
      names(state$stuff)[stuff_compressed] = names(state$stuff)[rev(stuff_compressed)]
      stuff_sparse <- match(c("sparse_chol", "proposed_sparse_chol"), names(state$stuff))
      stuff_sparse <- stuff_sparse[!is.na(stuff_sparse)]
      names(state$stuff)[stuff_sparse] = names(state$stuff)[rev(stuff_sparse)]
      state$params$range_beta[] <- new_range_beta
      state$params$range_log_scale[] <- new_range_log_scale
    }
  }
  return(state)
}

