renewMomentum <- function(momentum, kept_momentum = .9) {
  if (!inBounds(c(0,1), kept_momentum, includeBounds = TRUE))
    stop("kept_momentum must be between 0 and 1")
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
  kernel_value <- max(kernel_value, -8)
  kernel_value <- min(kernel_value, 0)
  kernel_value
}

#' Title
#'
#' @param nlocs numeric, number of locations
#' @param nvars numeric, number of explicatives variables
#' @param nppr numeric, number of PP in PP range
#' @param nppn numeric, number of PP in PP noise
#'
#' @returns a list
#' @export
#'
#' @examples
#' # TODO
initStateParams <- function(stateParams){
  lapply(stateParams, function(x){
    # return(matrix(NA_real_, nrow = nrow(x), ncol = ncol(x)))
    tmp <- rep(NA_real_, length(x))
    if(is.matrix(x)){
      dim(tmp) <- dim(x)
    }
    return(tmp)
  })
}

updateBeta <- function(params, stuff, X, vecchia_approx, observed_field){
  # centered parametrization of latent field
  beta_covmat = solve(crossprod(X$X/stuff$noise_var, X$X))
  if(all(!is.infinite(beta_covmat) & !is.nan(beta_covmat))) {
    if(all(eigen(beta_covmat)$val >0)) {
      beta_mean = beta_covmat %*% crossprod(X$X, ((observed_field-params$field[vecchia_approx$locs_match]) / stuff$noise_var))
      params$beta[] = (beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))[]
    }
  }
  # un-centered parametrization of latent field
  centered_field = as.vector(params$field + X$X_locs%*%matrix(params$beta[X$which_locs], ncol = 1))
  sparse_chol_X = as.matrix(stuff$sparse_chol %*% (X$X_locs))/exp(.5 * params$field_log_var[1,1])
  beta_precision = crossprod(x = sparse_chol_X, y = sparse_chol_X)
  beta_covmat = solve(beta_precision, tol = min(rcond(beta_precision),.Machine$double.eps))
  if(all(!is.infinite(beta_covmat) & !is.nan(beta_covmat))) {
    if(all(eigen(beta_covmat)$d >0)) {
      beta_mean =  c(as.vector(stuff$sparse_chol %*% (centered_field/exp(.5 * params$field_log_var[1,1])))  %*% sparse_chol_X %*% beta_covmat)
      params$beta[X$which_locs]   = as.vector(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
      params$field = centered_field - as.vector(X$X_locs %*% matrix(params$beta[X$which_locs], ncol = 1))
    }
  }
  # updating state$stuff 
  stuff$lm_fit[]       = as.vector(X$X %*% params$beta)
  stuff$lm_residuals = observed_field - stuff$lm_fit
  return(list("params" = params, "stuff"=stuff))
}


updateLatentField <- function(params, stuff, vecchia_approx, observed_field, iter, num_threads){
  cluster_idx = 1
  chosen_locs_partition_idx = 1 + iter%%ncol(vecchia_approx$locs_partition)
  locs_partition = vecchia_approx$locs_partition[,chosen_locs_partition_idx]
  precision_from_obs = vecchia_approx$locs_match_matrix %*% (1/stuff$noise_var)
  mean_from_obs = as.vector(vecchia_approx$locs_match_matrix %*% ((observed_field - stuff$lm_fit)/stuff$noise_var))
  
  d <-  1/exp(.5 * params$field_log_var[1,1])
  unique_clusters <- unique(locs_partition)
  
  chol_list = parallel::mclapply(
    mc.cores = num_threads,
    unique( locs_partition),
    function(cluster_idx){ #posterior_precision_subset
      idx = which(locs_partition == cluster_idx)
      pbs = Matrix::crossprod(stuff$sparse_chol[,idx])
      pbs = Matrix::Diagonal(length(idx), d) %*% pbs %*% Matrix::Diagonal(length(idx), d)
      Matrix::diag(pbs) <- Matrix::diag(pbs) + as.vector(precision_from_obs[idx])
      pbs = Matrix::expand(Matrix::Cholesky(pbs))
      return(pbs)
    }
  )
  
  for(pass in seq_len(1)){
    for(cluster_idx in unique(locs_partition)){
      selected_idx = which(locs_partition==cluster_idx)
      additional_mean_from_field =
        as.vector(
          (1/exp(.5 * params$field_log_var[1,1])) * 
            Matrix::crossprod( 
              stuff$sparse_chol[,selected_idx], 
              (stuff$sparse_chol %*% 
                 ((params$field*(locs_partition!=cluster_idx))/exp(.5 * params$field_log_var[1,1]))
              ))
        )
      chololo = chol_list[[match(cluster_idx, unique(locs_partition))]]
      params$field[selected_idx] = as.vector(
        Matrix::t(chololo$P) %*% 
          Matrix::solve(Matrix::t(chololo$L), 
                        rnorm(nrow(chololo$L)) + 
                          Matrix::solve(chololo$L, # inverse of precision matrix...
                                        chololo$P %*%
                                          (-additional_mean_from_field + mean_from_obs[selected_idx]))
          ))
    }
  }
  return(params$field)
}

updateFieldLogVar <- function(state, scale, vecchia_approx, iter, iter_start){
  for(iii in seq(1)){
 # ancillary 
    for(field_log_var_idx in seq_len(3)){
      new_field_log_var = state$params$field_log_var[1,1] + exp(.5 * state$ker_var$field_log_var_ancillary) * rnorm(1)
      new_field = state$params$field * exp(.5 * (new_field_log_var - state$params$field_log_var[1,1]))
      current_U =
        (
          - beta_prior_log_dens(beta = as.matrix(state$params$field_log_var[1,1]), n_PP = 0,log_scale = 0,
                                beta0_mean = scale$beta0_mean, 
                                beta0_var =  scale$beta0_sd^2) # normal prior
          + .5 * sum((state$stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$stuff$noise_var) # observation ll
        )
      
      proposed_U =
        (
          - beta_prior_log_dens(beta = as.matrix(new_field_log_var), n_PP = 0,log_scale = 0,
                                beta0_mean = scale$beta0_mean, 
                                beta0_var =  scale$beta0_sd^2) # normal prior
          + .5 * sum((state$stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$stuff$noise_var) # observation ll
        )
      state$ker_var$field_log_var_ancillary = updateKernel(iter_start = iter_start, 
                                                            kernel_value = state$ker_var$field_log_var_ancillary, iter = iter, mult = -.25)
      if(current_U - proposed_U > log(runif(1))){
        state$ker_var$field_log_var_ancillary = updateKernel(iter_start = iter_start,  
                                                              kernel_value = state$ker_var$field_log_var_ancillary, iter = iter, mult = 1)
        state$params$field_log_var[1,1] = new_field_log_var
        state$params$field = new_field
      }
    }
    
    # Sufficient
    fieldT_cholT_chol_field = sum((state$stuff$sparse_chol %*% state$params$field)^2)
    for(field_log_var_idx in seq_len(10)){
      new_field_log_var = state$params$field_log_var[1,1] + rnorm(1) * exp(state$ker_var$field_log_var_sufficient)
      current_U =
        (
          - beta_prior_log_dens(beta = as.matrix(state$params$field_log_var[1,1]), n_PP = 0,log_scale = 0,
                                beta0_mean = scale$beta0_mean, 
                                beta0_var =  scale$beta0_sd^2) # normal prior
          + .5 * fieldT_cholT_chol_field/exp(state$params$field_log_var[1,1])  # observation ll
          +  vecchia_approx$n_locs * (.5 * state$params$field_log_var[1,1])  # observation ll
        )
      
      proposed_U =
        (
          - beta_prior_log_dens(beta = as.matrix(new_field_log_var), n_PP = 0,log_scale = 0,
                                beta0_mean = scale$beta0_mean, 
                                beta0_var =  scale$beta0_sd^2) # normal prior
          + .5 * fieldT_cholT_chol_field/exp(new_field_log_var)  # observation ll
          +  vecchia_approx$n_locs * (.5 * new_field_log_var)  # observation ll
        )
      state$ker_var$field_log_var_sufficient = updateKernel(iter_start = iter_start,  
                                                            kernel_value = state$ker_var$field_log_var_sufficient, iter = iter, mult = -.25)
      if(current_U - proposed_U > log(runif(1))){
      state$ker_var$field_log_var_sufficient = updateKernel(iter_start = iter_start, 
                                                            kernel_value = state$ker_var$field_log_var_sufficient, iter = iter, mult = 1)
        state$params$field_log_var[1,1] = new_field_log_var
      }
    }
    }
  return(state)
}

updateRangeBeta <- function(state, hierarchical_model, vecchia_approx, range_X, iter, iter_start, num_threads){
  
  range_reparam_mat = matrix(1)
  if(hierarchical_model$anisotropic){
    range_reparam_mat = matrix(c(2, 2,  0, 
                                 2, -2, 0,
                                 0, 0,  2*sqrt(2)), 3)
  }
  new_compressed_sparse_chol = state$stuff$compressed_chol
  
  ##########################
  # Range beta (ancillary) #
  ##########################
  hmc_stepsize =diag(rep(exp(state$ker_var$range_beta_ancillary[1]), 1+  2*hierarchical_model$anisotropic))
  # initializing position 
  q = 0*state$params$range_beta
  q[] = range_X$L_minus_one %*% (state$params$range_beta)
  # initializing momentum
  state$momenta$range_beta_ancillary = 
    renewMomentum(
      state$momenta$range_beta_ancillary, 
      kept_momentum = .9)
  p = state$momenta$range_beta_ancillary
  # Make a half step for momentum at the beginning
  dens_grad = 0*q
  dens_grad[] =  crossprod(range_X$L, (  
    -beta_prior_log_dens_derivative(
      beta = state$params$range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
      beta0_mean = hierarchical_model$range$beta0_mean,
      beta0_var =  hierarchical_model$range$beta0_sd^2, 
      state$params$range_log_scale) # normal prior
    + X_PP_crossprod(
      X = range_X$X_locs,  vecchia_approx = vecchia_approx, permutate_PP_to_obs = FALSE,
      PP = hierarchical_model$range$PP,
      Y = # Jacobian of range field wrt range_beta
        t(
          # natural gradient of obs likelihood wrt range field
          derivative_sandwiches_(
            vecchia = state$stuff$compressed_chol, 
            left_vector = as.vector(
              Matrix::solve(
                Matrix::t(state$stuff$sparse_chol), 
                - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                              ((state$params$field[vecchia_approx$locs_match] - state$stuff$lm_residuals) / state$stuff$noise_var))
                * exp(.5 * state$params$field_log_var[1,1]) # part of sparse chol
              )), 
            right_vector = state$params$field/exp(.5 * state$params$field_log_var[1,1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
            NNarray = vecchia_approx$NNarray, 
            sauce_determinant_chef = FALSE, 
            num_threads = num_threads
          )
        )) %*% range_reparam_mat
  ))
  
  
  #    #testing the gradient
  #     gradient_test_chol = 0 * state$stuff$compressed_chol[,,1,drop= F]
  #     derivative_test = 0*q
  #     i=1
  #     j=1
  #     for(i in seq(nrow(state$params$range_beta))){
  #      for(j in seq(ncol(state$params$range_beta))){
  #        q_ = q
  #        q_[i,j] = q_[i,j] +  0.0000001
  #        range_beta_ = 0*state$params$range_beta
  #        range_beta_[] = range_X$L %*% q_
  #        vecchia_(
  #          log_range = t(compute_log_range(
  #            range_beta = range_beta_, 
  #            PP = hierarchical_model$range$PP,
  #            vecchia_approx = vecchia_approx, 
  #            range_X = range_X
  #          )), locs = vecchia_approx$t_locs, 
  #          NNarray = vecchia_approx$NNarray, 
  #          smoothness = hierarchical_model$matern_smoothness, 
  #          compute_derivative = F, num_threads = num_threads, 
  #          result = gradient_test_chol
  #        )
  #        sparse_chol_ = decompress_chol(
  #          vecchia_approx = vecchia_approx, 
  #          compressed_sparse_chol = gradient_test_chol
  #        ) 
  #        field_ = as.vector(Matrix::solve(sparse_chol_, state$stuff$sparse_chol %*% (state$params$field)))
  #        
  #     derivative_test[i,j] = 
  #        (1/ 0.0000001)*(
  #          beta_prior_log_dens(beta = state$params$range_beta, 
  #                              n_PP = hierarchical_model$range$PP$n_knots, 
  #                              beta0_mean = hierarchical_model$range$beta0_mean, 
  #                              beta0_var =  hierarchical_model$range$beta0_sd^2, 
  #                              state$params$range_log_scale)
  #          - beta_prior_log_dens(beta = range_beta_, 
  #                                n_PP = hierarchical_model$range$PP$n_knots, 
  #                                beta0_mean = hierarchical_model$range$beta0_mean, 
  #                                beta0_var =  hierarchical_model$range$beta0_sd^2, 
  #                                state$params$range_log_scale)
  #         - .5 * sum( (state$stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$stuff$noise_var) # observation ll
  #         + .5 * sum( (state$stuff$lm_residuals -  field_      [vecchia_approx$locs_match])^2/state$stuff$noise_var) # observation ll
  #        ) 
  #      }
  #     }
  #     par(mfrow = c(1,3))
  #     boxplot(c(derivative_test)/c(dens_grad))
  #     boxplot(c(derivative_test)-c(dens_grad))
  #     plot(c(derivative_test), c(dens_grad))
  #     abline(a=0, b=1)
  
  p[] = p[] - (dens_grad) %*%hmc_stepsize/ 2
  n_hmc_steps = 3
  for(hmc_step in seq_len(n_hmc_steps)){
    # Make a full step for the positio
    q = q + p %*% hmc_stepsize
    new_range_beta = state$params$range_beta
    new_range_beta[] = range_X$L %*% q
    log_range = compute_log_range(
      range_beta = new_range_beta, 
      PP = hierarchical_model$range$PP,
      vecchia_approx = vecchia_approx, 
      range_X = range_X
    )
    vecchia_(
      log_range = t(log_range), locs = vecchia_approx$t_locs, 
      NNarray = vecchia_approx$NNarray, 
      smoothness = hierarchical_model$matern_smoothness, 
      compute_derivative = TRUE, num_threads = num_threads, result = new_compressed_sparse_chol
    )
    new_sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, new_compressed_sparse_chol)
    new_field = exp(.5 * state$params$field_log_var[1,1]) * 
      as.vector(Matrix::solve(new_sparse_chol, state$stuff$sparse_chol %*% (state$params$field/exp(.5 * state$params$field_log_var[1,1]))))
    # Make a half step for momentum at the end.
    dens_grad[] =   crossprod(range_X$L, (  
      -beta_prior_log_dens_derivative(
        beta = new_range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
        beta0_mean = hierarchical_model$range$beta0_mean,
        beta0_var =  hierarchical_model$range$beta0_sd^2, 
        state$params$range_log_scale) # normal prior
      + X_PP_crossprod(
        X = range_X$X_locs, vecchia_approx = vecchia_approx, 
        permutate_PP_to_obs = FALSE,
        PP = hierarchical_model$range$PP,
        Y = # Jacobian of range field wrt range_beta
          t(
            # natural gradient of obs likelihood wrt range field
            derivative_sandwiches_(
              vecchia = new_compressed_sparse_chol, 
              left_vector = as.vector(
                Matrix::solve(
                  Matrix::t(new_sparse_chol), 
                  - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                ((new_field[vecchia_approx$locs_match] - state$stuff$lm_residuals) / state$stuff$noise_var))
                  * exp(.5 * state$params$field_log_var[1,1]) # part of sparse chol
                )), 
              right_vector = new_field/exp(.5 * state$params$field_log_var[1,1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
              NNarray = vecchia_approx$NNarray, 
              sauce_determinant_chef = FALSE, 
              num_threads = num_threads
            )
          )) %*% range_reparam_mat 
    ))
    
    #    #testing the gradient
    #     gradient_test_chol = 0 * state$stuff$compressed_chol
    #     derivative_test = 0*q
    #     j=1
    #     j=1
    #     for(i in seq(nrow(state$params$range_beta))){
    #      for(j in seq(ncol(state$params$range_beta))){
    #        q_ = q
    #        q_[i,j] = q_[i,j] +  0.0000001
    #        range_beta_ = 0*state$params$range_beta
    #        range_beta_[] = range_X$L %*% q_
    #        vecchia_(
    #          log_range = t(compute_log_range(
    #            range_beta = range_beta_, 
    #            PP = hierarchical_model$range$PP,
    #            vecchia_approx = vecchia_approx, 
    #            range_X = range_X
    #          )), locs = vecchia_approx$t_locs, 
    #          NNarray = vecchia_approx$NNarray, 
    #          smoothness = hierarchical_model$matern_smoothness, 
    #          compute_derivative = TRUE, num_threads = num_threads, 
    #          result = gradient_test_chol
    #        )
    #        sparse_chol_ = decompress_chol(
    #          vecchia_approx = vecchia_approx, 
    #          compressed_sparse_chol = gradient_test_chol
    #        ) 
    #        field_ = as.vector(Matrix::solve(sparse_chol_, new_sparse_chol %*% (new_field)))
    #        
    #     derivative_test[i,j] = 
    #        (1/ 0.0000001)*(
    #          beta_prior_log_dens(beta = new_range_beta, 
    #                              n_PP = hierarchical_model$range$PP$n_knots, 
    #                              beta0_mean = hierarchical_model$range$beta0_mean, 
    #                              beta0_var =  hierarchical_model$range$beta0_sd^2, 
    #                              state$params$range_log_scale)
    #          - beta_prior_log_dens(beta = range_beta_, 
    #                                n_PP = hierarchical_model$range$PP$n_knots, 
    #                                beta0_mean = hierarchical_model$range$beta0_mean, 
    #                                beta0_var =  hierarchical_model$range$beta0_sd^2, 
    #                                state$params$range_log_scale)
    #          - .5 * sum( (state$stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$stuff$noise_var) # observation ll
    #          + .5 * sum( (state$stuff$lm_residuals -  field_      [vecchia_approx$locs_match])^2/state$stuff$noise_var) # observation ll
    #        ) 
    #      }
    #     }
    #     par(mfrow = c(1,3))
    #     boxplot(c(derivative_test)/c(dens_grad))
    #     boxplot(c(derivative_test)-c(dens_grad))
    #     plot(c(derivative_test), c(dens_grad))
    #     abline(a=0, b=1)
    
    p[] = p[] - (dens_grad)%*%hmc_stepsize/ (1+  (hmc_step==n_hmc_steps))
  }
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_K = sum (state$momenta$range_beta_ancillary ^2) / 2
  proposed_K = sum(p^2) / 2
  
  current_U =
    (
      - beta_prior_log_dens(beta = state$params$range_beta, 
                            n_PP = hierarchical_model$range$PP$n_knots, 
                            beta0_mean = hierarchical_model$range$beta0_mean, 
                            beta0_var =  hierarchical_model$range$beta0_sd^2, 
                            state$params$range_log_scale) # normal prior
      + .5 * sum((state$stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$stuff$noise_var) # observation ll
    )
  proposed_U =
    (
      - beta_prior_log_dens(beta = new_range_beta, 
                            n_PP = hierarchical_model$range$PP$n_knots, 
                            beta0_mean = hierarchical_model$range$beta0_mean, 
                            beta0_var =  hierarchical_model$range$beta0_sd^2, 
                            state$params$range_log_scale) # normal prior
      + .5 * sum((state$stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$stuff$noise_var) # observation ll
    )
  
  # current_U-proposed_U
  # current_K- proposed_K
  
  state$ker_var$range_beta_ancillary[1] = updateKernel(
    iter = iter, iter_start = iter_start, 
    kernel_value = state$ker_var$range_beta_ancillary[1], mult = -.8
  )
  
  if(!is.nan(current_U-proposed_U+current_K- proposed_K))
  {
    if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
    {
      state$ker_var$range_beta_ancillary[1] = updateKernel(
        iter = iter, iter_start = iter_start, 
        kernel_value = state$ker_var$range_beta_ancillary[1], mult = 1
      )
      state$momenta$range_beta_ancillary = p
      state$params$field = new_field
      state$stuff$sparse_chol= new_sparse_chol
      state$stuff$compressed_chol[] = new_compressed_sparse_chol[]
      state$params$range_beta[] = new_range_beta
    }
  }
  ###########################
  # Range beta (sufficient) #
  ###########################
  hmc_stepsize =diag(rep(exp(state$ker_var$range_beta_sufficient[1]), 1+  2*hierarchical_model$anisotropic))
  q[] = range_X$L_minus_one %*% (state$params$range_beta)
  # initializing momentum
  state$momenta$range_beta_sufficient = 
    renewMomentum(
      state$momenta$range_beta_sufficient, 
      kept_momentum = .9)
  p = state$momenta$range_beta_sufficient
  # Make a half step for momentum at the beginning
  dens_grad = 0*q
  dens_grad[] =  crossprod(range_X$L, (
    - beta_prior_log_dens_derivative(
      beta = state$params$range_beta, 
      n_PP = hierarchical_model$range$PP$n_knots, 
      beta0_mean = hierarchical_model$range$beta0_mean, 
      beta0_var =  hierarchical_model$range$beta0_sd^2, 
      log_scale = state$params$range_log_scale) # normal prior
    # normal prior derivative                
    + X_PP_crossprod(
      X = range_X$X_locs, PP = hierarchical_model$range$PP,
      permutate_PP_to_obs = FALSE, 
      vecchia_approx = vecchia_approx,
      Y = # Jacobian of range field wrt range_beta
        t(# natural gradient of obs likelihood wrt range field
          derivative_sandwiches_(
            vecchia = state$stuff$compressed_chol, # derivative of the (unscaled) NNGP factor
            left_vector = as.vector(state$stuff$sparse_chol %*% (state$params$field/exp(.5 * state$params$field_log_var[1,1]))), # left vector = whitened latent field
            right_vector = state$params$field/exp(.5 * state$params$field_log_var[1,1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
            NNarray = vecchia_approx$NNarray, 
            sauce_determinant_chef = TRUE, 
            num_threads = num_threads  
          )
        ))  
    %*% range_reparam_mat
  ))
  p[] = p[] - (dens_grad) %*% hmc_stepsize/ 2
  
  
  ####testing the gradient
  # gradient_test_chol = 0 * state$stuff$compressed_chol[,,1,drop=F]
  # derivative_test = 0*q
  # i=1
  # j=1
  # for(i in seq(nrow(state$params$range_beta))){
  #   for(j in seq(ncol(state$params$range_beta))){
  #     q_ = q
  #     q_[i,j] = q_[i,j] + .000001
  #     range_beta_ = 0*state$params$range_beta
  #     range_beta_[] = range_X$L %*% q_
  #     vecchia_(
  #     log_range = t(compute_log_range(
  #       range_beta = range_beta_, 
  #       PP = hierarchical_model$range$PP,
  #       vecchia_approx = vecchia_approx, 
  #       range_X = range_X
  #     )), locs = vecchia_approx$t_locs, 
  #     NNarray = vecchia_approx$NNarray, 
  #     smoothness = hierarchical_model$matern_smoothness, 
  #     compute_derivative = F, num_threads = num_threads, 
  #     result = gradient_test_chol
  #   )
  #   sparse_chol_ = decompress_chol(
  #     vecchia_approx = vecchia_approx, 
  #     compressed_sparse_chol = gradient_test_chol
  #   ) 
  #     derivative_test[i,j] = 
  #       1000000*(
  #         beta_prior_log_dens(beta = state$params$range_beta, 
  #                               n_PP = hierarchical_model$range$PP$n_knots, 
  #                               beta0_mean = hierarchical_model$range$beta0_mean, 
  #                               beta0_var =  hierarchical_model$range$beta0_sd^2, 
  #                               state$params$range_log_scale) - 
  #          beta_prior_log_dens(beta = range_beta_, 
  #                               n_PP = hierarchical_model$range$PP$n_knots, 
  #                               beta0_mean = hierarchical_model$range$beta0_mean, 
  #                               beta0_var =  hierarchical_model$range$beta0_sd^2, 
  #                               state$params$range_log_scale) +
  #         (
  #           0
  #           + .5* sum((sparse_chol_ %*% (state$params$field/exp(.5 * state$params$field_log_var)))^2)
  #           - sum(log(Matrix::diag(sparse_chol_)))
  #          )-
  #           (
  #             0
  #             + .5* sum((state$stuff$sparse_chol %*% (state$params$field/exp(.5 * state$params$field_log_var)))^2)
  #             - sum(log(Matrix::diag(state$stuff$sparse_chol)))
  #           )
  #       )
  #   }
  # }
  # par(mfrow = c(1,3))
  # boxplot(c(derivative_test)/c(dens_grad))
  # boxplot(c(derivative_test)-c(dens_grad))
  # plot(c(derivative_test), c(dens_grad))
  # abline(a=0, b=1)
  
  
  
  n_hmc_steps = 3
  for(hmc_step in seq_len(n_hmc_steps)){
    # Make a full step for the position
    q = q + p %*% hmc_stepsize
    new_range_beta = state$params$range_beta      
    new_range_beta[] = range_X$L %*% q
    log_range = compute_log_range(
      range_beta = new_range_beta, 
      PP = hierarchical_model$range$PP,
      vecchia_approx = vecchia_approx, 
      range_X = range_X
    )
    vecchia_(
      log_range = t(log_range), locs = vecchia_approx$t_locs, 
      NNarray = vecchia_approx$NNarray, 
      smoothness = hierarchical_model$matern_smoothness, 
      compute_derivative = TRUE, num_threads = num_threads, result = new_compressed_sparse_chol
    )
    new_sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, new_compressed_sparse_chol)
    
    # Make a half step for momentum at the end.
    dens_grad = 0*q
    dens_grad[] = crossprod(range_X$L, (
      - beta_prior_log_dens_derivative(
        beta = new_range_beta, 
        n_PP = hierarchical_model$range$PP$n_knots, 
        beta0_mean = hierarchical_model$range$beta0_mean, 
        beta0_var =  hierarchical_model$range$beta0_sd^2, 
        log_scale = state$params$range_log_scale) # normal prior
      # normal prior derivative                
      + X_PP_crossprod(
        X = range_X$X_locs, PP = hierarchical_model$range$PP,
        permutate_PP_to_obs = FALSE, 
        vecchia_approx = vecchia_approx,
        Y = # Jacobian of range field wrt range_beta
          t(# natural gradient of obs likelihood wrt range field
            derivative_sandwiches_(
              vecchia = new_compressed_sparse_chol, # derivative of the (unscaled) NNGP factor
              left_vector = as.vector(new_sparse_chol %*% (state$params$field/exp(.5 * state$params$field_log_var[1,1]))), # left vector = whitened latent field
              right_vector = state$params$field/exp(.5 * state$params$field_log_var[1,1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
              NNarray = vecchia_approx$NNarray, 
              sauce_determinant_chef = TRUE, 
              num_threads = num_threads  
            )
          )
      )  %*% range_reparam_mat
    ))
    
    
    ## ###testing the gradient
    ## gradient_test_chol = 0 * state$stuff$compressed_chol[,,1,drop=F]
    ## derivative_test = 0*q
    ## i=1
    ## j=1
    ## for(i in seq(nrow(state$params$range_beta))){
    ##   for(j in seq(ncol(state$params$range_beta))){
    ##     q_ = q
    ##     q_[i,j] = q_[i,j] + .000001
    ##     range_beta_ = 0*state$params$range_beta
    ##     range_beta_[] = range_X$L %*% q_
    ##     vecchia_(
    ##       log_range = t(compute_log_range(
    ##         range_beta = range_beta_, 
    ##         PP = hierarchical_model$range$PP,
    ##         vecchia_approx = vecchia_approx, 
    ##         range_X = range_X
    ##       )), locs = vecchia_approx$t_locs, 
    ##       NNarray = vecchia_approx$NNarray, 
    ##       smoothness = hierarchical_model$matern_smoothness, 
    ##       compute_derivative = F, num_threads = num_threads, 
    ##       result = gradient_test_chol
    ##     )
    ##     sparse_chol_ = decompress_chol(
    ##       vecchia_approx = vecchia_approx, 
    ##       compressed_sparse_chol = gradient_test_chol
    ##     ) 
    ##     derivative_test[i,j] = 
    ##       1000000*(
    ##         beta_prior_log_dens(beta = new_range_beta, 
    ##                             n_PP = hierarchical_model$range$PP$n_knots, 
    ##                             beta0_mean = hierarchical_model$range$beta0_mean, 
    ##                             beta0_var =  hierarchical_model$range$beta0_sd^2, 
    ##                             state$params$range_log_scale) - 
    ##           beta_prior_log_dens(beta = range_beta_, 
    ##                               n_PP = hierarchical_model$range$PP$n_knots, 
    ##                               beta0_mean = hierarchical_model$range$beta0_mean, 
    ##                               beta0_var =  hierarchical_model$range$beta0_sd^2, 
    ##                               state$params$range_log_scale) +
    ##           (
    ##             0
    ##             + .5* sum((sparse_chol_ %*% (state$params$field/exp(.5 * state$params$field_log_var)))^2)
    ##             - sum(log(Matrix::diag(sparse_chol_)))
    ##           )-
    ##           (
    ##             0
    ##             + .5* sum((new_sparse_chol %*% (state$params$field/exp(.5 * state$params$field_log_var)))^2)
    ##             - sum(log(Matrix::diag(new_sparse_chol)))
    ##           )
    ##       )
    ##    }
    ##  }
    ##  par(mfrow = c(1,3))
    ##  boxplot(c(derivative_test)/c(dens_grad))
    ##  boxplot(c(derivative_test)-c(dens_grad))
    ##  plot(c(derivative_test), c(dens_grad))
    ##  abline(a=0, b=1)
    
    # updating momentum
    p[] = p[] - (dens_grad) %*% hmc_stepsize/ (1 + (hmc_step == n_hmc_steps))
  }
  
  # metropolis step
  current_K = sum (state$momenta$range_beta_sufficient ^2) / 2
  proposed_K = sum(p^2) / 2
  current_U = (
    - beta_prior_log_dens(beta = state$params$range_beta, 
                          n_PP = hierarchical_model$range$PP$n_knots, 
                          beta0_mean = hierarchical_model$range$beta0_mean, 
                          beta0_var =  hierarchical_model$range$beta0_sd^2, 
                          log_scale = state$params$range_log_scale)
    # normal prior 
    + .5* sum((state$stuff$sparse_chol %*% (state$params$field/exp(.5 * state$params$field_log_var[1,1])))^2)
    - sum(log(Matrix::diag(state$stuff$sparse_chol)))
  )
  proposed_U = (
    - beta_prior_log_dens(beta = new_range_beta, 
                          n_PP = hierarchical_model$range$PP$n_knots, 
                          beta0_mean = hierarchical_model$range$beta0_mean, 
                          beta0_var =  hierarchical_model$range$beta0_sd^2, 
                          state$params$range_log_scale)
    # normal prior 
    + .5* sum((new_sparse_chol %*% (state$params$field/exp(.5 * state$params$field_log_var[1,1])))^2)
    - sum(log(Matrix::diag(new_sparse_chol)))
  )
  
  # current_U - proposed_U
  # current_K - proposed_K
  state$ker_var$range_beta_sufficient[1] = updateKernel(
    iter = iter, iter_start = iter_start,  
    kernel_value = state$ker_var$range_beta_sufficient[1], mult = -.8
  )
  
  if(!is.nan(current_U-proposed_U+current_K- proposed_K)) {
    if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K) {
      state$ker_var$range_beta_sufficient[1] = updateKernel(
        iter = iter, iter_start = iter_start, 
        kernel_value = state$ker_var$range_beta_sufficient[1], mult = 1
      )
      state$momenta$range_beta_sufficient = p
      state$stuff$sparse_chol= new_sparse_chol
      state$stuff$compressed_chol[] = new_compressed_sparse_chol[]
      state$params$range_beta[] = new_range_beta
    }
  }
  return(list(state = state))
}


UpdateVarPPSuff = function(hm4params, beta4params, current_range_log_scale){
  for(i in seq_len(10))
  {
    q = current_range_log_scale + rnorm(length(current_range_log_scale), 0, .05)
    if(all(inBounds(hm4params$log_scale_bounds, q)) &
       ( 
         + beta_prior_log_dens(
           beta = beta4params, n_PP = hm4params$PP$n_knots, 
           beta0_mean = hm4params$beta0_mean,
           beta0_var =  hm4params$beta0_sd^2, 
           log_scale = q)
         - beta_prior_log_dens(
           beta = beta4params, n_PP = hm4params$PP$n_knots, 
           beta0_mean = hm4params$beta0_mean,
           beta0_var =  hm4params$beta0_sd^2, 
           log_scale = current_range_log_scale) 
         > log(runif(1))
       )
    )
    {
      current_range_log_scale[] = q
    }
  }
  return(current_range_log_scale)
}
  
PotSuffRangeLogScale = function(sparse_chol, field, field_log_var){
  return(
    (
      + .5* sum((sparse_chol %*%field)^2)/exp(field_log_var[1,1])
      - sum(log(Matrix::diag(sparse_chol)))
    )
  )
}
PotAnciRangeLogScale = function(field, lm_residuals, noise_var, vecchia_approx){
  return(+ .5 * sum((lm_residuals -field[vecchia_approx$locs_match])^2/noise_var))
}

  
updateVarianceRangepp <- function(state, hierarchical_model, range_X, vecchia_approx, iter, iter_start, num_threads){
  new_compressed_sparse_chol = state$stuff$compressed_chol[,,1,drop = F]
  n_range_log_scale_update = 2
  # ancillary - sufficient ####
  for(i in seq_len(n_range_log_scale_update)){
    q = state$params$range_log_scale + rnorm(1 + hierarchical_model$anisotropic, 0, exp(state$ker_var$range_log_scale_sufficient))
    
    if (all(inBounds(hierarchical_model$range$log_scale_bounds, q))){
      new_range_beta = state$params$range_beta
      new_range_beta[-seq_len(range_X$n_regressors),] = 
        new_range_beta[-seq_len(range_X$n_regressors),] %*% 
        diag(exp(-.5 * state$params$range_log_scale[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic) %*%
        diag(exp(.5 * q[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic)
      log_range = compute_log_range(
        range_beta = new_range_beta, 
        PP = hierarchical_model$range$PP,
        vecchia_approx = vecchia_approx, 
        range_X = range_X
      )
      vecchia_(
        log_range = t(log_range), locs = vecchia_approx$t_locs, 
        NNarray = vecchia_approx$NNarray, 
        smoothness = hierarchical_model$matern_smoothness, 
        compute_derivative = FALSE, num_threads = num_threads, result = new_compressed_sparse_chol
      )
      new_sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, new_compressed_sparse_chol)
    
    current_U =  PotSuffRangeLogScale(sparse_chol = state$stuff$sparse_chol, field = state$params$field, field_log_var = state$params$field_log_var)
    proposed_U = PotSuffRangeLogScale(sparse_chol = new_sparse_chol,         field = state$params$field, field_log_var = state$params$field_log_var)
    
    state$ker_var$range_log_scale_sufficient = updateKernel(
      iter = iter, iter_start = iter_start, 
      kernel_value = state$ker_var$range_log_scale_sufficient, 
      mult = -.25
    )
      if(!is.nan(current_U-proposed_U)){
        if(log(runif(1)) < (current_U-proposed_U))
        {
          
          state$ker_var$range_log_scale_sufficient = updateKernel(
            iter = iter, iter_start = iter_start, 
            kernel_value = state$ker_var$range_log_scale_sufficient, 
            mult = 1
          )
          
          state$params$range_beta = new_range_beta
          state$params$range_log_scale[] = q
          
          state$stuff$sparse_chol= new_sparse_chol
          state$stuff$compressed_chol[,,1] = new_compressed_sparse_chol[,,1]
        }}}
  }
  
  # sufficient - sufficient ####
  
  state$params$range_log_scale = UpdateVarPPSuff(hm4params = hierarchical_model$range, beta4param = state$params$range_beta, current_range_log_scale = state$params$range_log_scale)
  
  # ancillary-ancillary ####
  for(i in seq_len(n_range_log_scale_update)){
    q = state$params$range_log_scale + rnorm(1 + hierarchical_model$anisotropic, 0, exp(state$ker_var$range_log_scale_ancillary))
    
    if (all(inBounds(hierarchical_model$range$log_scale_bounds, q))){
    new_range_beta = state$params$range_beta
    new_range_beta[-seq_len(range_X$n_regressors),] = 
      new_range_beta[-seq_len(range_X$n_regressors),] %*% 
      diag(exp(-.5 * state$params$range_log_scale[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic) %*%
      diag(exp(.5 * q[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic)
    log_range = compute_log_range(
      range_beta = new_range_beta, 
      PP = hierarchical_model$range$PP,
      vecchia_approx = vecchia_approx, 
      range_X = range_X
    )
    vecchia_(
      log_range = t(log_range), locs = vecchia_approx$t_locs, 
      NNarray = vecchia_approx$NNarray, 
      smoothness = hierarchical_model$matern_smoothness, 
      compute_derivative = FALSE, num_threads = num_threads, result = new_compressed_sparse_chol
    )
    new_sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, new_compressed_sparse_chol)
    new_field = as.vector(Matrix::solve(new_sparse_chol, state$stuff$sparse_chol %*% (state$params$field)))
    
    current_U = PotAnciRangeLogScale(field = state$params$field, lm_residuals = state$stuff$lm_residuals, noise_var = state$stuff$noise_var, vecchia_approx = vecchia_approx)
    proposed_U = PotAnciRangeLogScale(field = new_field,         lm_residuals = state$stuff$lm_residuals, noise_var = state$stuff$noise_var, vecchia_approx = vecchia_approx)
    state$ker_var$range_log_scale_ancillary = updateKernel(
      iter = iter, iter_start = iter_start, 
      kernel_value = state$ker_var$range_log_scale_ancillary, 
      mult = -.25
    )
    
      if(!is.nan(current_U-proposed_U)){
        if(log(runif(1)) < (current_U-proposed_U)){
          state$ker_var$range_log_scale_ancillary = updateKernel(
            iter = iter, iter_start = iter_start,
            kernel_value = state$ker_var$range_log_scale_ancillary, 
            mult = 1
          )
          
          state$params$range_beta = new_range_beta
          state$params$range_log_scale[] = q
          state$params$field = new_field
          
          state$stuff$sparse_chol= new_sparse_chol
          state$stuff$compressed_chol[,,1] = new_compressed_sparse_chol[,,1]
        }}}
  }
  # sufficient - sufficient ####
  state$params$range_log_scale = UpdateVarPPSuff(hm4params = hierarchical_model$range, beta4param = state$params$range_beta, current_range_log_scale = state$params$range_log_scale)
  # return ####
  return(state)
}

updateNoiseBeta <- function(state, noise, noise_X, vecchia_approx, iter, iter_start){
  # VEWY IMPOWTANT don't remove or comment
  squared_residuals = as.matrix(state$stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2
  # HMC update
  q = noise_X$L_minus_one %*% state$params$noise_beta
  state$momenta$noise_beta = renewMomentum(state$momenta$noise_beta)
  p = state$momenta$noise_beta
  dens_grad = (
    - beta_prior_log_dens_derivative(
      beta = state$params$noise_beta, n_PP = noise$PP$n_knots, 
      beta0_mean = noise$beta0_mean,
      beta0_var =  noise$beta0_sd^2, 
      log_scale = state$params$noise_log_scale) # normal prior
    + X_PP_crossprod(
      X = noise_X$X, noise$PP, vecchia_approx = vecchia_approx, permutate_PP_to_obs = TRUE, 
      Y = 
        (
          + .5 # determinant part of normal likelihood
          - (squared_residuals/state$stuff$noise_var)/2 # exponential part of normal likelihood
        ))
  )
  # Make a half step for momentum at the beginning
  exp_noise_mala <- exp(state$ker_var$noise_beta_mala)
  p = p - exp_noise_mala * crossprod(noise_X$L, dens_grad) / 2
  
  n_hmc_steps = 10
  for(hmc_step in seq_len(n_hmc_steps)){
    # Make a full step for the position
    q = q + exp_noise_mala * p
    new_noise_beta = noise_X$L %*% q
    new_noise_var = as.vector(exp(X_PP_mult_right(
      X = noise_X$X, PP = noise$PP,
      vecchia_approx = vecchia_approx, Y = new_noise_beta, 
      permutate_PP_to_obs = TRUE
    )))
    # Make a half step for momentum at the end
    dens_grad = (
      - beta_prior_log_dens_derivative(
        beta = new_noise_beta, n_PP = noise$PP$n_knots, 
        beta0_mean = noise$beta0_mean,
        beta0_var =  noise$beta0_sd^2, 
        log_scale = state$params$noise_log_scale
      ) # normal prior
      + X_PP_crossprod(
        X = noise_X$X, 
        noise$PP, 
        vecchia_approx = vecchia_approx, 
        permutate_PP_to_obs = TRUE, 
        Y = (+ .5 # determinant part of normal likelihood
             - (squared_residuals/new_noise_var)/2 # exponential part of normal likelihood
        )
      )
    )
    p = p - exp_noise_mala * crossprod(noise_X$L, dens_grad) / (1 + (hmc_step == n_hmc_steps))
  }
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = (
    - beta_prior_log_dens(
      beta = state$params$noise_beta, 
      n_PP = noise$PP$n_knots, 
      beta0_mean = noise$beta0_mean,
      beta0_var =  noise$beta0_sd^2, 
      log_scale = state$params$noise_log_scale) # normal prior 
    +.5* sum(log(state$stuff$noise_var)) # det
    +.5*sum(squared_residuals/state$stuff$noise_var) # observations
  )
  current_K = sum (state$momenta$noise_beta ^2) / 2
  proposed_U = (
    - beta_prior_log_dens(beta = new_noise_beta, 
                          n_PP = noise$PP$n_knots, 
                          beta0_mean = noise$beta0_mean,
                          beta0_var =  noise$beta0_sd^2, 
                          log_scale = state$params$noise_log_scale) # normal prior        
    +.5* sum(log(new_noise_var)) # det
    +.5*sum(squared_residuals/new_noise_var) # observations
  )
  proposed_K = sum(p^2) / 2
  
  state$ker_var$noise_beta_mala = updateKernel(iter = iter, iter_start = iter_start, kernel_value = state$ker_var$noise_beta_mala, mult = -.8)
  if(!is.nan(current_U-proposed_U+current_K- proposed_K)) {
    if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K) {
      state$ker_var$noise_beta_mala = updateKernel(iter = iter, iter_start = iter_start, kernel_value = state$ker_var$noise_beta_mala, mult = 1)
      state$momenta$noise_beta = p
      state$params$noise_beta[] = new_noise_beta
      state$stuff$noise_var = new_noise_var
    }
  }
  
  return(list("state" = state, "squared_residuals" = squared_residuals))
}

updateNoiseLogScale <- function(state, noise, noise_X, vecchia_approx, squared_residuals, iter, iter_start){
  # ancillary -- sufficient ####
  for(i in seq_len(4)) {
    new_noise_log_scale = state$params$noise_log_scale + rnorm(1, 0, exp(state$ker_var$noise_log_scale))
    new_noise_beta = state$params$noise_beta
    new_noise_beta[-seq_len(noise_X$n_regressors)] = 
      new_noise_beta[-seq_len(noise_X$n_regressors)] *
      c(exp((new_noise_log_scale - state$params$noise_log_scale)/2))
    new_noise_var = as.vector(
      exp(
        X_PP_mult_right(
          X = noise_X$X, 
          PP = noise$PP,
          vecchia_approx = vecchia_approx, 
          Y = new_noise_beta, 
          permutate_PP_to_obs = TRUE
        )
      )
    )
    dens_ratio = (
      -.5* sum(log(new_noise_var)) 
      -.5*sum(squared_residuals/new_noise_var)
      +.5* sum(log(state$stuff$noise_var)) 
      +.5*sum(squared_residuals/state$stuff$noise_var)
    )
    if(!is.nan(dens_ratio)) {
      if(dens_ratio > log(runif(1))) {
        if(inBounds(noise$log_scale_bounds, new_noise_log_scale)) {
          state$params$noise_log_scale[] = new_noise_log_scale
          state$params$noise_beta = new_noise_beta 
          state$stuff$noise_var = new_noise_var
          state$ker_var$noise_log_scale = updateKernel(
            iter = iter, iter_start = iter_start, 
            kernel_value = state$ker_var$noise_log_scale, mult = 1)
        }
      }
    }
    state$ker_var$noise_log_scale = updateKernel(
      iter = iter, iter_start = iter_start, 
      kernel_value = state$ker_var$noise_log_scale, mult = -.25)
  }
  
  # sufficient -- sufficient ####
  for(i in seq_len(10)) {
    new_noise_log_scale = state$params$noise_log_scale + rnorm(1, 0, .1)
    if(inBounds(noise$log_scale_bounds, new_noise_log_scale)) {
      if(
        + beta_prior_log_dens(beta = state$params$noise_beta, n_PP = noise$PP$n_knots, 
                              beta0_mean = noise$beta0_mean,
                              beta0_var =  noise$beta0_sd^2, 
                              log_scale = new_noise_log_scale) - 
        beta_prior_log_dens(beta = state$params$noise_beta, n_PP = noise$PP$n_knots, 
                            beta0_mean = noise$beta0_mean,
                            beta0_var =  noise$beta0_sd^2, 
                            log_scale = state$params$noise_log_scale)
        > log(runif(1))) {
        state$params$noise_log_scale[] = new_noise_log_scale
      }
    }
  }
  return(state)
}