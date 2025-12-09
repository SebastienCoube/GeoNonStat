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
    if(is.matrix(x)){
      return(matrix(NA_real_, nrow = nrow(x), ncol = ncol(x)))
    } else if(is.numeric(x)) {
      return(rep(NA_real_, length(x)))
    }
  })
}

update_beta <- function(params, stuff, X, vecchia_approx, observed_field){
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


update_latent_field <- function(params, stuff, vecchia_approx, observed_field, iter, num_threads){
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
  
  for(pass in seq_len(3)){
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
  return(params)
}

update_field_log_var <- function(state, scale, vecchia_approx, iter, iter_start){
  # ancillary 
  for(idx_for_field_logvar in seq_len(4)){
    for(field_log_var_idx in seq_len(2)){
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
      state$ker_var$field_log_var_ancillary = update_kernel(iter_start = iter_start, update_kernel_groupsize = 1, 
                                                      kernel_value = state$ker_var$field_log_var_ancillary, iter = iter, mult = -.25)
      if(current_U - proposed_U > log(runif(1))){
        state$ker_var$field_log_var_ancillary = update_kernel(iter_start = iter_start, update_kernel_groupsize = 1, 
                                                        kernel_value = state$ker_var$field_log_var_ancillary, iter = iter, mult = 1)
        state$params$field_log_var[1,1] = new_field_log_var
        state$params$field = new_field
      }
    }
    
    # Sufficient
    fieldT_cholT_chol_field = sum((state$stuff$sparse_chol %*% state$params$field)^2)
    for(field_log_var_idx in seq_len(10)){
      new_field_log_var = state$params$field_log_var[1,1] + .1 * rnorm(1)
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
      if(current_U - proposed_U > log(runif(1))){
        state$params$field_log_var[1,1] = new_field_log_var
      }
    }
  }
  return(state)
}

update_range_beta <- function(state, hierarchical_model, range_X, iter, iter_start, num_threads){
  
  range_reparam_mat = matrix(1)
  if(hierarchical_model$anisotropic){
    range_reparam_mat = matrix(c(2, 2,  0, 
                                 2, -2, 0,
                                 0, 0,  2*sqrt(2)), 3)
  }
  L_minus_one =  solve(t(chol(
    solve(range_X$crossprod_X_locs)/max(solve(range_X$crossprod_X_locs))
  )))
  
  ##########################
  # Range beta (ancillary) #
  ##########################
  for(regime in seq_len(1 + hierarchical_model$anisotropic)){
    if(!hierarchical_model$anisotropic) hmc_stepsize =matrix(exp(state$ker_var$range_beta_ancillary))
    if(hierarchical_model$anisotropic & regime ==1) hmc_stepsize = diag(c(exp(state$ker_var$range_beta_ancillary[1]),0,0))
    if(hierarchical_model$anisotropic & regime ==2) hmc_stepsize = diag(c(0, rep(exp(state$ker_var$range_beta_ancillary[2]), 2)))
    # initializing position 
    q = 0*state$params$range_beta
    q[] = L_minus_one %*% (state$params$range_beta)
    # initializing momentum
    state$momenta$range_beta_ancillary = 
      renew_momentum(
        state$momenta$range_beta_ancillary, 
        kept_momentum = .9)
    p = state$momenta$range_beta_ancillary
    # Make a half step for momentum at the beginning
    dens_grad = 0*q
    dens_grad[] =  t(solve(L_minus_one)) %*% (  
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
            derivative_sandwiches(
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
    )
    p[] = p[] - (dens_grad) %*%hmc_stepsize/ 2
    
    n_hmc_steps = min(5, ceiling(sqrt(iter + iter_start)/3))
    for(hmc_step in seq_len(n_hmc_steps)){
      # Make a full step for the position
      q = q + p %*% hmc_stepsize
      new_range_beta = state$params$range_beta
      new_range_beta[] = solve(L_minus_one, (q))
      new_compressed_sparse_chol = 
        compute_sparse_chol(
          PP = hierarchical_model$range$PP,
          range_beta = new_range_beta, 
          vecchia_approx = vecchia_approx, 
          range_X = range_X, 
          matern_smoothness = hierarchical_model$matern_smoothness, 
          compute_derivative = TRUE, num_threads = num_threads
        )
      new_sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, new_compressed_sparse_chol)
      new_field = exp(.5 * state$params$field_log_var[1,1]) * 
        as.vector(Matrix::solve(new_sparse_chol, state$stuff$sparse_chol %*% (state$params$field/exp(.5 * state$params$field_log_var[1,1]))))
      # Make a half step for momentum at the end.
      dens_grad[] =   t(solve(L_minus_one)) %*% (  
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
              derivative_sandwiches(
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
      )
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
    
    state$ker_var$range_beta_ancillary[regime] = update_kernel(
      iter = iter, iter_start = iter_start, update_kernel_groupsize = 1, 
      kernel_value = state$ker_var$range_beta_ancillary[regime], mult = -.8
    )
    
    if(!is.nan(current_U-proposed_U+current_K- proposed_K))
    {
      if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
      {
        state$ker_var$range_beta_ancillary[regime] = update_kernel(
          iter = iter, iter_start = iter_start, update_kernel_groupsize = 1, 
          kernel_value = state$ker_var$range_beta_ancillary[regime], mult = 1
        )
        state$momenta$range_beta_ancillary = p
        state$params$field = new_field
        state$stuff$sparse_chol= new_sparse_chol
        state$stuff$compressed_chol = new_compressed_sparse_chol
        state$params$range_beta[] = new_range_beta
      }
    }
  }
  ###########################
  # Range beta (sufficient) #
  ###########################
  for(regime in seq_len(1 + hierarchical_model$anisotropic)){
    if(!hierarchical_model$anisotropic)hmc_stepsize =matrix(exp(state$ker_var$range_beta_sufficient))
    if(hierarchical_model$anisotropic & regime ==1)hmc_stepsize = diag(c(exp(state$ker_var$range_beta_sufficient[1]),0,0))
    if(hierarchical_model$anisotropic & regime ==2)hmc_stepsize = diag(c(0, rep(exp(state$ker_var$range_beta_sufficient[2]), 2)))
    q[] = L_minus_one %*% (state$params$range_beta)
    # initializing momentum
    state$momenta$range_beta_sufficient = 
      renew_momentum(
        state$momenta$range_beta_sufficient, 
        kept_momentum = .9)
    p = state$momenta$range_beta_sufficient
    # Make a half step for momentum at the beginning
    dens_grad = 0*q
    dens_grad[] =  t(solve(L_minus_one)) %*% (
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
            derivative_sandwiches(
              vecchia = state$stuff$compressed_chol, # derivative of the (unscaled) NNGP factor
              left_vector = as.vector(state$stuff$sparse_chol %*% (state$params$field/exp(.5 * state$params$field_log_var[1,1]))), # left vector = whitened latent field
              right_vector = state$params$field/exp(.5 * state$params$field_log_var[1,1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
              NNarray = vecchia_approx$NNarray, 
              sauce_determinant_chef = TRUE, 
              num_threads = num_threads  
            )
          ))  
      %*% range_reparam_mat
    )
    p[] = p[] - (dens_grad) %*% hmc_stepsize/ 2
    
    n_hmc_steps = min(5, ceiling(sqrt(iter + iter_start)/3))
    for(hmc_step in seq_len(n_hmc_steps)){
      # Make a full step for the position
      q = q + p %*% hmc_stepsize
      new_range_beta = state$params$range_beta      
      new_range_beta[] = solve(L_minus_one, (q))
      new_compressed_sparse_chol = 
        compute_sparse_chol(
          PP = hierarchical_model$range$PP,
          range_beta = new_range_beta, 
          vecchia_approx = vecchia_approx, 
          range_X = range_X, 
          matern_smoothness = hierarchical_model$matern_smoothness, 
          compute_derivative = TRUE, num_threads = num_threads
        )
      new_sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, new_compressed_sparse_chol)
      
      # Make a half step for momentum at the end.
      dens_grad = 0*q
      dens_grad[] = t(solve(L_minus_one)) %*% (
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
              derivative_sandwiches(
                vecchia = new_compressed_sparse_chol, # derivative of the (unscaled) NNGP factor
                left_vector = as.vector(new_sparse_chol %*% (state$params$field/exp(.5 * state$params$field_log_var[1,1]))), # left vector = whitened latent field
                right_vector = state$params$field/exp(.5 * state$params$field_log_var[1,1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                NNarray = vecchia_approx$NNarray, 
                sauce_determinant_chef = TRUE, 
                num_threads = num_threads  
              )
            )
        )  %*% range_reparam_mat
      )
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
    state$ker_var$range_beta_sufficient[regime] = update_kernel(
      iter = iter, iter_start = iter_start, update_kernel_groupsize = 1, 
      kernel_value = state$ker_var$range_beta_sufficient[regime], mult = -.8
    )
    
    if(!is.nan(current_U-proposed_U+current_K- proposed_K)) {
      if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K) {
        state$ker_var$range_beta_sufficient[regime] = update_kernel(
          iter = iter, iter_start = iter_start, update_kernel_groupsize = 1, 
          kernel_value = state$ker_var$range_beta_sufficient[regime], mult = 1
        )
        state$momenta$range_beta_sufficient = p
        state$stuff$sparse_chol= new_sparse_chol
        state$stuff$compressed_chol = new_compressed_sparse_chol
        state$params$range_beta[] = new_range_beta
      }
    }
  }
  return(list("state"=state, "p"=p))
}

update_variance_rangepp <- function(state, hierarchical_model, range_X, vecchia_approx, p, iter, iter_start, num_threads){
  n_range_log_scale_update = 2
  for(i in seq_len(n_range_log_scale_update)){
    # ancillary - sufficient
    q = state$params$range_log_scale + rnorm(1 + hierarchical_model$anisotropic, 0, exp(state$ker_var$range_log_scale_sufficient))
    
    new_range_beta = state$params$range_beta
    new_range_beta[-seq_len(range_X$n_regressors),] = 
      new_range_beta[-seq_len(range_X$n_regressors),] %*% 
      diag(exp(-.5 * state$params$range_log_scale[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic) %*% 
      diag(exp(.5 * q[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic)
    new_compressed_sparse_chol = 
      compute_sparse_chol(
        hierarchical_model$range$PP,
        range_beta = new_range_beta, 
        vecchia_approx = vecchia_approx, 
        range_X = range_X, 
        matern_smoothness = hierarchical_model$matern_smoothness, 
        compute_derivative = FALSE, num_threads = num_threads
      )
    new_sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, new_compressed_sparse_chol)
    
    
    current_U =
      (
        sum((state$params$range_log_scale - hierarchical_model$range$log_scale_bounds[1])^2)
        + .5* sum((state$stuff$sparse_chol %*% (state$params$field/exp(state$params$field_log_var[1,1] / 2)))^2)
        - sum(log(state$stuff$compressed_chol[1,1,]))
      )
    proposed_U =
      (
        sum((q - hierarchical_model$range$log_scale_bounds[1])^2)
        + .5* sum((new_sparse_chol %*% (state$params$field/exp(state$params$field_log_var[1,1] / 2)))^2)
        - sum(log(new_compressed_sparse_chol[1,1,]))
      )
    current_K = sum (state$momenta$range_beta_sufficient ^2) / 2
    proposed_K = sum(p^2) / 2
    
    state$ker_var$range_log_scale_sufficient = update_kernel(
      iter = iter, iter_start = iter_start, 
      update_kernel_groupsize = 3, kernel_value = state$ker_var$range_log_scale_sufficient, 
      mult = -.25
    )
    
    if (
      (all(q < hierarchical_model$range$log_scale_bounds[2] )) &
      (all(q > hierarchical_model$range$log_scale_bounds[1] ))
    ){
      if(!is.nan(current_U-proposed_U+current_K-proposed_K)){
        if(log(runif(1)) < (current_U-proposed_U+current_K-proposed_K))
        {
          
          state$ker_var$range_log_scale_sufficient = update_kernel(
            iter = iter, iter_start = iter_start, 
            update_kernel_groupsize = 3, kernel_value = state$ker_var$range_log_scale_sufficient, 
            mult = 1
          )
          state$momenta$range_beta_sufficient = p
          
          state$params$range_beta = new_range_beta
          state$params$range_log_scale[] = q
          
          state$stuff$sparse_chol= new_sparse_chol
          state$stuff$compressed_chol = new_compressed_sparse_chol
        }}}
  }
  # sufficient - sufficient ####
  for(i in seq_len(10))
  {
    q = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, .1)
    if(
      (all(q < hierarchical_model$range$log_scale_bounds[2] )) &
      (all(q > hierarchical_model$range$log_scale_bounds[1] )) &
      (
        -sum((q - hierarchical_model$range$log_scale_bounds[1])^2)
        +sum((state$params$range_log_scale - hierarchical_model$range$log_scale_bounds[1])^2)
        + beta_prior_log_dens(
          beta = state$params$range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
          beta0_mean = hierarchical_model$range$beta0_mean,
          beta0_var =  hierarchical_model$range$beta0_sd^2, 
          log_scale = q)
        - beta_prior_log_dens(
          beta = state$params$range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
          beta0_mean = hierarchical_model$range$beta0_mean,
          beta0_var =  hierarchical_model$range$beta0_sd^2, 
          log_scale = state$params$range_log_scale) 
        > log(runif(1))
      )
    )
    {
      state$params$range_log_scale[] = q
    }
  }
  # ancillary-ancillary ####
  
  for(i in seq_len(n_range_log_scale_update)){
    q = state$params$range_log_scale + rnorm(1 + hierarchical_model$anisotropic, 0, exp(state$ker_var$range_log_scale_ancillary))
    
    new_range_beta = state$params$range_beta
    new_range_beta[-seq_len(range_X$n_regressors),] = 
      new_range_beta[-seq_len(range_X$n_regressors),] %*% 
      diag(exp(-.5 * state$params$range_log_scale[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic) %*%
      diag(exp(.5 * q[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic)
    new_compressed_sparse_chol = 
      compute_sparse_chol(
        hierarchical_model$range$PP,
        range_beta = new_range_beta, 
        vecchia_approx = vecchia_approx, 
        range_X = range_X, 
        matern_smoothness = hierarchical_model$matern_smoothness, 
        compute_derivative = FALSE, num_threads = num_threads
      )
    new_sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, new_compressed_sparse_chol)
    new_field = as.vector(Matrix::solve(new_sparse_chol, state$stuff$sparse_chol %*% (state$params$field)))
    
    current_U =
      (
        sum((state$params$range_log_scale - hierarchical_model$range$log_scale_bounds[1])^2)
        + .5 * sum((state$stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$stuff$noise_var) # observation ll
      )
    proposed_U =
      (
        sum((q - hierarchical_model$range$log_scale_bounds[1])^2)
        + .5 * sum((state$stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$stuff$noise_var) # observation ll
      )
    state$ker_var$range_log_scale_ancillary = update_kernel(
      iter = iter, iter_start = iter_start, 
      update_kernel_groupsize = 2, kernel_value = state$ker_var$range_log_scale_ancillary, 
      mult = -.25
    )
    
    if (
      (all(q < hierarchical_model$range$log_scale_bounds[2] )) &
      (all(q > hierarchical_model$range$log_scale_bounds[1] ))
    )
    {
      if(!is.nan(current_U-proposed_U)){
        if(log(runif(1)) < (current_U-proposed_U)){
          state$ker_var$range_log_scale_ancillary = update_kernel(
            iter = iter, iter_start = iter_start, 
            update_kernel_groupsize = 2, kernel_value = state$ker_var$range_log_scale_ancillary, 
            mult = 1
          )
          
          state$params$range_beta = new_range_beta
          state$params$range_log_scale[] = q
          state$params$field = new_field
          
          state$stuff$sparse_chol= new_sparse_chol
          state$stuff$compressed_chol = new_compressed_sparse_chol
        }}}
  }
  # sufficient - sufficient 
  for(i in seq_len(10))
  {
    q = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, .1)
    if(
      (all(q < hierarchical_model$range$log_scale_bounds[2] )) &
      (all(q > hierarchical_model$range$log_scale_bounds[1] )) &
      (
        -sum((q - hierarchical_model$range$log_scale_bounds[1])^2)
        +sum((state$params$range_log_scale - hierarchical_model$range$log_scale_bounds[1])^2)
        
        + beta_prior_log_dens(
          beta = state$params$range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
          beta0_mean = hierarchical_model$range$beta0_mean,
          beta0_var =  hierarchical_model$range$beta0_sd^2, 
          log_scale = q)
        - beta_prior_log_dens(
          beta = state$params$range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
          beta0_mean = hierarchical_model$range$beta0_mean,
          beta0_var =  hierarchical_model$range$beta0_sd^2, 
          log_scale = state$params$range_log_scale) 
        > log(runif(1))
      )
    )
    {
      state$params$range_log_scale[] = q
    }
  }
  state$stuff$compressed_chol = 
    compute_sparse_chol(
      hierarchical_model$range$PP,
      range_beta = state$params$range_beta, 
      vecchia_approx = vecchia_approx, 
      range_X = range_X, 
      matern_smoothness = hierarchical_model$matern_smoothness, 
      compute_derivative = TRUE, num_threads = num_threads
    )
  state$stuff$sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, state$stuff$compressed_chol)
  return(state)
}

update_noise_beta <- function(state, noise, noise_X, vecchia_approx, iter, iter_start){
  L_minus_one =  solve(t(chol(
    solve(noise_X$crossprod_X)/max(solve(noise_X$crossprod_X))
  ))) 
  # VEWY IMPOWTANT don't remove or comment
  squared_residuals = as.matrix(state$stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2
  # HMC update
  q = L_minus_one %*% state$params$noise_beta
  state$momenta$noise_beta = renew_momentum(state$momenta$noise_beta)
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
  p = p - exp_noise_mala * solve(t(L_minus_one), dens_grad) / 2
  
  n_hmc_steps = min(5, ceiling(sqrt(iter + iter_start)/3))
  for(hmc_step in seq_len(n_hmc_steps)){
    # Make a full step for the position
    q = q + exp_noise_mala * p
    new_noise_beta = solve(L_minus_one, q)
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
    p = p - exp_noise_mala * solve(t(L_minus_one), dens_grad) / (1 + (hmc_step == n_hmc_steps))
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
  
  state$ker_var$noise_beta_mala = update_kernel(iter = iter, iter_start = iter_start, update_kernel_groupsize = 1, kernel_value = state$ker_var$noise_beta_mala, mult = -.7)
  if(!is.nan(current_U-proposed_U+current_K- proposed_K)) {
    if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K) {
      state$ker_var$noise_beta_mala = update_kernel(iter = iter, iter_start = iter_start, update_kernel_groupsize = 1, kernel_value = state$ker_var$noise_beta_mala, mult = 1)
      state$momenta$noise_beta = p
      state$params$noise_beta[] = new_noise_beta
      state$stuff$noise_var = new_noise_var
    }
  }
  
  return(list("state" = state, "squared_residuals" = squared_residuals))
}

update_noise_log_scale <- function(state, noise, noise_X, vecchia_approx, squared_residuals, iter, iter_start){
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
        if(
          (new_noise_log_scale > noise$log_scale_bounds[1])&
          (new_noise_log_scale < noise$log_scale_bounds[2])
        ) {
          state$params$noise_log_scale[] = new_noise_log_scale
          state$params$noise_beta = new_noise_beta 
          state$stuff$noise_var = new_noise_var
          state$ker_var$noise_log_scale = update_kernel(
            iter = iter, iter_start = iter_start, update_kernel_groupsize = 1, 
            kernel_value = state$ker_var$noise_log_scale, mult = 1)
        }
      }
    }
    state$ker_var$noise_log_scale = update_kernel(
      iter = iter, iter_start = iter_start, update_kernel_groupsize = 1, 
      kernel_value = state$ker_var$noise_log_scale, mult = -.25)
  }
  
  # sufficient -- sufficient ####
  for(i in seq_len(10)) {
    new_noise_log_scale = state$params$noise_log_scale + rnorm(1, 0, .1)
    if(new_noise_log_scale > noise$log_scale_bounds[1]
      && new_noise_log_scale < noise$log_scale_bounds[2]) {
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