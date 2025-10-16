# TO DO : Range parametrization changed

# With only one HMC step in range 

# numerically stable log(exp(x) + exp(y)) with 0 << x << y or 0 << y << x
# let's do: log(exp(x) + exp(y)) with 0 << x << y 
# log(exp(x) + exp(y))  = log(exp(y)(1 + exp(x) / exp(y)))
# = log(exp(y)) + log((1 + exp(x) / exp(y)))
# = y + log((1 + exp(x - y)))
log_of_exp_sum = function(x,y)return(max(x,y) + log((1 + exp(-abs(x - y)))))



renew_momentum  = function(momentum, kept_momentum = .9){
  if(kept_momentum <0 | kept_momentum>1)stop("kept_momentum must be between 0 and 1")
  momentum[] = sqrt(kept_momentum)*momentum[] + sqrt(1-kept_momentum)*rnorm(length(momentum[])) 
  momentum
}

update_kernel = function(
    iter, iter_start, update_kernel_groupsize, 
    kernel_value, mult
){
  idx = (((iter + iter_start-1) %/% update_kernel_groupsize +1)%%length(kernel_value))+1
  kernel_value [idx] =
    kernel_value [idx] + 
    length(kernel_value) * mult/sqrt(10 + iter + iter_start)
  kernel_value [idx] = max(kernel_value [idx], -8)
  kernel_value [idx] = min(kernel_value [idx], 4)
  kernel_value
}



#list2env(mcmc_nngp_list, envir = environment())
#state = mcmc_nngp_list$states$chain_1; n_iterations_update  =100; num_threads = 10; iter_start = 0; seed = 1; iter=1
#
##' @export
#MCMC = function(
    #    covariates,
#    observed_field,
#    hierarchical_model,
#    vecchia_approx,
#    state,
#    n_iterations_update = 100, 
#    num_threads = 1, 
#    iter_start,
#    seed=123 
#)
#{
message(paste("Starting MCMC chain for Nonstat NNGP model at iteration", iter_start, "for", n_iterations_update, "iterations"))
set.seed(seed)
#########################################
# Initializing chain storage structures #
#########################################
# this part re-creates a small portion of the $records objects of each chain. It fills it with chain state during the run, and then updates each chain with the new values
#params_records = list()

#par(mfrow = c(2, 1))

t1 = Sys.time()
for(iter in seq(iter, n_iterations_update)){
  if((iter + iter_start)/50 ==(iter + iter_start) %/% 50)print(paste("iteration", (iter + iter_start)))
  
  ###########################
  # Regression coefficients #
  ###########################
  # centered parametrization of latent field
  
  
  
  
  beta_covmat = solve(crossprod(covariates$X$X/stuff$noise_var, covariates$X$X))
  if(all(!is.infinite(beta_covmat) & !is.nan(beta_covmat)))
  {
    if(all(eigen(beta_covmat)$val >0))
    {
      beta_mean = beta_covmat %*% crossprod(covariates$X$X, ((observed_field-params$field[vecchia_approx$locs_match]) / stuff$noise_var))
      params$beta[] = (beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))[]
    }}
  
  # un-centered parametrization of latent field
  centered_field = as.vector(params$field + covariates$X$X_locs%*%matrix(params$beta[covariates$X$which_locs], ncol = 1))
  sparse_chol_X = as.matrix(stuff$sparse_chol %*% (covariates$X$X_locs))/exp(.5 * params$field_log_var)
  beta_precision = crossprod(x = sparse_chol_X, y = sparse_chol_X)
  beta_covmat = solve(beta_precision, tol = min(rcond(beta_precision),.Machine$double.eps))
  if(all(!is.infinite(beta_covmat) & !is.nan(beta_covmat)))
  {
    if(all(eigen(beta_covmat)$d >0))
    {
      beta_mean =  c(as.vector(stuff$sparse_chol %*% (centered_field/exp(.5 * params$field_log_var)))  %*% sparse_chol_X %*% beta_covmat)
      params$beta[covariates$X$which_locs]   = as.vector(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
      params$field = centered_field - as.vector(covariates$X$X_locs %*% matrix(params$beta[covariates$X$which_locs], ncol = 1))
    }}
  
  # updating stuff 
  stuff$lm_fit[]       = as.vector(covariates$X$X %*% params$beta)
  stuff$lm_residuals = observed_field - stuff$lm_fit
  
  ################
  # Latent field #
  ################
  
  cluster_idx = 1
  chosen_locs_partition_idx = 1 + iter%%ncol(vecchia_approx$locs_partition)
  locs_partition = vecchia_approx$locs_partition[,chosen_locs_partition_idx]
  precision_from_obs = vecchia_approx$locs_match_matrix %*% (1/stuff$noise_var)
  mean_from_obs = as.vector(vecchia_approx$locs_match_matrix %*% ((observed_field - stuff$lm_fit)/stuff$noise_var))
  
  t1 = Sys.time()
  chol_list = parallel::mclapply(mc.cores = num_threads,
                                 unique( locs_partition), 
                                 function(cluster_idx){
                                   selected_idx = which(locs_partition == cluster_idx)
                                   posterior_precision_subset = Matrix::crossprod(stuff$sparse_chol[,selected_idx])
                                   posterior_precision_subset = 
                                     Matrix::Diagonal(length(selected_idx), 1/exp(.5 * params$field_log_var)) %*% 
                                     posterior_precision_subset %*% Matrix::Diagonal(length(selected_idx), 1/exp(.5 * params$field_log_var))
                                   Matrix::diag(posterior_precision_subset) = Matrix::diag(posterior_precision_subset) + as.vector(precision_from_obs[selected_idx])
                                   posterior_precision_subset = Matrix::expand(Matrix::Cholesky(posterior_precision_subset))
                                   return(posterior_precision_subset)
                                 }
  )
  Sys.time()-t1
  
  t1 = Sys.time()
  for(pass in seq(1)){
    for(cluster_idx in unique(locs_partition)){
      selected_idx = which(locs_partition==cluster_idx)
      additional_mean_from_field =
        as.vector(
          (1/exp(.5 * params$field_log_var)) * 
            Matrix::crossprod( 
              stuff$sparse_chol[,selected_idx], 
              (stuff$sparse_chol %*% 
                 ((params$field*(locs_partition!=cluster_idx))/exp(.5 * params$field_log_var))
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
  Sys.time()-t1
  
  
  if(iter%/%20 == iter /20){
    par(mfrow = c(1,2))
    plot_pointillist_painting(vecchia_approx$locs, params$field, main = "sampled")
    plot_pointillist_painting(vecchia_approx$locs, fake_data$latent_field, main = "true")
    par(mar = rep(2, 4))
    par(mfrow = c(3,3))
    for(i in seq(nrow(params$range_beta))){
      for(j in seq(ncol(params$range_beta))){
        plot(c(rbind(coeff_list$range_X_coeff, coeff_list$range_PP_coeff)[i,j], sapply(params_records, function(x)x$range_beta[i,j])[seq(iter/5, iter)]), ylab = "", main  = paste(i, j))
        abline(h = rbind(coeff_list$range_X_coeff, coeff_list$range_PP_coeff)[i,j])
      }}
    range_log_scale_samples = rbind(range_log_scale, t(sapply(params_records, function(x)x$range_log_scale))[seq(iter/5, iter-1),])
    plot(range_log_scale_samples[,1], main = "range log scale")
    plot(range_log_scale_samples[,2], main = "aniso log scale")
  }
  
  
  par(mfrow = c(1,1 ))
  print(iter)
  
  ###############
  # Range beta  #
  ###############
  
  range_reparam_mat = matrix(1)
  if(hierarchical_model$anisotropic){
    range_reparam_mat = matrix(c(2, 2,  0, 
                                 2, -2, 0,
                                 0, 0,  2*sqrt(2)), 3)
  }
  L_minus_one =  solve(t(chol(
    solve(covariates$range_X$crossprod_X_locs)/max(solve(covariates$range_X$crossprod_X_locs))
  )))
  ##########################
  # Range beta (ancillary) #
  ##########################
  
  hmc_stepsize = diag(exp(ker_var$range_beta_ancillary)[c(1, rep(2, 2*hierarchical_model$anisotropic))], 1 + 2*hierarchical_model$anisotropic)
  # initializing position 
  q = 0*params$range_beta
  q[] = L_minus_one %*% (params$range_beta)
  # initializing momentum
  momenta$range_beta_ancillary = 
    renew_momentum(
      momenta$range_beta_ancillary, 
      kept_momentum = .9)
  p = momenta$range_beta_ancillary
  # Make a half step for momentum at the beginning
  dens_grad = 0*q
  dens_grad[] =  t(solve(L_minus_one)) %*% (  
    -beta_prior_log_dens_derivative(
      beta = params$range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
      beta0_mean = hierarchical_model$range$beta0_mean,
      beta0_var =  hierarchical_model$range$beta0_sd^2, 
      params$range_log_scale) # normal prior
    + X_PP_crossprod(
      X = covariates$range_X$X_locs,  vecchia_approx = vecchia_approx, permutate_PP_to_obs = F,
      PP = hierarchical_model$range$PP,
      Y = # Jacobian of range field wrt range_beta
        t(
          # natural gradient of obs likelihood wrt range field
          derivative_sandwiches(
            vecchia = stuff$compressed_chol, 
            left_vector = as.vector(
              Matrix::solve(
                Matrix::t(stuff$sparse_chol), 
                - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                              ((params$field[vecchia_approx$locs_match] - stuff$lm_residuals) / stuff$noise))
                * exp(.5 * params$field_log_var) # part of sparse chol
              )), 
            right_vector = params$field/exp(.5 * params$field_log_var), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
            NNarray = vecchia_approx$NNarray, 
            sauce_determinant_chef = F, 
            num_threads = num_threads
          )
        )) %*% range_reparam_mat
  )
  p[] = p[] - (dens_grad) %*%hmc_stepsize/ 2
  
  
  ###    #testing the gradient
  ###    derivative_test = 0*q
  ###    for(i in seq(nrow(params$range_beta))){
  ###      for(j in seq(ncol(params$range_beta))){
  ###        q_ = q
  ###        q_[i,j] = q_[i,j] +  0.0000001
  ###        range_beta_ = 0*params$range_beta
  ###        range_beta_[] = solve(L_minus_one, (q_))
  ###        sparse_chol_ = decompress_chol(
  ###          vecchia_approx = vecchia_approx, 
  ###          compressed_sparse_chol = compute_sparse_chol(
  ###            range_beta = range_beta_, vecchia_approx = vecchia_approx, 
  ###            range_X = covariates$range_X, PP = hierarchical_model$range$PP,
  ###            matern_smoothness = hierarchical_model$matern_smoothness, compute_derivative = F, num_threads = 10
  ###          )
  ###        ) 
  ###        field_ = exp(.5 * params$field_log_var) * as.vector(Matrix::solve(sparse_chol_, stuff$sparse_chol %*% (params$field/exp(.5 * params$field_log_var))))
  ###        
  ###    derivative_test[i,j] = 
  ###        (1/ 0.0000001)*(
  ###          beta_prior_log_dens(beta = params$range_beta, 
  ###                              n_PP = hierarchical_model$range$PP$n_knots, 
  ###                              beta0_mean = hierarchical_model$range$beta0_mean, 
  ###                              beta0_var =  hierarchical_model$range$beta0_sd^2, 
  ###                              params$range_log_scale)
  ###          - beta_prior_log_dens(beta = range_beta_, 
  ###                                n_PP = hierarchical_model$range$PP$n_knots, 
  ###                                beta0_mean = hierarchical_model$range$beta0_mean, 
  ###                                beta0_var =  hierarchical_model$range$beta0_sd^2, 
  ###                                params$range_log_scale)
  ###          - .5 * sum( (stuff$lm_residuals -  params$field[vecchia_approx$locs_match])^2/stuff$noise) # observation ll
  ###          + .5 * sum( (stuff$lm_residuals -  field_      [vecchia_approx$locs_match])^2/stuff$noise) # observation ll
  ###        ) 
  ###      }
  ###    }
  ###    par(mfrow = c(1,3))
  ###    boxplot(c(derivative_test)/c(dens_grad))
  ###    boxplot(c(derivative_test)-c(dens_grad))
  ###    plot(c(derivative_test), c(dens_grad))
  ###    abline(a=0, b=1)
  
  n_hmc_steps = 1
  for(hmc_step in seq(n_hmc_steps)){
  # Make a full step for the position
  q = q + p %*% hmc_stepsize
  new_range_beta = params$range_beta
  new_range_beta[] = solve(L_minus_one, (q))
  new_compressed_sparse_chol = 
    compute_sparse_chol(
      PP = hierarchical_model$range$PP,
      range_beta = new_range_beta, 
      vecchia_approx = vecchia_approx, 
      range_X = covariates$range_X, 
      matern_smoothness = hierarchical_model$matern_smoothness, 
      compute_derivative = T, num_threads = num_threads
    )
  new_sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, new_compressed_sparse_chol)
  new_field = exp(.5 * params$field_log_var) * as.vector(Matrix::solve(new_sparse_chol, stuff$sparse_chol %*% (params$field/exp(.5 * params$field_log_var))))
  # Make a half step for momentum at the end.
  dens_grad[] =   t(solve(L_minus_one)) %*% (  
    -beta_prior_log_dens_derivative(
      beta = new_range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
      beta0_mean = hierarchical_model$range$beta0_mean,
      beta0_var =  hierarchical_model$range$beta0_sd^2, 
      params$range_log_scale) # normal prior
    + X_PP_crossprod(
      X = covariates$range_X$X_locs, vecchia_approx = vecchia_approx, 
      permutate_PP_to_obs = F,
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
                              ((new_field[vecchia_approx$locs_match] - stuff$lm_residuals) / stuff$noise))
                * exp(.5 * params$field_log_var) # part of sparse chol
              )), 
            right_vector = new_field/exp(.5 * params$field_log_var), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
            NNarray = vecchia_approx$NNarray, 
            sauce_determinant_chef = F, 
            num_threads = num_threads
          )
        )) %*% range_reparam_mat 
  )
  p[] = p[] - (dens_grad)%*%hmc_stepsize/ (1+  (hmc_step==n_hmc_steps))
  }
  ###       #testing the gradient
  ###       derivative_test = 0*q
  ###       for(i in seq(nrow(params$range_beta))){
  ###         for(j in seq(ncol(params$range_beta))){
  ###           q_ = q
  ###           q_[i,j] = q_[i,j] + .00001
  ###           range_beta_ = 0*params$range_beta
  ###           range_beta_[] = solve(L_minus_one, (q_))
  ###           sparse_chol_ = decompress_chol(
  ###             vecchia_approx = vecchia_approx, 
  ###             compressed_sparse_chol = compute_sparse_chol(
  ###               range_beta = range_beta_, vecchia_approx = vecchia_approx, 
  ###               range_X = covariates$range_X, PP = hierarchical_model$range$PP,
  ###               matern_smoothness = hierarchical_model$matern_smoothness, compute_derivative = F, num_threads = 10
  ###             )
  ###           ) 
  ###           field_ = exp(.5 * params$field_log_var) * as.vector(Matrix::solve(sparse_chol_, stuff$sparse_chol %*% (params$field/exp(.5 * params$field_log_var))))
  ###           
  ###       derivative_test[i,j] = 
  ###           100000*(
  ###             beta_prior_log_dens(beta = new_range_beta, 
  ###                                 n_PP = hierarchical_model$range$PP$n_knots, 
  ###                                 beta0_mean = hierarchical_model$range$beta0_mean, 
  ###                                 beta0_var =  hierarchical_model$range$beta0_sd^2, 
  ###                                 params$range_log_scale)
  ###             - beta_prior_log_dens(beta = range_beta_, 
  ###                                   n_PP = hierarchical_model$range$PP$n_knots, 
  ###                                   beta0_mean = hierarchical_model$range$beta0_mean, 
  ###                                   beta0_var =  hierarchical_model$range$beta0_sd^2, 
  ###                                   params$range_log_scale)
  ###             - .5 * sum( (stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/stuff$noise) # observation ll
  ###             + .5 * sum( (stuff$lm_residuals -  field_   [vecchia_approx$locs_match])^2/stuff$noise) # observation ll
  ###           ) 
  ###         }
  ###       }
  ###       par(mfrow = c(1,3))
  ###       boxplot(c(derivative_test)/c(dens_grad))
  ###       boxplot(c(derivative_test)-c(dens_grad))
  ###       plot(c(derivative_test), c(dens_grad))
  ###       abline(a=0, b=1)
  
  
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_K = sum (momenta$range_beta_ancillary ^2) / 2
  proposed_K = sum(p^2) / 2
  
  current_U =
    (
      - beta_prior_log_dens(beta = params$range_beta, 
                            n_PP = hierarchical_model$range$PP$n_knots, 
                            beta0_mean = hierarchical_model$range$beta0_mean, 
                            beta0_var =  hierarchical_model$range$beta0_sd^2, 
                            params$range_log_scale) # normal prior
      + .5 * sum((stuff$lm_residuals -  params$field[vecchia_approx$locs_match])^2/stuff$noise) # observation ll
    )
  proposed_U =
    (
      - beta_prior_log_dens(beta = new_range_beta, 
                            n_PP = hierarchical_model$range$PP$n_knots, 
                            beta0_mean = hierarchical_model$range$beta0_mean, 
                            beta0_var =  hierarchical_model$range$beta0_sd^2, 
                            params$range_log_scale) # normal prior
      + .5 * sum((stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/stuff$noise) # observation ll
    )
  
  current_U-proposed_U
  current_K- proposed_K
  
  ker_var$range_beta_ancillary = update_kernel(
    iter = iter, iter_start = iter_start, update_kernel_groupsize = 2, 
    kernel_value = ker_var$range_beta_ancillary, mult = -.7
  )
  
  if(!is.nan(current_U-proposed_U+current_K- proposed_K))
  {
    if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
    {
      ker_var$range_beta_ancillary = update_kernel(
        iter = iter, iter_start = iter_start, update_kernel_groupsize = 2, 
        kernel_value = ker_var$range_beta_ancillary, mult = 1
      )
      print("tatato ancillary !")
      momenta$range_beta_ancillary = p
      params$field = new_field
      stuff$sparse_chol= new_sparse_chol
      stuff$compressed_chol = new_compressed_sparse_chol
      params$range_beta[] = new_range_beta
    }
  }
  
  
  ###########################
  # Range beta (sufficient) #
  ###########################
  hmc_stepsize = diag(exp(ker_var$range_beta_sufficient)[c(1, rep(2, 2*hierarchical_model$anisotropic))], 1 + 2*hierarchical_model$anisotropic)
  # initializing position 
  q = 0*params$range_beta
  q[] = L_minus_one %*% (params$range_beta)
  # initializing momentum
  momenta$range_beta_sufficient = 
    renew_momentum(
      momenta$range_beta_sufficient, 
      kept_momentum = .9)
  p = momenta$range_beta_sufficient
  # Make a half step for momentum at the beginning
  dens_grad = 0*q
  dens_grad[] =  t(solve(L_minus_one)) %*% (
    - beta_prior_log_dens_derivative(
      beta = params$range_beta, 
      n_PP = hierarchical_model$range$PP$n_knots, 
      beta0_mean = hierarchical_model$range$beta0_mean, 
      beta0_var =  hierarchical_model$range$beta0_sd^2, 
      log_scale = params$range_log_scale) # normal prior
    # normal prior derivative                
    + X_PP_crossprod(
      X = covariates$range_X$X_locs, PP = hierarchical_model$range$PP,
      permutate_PP_to_obs = F, 
      vecchia_approx = vecchia_approx,
      Y = # Jacobian of range field wrt range_beta
        t(# natural gradient of obs likelihood wrt range field
          derivative_sandwiches(
            vecchia = stuff$compressed_chol, # derivative of the (unscaled) NNGP factor
            left_vector = as.vector(stuff$sparse_chol %*% (params$field/exp(.5 * params$field_log_var))), # left vector = whitened latent field
            right_vector = params$field/exp(.5 * params$field_log_var), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
            NNarray = vecchia_approx$NNarray, 
            sauce_determinant_chef = T, 
            num_threads = num_threads  
          )
        ))  
    %*% range_reparam_mat
  )
  p[] = p[] - (dens_grad) %*% hmc_stepsize/ 2
  
  ###    ###testing the gradient
  ###     derivative_test = 0*q
  ###     for(i in seq(nrow(params$range_beta))){
  ###       for(j in seq(ncol(params$range_beta))){
  ###         q_ = q
  ###         q_[i,j] = q_[i,j] + .00001
  ###         range_beta_ = 0*params$range_beta
  ###         range_beta_[] = solve(L_minus_one, (q_))
  ###         sparse_chol_ = decompress_chol(
  ###           vecchia_approx = vecchia_approx, 
  ###           compressed_sparse_chol = compute_sparse_chol(
  ###             range_beta = range_beta_, vecchia_approx = vecchia_approx, 
  ###             range_X = covariates$range_X, PP = hierarchical_model$range$PP,
  ###             matern_smoothness = hierarchical_model$matern_smoothness, compute_derivative = F, num_threads = 10
  ###           )
  ###         ) 
  ###         derivative_test[i,j] = 
  ###           100000*(
  ###             beta_prior_log_dens(beta = params$range_beta, 
  ###                                   n_PP = hierarchical_model$range$PP$n_knots, 
  ###                                   beta0_mean = hierarchical_model$range$beta0_mean, 
  ###                                   beta0_var =  hierarchical_model$range$beta0_sd^2, 
  ###                                   params$range_log_scale) - 
  ###              beta_prior_log_dens(beta = range_beta_, 
  ###                                   n_PP = hierarchical_model$range$PP$n_knots, 
  ###                                   beta0_mean = hierarchical_model$range$beta0_mean, 
  ###                                   beta0_var =  hierarchical_model$range$beta0_sd^2, 
  ###                                   params$range_log_scale) +
  ###             (+ .5* sum((sparse_chol_ %*% (params$field/exp(.5 * params$field_log_var)))^2)
  ###              - sum(log(Matrix::diag(sparse_chol_))))-
  ###               (
  ###                 + .5* sum((stuff$sparse_chol %*% (params$field/exp(.5 * params$field_log_var)))^2)
  ###                 - sum(log(Matrix::diag(stuff$sparse_chol)))
  ###               )
  ###           ) 
  ###       }
  ###     }
  ###     par(mfrow = c(1,3))
  ###     boxplot(c(derivative_test)/c(dens_grad))
  ###     boxplot(c(derivative_test)-c(dens_grad))
  ###     plot(c(derivative_test), c(dens_grad))
  ###     abline(a=0, b=1)
  
  n_hmc_steps = 1
  for(hmc_step in seq(n_hmc_steps)){
  # Make a full step for the position
  q = q + p %*% hmc_stepsize
  new_range_beta = params$range_beta      
  new_range_beta[] = solve(L_minus_one, (q))
  new_compressed_sparse_chol = 
    compute_sparse_chol(
      PP = hierarchical_model$range$PP,
      range_beta = new_range_beta, 
      vecchia_approx = vecchia_approx, 
      range_X = covariates$range_X, 
      matern_smoothness = hierarchical_model$matern_smoothness, 
      compute_derivative = T, num_threads = num_threads
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
      log_scale = params$range_log_scale) # normal prior
    # normal prior derivative                
    + X_PP_crossprod(
      X = covariates$range_X$X_locs, PP = hierarchical_model$range$PP,
      permutate_PP_to_obs = F, 
      vecchia_approx = vecchia_approx,
      Y = # Jacobian of range field wrt range_beta
        t(# natural gradient of obs likelihood wrt range field
          derivative_sandwiches(
            vecchia = new_compressed_sparse_chol, # derivative of the (unscaled) NNGP factor
            left_vector = as.vector(new_sparse_chol %*% (params$field/exp(.5 * params$field_log_var))), # left vector = whitened latent field
            right_vector = params$field/exp(.5 * params$field_log_var), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
            NNarray = vecchia_approx$NNarray, 
            sauce_determinant_chef = T, 
            num_threads = num_threads  
          )
        )
    )  %*% range_reparam_mat
  )
  # updating momentum
  p[] = p[] - (dens_grad) %*% hmc_stepsize/ (1 + (hmc_step == n_hmc_steps))
}
  ####testing the gradient
  #     derivative_test = 0*q
  #     for(i in seq(nrow(params$range_beta))){
  #       for(j in seq(ncol(params$range_beta))){
  #         q_ = q
  #         q_[i,j] = q_[i,j] + .0000001
  #         range_beta_ = 0*params$range_beta
  #         range_beta_[] = solve(L_minus_one, (q_))
  #         sparse_chol_ = decompress_chol(
  #           vecchia_approx = vecchia_approx, 
  #           compressed_sparse_chol = compute_sparse_chol(
  #             range_beta = range_beta_, vecchia_approx = vecchia_approx, 
  #             range_X = covariates$range_X, PP = hierarchical_model$range$PP,
  #             matern_smoothness = hierarchical_model$matern_smoothness, compute_derivative = F, num_threads = 10
  #           )
  #         ) 
  #         derivative_test[i,j] = 
  #           10000000*(
  #             beta_prior_log_dens(beta = new_range_beta, 
  #                                 n_PP = hierarchical_model$range$PP$n_knots, 
  #                                 beta0_mean = hierarchical_model$range$beta0_mean, 
  #                                 beta0_var =  hierarchical_model$range$beta0_sd^2, 
  #                                 params$range_log_scale) - 
  #               beta_prior_log_dens(beta = range_beta_, 
  #                                   n_PP = hierarchical_model$range$PP$n_knots, 
  #                                   beta0_mean = hierarchical_model$range$beta0_mean, 
  #                                   beta0_var =  hierarchical_model$range$beta0_sd^2, 
  #                                   params$range_log_scale) 
  #             +
  #               (+ .5* sum((sparse_chol_ %*% (params$field/exp(.5 * params$field_log_var)))^2)
  #                - sum(log(Matrix::diag(sparse_chol_))))-
  #               (
  #                 + .5* sum((new_sparse_chol %*% (params$field/exp(.5 * params$field_log_var)))^2)
  #                 - sum(log(Matrix::diag(new_sparse_chol)))
  #               )
  #           ) 
  #       }
  #     }
  #     par(mfrow = c(1,3))
  #     boxplot(c(derivative_test)/c(dens_grad))
  #     boxplot(c(derivative_test)-c(dens_grad))
  #     plot(c(derivative_test), c(dens_grad))
  #     abline(a=0, b=1)
  
  # metropolis step
  current_K = sum (momenta$range_beta_sufficient ^2) / 2
  proposed_K = sum(p^2) / 2
  current_U =
    (
      - beta_prior_log_dens(beta = params$range_beta, 
                            n_PP = hierarchical_model$range$PP$n_knots, 
                            beta0_mean = hierarchical_model$range$beta0_mean, 
                            beta0_var =  hierarchical_model$range$beta0_sd^2, 
                            log_scale = params$range_log_scale)
      # normal prior 
      + .5* sum((stuff$sparse_chol %*% (params$field/exp(.5 * params$field_log_var)))^2)
      - sum(log(Matrix::diag(stuff$sparse_chol)))
    )
  proposed_U =
    (
      - beta_prior_log_dens(beta = new_range_beta, 
                            n_PP = hierarchical_model$range$PP$n_knots, 
                            beta0_mean = hierarchical_model$range$beta0_mean, 
                            beta0_var =  hierarchical_model$range$beta0_sd^2, 
                            params$range_log_scale)
      # normal prior 
      + .5* sum((new_sparse_chol %*% (params$field/exp(.5 * params$field_log_var)))^2)
      - sum(log(Matrix::diag(new_sparse_chol)))
    )
  
  current_U-proposed_U
  current_K- proposed_K
  
  ker_var$range_beta_sufficient = update_kernel(
    iter = iter, iter_start = iter_start, update_kernel_groupsize = 2, 
    kernel_value = ker_var$range_beta_sufficient, mult = -.7
  )
  if(!is.nan(current_U-proposed_U+current_K- proposed_K))
  {
    if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
    {
      ker_var$range_beta_sufficient = update_kernel(
        iter = iter, iter_start = iter_start, update_kernel_groupsize = 2, 
        kernel_value = ker_var$range_beta_sufficient, mult = 1
      )
      print("tatato sufficient!")
      momenta$range_beta_sufficient = p
      stuff$sparse_chol= new_sparse_chol
      stuff$compressed_chol = new_compressed_sparse_chol
      params$range_beta[] = new_range_beta
    }
  }
  
  print(ker_var$range_beta_sufficient)
  print(ker_var$range_beta_ancillary)
  
  
  #############################
  # Variance of the  range PP #
  #############################
  
  if(!is.null(hierarchical_model$range$PP)){
    
    # ancillary - sufficient
    q = params$range_log_scale + rnorm(1, 0, exp(ker_var$range_log_scale_sufficient))
    
    new_range_beta = params$range_beta
    new_range_beta[-seq(covariates$range_X$n_regressors),] = new_range_beta[-seq(covariates$range_X$n_regressors),] %*% 
      diag(exp(-.5 * params$range_log_scale[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic) %*% 
      diag(exp(.5 * q[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic)
    new_compressed_sparse_chol = 
      compute_sparse_chol(
        hierarchical_model$range$PP,
        range_beta = new_range_beta, 
        vecchia_approx = vecchia_approx, 
        range_X = covariates$range_X, 
        matern_smoothness = hierarchical_model$matern_smoothness, 
        compute_derivative = T, num_threads = num_threads
      )
    new_sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, new_compressed_sparse_chol)
    
    
    current_U =
      (
        3 * sum(sqrt(max(0, params$range_log_scale - hierarchical_model$range$log_scale_bounds[1])))
        + .5* sum((stuff$sparse_chol %*% (params$field/exp(params$field_log_var / 2)))^2)
        - sum(log(stuff$compressed_chol[1,1,]))
      )
    proposed_U =
      (
        3 * sum(sqrt(max(0, q - hierarchical_model$range$log_scale_bounds[1])))
        + .5* sum((new_sparse_chol %*% (params$field/exp(params$field_log_var / 2)))^2)
        - sum(log(new_compressed_sparse_chol[1,1,]))
      )
    current_K = sum (momenta$range_beta_sufficient ^2) / 2
    proposed_K = sum(p^2) / 2
    
    ker_var$range_log_scale_sufficient = update_kernel(
      iter = iter, iter_start = iter_start, 
      update_kernel_groupsize = 1, kernel_value = ker_var$range_log_scale_sufficient, 
      mult = -.25
    )
    
    if (
      (all(q < hierarchical_model$range$log_scale_bounds[2] )) &
      (all(q > hierarchical_model$range$log_scale_bounds[1] ))
    ){
      if(!is.nan(current_U-proposed_U+current_K-proposed_K)){
        if(log(runif(1)) < (current_U-proposed_U+current_K-proposed_K))
        {
          print("turlututu sufficient !")
          
          ker_var$range_log_scale_sufficient = update_kernel(
            iter = iter, iter_start = iter_start, 
            update_kernel_groupsize = 1, kernel_value = ker_var$range_log_scale_sufficient, 
            mult = 1
          )
          momenta$range_beta_sufficient = p
          
          params$range_beta = new_range_beta
          params$range_log_scale = q
          
          stuff$sparse_chol= new_sparse_chol
          stuff$compressed_chol = new_compressed_sparse_chol
        }}}
    # sufficient - sufficient ####
    for(i in seq(10))
    {
      q = params$range_log_scale + rnorm(length(params$range_log_scale), 0, .1)
      if(
        (all(q < hierarchical_model$range$log_scale_bounds[2] )) &
        (all(q > hierarchical_model$range$log_scale_bounds[1] )) &
        (
          - 3 * sum(sqrt(max(0, q - hierarchical_model$range$log_scale_bounds[1])))
          + 3 * sum(sqrt(max(0, params$range_log_scale - hierarchical_model$range$log_scale_bounds[1])))
          
          + beta_prior_log_dens(
            beta = params$range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
            beta0_mean = hierarchical_model$range$beta0_mean,
            beta0_var =  hierarchical_model$range$beta0_sd^2, 
            log_scale = q)
          - beta_prior_log_dens(
            beta = params$range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
            beta0_mean = hierarchical_model$range$beta0_mean,
            beta0_var =  hierarchical_model$range$beta0_sd^2, 
            log_scale = params$range_log_scale) 
          > log(runif(1))
        )
      )
      {
        params$range_log_scale = q
      }
    }
    # ancillary-ancillary ####
    
    q = params$range_log_scale 
    q = q + exp(ker_var$range_log_scale_ancillary) * rnorm(length(q))
    new_range_beta = params$range_beta
    new_range_beta[-seq(covariates$range_X$n_regressors),] = new_range_beta[-seq(covariates$range_X$n_regressors),] %*% 
      diag(exp(-.5 * params$range_log_scale[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic) %*%
      diag(exp(.5 * q[c(1, rep(2, 2*hierarchical_model$anisotropic))]), 1 + 2*hierarchical_model$anisotropic)
    new_compressed_sparse_chol = 
      compute_sparse_chol(
        hierarchical_model$range$PP,
        range_beta = new_range_beta, 
        vecchia_approx = vecchia_approx, 
        range_X = covariates$range_X, 
        matern_smoothness = hierarchical_model$matern_smoothness, 
        compute_derivative = T, num_threads = num_threads
      )
    new_sparse_chol = decompress_chol(vecchia_approx = vecchia_approx, new_compressed_sparse_chol)
    new_field = as.vector(Matrix::solve(new_sparse_chol, stuff$sparse_chol %*% (params$field)))
    
    current_U =
      (
        3 * sum(sqrt(max(0, params$range_log_scale - hierarchical_model$range$log_scale_bounds[1])))
        + .5 * sum((stuff$lm_residuals -  params$field[vecchia_approx$locs_match])^2/stuff$noise) # observation ll
      )
    proposed_U =
      (
        3 * sum(sqrt(max(0, q - hierarchical_model$range$log_scale_bounds[1])))
        + .5 * sum((stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/stuff$noise) # observation ll
      )
    ker_var$range_log_scale_ancillary = update_kernel(
      iter = iter, iter_start = iter_start, 
      update_kernel_groupsize = 1, kernel_value = ker_var$range_log_scale_ancillary, 
      mult = -.25
    )
    
    if (
      (all(q < hierarchical_model$range$log_scale_bounds[2] )) &
      (all(q > hierarchical_model$range$log_scale_bounds[1] ))
    )
    {
      if(!is.nan(current_U-proposed_U)){
        if(log(runif(1)) < (current_U-proposed_U)){
          ker_var$range_log_scale_ancillary = update_kernel(
            iter = iter, iter_start = iter_start, 
            update_kernel_groupsize = 1, kernel_value = ker_var$range_log_scale_ancillary, 
            mult = 1
          )
          print("turlututu ancillary !")
          
          params$range_beta = new_range_beta
          params$range_log_scale = q
          params$field = new_field
          
          stuff$sparse_chol= new_sparse_chol
          stuff$compressed_chol = new_compressed_sparse_chol
        }}}
    
    # sufficient - sufficient 
    for(i in seq(10))
    {
      q = params$range_log_scale + rnorm(length(params$range_log_scale), 0, .1)
      if(
        (all(q < hierarchical_model$range$log_scale_bounds[2] )) &
        (all(q > hierarchical_model$range$log_scale_bounds[1] )) &
        (
          - 3 * sum(sqrt(max(0, q - hierarchical_model$range$log_scale_bounds[1])))
          + 3 * sum(sqrt(max(0, params$range_log_scale - hierarchical_model$range$log_scale_bounds[1])))
          
          + beta_prior_log_dens(
            beta = params$range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
            beta0_mean = hierarchical_model$range$beta0_mean,
            beta0_var =  hierarchical_model$range$beta0_sd^2, 
            log_scale = q)
          - beta_prior_log_dens(
            beta = params$range_beta, n_PP = hierarchical_model$range$PP$n_knots, 
            beta0_mean = hierarchical_model$range$beta0_mean,
            beta0_var =  hierarchical_model$range$beta0_sd^2, 
            log_scale = params$range_log_scale) 
          > log(runif(1))
        )
      )
      {
        params$range_log_scale = q
      }
    }
  }
  
  
  
  
  #########
  # Noise #
  #########
  ##############
  # Noise beta #
  ##############
  #    # VEWY IMPOWTANT don't remove or comment
  #    squared_residuals = as.vector(stuff$lm_residuals - params$field[vecchia_approx$locs_match])^2
  #    # HMC update
  #    q = noise_X$chol_crossprod_X %*% params$noise_beta
  #    current_U =
  #      (
  #        - beta_prior_log_dens(
  #          beta = params$noise_beta, n_PP = hierarchical_model$PP$n_knots*hierarchical_model$noise_PP, 
  #          chol_crossprod_X = noise_X$chol_crossprod_X,
  #          beta0_mean = hierarchical_model$noise_beta0_mean,
  #          beta0_var =  hierarchical_model$noise_beta0_var, 
  #          log_scale = params$noise_log_scale) # normal prior 
  #        +.5* sum(log(stuff$noise)) # det
  #        +.5*sum(squared_residuals/stuff$noise) # observations
  #      )
  #    # HMC whitened
  #    momenta$noise_beta = sqrt(.9) * momenta$noise_beta + sqrt(.1)*rnorm(length(momenta$noise_beta))
  #    p = momenta$noise_beta
  #    
  #    # Make a half step for momentum at the beginning
  #    p = p - exp(ker_var$noise_beta_mala) *
  #      (
  #        + solve(t(noise_X$chol_crossprod_X), # solving by prior chol because of whitening
  #                - beta_prior_log_dens_derivative(
  #                  beta = params$noise_beta, n_PP = hierarchical_model$PP$n_knots*hierarchical_model$noise_PP, 
  #                  chol_crossprod_X = noise_X$chol_crossprod_X,
  #                  beta0_mean = hierarchical_model$noise_beta0_mean,
  #                  beta0_var =  hierarchical_model$noise_beta0_var, 
  #                  log_scale = params$noise_log_scale) # normal prior
  #                + X_PP_crossprod(X = noise_X$X, PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP, 
  #                                 Y = 
  #                                   (
  #                                     + .5 # determinant part of normal likelihood
  #                                     - (squared_residuals/stuff$noise)/2 # exponential part of normal likelihood
  #                                   ))
  #        ))/ 2
  #    # checking gradient with finite differences, to update
  #    
  #    #### idx = 200
  #    #### noise_beta_ = params$noise_beta 
  #    #### noise_beta_[idx] = noise_beta_[idx] + 0.0000001
  #    #### noise_ = variance_field(
  #    ####   beta = noise_beta_, X = noise_X$X, 
  #    ####   PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP
  #    #### )
  #    #### U_ =
  #    ####   (
  #    ####     - beta_prior_log_dens(beta = noise_beta_, n_PP = hierarchical_model$PP$n_knots*hierarchical_model$noise_PP, 
  #    ####                                   beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
  #    ####                                   beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
  #    ####                                   log_scale = params$noise_log_scale) # normal prior 
  #    ####     +.5* sum(log(noise_)) # det
  #    ####     +.5*sum(squared_residuals/noise_) # observations
  #    ####   )
  #    #### print(
  #    ####   (
  #    ####             - beta_prior_log_dens_derivative(beta = params$noise_beta, n_PP = hierarchical_model$PP$n_knots*hierarchical_model$noise_PP, 
  #    ####                                                      beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
  #    ####                                                      beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
  #    ####                                                      log_scale = params$noise_log_scale) # normal prior
  #    ####             + X_PP_crossprod(X = noise_X$X, PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP, 
  #    ####                                      Y = 
  #    ####                                        (
  #    ####                                          + .5 # determinant part of normal likelihood
  #    ####                                          - (squared_residuals/stuff$noise)/2 # exponential part of normal likelihood
  #    ####                                        ))
  #    ####     )[idx]
  #    #### )
  #    #### print(10000000*(U_- current_U))
  #    
  #    
  #    # Make a full step for the position
  #    q = q + exp(ker_var$noise_beta_mala) * p
  #    new_noise_beta = solve(noise_X$chol_crossprod_X, q)
  #    new_noise = variance_field(beta = new_noise_beta, PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP, 
  #                               X = noise_X$X, locs_idx = NULL)
  #    # Make a half step for momentum at the end
  #    p = p - exp(ker_var$noise_beta_mala) *
  #      (
  #        + solve(t(noise_X$chol_crossprod_X), # solving by prior sparse chol because of whitening
  #                - beta_prior_log_dens_derivative(
  #                  beta = new_noise_beta, n_PP = hierarchical_model$PP$n_knots*hierarchical_model$noise_PP, 
  #                  chol_crossprod_X = noise_X$chol_crossprod_X,
  #                  beta0_mean = hierarchical_model$noise_beta0_mean,
  #                  beta0_var =  hierarchical_model$noise_beta0_var, 
  #                  log_scale = params$noise_log_scale) # normal prior  
  #                + X_PP_crossprod(X = noise_X$X, PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP, 
  #                                 (
  #                                   + .5 # determinant part of normal likelihood
  #                                   - (squared_residuals/new_noise)/2 # exponential part of normal likelihood
  #                                 ))
  #        ))/ 2
  #    
  #    # Evaluate potential and kinetic energies at start and end of trajectory
  #    current_K = sum (momenta$noise_beta ^2) / 2
  #    proposed_U = 
  #      (
  #        - beta_prior_log_dens(beta = new_noise_beta, n_PP = hierarchical_model$PP$n_knots*hierarchical_model$noise_PP, 
  #                              chol_crossprod_X = noise_X$chol_crossprod_X,
  #                              beta0_mean = hierarchical_model$noise_beta0_mean,
  #                              beta0_var =  hierarchical_model$noise_beta0_var, 
  #                              log_scale = params$noise_log_scale) # normal prior        
  #        +.5* sum(log(new_noise)) # det
  #        +.5*sum(squared_residuals/new_noise) # observations
  #      )
  #    proposed_K = sum(p^2) / 2
  #    
  #    
  #    ker_var$noise_beta_mala = ker_var$noise_beta_mala- 1/sqrt(iter_start + iter +100)
  #    if(!is.nan(current_U-proposed_U+current_K- proposed_K))
  #    {
  #      if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
  #      {
  #        ker_var$noise_beta_mala = ker_var$noise_beta_mala + 2/sqrt(iter_start + iter +100)
  #        momenta$noise_beta = p
  #        params$noise_beta[] = new_noise_beta
  #        stuff$noise = new_noise
  #      }
  #    }
  #    
  #    ###################
  #    # Noise log scale # 
  #    ###################
  #    if(hierarchical_model$noise_PP)
  #    {
  #      # ancillary -- sufficient ####
  #      new_noise_log_scale = params$noise_log_scale + rnorm(1, 0, exp(ker_var$noise_log_scale))
  #      new_noise_beta = params$noise_beta
  #      new_noise_beta[-seq(noise_X$n_regressors)] = new_noise_beta[-seq(noise_X$n_regressors)] *
  #        exp((new_noise_log_scale - params$noise_log_scale)/2)
  #      new_noise =variance_field(beta = new_noise_beta, PP = hierarchical_model$PP, 
  #                                use_PP = hierarchical_model$noise_PP, X = noise_X$X, 
  #                                locs_idx = NULL)
  #      ll_ratio = (
  #        -.5* sum(log(new_noise)) 
  #        -.5*sum(squared_residuals/new_noise)
  #        +.5* sum(log(stuff$noise)) 
  #        +.5*sum(squared_residuals/stuff$noise)
  #      )
  #      if(!is.nan(ll_ratio))
  #      {
  #        if(ll_ratio > log(runif(1)))
  #        {
  #          if(
  #            (new_noise_log_scale > hierarchical_model$noise_log_scale_prior[1])&
  #            (new_noise_log_scale < hierarchical_model$noise_log_scale_prior[2])
  #          )
  #          {
  #            params$noise_log_scale = new_noise_log_scale
  #            params$noise_beta = new_noise_beta 
  #            stuff$noise = new_noise
  #            ker_var$noise_log_scale = ker_var$noise_log_scale + 4/sqrt(iter_start + iter +100)
  #          }
  #        }
  #      }
  #      ker_var$noise_log_scale = ker_var$noise_log_scale - 1/sqrt(iter_start + iter +100)
  #      
  #      # sufficient -- sufficient ####
  #      for(i in seq(10))
  #      {
  #        new_noise_log_scale = params$noise_log_scale + rnorm(1, 0, .1)
  #        if(
  #          (
  #            + beta_prior_log_dens(beta = params$noise_beta, n_PP = hierarchical_model$PP$n_knots*hierarchical_model$noise_PP, 
  #                                  chol_crossprod_X = noise_X$chol_crossprod_X,
  #                                  beta0_mean = hierarchical_model$noise_beta0_mean,
  #                                  beta0_var =  hierarchical_model$noise_beta0_var, 
  #                                  log_scale = new_noise_log_scale) - 
  #            beta_prior_log_dens(beta = params$noise_beta, n_PP = hierarchical_model$PP$n_knots*hierarchical_model$noise_PP, 
  #                                chol_crossprod_X = noise_X$chol_crossprod_X,
  #                                beta0_mean = hierarchical_model$noise_beta0_mean,
  #                                beta0_var =  hierarchical_model$noise_beta0_var, 
  #                                log_scale = params$noise_log_scale)
  #          ) > log(runif(1))
  #          &(new_noise_log_scale > hierarchical_model$noise_log_scale_prior[1])
  #          &(new_noise_log_scale < hierarchical_model$noise_log_scale_prior[2])
  #        )
  #        {
  #          params$noise_log_scale = new_noise_log_scale
  #        }
  #      }
  #    }
  #    
  #    
  #    
  #######################
  # Storing the samples #
  #######################
  params_records[[iter]] = params
  
  #     if(iter>10){
  #     HMC_range_cond = solve(chol(var(t(sapply(params_records[seq(iter/5, iter)], function(x)c(x$range_beta)))) + diag(.05, length(params$range_beta))))}
}
print(Sys.time()-t1)
#return(list("state" = state, "params_records" = params_records))
#}



