# With only one HMC step in range 

# numerically stable log(exp(x) + exp(y)) with 0 << x << y or 0 << y << x
# let's do: log(exp(x) + exp(y)) with 0 << x << y 
# log(exp(x) + exp(y))  = log(exp(y)(1 + exp(x) / exp(y)))
# = log(exp(y)) + log((1 + exp(x) / exp(y)))
# = y + log((1 + exp(x - y)))
log_of_exp_sum = function(x,y)return(max(x,y) + log((1 + exp(-abs(x - y)))))


#list2env(mcmc_nngp_list, envir = environment())
#MCMC_NNGP$state = mcmc_nngp_list$states$chain_1; n_iterations_update  =100; num_threads = 10; iter_start = 0; seed = 1; iter=1

#' @export
mcmc_nngp_update_Gaussian = function(MCMC_NNGP, # model MCMC_NNGP$state
                                     n_iterations_update = 100, 
                                     num_threads = 1 # practical settings
)
{
  iter_start = MCMC_NNGP$iterations$checkpoints[nrow(MCMC_NNGP$iterations$checkpoints), 1]
  message(paste("Starting MCMC chain for Nonstat NNGP model at iteration", iter_start, "for", n_iterations_update, "iterations"))
  #################
  # Sanity checks #
  #################
  # iterations
  if((floor(n_iterations_update)!=n_iterations_update) |  n_iterations_update<1)stop("n_iterations_update must be a positive round number")
  if((floor(iter_start)!=iter_start)|  iter_start<0)stop("iter_start must be a positive round number")
  # set seed 
  set.seed(seed)
  #########################################
  # Initializing chain storage structures #
  #########################################
  # this part re-creates a small portion of the $records objects of each chain. It fills it with chain MCMC_NNGP$state during the run, and then updates each chain with the new values
  params_records = lapply(MCMC_NNGP$state$params, function(x)array(0, c(dim(as.matrix(x)), n_iterations_update)))
  #par(mfrow = c(2, 1))
  
  #################
  # Gibbs sampler #
  #################
  for(iter in seq(1, n_iterations_update))
  {
    gc()
    if((iter + iter_start)/10 ==(iter + iter_start) %/% 10)print(paste("iteration", (iter + iter_start)))
    
    ###########################
    # Regression coefficients #
    ###########################
    t1 = Sys.time()
    # sufficient parametrization of latent field
    beta_covmat = solve(as.matrix(t(MCMC_NNGP$data$covariates$X$X) %*% Matrix::Diagonal(x = 1/MCMC_NNGP$state$sparse_chol_and_stuff$noise) %*% MCMC_NNGP$data$covariates$X$X))
    if(all(!is.infinite(beta_covmat) & !is.nan(beta_covmat)))
    {
      if(all(eigen(beta_covmat)$d >0))
      {
        beta_mean = c((((MCMC_NNGP$data$observed_field-MCMC_NNGP$state$params$field[MCMC_NNGP$vecchia_approx$locs_match]) / MCMC_NNGP$state$sparse_chol_and_stuff$noise) %*% MCMC_NNGP$data$covariates$X$X) %*% beta_covmat)
        MCMC_NNGP$state$params$beta[]   = c(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
      }}
    # interweaving centered sampling in case of location-wise MCMC_NNGP$data to improve beta sampling
    centered_field = as.vector(MCMC_NNGP$state$params$field + MCMC_NNGP$data$covariates$X$X_locs%*%matrix(MCMC_NNGP$state$params$beta[MCMC_NNGP$data$covariates$X$which_locs], ncol = 1))
    sparse_chol_X = as.matrix(MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% Matrix::Diagonal(x = 1/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale)) %*% MCMC_NNGP$data$covariates$X$X_locs)
    beta_precision = crossprod(x = sparse_chol_X, y = sparse_chol_X)
    beta_covmat = solve(beta_precision, tol = min(rcond(beta_precision),.Machine$double.eps))
    if(all(!is.infinite(beta_covmat) & !is.nan(beta_covmat)))
    {
      if(all(eigen(beta_covmat)$d >0))
      {
        beta_mean =  c(as.vector(MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (centered_field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale)))  %*% sparse_chol_X %*% beta_covmat)
        MCMC_NNGP$state$params$beta[MCMC_NNGP$data$covariates$X$which_locs]   = as.vector(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
        MCMC_NNGP$state$params$field = centered_field - as.vector(MCMC_NNGP$data$covariates$X$X_locs %*% matrix(MCMC_NNGP$state$params$beta[MCMC_NNGP$data$covariates$X$which_locs], ncol = 1))
      }}
    # updating stuff 
    MCMC_NNGP$state$sparse_chol_and_stuff$lm_fit       = as.vector(MCMC_NNGP$data$covariates$X$X%*%     matrix(MCMC_NNGP$state$params$beta, ncol = 1))
    MCMC_NNGP$state$sparse_chol_and_stuff$lm_fit_locs  = as.vector(MCMC_NNGP$data$covariates$X$X_locs%*%matrix(MCMC_NNGP$state$params$beta[MCMC_NNGP$data$covariates$X$which_locs], ncol = 1))
    MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals = as.vector(MCMC_NNGP$data$observed_field-              MCMC_NNGP$state$sparse_chol_and_stuff$lm_fit)
    #print(Sys.time()-t1)
    ################
    # Latent field #
    ################
    locs_partition = MCMC_NNGP$vecchia_approx$locs_partition[,runif(1, 1, ncol(MCMC_NNGP$vecchia_approx$locs_partition))]
    additional_precision = MCMC_NNGP$vecchia_approx$locs_match_matrix %*% (1/MCMC_NNGP$state$sparse_chol_and_stuff$noise)
    additional_mean = MCMC_NNGP$vecchia_approx$locs_match_matrix %*% ((MCMC_NNGP$data$observed_field - MCMC_NNGP$state$sparse_chol_and_stuff$lm_fit)/MCMC_NNGP$state$sparse_chol_and_stuff$noise)
    chol_sub_list = lapply(
      seq(length(unique(locs_partition))), function(idx)
      {
        selected_locs = which(locs_partition==idx)
        MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol[,selected_locs] %*% Matrix::Diagonal(x = 1/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale[selected_locs]))
      }
    )
    chol_Q_AA_sub_list = 
      lapply(
        seq(length(unique(locs_partition))), function(idx)
        {
          selected_locs = which(locs_partition==idx)
          Q_AA_sub =  Matrix::crossprod(chol_sub_list[[idx]])
          Matrix::diag(Q_AA_sub) = Matrix::diag(Q_AA_sub) + additional_precision[selected_locs]
          expanded = Matrix::expand(Matrix::Cholesky(Q_AA_sub, perm = T))
        }
      )
    
    for(idx in seq(length(unique(locs_partition))))
    {
      selected_locs = which(locs_partition==idx)
      cond_mean = - 
        Matrix::t(chol_Q_AA_sub_list[[idx]]$P) %*%
        Matrix::solve(Matrix::t(chol_Q_AA_sub_list[[idx]]$L),
                      Matrix::solve(chol_Q_AA_sub_list[[idx]]$L, # inverse of precision matrix...
                                    chol_Q_AA_sub_list[[idx]]$P %*%
                                      (as.vector(Matrix::crossprod(chol_sub_list[[idx]], MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol%*%((MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale))*(locs_partition!=idx)))) # rest of the latent field
                                       - additional_mean[selected_locs]))) # Gaussian observations 
      MCMC_NNGP$state$params$field [selected_locs] = as.vector(as.vector(cond_mean) +  Matrix::t(chol_Q_AA_sub_list[[idx]]$P) %*% Matrix::solve(Matrix::t(chol_Q_AA_sub_list[[idx]]$L), rnorm(length(selected_locs))) )
    }
    
    
    
    
    
    ###############
    # Range beta  #
    ###############
    ##########################
    # Range beta (ancillary) #
    ##########################
    if(!MCMC_NNGP$hierarchical_model$anisotropic) eps_mat = matrix(exp(MCMC_NNGP$state$transition_kernels$range_beta_ancillary[1]))
    if(MCMC_NNGP$hierarchical_model$anisotropic) eps_mat = diag(exp(MCMC_NNGP$state$transition_kernels$range_beta_ancillary[c(1,2,2)]))
    q = MCMC_NNGP$data$covariates$range$chol_crossprod_X_locs %*% MCMC_NNGP$state$params$range_beta # whitening wrt covariates of the range
    current_U =
      (
        - beta_prior_log_dens(beta = MCMC_NNGP$state$params$range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
                                      chol_crossprod_X = MCMC_NNGP$data$covariates$range_X$chol_crossprod_X,
                                      beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean, 
                                      beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
                                      log_scale = MCMC_NNGP$state$params$range_log_scale) # normal prior
        + .5 * sum((MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals -  MCMC_NNGP$state$params$field[MCMC_NNGP$vecchia_approx$locs_match])^2/MCMC_NNGP$state$sparse_chol_and_stuff$noise) # observation ll
      )
    # MALA whitened
    MCMC_NNGP$state$momenta$range_beta_ancillary = sqrt(.9) * MCMC_NNGP$state$momenta$range_beta_ancillary + sqrt(.1)*matrix(rnorm(length(q)), nrow = nrow(q))
    p = MCMC_NNGP$state$momenta$range_beta_ancillary
    # Make a half step for momentum at the beginning
    p = p - 
      as.matrix(
        solve(t(MCMC_NNGP$data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
              - beta_prior_log_dens_derivative(
                beta = MCMC_NNGP$state$params$range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
                chol_crossprod_X = MCMC_NNGP$data$covariates$range_X$chol_crossprod_X,
                beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean,
                beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
                log_scale = MCMC_NNGP$state$params$range_log_scale) # normal prior
              # normal prior derivative                
              + X_PP_crossprod(
                X = MCMC_NNGP$data$covariates$range_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
                Y = # Jacobian of range field wrt range_beta
                  (
                    # natural gradient of obs likelihood wrt range field
                    derivative_sandwiches(derivatives = MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                  left_vector = as.vector(
                                                    Matrix::solve(
                                                      Matrix::t(MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol), 
                                                      - as.vector(MCMC_NNGP$vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                                    ((MCMC_NNGP$state$params$field[MCMC_NNGP$vecchia_approx$locs_match] - MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals) / MCMC_NNGP$state$sparse_chol_and_stuff$noise))
                                                      * sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale) # part of sparse chol
                                                    )), 
                                                  right_vector = MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                  NNarray = MCMC_NNGP$vecchia_approx$NNarray  
                    )
                  ))
        ))%*%eps_mat/ 2
    
    
    #######testing the gradient
    ##source("Bidart/R/Useful_stuff.R")
    ##d1 =         - beta_prior_log_dens_derivative(
    ##  beta = MCMC_NNGP$state$params$range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
    ##  beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$range_beta_mean, 
    ##  beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$range_beta_precision, 
    ##  log_scale = MCMC_NNGP$state$params$range_log_scale) # normal prior
    ##beta_ = MCMC_NNGP$state$params$range_beta 
    ##beta_[2, 2] = beta_[2, 2] + 1/10000
    ##d2 = 10000*(
    ##  -
    ##  beta_prior_log_dens(beta = beta_, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
    ##                        beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$range_beta_mean, 
    ##                        beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$range_beta_precision, 
    ##                        log_scale = MCMC_NNGP$state$params$range_log_scale) # normal prior
    ##  + 
    ##  beta_prior_log_dens(beta = MCMC_NNGP$state$params$range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
    ##                        beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$range_beta_mean, 
    ##                        beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$range_beta_precision, 
    ##                        log_scale = MCMC_NNGP$state$params$range_log_scale) # normal prior
    ##)
    ##d1/d2
    #### i = 4
    #### j = 2
    #### 
    #### range_beta_ = MCMC_NNGP$state$params$range_beta
    #### range_beta_[i,j] = MCMC_NNGP$state$params$range_beta[i,j] + .0001
    #### 
    #### compressed_sparse_chol_and_grad_ = 
    ####   compute_sparse_chol(
    ####     num_threads = num_threads, 
    ####     anisotropic = MCMC_NNGP$hierarchical_model$anisotropic, 
    ####     sphere = MCMC_NNGP$hierarchical_model$sphere, 
    ####     range_beta = range_beta_, NNarray = MCMC_NNGP$vecchia_approx$NNarray, 
    ####     locs = MCMC_NNGP$data$locs, 
    ####     range_X = MCMC_NNGP$data$covariates$range_X$X_locs, 
    ####     nu = MCMC_NNGP$hierarchical_model$nu, 
    ####     PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP,
    ####     compute_derivative = T, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1
    ####   )
    #### sparse_chol_ = Matrix::sparseMatrix(i = MCMC_NNGP$vecchia_approx$sparse_chol_row_idx, j = MCMC_NNGP$vecchia_approx$sparse_chol_column_idx, x = compressed_sparse_chol_and_grad_[[1]][MCMC_NNGP$vecchia_approx$NNarray_non_NA], triangular = T)
    #### field_ = sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(sparse_chol_, MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale))))
    #### 
    #### 10000*(
    ####   + .5 * sum((MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals -  field_[MCMC_NNGP$vecchia_approx$locs_match])^2/MCMC_NNGP$state$sparse_chol_and_stuff$noise) - # observation ll
    ####     + .5 * sum((MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals -  MCMC_NNGP$state$params$field[MCMC_NNGP$vecchia_approx$locs_match])^2/MCMC_NNGP$state$sparse_chol_and_stuff$noise) # observation ll
    #### )
    #### 
    #### 
    #### + X_PP_crossprod(
    ####   X = MCMC_NNGP$data$covariates$range_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
    ####   Y = # Jacobian of range field wrt range_beta
    ####     (
    ####       # natural gradient of obs likelihood wrt range field
    ####       derivative_sandwiches(derivatives = MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
    ####                                     left_vector = as.vector(
    ####                                       Matrix::solve(
    ####                                         Matrix::t(MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol), 
    ####                                         - as.vector(MCMC_NNGP$vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
    ####                                                       ((MCMC_NNGP$state$params$field[MCMC_NNGP$vecchia_approx$locs_match] - MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals) / MCMC_NNGP$state$sparse_chol_and_stuff$noise))
    ####                                         * sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale) # part of sparse chol
    ####                                       )), 
    ####                                     right_vector = MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
    ####                                     NNarray = MCMC_NNGP$vecchia_approx$NNarray  
    ####       )
    ####     ))[i,j]
    
    
    # Make a full step for the position
    q = q + p %*%eps_mat
    new_range_beta = MCMC_NNGP$state$params$range_beta
    new_range_beta = solve(MCMC_NNGP$data$covariates$range_X$chol_crossprod_X_locs, q )
    new_compressed_sparse_chol_and_grad = 
      compute_sparse_chol(
        num_threads = num_threads, 
        anisotropic = MCMC_NNGP$hierarchical_model$anisotropic, 
        sphere = MCMC_NNGP$hierarchical_model$sphere, 
        range_beta = new_range_beta, NNarray = MCMC_NNGP$vecchia_approx$NNarray, 
        locs = MCMC_NNGP$data$locs, 
        range_X = MCMC_NNGP$data$covariates$range_X$X_locs, 
        nu = MCMC_NNGP$hierarchical_model$nu, 
        PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP,
        compute_derivative = T, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1
      )
    new_sparse_chol = Matrix::sparseMatrix(i = MCMC_NNGP$vecchia_approx$sparse_chol_row_idx, j = MCMC_NNGP$vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][MCMC_NNGP$vecchia_approx$NNarray_non_NA], triangular = T)
    new_field = sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale))))
    # Make a half step for momentum at the end.
    p = p - 
      as.matrix(
        solve(t(MCMC_NNGP$data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
              - beta_prior_log_dens_derivative
              (beta = new_range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
                chol_crossprod_X = MCMC_NNGP$data$covariates$range_X$chol_crossprod_X,
                beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean,
                beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
                log_scale = MCMC_NNGP$state$params$range_log_scale) # normal prior
              #normal prior derivative                
              + X_PP_crossprod(
                X = MCMC_NNGP$data$covariates$range_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
                Y = # Jacobian of range field wrt range_beta
                  (
                    # natural gradient of obs likelihood wrt range field
                    derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                  left_vector = as.vector(
                                                    Matrix::solve(
                                                      Matrix::t(new_sparse_chol), 
                                                      - as.vector(MCMC_NNGP$vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                                    ((new_field[MCMC_NNGP$vecchia_approx$locs_match] - MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals) / MCMC_NNGP$state$sparse_chol_and_stuff$noise))
                                                      * sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale) # part of sparse chol
                                                    )), 
                                                  right_vector = new_field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                  NNarray = MCMC_NNGP$vecchia_approx$NNarray  
                    )
                  )
              )
        ))%*%eps_mat/ 2
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (MCMC_NNGP$state$momenta$range_beta_ancillary ^2) / 2
    proposed_U =
      (
        - beta_prior_log_dens(
          beta = new_range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
          chol_crossprod_X = MCMC_NNGP$data$covariates$range_X$chol_crossprod_X,
          beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean,
          beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
          log_scale = MCMC_NNGP$state$params$range_log_scale) # normal prior
        # normal prior 
        + .5 * sum((MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals -  new_field[MCMC_NNGP$vecchia_approx$locs_match])^2/MCMC_NNGP$state$sparse_chol_and_stuff$noise) # observation ll
      )
    proposed_K = sum(p^2) / 2
    
    #print(MCMC_NNGP$state$transition_kernels$range_beta_ancillary)
    
    MCMC_NNGP$state$transition_kernels$range_beta_ancillary[iter%%2+1] = MCMC_NNGP$state$transition_kernels$range_beta_ancillary[iter%%2+1]- 15/(iter_start + iter + 100)
    if(!is.nan(current_U-proposed_U+current_K- proposed_K))
    {
      if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
      {
        MCMC_NNGP$state$transition_kernels$range_beta_ancillary[iter%%2+1] = MCMC_NNGP$state$transition_kernels$range_beta_ancillary[iter%%2+1] +20/(iter_start + iter + 100)
        #print("tatato!")
        MCMC_NNGP$state$momenta$range_beta_ancillary = p
        MCMC_NNGP$state$params$field = new_field
        MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
        MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
        MCMC_NNGP$state$params$range_beta[] = new_range_beta
        MCMC_NNGP$state$sparse_chol_and_stuff$precision_diag = as.vector((MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][MCMC_NNGP$vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(MCMC_NNGP$vecchia_approx$sparse_chol_column_idx)), j = MCMC_NNGP$vecchia_approx$sparse_chol_column_idx, x = rep(1, length(MCMC_NNGP$vecchia_approx$sparse_chol_row_idx))))
      }
    }
    
    ###########################
    # Range beta (sufficient) #
    ###########################
    if(!MCMC_NNGP$hierarchical_model$anisotropic) eps_mat = matrix(exp(MCMC_NNGP$state$transition_kernels$range_beta_sufficient[1]))
    if(MCMC_NNGP$hierarchical_model$anisotropic) eps_mat = diag(exp(MCMC_NNGP$state$transition_kernels$range_beta_sufficient[c(1,2,2)]))
    q = MCMC_NNGP$data$covariates$range$chol_crossprod_X_locs %*% MCMC_NNGP$state$params$range_beta # whitening wrt covariates of the range
    current_U =
      (
        - beta_prior_log_dens(
          beta = MCMC_NNGP$state$params$range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
          chol_crossprod_X = MCMC_NNGP$data$covariates$range_X$chol_crossprod_X,
          beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean,
          beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
          log_scale = MCMC_NNGP$state$params$range_log_scale) # normal prior
        # normal prior 
        + .5* sum((MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale)))^2)
        - sum(log(MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
      )
    
    # MALA whitened
    MCMC_NNGP$state$momenta$range_beta_sufficient = sqrt(.9) * MCMC_NNGP$state$momenta$range_beta_sufficient + sqrt(.1)*matrix(rnorm(length(q)), nrow = nrow(q))
    p = MCMC_NNGP$state$momenta$range_beta_sufficient
    # Make a half step for momentum at the beginning
    p = p - 
      as.matrix(
        solve(t(MCMC_NNGP$data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
              - beta_prior_log_dens_derivative
              (beta = MCMC_NNGP$state$params$range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
                chol_crossprod_X = MCMC_NNGP$data$covariates$range_X$chol_crossprod_X,
                beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean,
                beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
                log_scale = MCMC_NNGP$state$params$range_log_scale) # normal prior
              # normal prior derivative                
              + X_PP_crossprod(
                X = MCMC_NNGP$data$covariates$range_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
                Y = # Jacobian of range field wrt range_beta
                  (# natural gradient of obs likelihood wrt range field
                    derivative_sandwiches(derivatives = MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                  left_vector = as.vector(MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                                  right_vector = MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                  NNarray = MCMC_NNGP$vecchia_approx$NNarray  
                    )
                    - log_determinant_derivatives(sparse_chol_and_grad = MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad, NNarray = MCMC_NNGP$vecchia_approx$NNarray)# derivative of determinant
                  ))))%*% eps_mat/ 2
    #### Checking the gradient
    ##      source("Bidart/R/Useful_stuff.R")
    ##       # recomputing current sparse chol
    ##       MCMC_NNGP$state$params$range_beta[] = rnorm(length(MCMC_NNGP$state$params$range_beta[]))
    ##       MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = 
    ##         compute_sparse_chol(covfun_name = MCMC_NNGP$hierarchical_model$covfun, 
    ##                             range_beta = MCMC_NNGP$state$params$range_beta, NNarray = MCMC_NNGP$vecchia_approx$NNarray, 
    ##                             locs = MCMC_NNGP$data$locs, 
    ##                             range_X = MCMC_NNGP$data$covariates$range_X$X_locs, 
    ##                             nu = MCMC_NNGP$hierarchical_model$nu, 
    ##                             PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP,
    ##                             compute_derivative = T, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1
    ##         )
    ##       MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol = Matrix::sparseMatrix(i = MCMC_NNGP$vecchia_approx$sparse_chol_row_idx, j = MCMC_NNGP$vecchia_approx$sparse_chol_column_idx, x = MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][MCMC_NNGP$vecchia_approx$NNarray_non_NA], triangular = T)
    ##       
    ##       # compute gradient using derivative of sparse chol
    ##       d1 = X_PP_crossprod(
    ##         X = MCMC_NNGP$data$covariates$range_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
    ##         Y = # Jacobian of range field wrt range_beta
    ##           (
    ##             # natural gradient of obs likelihood wrt range field
    ##             derivative_sandwiches(derivatives = MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
    ##                                           left_vector = as.vector(MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
    ##                                           right_vector = MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
    ##                                           NNarray = MCMC_NNGP$vecchia_approx$NNarray  
    ##             )
    ##             - log_determinant_derivatives(sparse_chol_and_grad = MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad, NNarray = MCMC_NNGP$vecchia_approx$NNarray)# derivative of determinant
    ##             ))
    ##       # compute gradient using finite diff
    ##       new_range_beta = MCMC_NNGP$state$params$range_beta
    ##       new_range_beta[10,1] = new_range_beta[10,1] + .0001
    ##       new_compressed_sparse_chol_and_grad = 
    ##         compute_sparse_chol(covfun_name = MCMC_NNGP$hierarchical_model$covfun, 
    ##                                     range_beta = new_range_beta, NNarray = MCMC_NNGP$vecchia_approx$NNarray, 
    ##                                     locs = MCMC_NNGP$data$locs, 
    ##                                     range_X = MCMC_NNGP$data$covariates$range_X$X_locs, 
    ##                                     nu = MCMC_NNGP$hierarchical_model$nu, 
    ##                                     PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP,
    ##                                     compute_derivative = T, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1
    ##         )
    ##       new_sparse_chol = Matrix::sparseMatrix(i = MCMC_NNGP$vecchia_approx$sparse_chol_row_idx, j = MCMC_NNGP$vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][MCMC_NNGP$vecchia_approx$NNarray_non_NA], triangular = T)
    ##       new_field = sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale))))
    ##       d2 = 10000 * (
    ##         
    ##         (
    ##           + .5* sum((new_sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale)))^2)
    ##           - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
    ##         )
    ##         -
    ##           (
    ##             + .5* sum((MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale)))^2)
    ##             - sum(log(MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
    ##           )
    ##       )
    ##       d1/d2
    # Make a full step for the position
    q = q + p %*% eps_mat
    new_range_beta = MCMC_NNGP$state$params$range_beta
    new_range_beta = solve(MCMC_NNGP$data$covariates$range_X$chol_crossprod_X_locs, q)
    new_compressed_sparse_chol_and_grad =
      compute_sparse_chol(
        num_threads = num_threads, 
        anisotropic = MCMC_NNGP$hierarchical_model$anisotropic, 
        sphere = MCMC_NNGP$hierarchical_model$sphere, 
        range_beta = new_range_beta, NNarray = MCMC_NNGP$vecchia_approx$NNarray, 
        locs = MCMC_NNGP$data$locs, 
        range_X = MCMC_NNGP$data$covariates$range_X$X_locs, 
        nu = MCMC_NNGP$hierarchical_model$nu, 
        PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP,
        compute_derivative = T, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1
      )
    new_sparse_chol = Matrix::sparseMatrix(i = MCMC_NNGP$vecchia_approx$sparse_chol_row_idx, j = MCMC_NNGP$vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][MCMC_NNGP$vecchia_approx$NNarray_non_NA], triangular = T)
    # Make a half step for momentum at the end.
    p = p - 
      as.matrix(
        solve(t(MCMC_NNGP$data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
              - beta_prior_log_dens_derivative
              (beta = new_range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
                chol_crossprod_X = MCMC_NNGP$data$covariates$range_X$chol_crossprod_X,
                beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean,
                beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
                log_scale = MCMC_NNGP$state$params$range_log_scale) # normal prior
              # normal prior derivative                
              + X_PP_crossprod(
                X = MCMC_NNGP$data$covariates$range_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
                Y = # Jacobian of range field wrt range_beta
                  (# natural gradient of obs likelihood wrt range field
                    derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                  left_vector = as.vector(new_sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                                  right_vector = MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                  NNarray = MCMC_NNGP$vecchia_approx$NNarray  
                    )
                    - log_determinant_derivatives(sparse_chol_and_grad = new_compressed_sparse_chol_and_grad, NNarray = MCMC_NNGP$vecchia_approx$NNarray)# derivative of determinant
                  )))) %*% eps_mat / 2
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (MCMC_NNGP$state$momenta$range_beta_sufficient ^2) / 2
    proposed_U =
      (
        - beta_prior_log_dens(
          beta = new_range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
          chol_crossprod_X = MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X,
          beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean,
          beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
          log_scale = MCMC_NNGP$state$params$range_log_scale) # normal prior
        + .5* sum((new_sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale)))^2)
        - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
      )
    
    MCMC_NNGP$state$transition_kernels$range_beta_sufficient[iter%%2+1] = MCMC_NNGP$state$transition_kernels$range_beta_sufficient[iter%%2+1]- 15/(iter_start + iter + 100)
    proposed_K = sum(p^2) / 2
    if(!is.nan(current_U-proposed_U+current_K- proposed_K))
    {
      if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
      {
        MCMC_NNGP$state$transition_kernels$range_beta_sufficient[iter%%2+1] = MCMC_NNGP$state$transition_kernels$range_beta_sufficient[iter%%2+1]+ 20/(iter_start + iter + 100)
        MCMC_NNGP$state$momenta$range_beta_sufficient = p
        MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
        MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
        MCMC_NNGP$state$params$range_beta[] = new_range_beta
        MCMC_NNGP$state$sparse_chol_and_stuff$precision_diag = as.vector((MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][MCMC_NNGP$vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(MCMC_NNGP$vecchia_approx$sparse_chol_column_idx)), j = MCMC_NNGP$vecchia_approx$sparse_chol_column_idx, x = rep(1, length(MCMC_NNGP$vecchia_approx$sparse_chol_row_idx))))
      }
    }
    
    #############################
    # Variance of the  range PP #
    #############################
    
    if(MCMC_NNGP$hierarchical_model$range_PP){
      # ancillary-sufficient ####
      q = MCMC_NNGP$state$params$range_log_scale
      current_U =
        (
          + .5* sum((MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale)))^2)
          - sum(log(MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
        )
      
      # MALA whitened
      MCMC_NNGP$state$momenta$range_log_scale_sufficient = sqrt(.9) * MCMC_NNGP$state$momenta$range_log_scale_sufficient + sqrt(.1)*rnorm(length(q))
      p = MCMC_NNGP$state$momenta$range_log_scale_sufficient
      # Make a half step for momentum at the beginning
      d_beta_d_scale = 
        derivative_field_wrt_scale(
          MCMC_NNGP$state$params$range_beta[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),,drop = F], 
          MCMC_NNGP$state$params$range_log_scale
        )
      # derivative of potential wrt range beta
      d_potential_d_beta = as.matrix(
        + X_PP_crossprod(
          X = MCMC_NNGP$data$covariates$range_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
          Y = # Jacobian of range field wrt range_beta
            (# natural gradient of obs likelihood wrt range field
              derivative_sandwiches(derivatives = MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                            left_vector = as.vector(MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                            right_vector = MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                            NNarray = MCMC_NNGP$vecchia_approx$NNarray  
              )
              - log_determinant_derivatives(sparse_chol_and_grad = MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad, NNarray = MCMC_NNGP$vecchia_approx$NNarray)# derivative of determinant
            ))[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),,drop = F]
      )
      p = p - exp(MCMC_NNGP$state$transition_kernels$range_log_scale_sufficient) *
        (
          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
        )/2
      
      # Make a full step for the position
      q = q + exp(MCMC_NNGP$state$transition_kernels$range_log_scale_sufficient) * p
      new_range_beta = MCMC_NNGP$state$params$range_beta
      new_range_beta[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),] = new_range_beta[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),] %*% 
        solve(chol(expmat(MCMC_NNGP$state$params$range_log_scale))) %*% chol(expmat(q))
      
      new_compressed_sparse_chol_and_grad =
        compute_sparse_chol(
          num_threads = num_threads, 
          anisotropic = MCMC_NNGP$hierarchical_model$anisotropic, 
          sphere = MCMC_NNGP$hierarchical_model$sphere, 
          range_beta = new_range_beta, NNarray = MCMC_NNGP$vecchia_approx$NNarray, 
          locs = MCMC_NNGP$data$locs, 
          range_X = MCMC_NNGP$data$covariates$range_X$X_locs, 
          nu = MCMC_NNGP$hierarchical_model$nu, 
          PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP,
          compute_derivative = T, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1
        )
      new_sparse_chol = Matrix::sparseMatrix(i = MCMC_NNGP$vecchia_approx$sparse_chol_row_idx, j = MCMC_NNGP$vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][MCMC_NNGP$vecchia_approx$NNarray_non_NA], triangular = T)
      # Make a half step for momentum at the end.
      d_beta_d_scale = 
        derivative_field_wrt_scale(
          new_range_beta[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),,drop = F], 
          q
        )
      # derivative of potential wrt range beta
      d_potential_d_beta = as.matrix(
        + X_PP_crossprod(
          X = MCMC_NNGP$data$covariates$range_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
          Y = # Jacobian of range field wrt range_beta
            (# natural gradient of obs likelihood wrt range field
              derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                            left_vector = as.vector(new_sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                            right_vector = MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                            NNarray = MCMC_NNGP$vecchia_approx$NNarray  
              )
              - log_determinant_derivatives(sparse_chol_and_grad = new_compressed_sparse_chol_and_grad, NNarray = MCMC_NNGP$vecchia_approx$NNarray)# derivative of determinant
            ))[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),,drop = F]
      )
      p = p - exp(MCMC_NNGP$state$transition_kernels$range_log_scale_sufficient) *
        (
          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
        )/2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (MCMC_NNGP$state$momenta$range_log_scale_sufficient ^2) / 2
      proposed_U =
        (
          + .5* sum((new_sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale)))^2)
          - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
        )
      proposed_K = sum(p^2) / 2
      MCMC_NNGP$state$transition_kernels$range_log_scale_sufficient  = MCMC_NNGP$state$transition_kernels$range_log_scale_sufficient - 15/(iter_start + iter + 100)
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
        {
          MCMC_NNGP$state$transition_kernels$range_log_scale_sufficient  = MCMC_NNGP$state$transition_kernels$range_log_scale_sufficient + 20/(iter_start + iter + 100)
          new_eigen = eigen(expmat(q))$values
          if(
            all(min(new_eigen)>exp(MCMC_NNGP$hierarchical_model$range_log_scale_prior[1]))&
            all(max(new_eigen)<exp(MCMC_NNGP$hierarchical_model$range_log_scale_prior[2])))
          {
            MCMC_NNGP$state$momenta$range_log_scale_sufficient = p
            
            MCMC_NNGP$state$params$range_beta = new_range_beta
            MCMC_NNGP$state$params$range_log_scale[] = q
            
            MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
            MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
            MCMC_NNGP$state$sparse_chol_and_stuff$precision_diag = as.vector((MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][MCMC_NNGP$vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(MCMC_NNGP$vecchia_approx$sparse_chol_column_idx)), j = MCMC_NNGP$vecchia_approx$sparse_chol_column_idx, x = rep(1, length(MCMC_NNGP$vecchia_approx$sparse_chol_row_idx))))
          }
        }
      }
      
      
      # sufficient - sufficient ####
      for(i in seq(10))
      {
        new_range_log_scale = MCMC_NNGP$state$params$range_log_scale + rnorm(length(MCMC_NNGP$state$params$range_log_scale), 0, .1)
        old_eigen = eigen(expmat(MCMC_NNGP$state$params$range_log_scale))$values
        new_eigen = eigen(expmat(new_range_log_scale))$values
        if(
          all(min(new_eigen)>exp(MCMC_NNGP$hierarchical_model$range_log_scale_prior[1]))&
          all(max(new_eigen)<exp(MCMC_NNGP$hierarchical_model$range_log_scale_prior[2]))&
          (
            + beta_prior_log_dens(
              beta = MCMC_NNGP$state$params$range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
              chol_crossprod_X = MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X,
              beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean,
              beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
              log_scale = new_range_log_scale)
            - beta_prior_log_dens(
              beta = MCMC_NNGP$state$params$range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
              chol_crossprod_X = MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X,
              beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean,
              beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
              log_scale = MCMC_NNGP$state$params$range_log_scale) 
            > log(runif(1))
          )
          
        )
        {
          MCMC_NNGP$state$params$range_log_scale = new_range_log_scale
        }
      }
      # ancillary - ancillary ####
      q = MCMC_NNGP$state$params$range_log_scale
      current_U =
        (
          + .5 * sum((MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals -  MCMC_NNGP$state$params$field[MCMC_NNGP$vecchia_approx$locs_match])^2/MCMC_NNGP$state$sparse_chol_and_stuff$noise) # observation ll
        )
      # MALA whitened
      MCMC_NNGP$state$momenta$range_log_scale_ancillary = sqrt(.9) * MCMC_NNGP$state$momenta$range_log_scale_ancillary + sqrt(.1)*rnorm(length(q))
      p = MCMC_NNGP$state$momenta$range_log_scale_ancillary
      # Make a half step for momentum at the beginning
      d_beta_d_scale = 
        derivative_field_wrt_scale(
          MCMC_NNGP$state$params$range_beta[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),,drop = F], 
          MCMC_NNGP$state$params$range_log_scale
        )
      # derivative of potential wrt range beta
      d_potential_d_beta = X_PP_crossprod(
        X = MCMC_NNGP$data$covariates$range_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
        Y = # Jacobian of range field wrt range_beta
          (
            # natural gradient of obs likelihood wrt range field
            derivative_sandwiches(derivatives = MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                          left_vector = as.vector(
                                            Matrix::solve(
                                              Matrix::t(MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol), 
                                              - as.vector(MCMC_NNGP$vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                            ((MCMC_NNGP$state$params$field[MCMC_NNGP$vecchia_approx$locs_match] - MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals) / MCMC_NNGP$state$sparse_chol_and_stuff$noise))
                                              * sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale) # part of sparse chol
                                            )), 
                                          right_vector = MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                          NNarray = MCMC_NNGP$vecchia_approx$NNarray  
            )
          ))[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),,drop = F]
      p = p - exp(MCMC_NNGP$state$transition_kernels$range_log_scale_ancillary) *
        (
          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
        )/2
      
      q = q + exp(MCMC_NNGP$state$transition_kernels$range_log_scale_ancillary) * p
      new_range_beta = MCMC_NNGP$state$params$range_beta
      new_range_beta[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),] = new_range_beta[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),] %*% 
        solve(chol(expmat(MCMC_NNGP$state$params$range_log_scale))) %*% chol(expmat(q))
      
      new_compressed_sparse_chol_and_grad =
        compute_sparse_chol(
          num_threads = num_threads,
          anisotropic = MCMC_NNGP$hierarchical_model$anisotropic, 
          sphere = MCMC_NNGP$hierarchical_model$sphere, 
          range_beta = new_range_beta, NNarray = MCMC_NNGP$vecchia_approx$NNarray, 
          locs = MCMC_NNGP$data$locs, 
          range_X = MCMC_NNGP$data$covariates$range_X$X_locs, 
          nu = MCMC_NNGP$hierarchical_model$nu, 
          PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP,
          compute_derivative = T, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1
        )
      new_sparse_chol = Matrix::sparseMatrix(i = MCMC_NNGP$vecchia_approx$sparse_chol_row_idx, j = MCMC_NNGP$vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][MCMC_NNGP$vecchia_approx$NNarray_non_NA], triangular = T)
      new_field = sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale))))
      # Make a half step for momentum at the end.
      d_beta_d_scale = 
        derivative_field_wrt_scale(
          new_range_beta[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),,drop = F], 
          q
        )
      # derivative of potential wrt range beta
      d_potential_d_beta = X_PP_crossprod(
        X = MCMC_NNGP$data$covariates$range_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$range_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
        Y = # Jacobian of range field wrt range_beta
          (
            # natural gradient of obs likelihood wrt range field
            derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                          left_vector = as.vector(
                                            Matrix::solve(
                                              Matrix::t(new_sparse_chol), 
                                              - as.vector(MCMC_NNGP$vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                            ((new_field[MCMC_NNGP$vecchia_approx$locs_match] - MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals) / MCMC_NNGP$state$sparse_chol_and_stuff$noise))
                                              * sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale) # part of sparse chol
                                            )), 
                                          right_vector = new_field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                          NNarray = MCMC_NNGP$vecchia_approx$NNarray  
            )
          )
      )[-seq(MCMC_NNGP$data$covariates$range_X$n_regressors),,drop = F]
      p = p - exp(MCMC_NNGP$state$transition_kernels$range_log_scale_ancillary) *
        (
          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
        )/2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (MCMC_NNGP$state$momenta$range_log_scale_ancillary ^2) / 2
      proposed_U =
        (
          + .5 * sum((MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals -  new_field[MCMC_NNGP$vecchia_approx$locs_match])^2/MCMC_NNGP$state$sparse_chol_and_stuff$noise) # observation ll
        )
      proposed_K = sum(p^2) / 2
      MCMC_NNGP$state$transition_kernels$range_log_scale_ancillary = MCMC_NNGP$state$transition_kernels$range_log_scale_ancillary -15/(iter_start + iter + 100)
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
        {
          new_eigen = eigen(expmat(q))$values
          if(
            all(min(new_eigen)>exp(MCMC_NNGP$hierarchical_model$range_log_scale_prior[1]))&
            all(max(new_eigen)<exp(MCMC_NNGP$hierarchical_model$range_log_scale_prior[2])))
          {
            
            MCMC_NNGP$state$transition_kernels$range_log_scale_ancillary = MCMC_NNGP$state$transition_kernels$range_log_scale_ancillary + 20/(iter_start + iter + 100)
            MCMC_NNGP$state$momenta$range_log_scale_ancillary = p
            
            MCMC_NNGP$state$params$range_beta[] = new_range_beta
            MCMC_NNGP$state$params$range_log_scale[] = q
            MCMC_NNGP$state$params$field = new_field
            
            MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
            MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
            MCMC_NNGP$state$sparse_chol_and_stuff$precision_diag = as.vector((MCMC_NNGP$state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][MCMC_NNGP$vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(MCMC_NNGP$vecchia_approx$sparse_chol_column_idx)), j = MCMC_NNGP$vecchia_approx$sparse_chol_column_idx, x = rep(1, length(MCMC_NNGP$vecchia_approx$sparse_chol_row_idx))))
          }
        }
      }
      
      # sufficient - ancillary (equivalent to sufficient-sufficient) ####
      for(i in seq(10))
      {
        new_range_log_scale = MCMC_NNGP$state$params$range_log_scale + rnorm(length(MCMC_NNGP$state$params$range_log_scale), 0, .1)
        old_eigen = eigen(expmat(MCMC_NNGP$state$params$range_log_scale))$values
        new_eigen = eigen(expmat(new_range_log_scale))$values
        if(
          all(min(new_eigen)>exp(MCMC_NNGP$hierarchical_model$range_log_scale_prior[1]))&
          all(max(new_eigen)<exp(MCMC_NNGP$hierarchical_model$range_log_scale_prior[2]))&
          (
            + beta_prior_log_dens(
              beta = MCMC_NNGP$state$params$range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
              chol_crossprod_X = MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X,
              beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean,
              beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
              log_scale = new_range_log_scale)
            - beta_prior_log_dens(
              beta = MCMC_NNGP$state$params$range_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$range_PP, 
              chol_crossprod_X = MCMC_NNGP$data$covariates$range_X$chol_crossprod_X,
              beta0_mean = MCMC_NNGP$hierarchical_model$range_beta0_mean,
              beta0_var =  MCMC_NNGP$hierarchical_model$range_beta0_var, 
              log_scale = MCMC_NNGP$state$params$range_log_scale) 
            > log(runif(1))
          )
          
        )
        {
          MCMC_NNGP$state$params$range_log_scale = new_range_log_scale
        }
      }
    }
    
    
    
    #########
    # Noise #
    #########
    ##############
    # Noise beta #
    ##############
    # recomputation in order to avoid errors
    MCMC_NNGP$state$sparse_chol_and_stuff$noise = variance_field(
      beta = MCMC_NNGP$state$params$noise_beta, X = MCMC_NNGP$data$covariates$noise_X$X, 
      PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$noise_PP, 
      locs_idx = NULL
    )
    # VEWY IMPOWTANT don't remove or comment
    squared_residuals = as.vector(MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals - MCMC_NNGP$state$params$field[MCMC_NNGP$vecchia_approx$locs_match])^2
    # HMC update
    q = MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X %*% MCMC_NNGP$state$params$noise_beta
    current_U =
      (
        - beta_prior_log_dens(
          beta = MCMC_NNGP$state$params$noise_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$noise_PP, 
          chol_crossprod_X = MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X,
          beta0_mean = MCMC_NNGP$hierarchical_model$noise_beta0_mean,
          beta0_var =  MCMC_NNGP$hierarchical_model$noise_beta0_var, 
          log_scale = MCMC_NNGP$state$params$noise_log_scale) # normal prior 
        +.5* sum(log(MCMC_NNGP$state$sparse_chol_and_stuff$noise)) # det
        +.5*sum(squared_residuals/MCMC_NNGP$state$sparse_chol_and_stuff$noise) # observations
      )
    # HMC whitened
    MCMC_NNGP$state$momenta$noise_beta = sqrt(.9) * MCMC_NNGP$state$momenta$noise_beta + sqrt(.1)*rnorm(length(MCMC_NNGP$state$momenta$noise_beta))
    p = MCMC_NNGP$state$momenta$noise_beta
    
    # Make a half step for momentum at the beginning
    p = p - exp(MCMC_NNGP$state$transition_kernels$noise_beta_mala) *
      (
        + solve(t(MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X), # solving by prior chol because of whitening
                - beta_prior_log_dens_derivative(
                  beta = MCMC_NNGP$state$params$noise_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$noise_PP, 
                  chol_crossprod_X = MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X,
                  beta0_mean = MCMC_NNGP$hierarchical_model$noise_beta0_mean,
                  beta0_var =  MCMC_NNGP$hierarchical_model$noise_beta0_var, 
                  log_scale = MCMC_NNGP$state$params$noise_log_scale) # normal prior
                + X_PP_crossprod(X = MCMC_NNGP$data$covariates$noise_X$X, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$noise_PP, 
                                         Y = 
                                           (
                                             + .5 # determinant part of normal likelihood
                                             - (squared_residuals/MCMC_NNGP$state$sparse_chol_and_stuff$noise)/2 # exponential part of normal likelihood
                                           ))
        ))/ 2
    # checking gradient with finite differences, to update
    
    #### idx = 200
    #### noise_beta_ = MCMC_NNGP$state$params$noise_beta 
    #### noise_beta_[idx] = noise_beta_[idx] + 0.0000001
    #### noise_ = variance_field(
    ####   beta = noise_beta_, X = MCMC_NNGP$data$covariates$noise_X$X, 
    ####   PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$noise_PP
    #### )
    #### U_ =
    ####   (
    ####     - beta_prior_log_dens(beta = noise_beta_, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$noise_PP, 
    ####                                   beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$noise_beta_mean, 
    ####                                   beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$noise_beta_precision, 
    ####                                   log_scale = MCMC_NNGP$state$params$noise_log_scale) # normal prior 
    ####     +.5* sum(log(noise_)) # det
    ####     +.5*sum(squared_residuals/noise_) # observations
    ####   )
    #### print(
    ####   (
    ####             - beta_prior_log_dens_derivative(beta = MCMC_NNGP$state$params$noise_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$noise_PP, 
    ####                                                      beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$noise_beta_mean, 
    ####                                                      beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$noise_beta_precision, 
    ####                                                      log_scale = MCMC_NNGP$state$params$noise_log_scale) # normal prior
    ####             + X_PP_crossprod(X = MCMC_NNGP$data$covariates$noise_X$X, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$noise_PP, 
    ####                                      Y = 
    ####                                        (
    ####                                          + .5 # determinant part of normal likelihood
    ####                                          - (squared_residuals/MCMC_NNGP$state$sparse_chol_and_stuff$noise)/2 # exponential part of normal likelihood
    ####                                        ))
    ####     )[idx]
    #### )
    #### print(10000000*(U_- current_U))
    
    
    # Make a full step for the position
    q = q + exp(MCMC_NNGP$state$transition_kernels$noise_beta_mala) * p
    new_noise_beta = solve(MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X, q)
    new_noise = variance_field(beta = new_noise_beta, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$noise_PP, 
                                       X = MCMC_NNGP$data$covariates$noise_X$X, locs_idx = NULL)
    # Make a half step for momentum at the end
    p = p - exp(MCMC_NNGP$state$transition_kernels$noise_beta_mala) *
      (
        + solve(t(MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X), # solving by prior sparse chol because of whitening
                - beta_prior_log_dens_derivative(
                  beta = new_noise_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$noise_PP, 
                  chol_crossprod_X = MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X,
                  beta0_mean = MCMC_NNGP$hierarchical_model$noise_beta0_mean,
                  beta0_var =  MCMC_NNGP$hierarchical_model$noise_beta0_var, 
                  log_scale = MCMC_NNGP$state$params$noise_log_scale) # normal prior  
                + X_PP_crossprod(X = MCMC_NNGP$data$covariates$noise_X$X, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$noise_PP, 
                                         (
                                           + .5 # determinant part of normal likelihood
                                           - (squared_residuals/new_noise)/2 # exponential part of normal likelihood
                                         ))
        ))/ 2
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (MCMC_NNGP$state$momenta$noise_beta ^2) / 2
    proposed_U = 
      (
        - beta_prior_log_dens(beta = new_noise_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$noise_PP, 
                                      chol_crossprod_X = MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X,
                                      beta0_mean = MCMC_NNGP$hierarchical_model$noise_beta0_mean,
                                      beta0_var =  MCMC_NNGP$hierarchical_model$noise_beta0_var, 
                                      log_scale = MCMC_NNGP$state$params$noise_log_scale) # normal prior        
        +.5* sum(log(new_noise)) # det
        +.5*sum(squared_residuals/new_noise) # observations
      )
    proposed_K = sum(p^2) / 2
    
    
    MCMC_NNGP$state$transition_kernels$noise_beta_mala = MCMC_NNGP$state$transition_kernels$noise_beta_mala- 1/sqrt(iter_start + iter +100)
    if(!is.nan(current_U-proposed_U+current_K- proposed_K))
    {
      if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
      {
        MCMC_NNGP$state$transition_kernels$noise_beta_mala = MCMC_NNGP$state$transition_kernels$noise_beta_mala + 2/sqrt(iter_start + iter +100)
        MCMC_NNGP$state$momenta$noise_beta = p
        MCMC_NNGP$state$params$noise_beta[] = new_noise_beta
        MCMC_NNGP$state$sparse_chol_and_stuff$noise = new_noise
      }
    }
    
    ###################
    # Noise log scale # 
    ###################
    if(MCMC_NNGP$hierarchical_model$noise_PP)
    {
      # ancillary -- sufficient ####
      new_noise_log_scale = MCMC_NNGP$state$params$noise_log_scale + rnorm(1, 0, exp(MCMC_NNGP$state$transition_kernels$noise_log_scale))
      new_noise_beta = MCMC_NNGP$state$params$noise_beta
      new_noise_beta[-seq(MCMC_NNGP$data$covariates$noise_X$n_regressors)] = new_noise_beta[-seq(MCMC_NNGP$data$covariates$noise_X$n_regressors)] *
        exp((new_noise_log_scale - MCMC_NNGP$state$params$noise_log_scale)/2)
      new_noise =variance_field(beta = new_noise_beta, PP = MCMC_NNGP$hierarchical_model$PP, 
                                        use_PP = MCMC_NNGP$hierarchical_model$noise_PP, X = MCMC_NNGP$data$covariates$noise_X$X, 
                                        locs_idx = NULL)
      ll_ratio = (
        -.5* sum(log(new_noise)) 
        -.5*sum(squared_residuals/new_noise)
        +.5* sum(log(MCMC_NNGP$state$sparse_chol_and_stuff$noise)) 
        +.5*sum(squared_residuals/MCMC_NNGP$state$sparse_chol_and_stuff$noise)
      )
      if(!is.nan(ll_ratio))
      {
        if(ll_ratio > log(runif(1)))
        {
          if(
            (new_noise_log_scale > MCMC_NNGP$hierarchical_model$noise_log_scale_prior[1])&
            (new_noise_log_scale < MCMC_NNGP$hierarchical_model$noise_log_scale_prior[2])
          )
          {
            MCMC_NNGP$state$params$noise_log_scale = new_noise_log_scale
            MCMC_NNGP$state$params$noise_beta = new_noise_beta 
            MCMC_NNGP$state$sparse_chol_and_stuff$noise = new_noise
            MCMC_NNGP$state$transition_kernels$noise_log_scale = MCMC_NNGP$state$transition_kernels$noise_log_scale + 4/sqrt(iter_start + iter +100)
          }
        }
      }
      MCMC_NNGP$state$transition_kernels$noise_log_scale = MCMC_NNGP$state$transition_kernels$noise_log_scale - 1/sqrt(iter_start + iter +100)
      
      # sufficient -- sufficient ####
      for(i in seq(10))
      {
        new_noise_log_scale = MCMC_NNGP$state$params$noise_log_scale + rnorm(1, 0, .1)
        if(
          (
            + beta_prior_log_dens(beta = MCMC_NNGP$state$params$noise_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$noise_PP, 
                                          chol_crossprod_X = MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X,
                                          beta0_mean = MCMC_NNGP$hierarchical_model$noise_beta0_mean,
                                          beta0_var =  MCMC_NNGP$hierarchical_model$noise_beta0_var, 
                                          log_scale = new_noise_log_scale) - 
            beta_prior_log_dens(beta = MCMC_NNGP$state$params$noise_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$noise_PP, 
                                        chol_crossprod_X = MCMC_NNGP$data$covariates$noise_X$chol_crossprod_X,
                                        beta0_mean = MCMC_NNGP$hierarchical_model$noise_beta0_mean,
                                        beta0_var =  MCMC_NNGP$hierarchical_model$noise_beta0_var, 
                                        log_scale = MCMC_NNGP$state$params$noise_log_scale)
          ) > log(runif(1))
          &(new_noise_log_scale > MCMC_NNGP$hierarchical_model$noise_log_scale_prior[1])
          &(new_noise_log_scale < MCMC_NNGP$hierarchical_model$noise_log_scale_prior[2])
        )
        {
          MCMC_NNGP$state$params$noise_log_scale = new_noise_log_scale
        }
      }
    }
    
    
    
    #########
    # Scale #
    #########
    {
      ##############################
      # scale beta sufficient MALA #
      ##############################
      sparse_chol_diag_field = MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% Matrix::Diagonal(x = MCMC_NNGP$state$params$field)
      q = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X_locs %*% MCMC_NNGP$state$params$scale_beta
      current_U =
        (
          - beta_prior_log_dens(beta = MCMC_NNGP$state$params$scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                        chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                        beta0_mean = MCMC_NNGP$hierarchical_model$scale_beta0_mean,
                                        beta0_var =  MCMC_NNGP$hierarchical_model$scale_beta0_var, 
                                        log_scale = MCMC_NNGP$state$params$scale_log_scale) # normal prior 
          +0.5*sum(log(MCMC_NNGP$state$sparse_chol_and_stuff$scale))# determinant part
          +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/MCMC_NNGP$state$sparse_chol_and_stuff$scale))^2)# covmat product part
        )
      # HMC whitened
      MCMC_NNGP$state$momenta$scale_beta_sufficient = sqrt(.9) * MCMC_NNGP$state$momenta$scale_beta_sufficient + sqrt(.1)*rnorm(length(MCMC_NNGP$state$momenta$scale_beta_sufficient))
      p = MCMC_NNGP$state$momenta$scale_beta_sufficient
      # Make a half step for momentum at the beginning
      p = p - exp(MCMC_NNGP$state$transition_kernels$scale_beta_sufficient_mala) *
        (
          + solve(t(MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X_locs), # solving by t chol because of whitening
                  - beta_prior_log_dens_derivative(beta = MCMC_NNGP$state$params$scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                                           chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                                           beta0_mean = MCMC_NNGP$hierarchical_model$scale_beta0_mean,
                                                           beta0_var =  MCMC_NNGP$hierarchical_model$scale_beta0_var, 
                                                           log_scale = MCMC_NNGP$state$params$scale_log_scale) # normal prior
                  +  X_PP_crossprod(X = MCMC_NNGP$data$covariates$scale_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$scale_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
                                            Y = (
                                              .5  # determinant part 
                                              -.5 * sqrt(1/MCMC_NNGP$state$sparse_chol_and_stuff$scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/MCMC_NNGP$state$sparse_chol_and_stuff$scale)))# natural derivative
                                            ))
          ))/ 2
      
      # Make a full step for the position
      q = q + exp(MCMC_NNGP$state$transition_kernels$scale_beta_sufficient_mala) * p
      new_scale_beta = solve(MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X_locs, q)
      new_scale =variance_field(beta = new_scale_beta, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$scale_PP, 
                                        X = MCMC_NNGP$data$covariates$scale_X$X_locs, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1)
      
      # Make a half step for momentum at the end.
      p = p - exp(MCMC_NNGP$state$transition_kernels$scale_beta_sufficient_mala) *
        (
          + solve(t(MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                  - beta_prior_log_dens_derivative(beta = new_scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                                           chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                                           beta0_mean = MCMC_NNGP$hierarchical_model$scale_beta0_mean,
                                                           beta0_var =  MCMC_NNGP$hierarchical_model$scale_beta0_var, 
                                                           log_scale = MCMC_NNGP$state$params$scale_log_scale) # normal prior  
                  + X_PP_crossprod(X = MCMC_NNGP$data$covariates$scale_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$scale_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
                                           (
                                             .5  # determinant part 
                                             -.5 * sqrt(1/new_scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/new_scale)))# natural derivative
                                           ))
          ))/ 2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (MCMC_NNGP$state$momenta$scale_beta_sufficient ^2) / 2
      proposed_U = 
        (
          - beta_prior_log_dens(beta = new_scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                        chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                        beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_mean, 
                                        beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_precision, 
                                        log_scale = MCMC_NNGP$state$params$scale_log_scale) # normal prior 
          +0.5*sum(log(new_scale))# determinant part
          +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/new_scale))^2)# covmat product part
        )
      proposed_K = sum(p^2) / 2
      MCMC_NNGP$state$transition_kernels$scale_beta_sufficient_mala = MCMC_NNGP$state$transition_kernels$scale_beta_sufficient_mala - 1/sqrt(iter_start + iter +100)
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
        {
          MCMC_NNGP$state$transition_kernels$scale_beta_sufficient_mala = MCMC_NNGP$state$transition_kernels$scale_beta_sufficient_mala + 2/sqrt(iter_start + iter +100)
          MCMC_NNGP$state$momenta$scale_beta_sufficient = p
          MCMC_NNGP$state$params$scale_beta[] = new_scale_beta
          MCMC_NNGP$state$sparse_chol_and_stuff$scale = new_scale
        }
      }
      
      #############################
      # scale beta ancillary MALA #
      #############################
      q = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X_locs %*% MCMC_NNGP$state$params$scale_beta
      current_U =
        (
          - beta_prior_log_dens(beta = MCMC_NNGP$state$params$scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                        chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                        beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_mean, 
                                        beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_precision, 
                                        log_scale = MCMC_NNGP$state$params$scale_log_scale) # normal prior 
          +.5 * sum((MCMC_NNGP$data$observed_field - MCMC_NNGP$state$sparse_chol_and_stuff$lm_fit - MCMC_NNGP$state$params$field[MCMC_NNGP$vecchia_approx$locs_match] )^2/MCMC_NNGP$state$sparse_chol_and_stuff$noise)
        )
      # MALA whitened
      MCMC_NNGP$state$momenta$scale_beta_ancillary = sqrt(.9) * MCMC_NNGP$state$momenta$scale_beta_ancillary   + sqrt(.1)*rnorm(length(MCMC_NNGP$state$momenta$scale_beta_ancillary))
      p = MCMC_NNGP$state$momenta$scale_beta_ancillary
      # Make a. half step for momentum at the beginning
      p = p - exp(MCMC_NNGP$state$transition_kernels$scale_beta_ancillary_mala) *
        (
          + solve(t(MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X_locs), # solving by prior chol because of whitening
                  - beta_prior_log_dens_derivative(beta = MCMC_NNGP$state$params$scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                                           chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                                           beta0_mean = MCMC_NNGP$hierarchical_model$scale_beta0_mean,
                                                           beta0_var =  MCMC_NNGP$hierarchical_model$scale_beta0_var, 
                                                           log_scale = MCMC_NNGP$state$params$scale_log_scale) # normal prior 
                  + X_PP_crossprod(X = MCMC_NNGP$data$covariates$scale_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$scale_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
                                           (.5 * MCMC_NNGP$state$params$field * 
                                              as.vector(MCMC_NNGP$vecchia_approx$locs_match_matrix %*% 
                                                          ((MCMC_NNGP$state$params$field[MCMC_NNGP$vecchia_approx$locs_match] - MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals)/
                                                             MCMC_NNGP$state$sparse_chol_and_stuff$noise)
                                              )))
          ))/ 2
      # Make a full step for the position
      q = q + exp(MCMC_NNGP$state$transition_kernels$scale_beta_ancillary_mala) * p
      new_scale_beta = solve(MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X_locs, q)
      new_scale =variance_field(beta = new_scale_beta, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$scale_PP, 
                                        X = MCMC_NNGP$data$covariates$scale_X$X_locs, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1)
      new_field  = MCMC_NNGP$state$params$field*sqrt(new_scale)/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale)
      # Make a half step for momentum at the end.
      p = p - exp(MCMC_NNGP$state$transition_kernels$scale_beta_ancillary_mala) *
        (
          + solve(t(MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                  - beta_prior_log_dens_derivative(beta = new_scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                                           chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                                           beta0_mean = MCMC_NNGP$hierarchical_model$scale_beta0_mean,
                                                           beta0_var =  MCMC_NNGP$hierarchical_model$scale_beta0_var, 
                                                           log_scale = MCMC_NNGP$state$params$scale_log_scale) # normal prior  
                  + X_PP_crossprod(X = MCMC_NNGP$data$covariates$scale_X$X_locs, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$scale_PP, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1,
                                           (.5 * new_field * 
                                              as.vector(MCMC_NNGP$vecchia_approx$locs_match_matrix %*% 
                                                          ((new_field[MCMC_NNGP$vecchia_approx$locs_match] - MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals)/
                                                             MCMC_NNGP$state$sparse_chol_and_stuff$noise)
                                              )))
          ))/ 2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (MCMC_NNGP$state$momenta$scale_beta_ancillary ^2) / 2
      proposed_U = 
        (
          - beta_prior_log_dens(beta = new_scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                        chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                        beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_mean, 
                                        beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_precision, 
                                        log_scale = MCMC_NNGP$state$params$scale_log_scale) # normal prior  
          +.5 * sum((MCMC_NNGP$data$observed_field - MCMC_NNGP$state$sparse_chol_and_stuff$lm_fit - new_field[MCMC_NNGP$vecchia_approx$locs_match] )^2/MCMC_NNGP$state$sparse_chol_and_stuff$noise)
        )
      proposed_K = sum(p^2) / 2
      MCMC_NNGP$state$transition_kernels$scale_beta_ancillary_mala = MCMC_NNGP$state$transition_kernels$scale_beta_ancillary_mala - 1/sqrt(iter_start + iter +100)
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
        {
          MCMC_NNGP$state$transition_kernels$scale_beta_ancillary_mala = MCMC_NNGP$state$transition_kernels$scale_beta_ancillary_mala + 2/sqrt(iter_start + iter +100)
          MCMC_NNGP$state$momenta$scale_beta_ancillary = p
          MCMC_NNGP$state$params$scale_beta[] = new_scale_beta
          MCMC_NNGP$state$params$field = new_field
          MCMC_NNGP$state$sparse_chol_and_stuff$scale = new_scale
        }
      }
      
      ###################
      # scale log scale #
      ###################
      if(MCMC_NNGP$hierarchical_model$scale_PP){
        # ancillary -- sufficient  ####
        # a change in hyperprior scale changes (rescales) the scale, which is then compared with the latent field
        new_scale_log_scale = MCMC_NNGP$state$params$scale_log_scale + rnorm(1, 0, exp(MCMC_NNGP$state$transition_kernels$scale_log_scale_sufficient)) 
        new_scale_beta = MCMC_NNGP$state$params$scale_beta
        new_scale_beta[-seq(MCMC_NNGP$data$covariates$scale_X$n_regressors)] = new_scale_beta[-seq(MCMC_NNGP$data$covariates$scale_X$n_regressors)] *
          exp((new_scale_log_scale - MCMC_NNGP$state$params$scale_log_scale)/2)
        new_scale = variance_field(beta = new_scale_beta, PP = MCMC_NNGP$hierarchical_model$PP, 
                                           use_PP = MCMC_NNGP$hierarchical_model$scale_PP, X = MCMC_NNGP$data$covariates$scale_X$X_locs, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1)
        MCMC_NNGP$state$transition_kernels$scale_log_scale_sufficient = MCMC_NNGP$state$transition_kernels$scale_log_scale_sufficient - .25/sqrt(iter_start + iter +100)
        if(
          (
            (
              +.5* sum(log(MCMC_NNGP$state$sparse_chol_and_stuff$scale)) 
              -.5* sum(log(new_scale)) # log determinant
              +.5*sum((MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale)))^2)
              -.5*sum((MCMC_NNGP$state$sparse_chol_and_stuff$sparse_chol %*% (MCMC_NNGP$state$params$field/sqrt(new_scale)))^2) # Gaussian density of the latent field
            )
            > log(runif(1))
          )
        )
        {
          if(
            (new_scale_log_scale > MCMC_NNGP$hierarchical_model$scale_log_scale_prior[1])&
            (new_scale_log_scale < MCMC_NNGP$hierarchical_model$scale_log_scale_prior[2])
          )
          {
            MCMC_NNGP$state$transition_kernels$scale_log_scale_sufficient = MCMC_NNGP$state$transition_kernels$scale_log_scale_sufficient + 1/sqrt(iter_start + iter +100)
            MCMC_NNGP$state$params$scale_log_scale = new_scale_log_scale
            MCMC_NNGP$state$params$scale_beta = new_scale_beta 
            MCMC_NNGP$state$sparse_chol_and_stuff$scale = new_scale
          }
        }
        
        
        # sufficient -- sufficient ####
        for(i in seq(4))
        {
          new_scale_log_scale = MCMC_NNGP$state$params$scale_log_scale + rnorm(1, 0, .1)
          if(
            (
              + beta_prior_log_dens(beta = MCMC_NNGP$state$params$scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                            chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                            beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_mean, 
                                            beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_precision, 
                                            log_scale = new_scale_log_scale)     
              - beta_prior_log_dens(beta = MCMC_NNGP$state$params$scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                            chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                            beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_mean, 
                                            beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_precision, 
                                            log_scale = MCMC_NNGP$state$params$scale_log_scale)     
              > log(runif(1))
            )
            &(new_scale_log_scale > MCMC_NNGP$hierarchical_model$scale_log_scale_prior[1])
            &(new_scale_log_scale < MCMC_NNGP$hierarchical_model$scale_log_scale_prior[2])
          )
          {
            MCMC_NNGP$state$params$scale_log_scale = new_scale_log_scale
          }
        }
        # ancillary -- ancillary ####
        new_scale_log_scale = MCMC_NNGP$state$params$scale_log_scale + rnorm(1, 0, exp(MCMC_NNGP$state$transition_kernels$scale_log_scale_ancillary)) 
        new_scale_beta = MCMC_NNGP$state$params$scale_beta
        new_scale_beta[-seq(MCMC_NNGP$data$covariates$scale_X$n_regressors)] = new_scale_beta[-seq(MCMC_NNGP$data$covariates$scale_X$n_regressors)] *
          exp((new_scale_log_scale - MCMC_NNGP$state$params$scale_log_scale)/2)
        new_scale =variance_field(new_scale_beta, PP = MCMC_NNGP$hierarchical_model$PP, use_PP = MCMC_NNGP$hierarchical_model$scale_PP, 
                                          X = MCMC_NNGP$data$covariates$scale_X$X_locs, locs_idx = MCMC_NNGP$vecchia_approx$hctam_scol_1)
        new_field = MCMC_NNGP$state$params$field * sqrt(new_scale)/sqrt(MCMC_NNGP$state$sparse_chol_and_stuff$scale)
        MCMC_NNGP$state$transition_kernels$scale_log_scale_ancillary = MCMC_NNGP$state$transition_kernels$scale_log_scale_ancillary - .25/sqrt(iter_start + iter +100)
        if(
          (
            -.5* sum((MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals -          new_field[MCMC_NNGP$vecchia_approx$locs_match])^2/MCMC_NNGP$state$sparse_chol_and_stuff$noise) 
            +.5* sum((MCMC_NNGP$state$sparse_chol_and_stuff$lm_residuals - MCMC_NNGP$state$params$field[MCMC_NNGP$vecchia_approx$locs_match])^2/MCMC_NNGP$state$sparse_chol_and_stuff$noise) 
            > log(runif(1))
          )
        )
        {
          if(
            (new_scale_log_scale >  MCMC_NNGP$hierarchical_model$scale_log_scale_prior[1])
            &(new_scale_log_scale < MCMC_NNGP$hierarchical_model$scale_log_scale_prior[2])
          )
          {  
            MCMC_NNGP$state$transition_kernels$scale_log_scale_ancillary = MCMC_NNGP$state$transition_kernels$scale_log_scale_ancillary + 1/sqrt(iter_start + iter +100)
            MCMC_NNGP$state$params$field = new_field
            MCMC_NNGP$state$params$scale_log_scale = new_scale_log_scale
            MCMC_NNGP$state$params$scale_beta = new_scale_beta 
            MCMC_NNGP$state$sparse_chol_and_stuff$scale = new_scale
          }
        }
        # sufficient -- ancillary ####
        for(i in seq(4))
        {
          new_scale_log_scale = MCMC_NNGP$state$params$scale_log_scale + rnorm(1, 0, .1)
          if(
            (
              + beta_prior_log_dens(beta = MCMC_NNGP$state$params$scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                            chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                            beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_mean, 
                                            beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_precision, 
                                            log_scale = new_scale_log_scale)     
              - beta_prior_log_dens(beta = MCMC_NNGP$state$params$scale_beta, n_PP = MCMC_NNGP$hierarchical_model$PP$n_PP*MCMC_NNGP$hierarchical_model$scale_PP, 
                                            chol_crossprod_X = MCMC_NNGP$data$covariates$scale_X$chol_crossprod_X,
                                            beta_mean = MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_mean, 
                                            beta_precision =  MCMC_NNGP$hierarchical_model$beta_priors$scale_beta_precision, 
                                            log_scale = MCMC_NNGP$state$params$scale_log_scale)     
              > log(runif(1))
            )
            &(new_scale_log_scale > MCMC_NNGP$hierarchical_model$scale_log_scale_prior[1])
            &(new_scale_log_scale < MCMC_NNGP$hierarchical_model$scale_log_scale_prior[2])
          )
          {
            MCMC_NNGP$state$params$scale_log_scale = new_scale_log_scale
          }
        }
        
      }
      
    }
    
    
    #######################
    # Storing the samples #
    #######################
    for(name in names(MCMC_NNGP$state$params))
    {
      params_records[[name]][,,iter] = MCMC_NNGP$state$params[[name]]
    }
  }
  return(list("state" = MCMC_NNGP$state, "params_records" = params_records))
}
