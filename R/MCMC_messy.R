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



if(F){
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
  list2env(state, environment())
  #########################################
  # Initializing chain storage structures #
  #########################################
  # this part re-creates a small portion of the $records objects of each chain. It fills it with chain state during the run, and then updates each chain with the new values
  params_records = list()
  
  #par(mfrow = c(2, 1))
  
  t1 = Sys.time()
  for(iter in seq(1, n_iterations_update)){
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
       beta_mean = c(( ((observed_field-params$field[vecchia_approx$locs_match]) / stuff$noise_var) %*% covariates$X$X) %*% beta_covmat)
       params$beta[]   = c(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
     }}
    
    # un-centered parametrization of latent field
    centered_field = as.vector(params$field + covariates$X$X_locs%*%matrix(params$beta[covariates$X$which_locs], ncol = 1))
    sparse_chol_X = as.matrix(stuff$sparse_chol %*% (covariates$X$X_locs/stuff$field_sd))
    beta_precision = crossprod(x = sparse_chol_X, y = sparse_chol_X)
    beta_covmat = solve(beta_precision, tol = min(rcond(beta_precision),.Machine$double.eps))
    if(all(!is.infinite(beta_covmat) & !is.nan(beta_covmat)))
    {
      if(all(eigen(beta_covmat)$d >0))
      {
        beta_mean =  c(as.vector(stuff$sparse_chol %*% (centered_field/stuff$field_sd))  %*% sparse_chol_X %*% beta_covmat)
        params$beta[covariates$X$which_locs]   = as.vector(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
        params$field = centered_field - as.vector(covariates$X$X_locs %*% matrix(params$beta[covariates$X$which_locs], ncol = 1))
      }}
    
    # updating stuff 
    stuff$lm_fit[]       = as.vector(covariates$X$X%*%     params$beta)
    stuff$lm_residuals = observed_field - stuff$lm_fit
    
    ################
    # Latent field #
    ################
    t1 = Sys.time()
    additional_precision = vecchia_approx$locs_match_matrix %*% (1/stuff$noise_var)
    additional_mean = as.vector(vecchia_approx$locs_match_matrix %*% ((observed_field - stuff$lm_fit)/stuff$noise_var))
    posterior_precision = Matrix::crossprod(stuff$sparse_chol)
    dp = diff(posterior_precision@p)
    j = rep(seq_along(dp), dp)
    posterior_precision@x = posterior_precision@x / (stuff$field_sd[posterior_precision@i+1] * stuff$field_sd[j])
    Matrix::diag(posterior_precision) = Matrix::diag(posterior_precision) + as.vector(additional_precision)
    posterior_precision_chol = Matrix::expand(Matrix::Cholesky(posterior_precision))
      params$field = as.vector(
        Matrix::t(posterior_precision_chol$P) %*% 
          Matrix::solve(Matrix::t(posterior_precision_chol$L), 
                        rnorm(vecchia_approx$n_locs) + 
                          Matrix::solve(posterior_precision_chol$L, # inverse of precision matrix...
                                        posterior_precision_chol$P %*%
                                          (additional_mean))
          ))
    print(Sys.time()-t1)
    
    ###############
    # Range beta  #
    ###############
    ##########################
    # Range beta (ancillary) #
    ##########################
##    if(is.null(hierarchical_model$range_PP))n_PP = 0
##    if(!is.null(hierarchical_model$range_PP))n_PP = hierarchical_model$range_PP$n_knots
##    hmc_stepsize = matrix(c(
##      rep(ker_var$range_beta_ancillary[1], covariates$range_X$n_regressors),
##      rep(ker_var$range_beta_ancillary[1], n_PP),
##      rep(ker_var$range_beta_ancillary[3], covariates$range_X$n_regressors),
##      rep(ker_var$range_beta_ancillary[4], n_PP),
##      rep(ker_var$range_beta_ancillary[3], covariates$range_X$n_regressors),
##      rep(ker_var$range_beta_ancillary[4], n_PP)
##    ), ncol=3)[,seq(1+2*hierarchical_model$anisotropic)]
##    
##    q = covariates$range_X$chol_crossprod_X_locs %*% params$range_beta # whitening wrt covariates of the range
##    current_U =
##      (
##        - beta_prior_log_dens(beta = params$range_beta, 
##                              n_PP = n_PP, 
##                              chol_crossprod_X = covariates$range_X$chol_crossprod_X,
##                              beta0_mean = hierarchical_model$range_beta0_mean, 
##                              beta0_var =  hierarchical_model$range_beta0_var, 
##                              log_scale = params$range_log_scale) # normal prior
##        + .5 * sum((stuff$lm_residuals -  params$field[vecchia_approx$locs_match])^2/stuff$noise) # observation ll
##      )
##    # MALA whitened
##    momenta$range_beta_ancillary = renew_momentum(momenta$range_beta_ancillary)
##    p = momenta$range_beta_ancillary
##    # Make a half step for momentum at the beginning
##    p = p - 
##      as.matrix(
##        solve(t(range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
##              - beta_prior_log_dens_derivative(
##                beta = params$range_beta, n_PP = n_PP, 
##                chol_crossprod_X = covariates$range_X$chol_crossprod_X,
##                beta0_mean = hierarchical_model$range_beta0_mean,
##                beta0_var =  hierarchical_model$range_beta0_var, 
##                log_scale = params$range_log_scale) # normal prior
##              # normal prior derivative                
##              + X_PP_crossprod(
##                X = range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
##                Y = # Jacobian of range field wrt range_beta
##                  (
##                    # natural gradient of obs likelihood wrt range field
##                    derivative_sandwiches(derivatives = stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
##                                          left_vector = as.vector(
##                                            Matrix::solve(
##                                              Matrix::t(stuff$sparse_chol), 
##                                              - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
##                                                            ((params$field[vecchia_approx$locs_match] - stuff$lm_residuals) / stuff$noise))
##                                              * sqrt(stuff$scale) # part of sparse chol
##                                            )), 
##                                          right_vector = params$field/sqrt(stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
##                                          NNarray = vecchia_approx$NNarray  
##                    )
##                  ))
##        ))%*%eps_mat/ 2
##    
##    
##    #######testing the gradient
##    ##source("Bidart/R/Useful_stuff.R")
##    ##d1 =         - beta_prior_log_dens_derivative(
##    ##  beta = params$range_beta, n_PP = n_PP, 
##    ##  beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
##    ##  beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
##    ##  log_scale = params$range_log_scale) # normal prior
    ##beta_ = params$range_beta 
    ##beta_[2, 2] = beta_[2, 2] + 1/10000
    ##d2 = 10000*(
    ##  -
    ##  beta_prior_log_dens(beta = beta_, n_PP = n_PP, 
    ##                        beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
    ##                        beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
    ##                        log_scale = params$range_log_scale) # normal prior
    ##  + 
    ##  beta_prior_log_dens(beta = params$range_beta, n_PP = n_PP, 
    ##                        beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
    ##                        beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
    ##                        log_scale = params$range_log_scale) # normal prior
    ##)
    ##d1/d2
    #### i = 4
    #### j = 2
    #### 
    #### range_beta_ = params$range_beta
    #### range_beta_[i,j] = params$range_beta[i,j] + .0001
    #### 
    #### compressed_sparse_chol_and_grad_ = 
    ####   compute_sparse_chol(
    ####     num_threads = num_threads, 
    ####     anisotropic = hierarchical_model$anisotropic, 
    ####     sphere = hierarchical_model$sphere, 
    ####     range_beta = range_beta_, NNarray = vecchia_approx$NNarray, 
    ####     locs = locs, 
    ####     range_X = range_X$X_locs, 
    ####     nu = hierarchical_model$nu, 
    ####     PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
    ####     compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
    ####   )
    #### sparse_chol_ = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = compressed_sparse_chol_and_grad_[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
    #### field_ = sqrt(stuff$scale) * as.vector(Matrix::solve(sparse_chol_, stuff$sparse_chol %*% (params$field/sqrt(stuff$scale))))
    #### 
    #### 10000*(
    ####   + .5 * sum((stuff$lm_residuals -  field_[vecchia_approx$locs_match])^2/stuff$noise) - # observation ll
    ####     + .5 * sum((stuff$lm_residuals -  params$field[vecchia_approx$locs_match])^2/stuff$noise) # observation ll
    #### )
    #### 
    #### 
    #### + X_PP_crossprod(
    ####   X = range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
    ####   Y = # Jacobian of range field wrt range_beta
    ####     (
    ####       # natural gradient of obs likelihood wrt range field
    ####       derivative_sandwiches(derivatives = stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
    ####                                     left_vector = as.vector(
    ####                                       Matrix::solve(
    ####                                         Matrix::t(stuff$sparse_chol), 
    ####                                         - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
    ####                                                       ((params$field[vecchia_approx$locs_match] - stuff$lm_residuals) / stuff$noise))
    ####                                         * sqrt(stuff$scale) # part of sparse chol
    ####                                       )), 
    ####                                     right_vector = params$field/sqrt(stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
    ####                                     NNarray = vecchia_approx$NNarray  
    ####       )
####    ####     ))[i,j]
####    
####    
####    # Make a full step for the position
####    q = q + p %*%eps_mat
####    new_range_beta = params$range_beta
####    new_range_beta = solve(covariates$range_X$chol_crossprod_X_locs, q )
####    new_compressed_sparse_chol_and_grad = 
####      compute_sparse_chol(
####        num_threads = num_threads, 
####        anisotropic = hierarchical_model$anisotropic, 
####        sphere = hierarchical_model$sphere, 
####        range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
####        locs = locs, 
####        range_X = range_X$X_locs, 
####        nu = hierarchical_model$nu, 
####        PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
####        compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
####      )
####    new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
####    new_field = sqrt(stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, stuff$sparse_chol %*% (params$field/sqrt(stuff$scale))))
####    # Make a half step for momentum at the end.
####    p = p - 
####      as.matrix(
####        solve(t(range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
####              - beta_prior_log_dens_derivative
####              (beta = new_range_beta, n_PP = n_PP, 
####                chol_crossprod_X = covariates$range_X$chol_crossprod_X,
####                beta0_mean = hierarchical_model$range_beta0_mean,
####                beta0_var =  hierarchical_model$range_beta0_var, 
####                log_scale = params$range_log_scale) # normal prior
####              #normal prior derivative                
####              + X_PP_crossprod(
####                X = range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
####                Y = # Jacobian of range field wrt range_beta
####                  (
####                    # natural gradient of obs likelihood wrt range field
####                    derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
####                                          left_vector = as.vector(
####                                            Matrix::solve(
####                                              Matrix::t(new_sparse_chol), 
####                                              - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
####                                                            ((new_field[vecchia_approx$locs_match] - stuff$lm_residuals) / stuff$noise))
####                                              * sqrt(stuff$scale) # part of sparse chol
####                                            )), 
####                                          right_vector = new_field/sqrt(stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
####                                          NNarray = vecchia_approx$NNarray  
####                    )
####                  )
####              )
####        ))%*%eps_mat/ 2
####    # Evaluate potential and kinetic energies at start and end of trajectory
####    current_K = sum (momenta$range_beta_ancillary ^2) / 2
####    proposed_U =
####      (
####        - beta_prior_log_dens(
####          beta = new_range_beta, n_PP = n_PP, 
####          chol_crossprod_X = covariates$range_X$chol_crossprod_X,
####          beta0_mean = hierarchical_model$range_beta0_mean,
####          beta0_var =  hierarchical_model$range_beta0_var, 
####          log_scale = params$range_log_scale) # normal prior
####        # normal prior 
####        + .5 * sum((stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/stuff$noise) # observation ll
####      )
####    proposed_K = sum(p^2) / 2
####    
####    #print(ker_var$range_beta_ancillary)
####    
####    ker_var$range_beta_ancillary[iter%%2+1] = ker_var$range_beta_ancillary[iter%%2+1]- 15/(iter_start + iter + 100)
####    if(!is.nan(current_U-proposed_U+current_K- proposed_K))
####    {
####      if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
####      {
####        ker_var$range_beta_ancillary[iter%%2+1] = ker_var$range_beta_ancillary[iter%%2+1] +20/(iter_start + iter + 100)
####        #print("tatato!")
####        momenta$range_beta_ancillary = p
####        params$field = new_field
####        stuff$sparse_chol= new_sparse_chol
####        stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
####        params$range_beta[] = new_range_beta
####      }
####    }
    
#    ###########################
#    # Range beta (sufficient) #
#    ###########################
#    if(!hierarchical_model$anisotropic) eps_mat = matrix(exp(ker_var$range_beta_sufficient[1]))
#    if(hierarchical_model$anisotropic) eps_mat = diag(exp(ker_var$range_beta_sufficient[c(1,2,2)]))
#    q = range$chol_crossprod_X_locs %*% params$range_beta # whitening wrt covariates of the range
#    current_U =
#      (
#        - beta_prior_log_dens(
#          beta = params$range_beta, n_PP = n_PP, 
#          chol_crossprod_X = covariates$range_X$chol_crossprod_X,
#          beta0_mean = hierarchical_model$range_beta0_mean,
#          beta0_var =  hierarchical_model$range_beta0_var, 
#          log_scale = params$range_log_scale) # normal prior
#        # normal prior 
#        + .5* sum((stuff$sparse_chol %*% (params$field/sqrt(stuff$scale)))^2)
#        - sum(log(stuff$compressed_sparse_chol_and_grad[[1]][,1]))
#      )
#    
#    # MALA whitened
#    momenta$range_beta_sufficient = sqrt(.9) * momenta$range_beta_sufficient + sqrt(.1)*matrix(rnorm(length(q)), nrow = nrow(q))
#    p = momenta$range_beta_sufficient
#    # Make a half step for momentum at the beginning
#    p = p - 
#      as.matrix(
#        solve(t(range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
#              - beta_prior_log_dens_derivative
#              (beta = params$range_beta, n_PP = n_PP, 
#                chol_crossprod_X = covariates$range_X$chol_crossprod_X,
#                beta0_mean = hierarchical_model$range_beta0_mean,
#                beta0_var =  hierarchical_model$range_beta0_var, 
#                log_scale = params$range_log_scale) # normal prior
#              # normal prior derivative                
#              + X_PP_crossprod(
#                X = range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
#                Y = # Jacobian of range field wrt range_beta
#                  (# natural gradient of obs likelihood wrt range field
#                    derivative_sandwiches(derivatives = stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
#                                          left_vector = as.vector(stuff$sparse_chol %*% (params$field/sqrt(stuff$scale))), # left vector = whitened latent field
#                                          right_vector = params$field/sqrt(stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
#                                          NNarray = vecchia_approx$NNarray  
#                    )
#                    - log_determinant_derivatives(sparse_chol_and_grad = stuff$compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
#                  ))))%*% eps_mat/ 2
#    #### Checking the gradient
#    ##      source("Bidart/R/Useful_stuff.R")
#    ##       # recomputing current sparse chol
#    ##       params$range_beta[] = rnorm(length(params$range_beta[]))
#    ##       stuff$compressed_sparse_chol_and_grad = 
#    ##         compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
#    ##                             range_beta = params$range_beta, NNarray = vecchia_approx$NNarray, 
#    ##                             locs = locs, 
#    ##                             range_X = range_X$X_locs, 
#    ##                             nu = hierarchical_model$nu, 
#    ##                             PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
#    ##                             compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
#    ##         )
#    ##       stuff$sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
#    ##       
#    ##       # compute gradient using derivative of sparse chol
#    ##       d1 = X_PP_crossprod(
#    ##         X = range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
#    ##         Y = # Jacobian of range field wrt range_beta
#    ##           (
#    ##             # natural gradient of obs likelihood wrt range field
#    ##             derivative_sandwiches(derivatives = stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
#    ##                                           left_vector = as.vector(stuff$sparse_chol %*% (params$field/sqrt(stuff$scale))), # left vector = whitened latent field
#    ##                                           right_vector = params$field/sqrt(stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
#    ##                                           NNarray = vecchia_approx$NNarray  
#    ##             )
#    ##             - log_determinant_derivatives(sparse_chol_and_grad = stuff$compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
#    ##             ))
#    ##       # compute gradient using finite diff
#    ##       new_range_beta = params$range_beta
#    ##       new_range_beta[10,1] = new_range_beta[10,1] + .0001
#    ##       new_compressed_sparse_chol_and_grad = 
#    ##         compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
#    ##                                     range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
#    ##                                     locs = locs, 
#    ##                                     range_X = range_X$X_locs, 
#    ##                                     nu = hierarchical_model$nu, 
#    ##                                     PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
#    ##                                     compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
#    ##         )
#    ##       new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
#    ##       new_field = sqrt(stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, stuff$sparse_chol %*% (params$field/sqrt(stuff$scale))))
#    ##       d2 = 10000 * (
#    ##         
#    ##         (
#    ##           + .5* sum((new_sparse_chol %*% (params$field/sqrt(stuff$scale)))^2)
#    ##           - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
#    ##         )
#    ##         -
#    ##           (
#    ##             + .5* sum((stuff$sparse_chol %*% (params$field/sqrt(stuff$scale)))^2)
#    ##             - sum(log(stuff$compressed_sparse_chol_and_grad[[1]][,1]))
#    ##           )
#    ##       )
#    ##       d1/d2
#    # Make a full step for the position
#    q = q + p %*% eps_mat
#    new_range_beta = params$range_beta
#    new_range_beta = solve(covariates$range_X$chol_crossprod_X_locs, q)
#    new_compressed_sparse_chol_and_grad =
#      compute_sparse_chol(
#        num_threads = num_threads, 
#        anisotropic = hierarchical_model$anisotropic, 
#        sphere = hierarchical_model$sphere, 
#        range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
#        locs = locs, 
#        range_X = range_X$X_locs, 
#        nu = hierarchical_model$nu, 
#        PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
#        compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
#      )
#    new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
#    # Make a half step for momentum at the end.
#    p = p - 
#      as.matrix(
#        solve(t(range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
#              - beta_prior_log_dens_derivative
#              (beta = new_range_beta, n_PP = n_PP, 
#                chol_crossprod_X = covariates$range_X$chol_crossprod_X,
#                beta0_mean = hierarchical_model$range_beta0_mean,
#                beta0_var =  hierarchical_model$range_beta0_var, 
#                log_scale = params$range_log_scale) # normal prior
#              # normal prior derivative                
#              + X_PP_crossprod(
#                X = range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
#                Y = # Jacobian of range field wrt range_beta
#                  (# natural gradient of obs likelihood wrt range field
#                    derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
#                                          left_vector = as.vector(new_sparse_chol %*% (params$field/sqrt(stuff$scale))), # left vector = whitened latent field
#                                          right_vector = params$field/sqrt(stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
#                                          NNarray = vecchia_approx$NNarray  
#                    )
#                    - log_determinant_derivatives(sparse_chol_and_grad = new_compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
#                  )))) %*% eps_mat / 2
#    # Evaluate potential and kinetic energies at start and end of trajectory
#    current_K = sum (momenta$range_beta_sufficient ^2) / 2
#    proposed_U =
#      (
#        - beta_prior_log_dens(
#          beta = new_range_beta, n_PP = n_PP, 
#          chol_crossprod_X = noise_X$chol_crossprod_X,
#          beta0_mean = hierarchical_model$range_beta0_mean,
#          beta0_var =  hierarchical_model$range_beta0_var, 
#          log_scale = params$range_log_scale) # normal prior
#        + .5* sum((new_sparse_chol %*% (params$field/sqrt(stuff$scale)))^2)
#        - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
#      )
#    
#    ker_var$range_beta_sufficient[iter%%2+1] = ker_var$range_beta_sufficient[iter%%2+1]- 15/(iter_start + iter + 100)
#    proposed_K = sum(p^2) / 2
#    if(!is.nan(current_U-proposed_U+current_K- proposed_K))
#    {
#      if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
#      {
#        ker_var$range_beta_sufficient[iter%%2+1] = ker_var$range_beta_sufficient[iter%%2+1]+ 20/(iter_start + iter + 100)
#        momenta$range_beta_sufficient = p
#        stuff$sparse_chol= new_sparse_chol
#        stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
#        params$range_beta[] = new_range_beta
#      }
#    }
#    
#    #############################
#    # Variance of the  range PP #
#    #############################
#    
#    if(hierarchical_model$range_PP){
#      # ancillary-sufficient ####
#      q = params$range_log_scale
#      current_U =
#        (
#          + .5* sum((stuff$sparse_chol %*% (params$field/sqrt(stuff$scale)))^2)
#          - sum(log(stuff$compressed_sparse_chol_and_grad[[1]][,1]))
#        )
#      
#      # MALA whitened
#      momenta$range_log_scale_sufficient = sqrt(.9) * momenta$range_log_scale_sufficient + sqrt(.1)*rnorm(length(q))
#      p = momenta$range_log_scale_sufficient
#      # Make a half step for momentum at the beginning
#      d_beta_d_scale = 
#        derivative_field_wrt_scale(
#          params$range_beta[-seq(range_X$n_regressors),,drop = F], 
#          params$range_log_scale
#        )
#      # derivative of potential wrt range beta
#      d_potential_d_beta = as.matrix(
#        + X_PP_crossprod(
#          X = range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
#          Y = # Jacobian of range field wrt range_beta
#            (# natural gradient of obs likelihood wrt range field
#              derivative_sandwiches(derivatives = stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
#                                    left_vector = as.vector(stuff$sparse_chol %*% (params$field/sqrt(stuff$scale))), # left vector = whitened latent field
#                                    right_vector = params$field/sqrt(stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
#                                    NNarray = vecchia_approx$NNarray  
#              )
#              - log_determinant_derivatives(sparse_chol_and_grad = stuff$compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
#            ))[-seq(range_X$n_regressors),,drop = F]
#      )
#      p = p - exp(ker_var$range_log_scale_sufficient) *
#        (
#          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
#        )/2
#      
#      # Make a full step for the position
#      q = q + exp(ker_var$range_log_scale_sufficient) * p
#      new_range_beta = params$range_beta
#      new_range_beta[-seq(range_X$n_regressors),] = new_range_beta[-seq(range_X$n_regressors),] %*% 
#        solve(chol(expmat(params$range_log_scale))) %*% chol(expmat(q))
#      
#      new_compressed_sparse_chol_and_grad =
#        compute_sparse_chol(
#          num_threads = num_threads, 
#          anisotropic = hierarchical_model$anisotropic, 
#          sphere = hierarchical_model$sphere, 
#          range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
#          locs = locs, 
#          range_X = range_X$X_locs, 
#          nu = hierarchical_model$nu, 
#          PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
#          compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
#        )
#      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
#      # Make a half step for momentum at the end.
#      d_beta_d_scale = 
#        derivative_field_wrt_scale(
#          new_range_beta[-seq(range_X$n_regressors),,drop = F], 
#          q
#        )
#      # derivative of potential wrt range beta
#      d_potential_d_beta = as.matrix(
#        + X_PP_crossprod(
#          X = range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
#          Y = # Jacobian of range field wrt range_beta
#            (# natural gradient of obs likelihood wrt range field
#              derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
#                                    left_vector = as.vector(new_sparse_chol %*% (params$field/sqrt(stuff$scale))), # left vector = whitened latent field
#                                    right_vector = params$field/sqrt(stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
#                                    NNarray = vecchia_approx$NNarray  
#              )
#              - log_determinant_derivatives(sparse_chol_and_grad = new_compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
#            ))[-seq(range_X$n_regressors),,drop = F]
#      )
#      p = p - exp(ker_var$range_log_scale_sufficient) *
#        (
#          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
#        )/2
#      # Evaluate potential and kinetic energies at start and end of trajectory
#      current_K = sum (momenta$range_log_scale_sufficient ^2) / 2
#      proposed_U =
#        (
#          + .5* sum((new_sparse_chol %*% (params$field/sqrt(stuff$scale)))^2)
#          - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
#        )
#      proposed_K = sum(p^2) / 2
#      ker_var$range_log_scale_sufficient  = ker_var$range_log_scale_sufficient - 15/(iter_start + iter + 100)
#      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
#      {
#        if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
#        {
#          ker_var$range_log_scale_sufficient  = ker_var$range_log_scale_sufficient + 20/(iter_start + iter + 100)
#          new_eigen = eigen(expmat(q))$values
#          if(
#            all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
#            all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2])))
#          {
#            momenta$range_log_scale_sufficient = p
#            
#            params$range_beta = new_range_beta
#            params$range_log_scale[] = q
#            
#            stuff$sparse_chol= new_sparse_chol
#            stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
#          }
#        }
#      }
#      
#      
#      # sufficient - sufficient ####
#      for(i in seq(10))
#      {
#        new_range_log_scale = params$range_log_scale + rnorm(length(params$range_log_scale), 0, .1)
#        old_eigen = eigen(expmat(params$range_log_scale))$values
#        new_eigen = eigen(expmat(new_range_log_scale))$values
#        if(
#          all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
#          all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2]))&
#          (
#            + beta_prior_log_dens(
#              beta = params$range_beta, n_PP = n_PP, 
#              chol_crossprod_X = noise_X$chol_crossprod_X,
#              beta0_mean = hierarchical_model$range_beta0_mean,
#              beta0_var =  hierarchical_model$range_beta0_var, 
#              log_scale = new_range_log_scale)
#            - beta_prior_log_dens(
#              beta = params$range_beta, n_PP = n_PP, 
#              chol_crossprod_X = noise_X$chol_crossprod_X,
#              beta0_mean = hierarchical_model$range_beta0_mean,
#              beta0_var =  hierarchical_model$range_beta0_var, 
#              log_scale = params$range_log_scale) 
#            > log(runif(1))
#          )
#          
#        )
#        {
#          params$range_log_scale = new_range_log_scale
#        }
#      }
#      # ancillary - ancillary ####
#      q = params$range_log_scale
#      current_U =
#        (
#          + .5 * sum((stuff$lm_residuals -  params$field[vecchia_approx$locs_match])^2/stuff$noise) # observation ll
#        )
#      # MALA whitened
#      momenta$range_log_scale_ancillary = sqrt(.9) * momenta$range_log_scale_ancillary + sqrt(.1)*rnorm(length(q))
#      p = momenta$range_log_scale_ancillary
#      # Make a half step for momentum at the beginning
#      d_beta_d_scale = 
#        derivative_field_wrt_scale(
#          params$range_beta[-seq(range_X$n_regressors),,drop = F], 
#          params$range_log_scale
#        )
#      # derivative of potential wrt range beta
#      d_potential_d_beta = X_PP_crossprod(
#        X = range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
#        Y = # Jacobian of range field wrt range_beta
#          (
#            # natural gradient of obs likelihood wrt range field
#            derivative_sandwiches(derivatives = stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
#                                  left_vector = as.vector(
#                                    Matrix::solve(
#                                      Matrix::t(stuff$sparse_chol), 
#                                      - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
#                                                    ((params$field[vecchia_approx$locs_match] - stuff$lm_residuals) / stuff$noise))
#                                      * sqrt(stuff$scale) # part of sparse chol
#                                    )), 
#                                  right_vector = params$field/sqrt(stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
#                                  NNarray = vecchia_approx$NNarray  
#            )
#          ))[-seq(range_X$n_regressors),,drop = F]
#      p = p - exp(ker_var$range_log_scale_ancillary) *
#        (
#          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
#        )/2
#      
#      q = q + exp(ker_var$range_log_scale_ancillary) * p
#      new_range_beta = params$range_beta
#      new_range_beta[-seq(range_X$n_regressors),] = new_range_beta[-seq(range_X$n_regressors),] %*% 
#        solve(chol(expmat(params$range_log_scale))) %*% chol(expmat(q))
#      
#      new_compressed_sparse_chol_and_grad =
#        compute_sparse_chol(
#          num_threads = num_threads,
#          anisotropic = hierarchical_model$anisotropic, 
#          sphere = hierarchical_model$sphere, 
#          range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
#          locs = locs, 
#          range_X = range_X$X_locs, 
#          nu = hierarchical_model$nu, 
#          PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
#          compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
#        )
#      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
#      new_field = sqrt(stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, stuff$sparse_chol %*% (params$field/sqrt(stuff$scale))))
#      # Make a half step for momentum at the end.
#      d_beta_d_scale = 
#        derivative_field_wrt_scale(
#          new_range_beta[-seq(range_X$n_regressors),,drop = F], 
#          q
#        )
#      # derivative of potential wrt range beta
#      d_potential_d_beta = X_PP_crossprod(
#        X = range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
#        Y = # Jacobian of range field wrt range_beta
#          (
#            # natural gradient of obs likelihood wrt range field
#            derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
#                                  left_vector = as.vector(
#                                    Matrix::solve(
#                                      Matrix::t(new_sparse_chol), 
#                                      - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
#                                                    ((new_field[vecchia_approx$locs_match] - stuff$lm_residuals) / stuff$noise))
#                                      * sqrt(stuff$scale) # part of sparse chol
#                                    )), 
#                                  right_vector = new_field/sqrt(stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
#                                  NNarray = vecchia_approx$NNarray  
#            )
#          )
#      )[-seq(range_X$n_regressors),,drop = F]
#      p = p - exp(ker_var$range_log_scale_ancillary) *
#        (
#          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
#        )/2
#      # Evaluate potential and kinetic energies at start and end of trajectory
#      current_K = sum (momenta$range_log_scale_ancillary ^2) / 2
#      proposed_U =
#        (
#          + .5 * sum((stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/stuff$noise) # observation ll
#        )
#      proposed_K = sum(p^2) / 2
#      ker_var$range_log_scale_ancillary = ker_var$range_log_scale_ancillary -15/(iter_start + iter + 100)
#      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
#      {
#        if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
#        {
#          new_eigen = eigen(expmat(q))$values
#          if(
#            all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
#            all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2])))
#          {
#            
#            ker_var$range_log_scale_ancillary = ker_var$range_log_scale_ancillary + 20/(iter_start + iter + 100)
#            momenta$range_log_scale_ancillary = p
#            
#            params$range_beta[] = new_range_beta
#            params$range_log_scale[] = q
#            params$field = new_field
#            
#            stuff$sparse_chol= new_sparse_chol
#            stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
#          }
#        }
#      }
#      
#      # sufficient - ancillary (equivalent to sufficient-sufficient) ####
#      for(i in seq(10))
#      {
#        new_range_log_scale = params$range_log_scale + rnorm(length(params$range_log_scale), 0, .1)
#        old_eigen = eigen(expmat(params$range_log_scale))$values
#        new_eigen = eigen(expmat(new_range_log_scale))$values
#        if(
#          all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
#          all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2]))&
#          (
#            + beta_prior_log_dens(
#              beta = params$range_beta, n_PP = n_PP, 
#              chol_crossprod_X = noise_X$chol_crossprod_X,
#              beta0_mean = hierarchical_model$range_beta0_mean,
#              beta0_var =  hierarchical_model$range_beta0_var, 
#              log_scale = new_range_log_scale)
#            - beta_prior_log_dens(
#              beta = params$range_beta, n_PP = n_PP, 
#              chol_crossprod_X = covariates$range_X$chol_crossprod_X,
#              beta0_mean = hierarchical_model$range_beta0_mean,
#              beta0_var =  hierarchical_model$range_beta0_var, 
#              log_scale = params$range_log_scale) 
#            > log(runif(1))
#          )
#          
#        )
#        {
#          params$range_log_scale = new_range_log_scale
#        }
#      }
#    }
#    
#    
#    
#    #########
#    # Noise #
#    #########
#    ##############
#    # Noise beta #
#    ##############
#    # recomputation in order to avoid errors
#    stuff$noise = variance_field(
#      beta = params$noise_beta, X = noise_X$X, 
#      PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP, 
#      locs_idx = NULL
#    )
#    # VEWY IMPOWTANT don't remove or comment
#    squared_residuals = as.vector(stuff$lm_residuals - params$field[vecchia_approx$locs_match])^2
#    # HMC update
#    q = noise_X$chol_crossprod_X %*% params$noise_beta
#    current_U =
#      (
#        - beta_prior_log_dens(
#          beta = params$noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
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
#                  beta = params$noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
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
#    ####     - beta_prior_log_dens(beta = noise_beta_, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
#    ####                                   beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
#    ####                                   beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
#    ####                                   log_scale = params$noise_log_scale) # normal prior 
#    ####     +.5* sum(log(noise_)) # det
#    ####     +.5*sum(squared_residuals/noise_) # observations
#    ####   )
#    #### print(
#    ####   (
#    ####             - beta_prior_log_dens_derivative(beta = params$noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
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
#                  beta = new_noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
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
#        - beta_prior_log_dens(beta = new_noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
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
#            + beta_prior_log_dens(beta = params$noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
#                                  chol_crossprod_X = noise_X$chol_crossprod_X,
#                                  beta0_mean = hierarchical_model$noise_beta0_mean,
#                                  beta0_var =  hierarchical_model$noise_beta0_var, 
#                                  log_scale = new_noise_log_scale) - 
#            beta_prior_log_dens(beta = params$noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
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
#    #########
#    # Scale #
#    #########
#    {
#      ##############################
#      # scale beta sufficient MALA #
#      ##############################
#      sparse_chol_diag_field = stuff$sparse_chol %*% Matrix::Diagonal(x = params$field)
#      q = scale_X$chol_crossprod_X_locs %*% params$scale_beta
#      current_U =
#        (
#          - beta_prior_log_dens(beta = params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                chol_crossprod_X = scale_X$chol_crossprod_X,
#                                beta0_mean = hierarchical_model$scale_beta0_mean,
#                                beta0_var =  hierarchical_model$scale_beta0_var, 
#                                log_scale = params$scale_log_scale) # normal prior 
#          +0.5*sum(log(stuff$scale))# determinant part
#          +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/stuff$scale))^2)# covmat product part
#        )
#      # HMC whitened
#      momenta$scale_beta_sufficient = sqrt(.9) * momenta$scale_beta_sufficient + sqrt(.1)*rnorm(length(momenta$scale_beta_sufficient))
#      p = momenta$scale_beta_sufficient
#      # Make a half step for momentum at the beginning
#      p = p - exp(ker_var$scale_beta_sufficient_mala) *
#        (
#          + solve(t(scale_X$chol_crossprod_X_locs), # solving by t chol because of whitening
#                  - beta_prior_log_dens_derivative(beta = params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                                   chol_crossprod_X = scale_X$chol_crossprod_X,
#                                                   beta0_mean = hierarchical_model$scale_beta0_mean,
#                                                   beta0_var =  hierarchical_model$scale_beta0_var, 
#                                                   log_scale = params$scale_log_scale) # normal prior
#                  +  X_PP_crossprod(X = scale_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, locs_idx = vecchia_approx$hctam_scol_1,
#                                    Y = (
#                                      .5  # determinant part 
#                                      -.5 * sqrt(1/stuff$scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/stuff$scale)))# natural derivative
#                                    ))
#          ))/ 2
#      
#      # Make a full step for the position
#      q = q + exp(ker_var$scale_beta_sufficient_mala) * p
#      new_scale_beta = solve(scale_X$chol_crossprod_X_locs, q)
#      new_scale =variance_field(beta = new_scale_beta, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, 
#                                X = scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
#      
#      # Make a half step for momentum at the end.
#      p = p - exp(ker_var$scale_beta_sufficient_mala) *
#        (
#          + solve(t(scale_X$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
#                  - beta_prior_log_dens_derivative(beta = new_scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                                   chol_crossprod_X = scale_X$chol_crossprod_X,
#                                                   beta0_mean = hierarchical_model$scale_beta0_mean,
#                                                   beta0_var =  hierarchical_model$scale_beta0_var, 
#                                                   log_scale = params$scale_log_scale) # normal prior  
#                  + X_PP_crossprod(X = scale_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, locs_idx = vecchia_approx$hctam_scol_1,
#                                   (
#                                     .5  # determinant part 
#                                     -.5 * sqrt(1/new_scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/new_scale)))# natural derivative
#                                   ))
#          ))/ 2
#      # Evaluate potential and kinetic energies at start and end of trajectory
#      current_K = sum (momenta$scale_beta_sufficient ^2) / 2
#      proposed_U = 
#        (
#          - beta_prior_log_dens(beta = new_scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                chol_crossprod_X = scale_X$chol_crossprod_X,
#                                beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
#                                beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
#                                log_scale = params$scale_log_scale) # normal prior 
#          +0.5*sum(log(new_scale))# determinant part
#          +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/new_scale))^2)# covmat product part
#        )
#      proposed_K = sum(p^2) / 2
#      ker_var$scale_beta_sufficient_mala = ker_var$scale_beta_sufficient_mala - 1/sqrt(iter_start + iter +100)
#      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
#      {
#        if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
#        {
#          ker_var$scale_beta_sufficient_mala = ker_var$scale_beta_sufficient_mala + 2/sqrt(iter_start + iter +100)
#          momenta$scale_beta_sufficient = p
#          params$scale_beta[] = new_scale_beta
#          stuff$scale = new_scale
#        }
#      }
#      
#      #############################
#      # scale beta ancillary MALA #
#      #############################
#      q = scale_X$chol_crossprod_X_locs %*% params$scale_beta
#      current_U =
#        (
#          - beta_prior_log_dens(beta = params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                chol_crossprod_X = scale_X$chol_crossprod_X,
#                                beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
#                                beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
#                                log_scale = params$scale_log_scale) # normal prior 
#          +.5 * sum((observed_field - stuff$lm_fit - params$field[vecchia_approx$locs_match] )^2/stuff$noise)
#        )
#      # MALA whitened
#      momenta$scale_beta_ancillary = sqrt(.9) * momenta$scale_beta_ancillary   + sqrt(.1)*rnorm(length(momenta$scale_beta_ancillary))
#      p = momenta$scale_beta_ancillary
#      # Make a. half step for momentum at the beginning
#      p = p - exp(ker_var$scale_beta_ancillary_mala) *
#        (
#          + solve(t(scale_X$chol_crossprod_X_locs), # solving by prior chol because of whitening
#                  - beta_prior_log_dens_derivative(beta = params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                                   chol_crossprod_X = scale_X$chol_crossprod_X,
#                                                   beta0_mean = hierarchical_model$scale_beta0_mean,
#                                                   beta0_var =  hierarchical_model$scale_beta0_var, 
#                                                   log_scale = params$scale_log_scale) # normal prior 
#                  + X_PP_crossprod(X = scale_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, locs_idx = vecchia_approx$hctam_scol_1,
#                                   (.5 * params$field * 
#                                      as.vector(vecchia_approx$locs_match_matrix %*% 
#                                                  ((params$field[vecchia_approx$locs_match] - stuff$lm_residuals)/
#                                                     stuff$noise)
#                                      )))
#          ))/ 2
#      # Make a full step for the position
#      q = q + exp(ker_var$scale_beta_ancillary_mala) * p
#      new_scale_beta = solve(scale_X$chol_crossprod_X_locs, q)
#      new_scale =variance_field(beta = new_scale_beta, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, 
#                                X = scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
#      new_field  = params$field*sqrt(new_scale)/sqrt(stuff$scale)
#      # Make a half step for momentum at the end.
#      p = p - exp(ker_var$scale_beta_ancillary_mala) *
#        (
#          + solve(t(scale_X$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
#                  - beta_prior_log_dens_derivative(beta = new_scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                                   chol_crossprod_X = scale_X$chol_crossprod_X,
#                                                   beta0_mean = hierarchical_model$scale_beta0_mean,
#                                                   beta0_var =  hierarchical_model$scale_beta0_var, 
#                                                   log_scale = params$scale_log_scale) # normal prior  
#                  + X_PP_crossprod(X = scale_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, locs_idx = vecchia_approx$hctam_scol_1,
#                                   (.5 * new_field * 
#                                      as.vector(vecchia_approx$locs_match_matrix %*% 
#                                                  ((new_field[vecchia_approx$locs_match] - stuff$lm_residuals)/
#                                                     stuff$noise)
#                                      )))
#          ))/ 2
#      # Evaluate potential and kinetic energies at start and end of trajectory
#      current_K = sum (momenta$scale_beta_ancillary ^2) / 2
#      proposed_U = 
#        (
#          - beta_prior_log_dens(beta = new_scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                chol_crossprod_X = scale_X$chol_crossprod_X,
#                                beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
#                                beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
#                                log_scale = params$scale_log_scale) # normal prior  
#          +.5 * sum((observed_field - stuff$lm_fit - new_field[vecchia_approx$locs_match] )^2/stuff$noise)
#        )
#      proposed_K = sum(p^2) / 2
#      ker_var$scale_beta_ancillary_mala = ker_var$scale_beta_ancillary_mala - 1/sqrt(iter_start + iter +100)
#      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
#      {
#        if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
#        {
#          ker_var$scale_beta_ancillary_mala = ker_var$scale_beta_ancillary_mala + 2/sqrt(iter_start + iter +100)
#          momenta$scale_beta_ancillary = p
#          params$scale_beta[] = new_scale_beta
#          params$field = new_field
#          stuff$scale = new_scale
#        }
#      }
#      
#      ###################
#      # scale log scale #
#      ###################
#      if(hierarchical_model$scale_PP){
#        # ancillary -- sufficient  ####
#        # a change in hyperprior scale changes (rescales) the scale, which is then compared with the latent field
#        new_scale_log_scale = params$scale_log_scale + rnorm(1, 0, exp(ker_var$scale_log_scale_sufficient)) 
#        new_scale_beta = params$scale_beta
#        new_scale_beta[-seq(scale_X$n_regressors)] = new_scale_beta[-seq(scale_X$n_regressors)] *
#          exp((new_scale_log_scale - params$scale_log_scale)/2)
#        new_scale = variance_field(beta = new_scale_beta, PP = hierarchical_model$PP, 
#                                   use_PP = hierarchical_model$scale_PP, X = scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
#        ker_var$scale_log_scale_sufficient = ker_var$scale_log_scale_sufficient - .25/sqrt(iter_start + iter +100)
#        if(
#          (
#            (
#              +.5* sum(log(stuff$scale)) 
#              -.5* sum(log(new_scale)) # log determinant
#              +.5*sum((stuff$sparse_chol %*% (params$field/sqrt(stuff$scale)))^2)
#              -.5*sum((stuff$sparse_chol %*% (params$field/sqrt(new_scale)))^2) # Gaussian density of the latent field
#            )
#            > log(runif(1))
#          )
#        )
#        {
#          if(
#            (new_scale_log_scale > hierarchical_model$scale_log_scale_prior[1])&
#            (new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
#          )
#          {
#            ker_var$scale_log_scale_sufficient = ker_var$scale_log_scale_sufficient + 1/sqrt(iter_start + iter +100)
#            params$scale_log_scale = new_scale_log_scale
#            params$scale_beta = new_scale_beta 
#            stuff$scale = new_scale
#          }
#        }
#        
#        
#        # sufficient -- sufficient ####
#        for(i in seq(4))
#        {
#          new_scale_log_scale = params$scale_log_scale + rnorm(1, 0, .1)
#          if(
#            (
#              + beta_prior_log_dens(beta = params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                    chol_crossprod_X = scale_X$chol_crossprod_X,
#                                    beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
#                                    beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
#                                    log_scale = new_scale_log_scale)     
#              - beta_prior_log_dens(beta = params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                    chol_crossprod_X = scale_X$chol_crossprod_X,
#                                    beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
#                                    beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
#                                    log_scale = params$scale_log_scale)     
#              > log(runif(1))
#            )
#            &(new_scale_log_scale > hierarchical_model$scale_log_scale_prior[1])
#            &(new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
#          )
#          {
#            params$scale_log_scale = new_scale_log_scale
#          }
#        }
#        # ancillary -- ancillary ####
#        new_scale_log_scale = params$scale_log_scale + rnorm(1, 0, exp(ker_var$scale_log_scale_ancillary)) 
#        new_scale_beta = params$scale_beta
#        new_scale_beta[-seq(scale_X$n_regressors)] = new_scale_beta[-seq(scale_X$n_regressors)] *
#          exp((new_scale_log_scale - params$scale_log_scale)/2)
#        new_scale =variance_field(new_scale_beta, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, 
#                                  X = scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
#        new_field = params$field * sqrt(new_scale)/sqrt(stuff$scale)
#        ker_var$scale_log_scale_ancillary = ker_var$scale_log_scale_ancillary - .25/sqrt(iter_start + iter +100)
#        if(
#          (
#            -.5* sum((stuff$lm_residuals -          new_field[vecchia_approx$locs_match])^2/stuff$noise) 
#            +.5* sum((stuff$lm_residuals - params$field[vecchia_approx$locs_match])^2/stuff$noise) 
#            > log(runif(1))
#          )
#        )
#        {
#          if(
#            (new_scale_log_scale >  hierarchical_model$scale_log_scale_prior[1])
#            &(new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
#          )
#          {  
#            ker_var$scale_log_scale_ancillary = ker_var$scale_log_scale_ancillary + 1/sqrt(iter_start + iter +100)
#            params$field = new_field
#            params$scale_log_scale = new_scale_log_scale
#            params$scale_beta = new_scale_beta 
#            stuff$scale = new_scale
#          }
#        }
#        # sufficient -- ancillary ####
#        for(i in seq(4))
#        {
#          new_scale_log_scale = params$scale_log_scale + rnorm(1, 0, .1)
#          if(
#            (
#              + beta_prior_log_dens(beta = params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                    chol_crossprod_X = scale_X$chol_crossprod_X,
#                                    beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
#                                    beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
#                                    log_scale = new_scale_log_scale)     
#              - beta_prior_log_dens(beta = params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
#                                    chol_crossprod_X = scale_X$chol_crossprod_X,
#                                    beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
#                                    beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
#                                    log_scale = params$scale_log_scale)     
#              > log(runif(1))
#            )
#            &(new_scale_log_scale > hierarchical_model$scale_log_scale_prior[1])
#            &(new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
#          )
#          {
#            params$scale_log_scale = new_scale_log_scale
#          }
#        }
#        
#      }
#      
#    }
#    
#    
    #######################
    # Storing the samples #
    #######################
    params_records[[iter]] = params
  }
  print(Sys.time()-t1)
  return(list("state" = state, "params_records" = params_records))
#}

}

