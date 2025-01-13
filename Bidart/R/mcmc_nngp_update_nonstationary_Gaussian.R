# With only one HMC step in range 

# numerically stable log(exp(x) + exp(y)) with 0 << x << y or 0 << y << x
# let's do: log(exp(x) + exp(y)) with 0 << x << y 
# log(exp(x) + exp(y))  = log(exp(y)(1 + exp(x) / exp(y)))
# = log(exp(y)) + log((1 + exp(x) / exp(y)))
# = y + log((1 + exp(x - y)))
log_of_exp_sum = function(x,y)return(max(x,y) + log((1 + exp(-abs(x - y)))))


#list2env(mcmc_nngp_list, envir = environment())
#state = mcmc_nngp_list$states$chain_1; n_iterations_update  =100; num_threads = 10; thinning = .1; iter_start = 0; seed = 1; iter=1

#' @export
mcmc_nngp_update_Gaussian = function(data,
                                     hierarchical_model, vecchia_approx, # model architecture
                                     state, # model state
                                     n_iterations_update = 400, num_threads = 1, thinning = .1, iter_start = 0, seed = 1,# practical settings
                                     lib.loc = NULL
)
{
  #################
  # Sanity checks #
  #################
  # iterations and thinning
  if((floor(n_iterations_update)!=n_iterations_update) |  n_iterations_update<1)stop("n_iterations_update must be a positive round number")
  if((thinning<0)|(thinning>1))stop("thinning is a proportion and must be between 0 and 1")
  if((floor(iter_start)!=iter_start)|  iter_start<0)stop("iter_start must be a positive round number")
  # set seed 
  set.seed(seed)
  #########################################
  # Initializing chain storage structures #
  #########################################
  # this part re-creates a small portion of the $records objects of each chain. It fills it with chain state during the run, and then updates each chain with the new values
  params_records = lapply(state$params, function(x)array(0, c(dim(as.matrix(x)), round(thinning * n_iterations_update))))
  #par(mfrow = c(2, 1))
  
  library(Bidart, lib.loc = lib.loc)
  library(expm, lib.loc = lib.loc)
  library(Matrix, lib.loc = lib.loc)
  # do NOT remove or the code will bug in parallel. This is magic
  if(!is.null(state$params$range_log_scale))Bidart::expmat(state$params$range_log_scale)
  Matrix::t(state$sparse_chol_and_stuff$sparse_chol) %*% matrix(rep(1, vecchia_approx$n_locs))
  #################
  # Gibbs sampler #
  #################
  for(iter in seq(1, n_iterations_update))
  {
    gc()
    if(iter/10 ==iter %/% 10)print(paste("iteration", iter))
    
    ###########################
    # Regression coefficients #
    ###########################
    t1 = Sys.time()
    # sufficient parametrization of latent field
    beta_covmat = solve(as.matrix(t(data$covariates$X$X) %*% Matrix::Diagonal(x = 1/state$sparse_chol_and_stuff$noise) %*% data$covariates$X$X))
    if(all(!is.infinite(beta_covmat) & !is.nan(beta_covmat)))
    {
      if(all(eigen(beta_covmat)$d >0))
      {
        beta_mean = c((((data$observed_field-state$params$field[vecchia_approx$locs_match]) / state$sparse_chol_and_stuff$noise) %*% data$covariates$X$X) %*% beta_covmat)
        state$params$beta[]   = c(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
      }}
    # interweaving centered sampling in case of location-wise data to improve beta sampling
    centered_field = as.vector(state$params$field + data$covariates$X$X_locs%*%matrix(state$params$beta[data$covariates$X$which_locs], ncol = 1))
    sparse_chol_X = as.matrix(state$sparse_chol_and_stuff$sparse_chol %*% Matrix::Diagonal(x = 1/sqrt(state$sparse_chol_and_stuff$scale)) %*% data$covariates$X$X_locs)
    beta_precision = crossprod(x = sparse_chol_X, y = sparse_chol_X)
    beta_covmat = solve(beta_precision, tol = min(rcond(beta_precision),.Machine$double.eps))
    if(all(!is.infinite(beta_covmat) & !is.nan(beta_covmat)))
    {
      if(all(eigen(beta_covmat)$d >0))
      {
        beta_mean =  c(as.vector(state$sparse_chol_and_stuff$sparse_chol %*% (centered_field/sqrt(state$sparse_chol_and_stuff$scale)))  %*% sparse_chol_X %*% beta_covmat)
        state$params$beta[data$covariates$X$which_locs]   = as.vector(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
        state$params$field = centered_field - as.vector(data$covariates$X$X_locs %*% matrix(state$params$beta[data$covariates$X$which_locs], ncol = 1))
      }}
    # updating stuff 
    state$sparse_chol_and_stuff$lm_fit       = as.vector(data$covariates$X$X%*%     matrix(state$params$beta, ncol = 1))
    state$sparse_chol_and_stuff$lm_fit_locs  = as.vector(data$covariates$X$X_locs%*%matrix(state$params$beta[data$covariates$X$which_locs], ncol = 1))
    state$sparse_chol_and_stuff$lm_residuals = as.vector(data$observed_field-              state$sparse_chol_and_stuff$lm_fit)
    #print(Sys.time()-t1)
    ################
    # Latent field #
    ################
    locs_partition = vecchia_approx$locs_partition[,runif(1, 1, ncol(vecchia_approx$locs_partition))]
    additional_precision = vecchia_approx$locs_match_matrix %*% (1/state$sparse_chol_and_stuff$noise)
    additional_mean = vecchia_approx$locs_match_matrix %*% ((data$observed_field - state$sparse_chol_and_stuff$lm_fit)/state$sparse_chol_and_stuff$noise)
    chol_sub_list = lapply(
      seq(length(unique(locs_partition))), function(idx)
      {
        selected_locs = which(locs_partition==idx)
        state$sparse_chol_and_stuff$sparse_chol[,selected_locs] %*% Matrix::Diagonal(x = 1/sqrt(state$sparse_chol_and_stuff$scale[selected_locs]))
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
                                      (as.vector(Matrix::crossprod(chol_sub_list[[idx]], state$sparse_chol_and_stuff$sparse_chol%*%((state$params$field/sqrt(state$sparse_chol_and_stuff$scale))*(locs_partition!=idx)))) # rest of the latent field
                                       - additional_mean[selected_locs]))) # Gaussian observations 
      state$params$field [selected_locs] = as.vector(as.vector(cond_mean) +  Matrix::t(chol_Q_AA_sub_list[[idx]]$P) %*% Matrix::solve(Matrix::t(chol_Q_AA_sub_list[[idx]]$L), rnorm(length(selected_locs))) )
    }
    
    
    
    ####################
    # Range scale swap #
    ####################
    
    innovation = exp(state$transition_kernels$range_scale_blocked_KHR) *  rnorm(1)
    new_range_beta = state$params$range_beta
    new_range_beta[1,1] = new_range_beta[1,1] + innovation / hierarchical_model$nu
    new_compressed_sparse_chol_and_grad =
      Bidart::compute_sparse_chol(
        num_threads = num_threads, 
        anisotropic = hierarchical_model$anisotropic, 
        sphere = hierarchical_model$sphere, 
        range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
        locs = data$locs, 
        range_X = data$covariates$range_X$X_locs, 
        nu = hierarchical_model$nu, 
        PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
        compute_derivative = nrow(state$params$range_beta)>1
        , locs_idx = vecchia_approx$hctam_scol_1
      )
    new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
    
    new_scale_beta = state$params$scale_beta
    new_scale_beta[1,1] = new_scale_beta[1,1] + innovation
    new_scale =Bidart::variance_field(beta = new_scale_beta, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, 
                                      X = data$covariates$scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
    
    
    current_U =
      (
        - Bidart::beta_prior_log_dens(beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                      beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                      beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                      log_scale = state$params$range_log_scale) # normal prior
        - Bidart::beta_prior_log_dens(beta = state$params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                      beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                      beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                      log_scale = state$params$scale_log_scale) # normal prior 
        + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2) # covmat product part
        - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
        +0.5*sum(log(state$sparse_chol_and_stuff$scale))# determinant part
      )
    proposed_U =
      (
        - Bidart::beta_prior_log_dens(beta = new_range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                      beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                      beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                      log_scale = state$params$range_log_scale) # normal prior
        - Bidart::beta_prior_log_dens(beta = new_scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                      beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                      beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                      log_scale = state$params$scale_log_scale) # normal prior 
        + .5* sum((new_sparse_chol %*% (state$params$field/sqrt(new_scale)))^2)
        - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
        +0.5*sum(log(new_scale))# determinant part
      )
    state$transition_kernels$range_scale_blocked_KHR = state$transition_kernels$range_scale_blocked_KHR -1/sqrt(iter_start + iter +100)
    if(!is.nan(current_U-proposed_U))
    {
      if (log(runif(1)) < current_U-proposed_U)
      {
        #print("range scale swap !")
        state$transition_kernels$range_scale_blocked_KHR = state$transition_kernels$range_scale_blocked_KHR + 4/sqrt(iter_start + iter +100)
        
        state$params$scale_beta[] = new_scale_beta
        state$sparse_chol_and_stuff$scale = new_scale
        
        state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
        state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
        state$params$range_beta[] = new_range_beta
        state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
      }
    }
    
    ###############
    # Range beta  #
    ###############
    if(nrow(state$params$range_beta)==1){
      n_regimes = 1+(ncol(state$params$range_beta)>1)
      updated_comps = list(c(1), c(2,3))
      for (regime in seq(n_regimes)){
        ##########################
        # Range beta (ancillary) #
        ##########################
        q = data$covariates$range$chol_crossprod_X_locs %*% state$params$range_beta[,updated_comps[[regime]]] # whitening wrt covariates of the range
        current_U =
          (
            - Bidart::beta_prior_log_dens(beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = state$params$range_log_scale) # normal prior
            + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
          )
        # Make a full step for the position
        q = q + exp(state$transition_kernels$range_beta_ancillary[regime]) * rnorm(length(updated_comps[[regime]]))
        new_range_beta = state$params$range_beta
        new_range_beta[,updated_comps[[regime]]] = solve(data$covariates$range_X$chol_crossprod_X_locs, q )
        new_compressed_sparse_chol_and_grad = 
          Bidart::compute_sparse_chol(
            num_threads = num_threads, 
            anisotropic = hierarchical_model$anisotropic, 
            sphere = hierarchical_model$sphere, 
            range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
            locs = data$locs, 
            range_X = data$covariates$range_X$X_locs, 
            nu = hierarchical_model$nu, 
            PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
            compute_derivative = F, locs_idx = vecchia_approx$hctam_scol_1
          )
        new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
        proposed_U =
          (
            - Bidart::beta_prior_log_dens(beta = new_range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = state$params$range_log_scale) # normal prior
            # normal prior 
            + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
          )
        state$transition_kernels$range_beta_ancillary[regime] = state$transition_kernels$range_beta_ancillary[regime]- 1/sqrt(iter_start + iter +100)
        if(!is.nan(current_U-proposed_U))
        {
          if (log(runif(1)) < current_U-proposed_U)
          {
            state$transition_kernels$range_beta_ancillary[regime] = state$transition_kernels$range_beta_ancillary[regime] + 4/sqrt(iter_start + iter +100)
            # print("ancillary!")
            # print(regime)
            state$params$field = new_field
            state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
            state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
            state$params$range_beta[] = new_range_beta
            state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
          }
        }
        
        ###########################
        # Range beta (sufficient) #
        ###########################
        q = data$covariates$range$chol_crossprod_X_locs %*% state$params$range_beta[,updated_comps[[regime]]] # whitening wrt covariates of the range
        current_U =
          (
            - Bidart::beta_prior_log_dens(beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = state$params$range_log_scale) # normal prior
            # normal prior 
            + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
            - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
          )
        # Make a full step for the position
        q = q + exp(state$transition_kernels$range_beta_sufficient[regime])  * rnorm(length(updated_comps[[regime]]))
        new_range_beta = state$params$range_beta
        new_range_beta[,updated_comps[[regime]]] = solve(data$covariates$range_X$chol_crossprod_X_locs, q)
        new_compressed_sparse_chol_and_grad =
          Bidart::compute_sparse_chol(
            num_threads = num_threads, 
            anisotropic = hierarchical_model$anisotropic, 
            sphere = hierarchical_model$sphere, 
            range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
            locs = data$locs, 
            range_X = data$covariates$range_X$X_locs, 
            nu = hierarchical_model$nu, 
            PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
            compute_derivative = F, locs_idx = vecchia_approx$hctam_scol_1
          )
        new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        proposed_U =
          (
            - Bidart::beta_prior_log_dens(beta = new_range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = state$params$range_log_scale) # normal prior
            + .5* sum((new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
            - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
          )
        
        state$transition_kernels$range_beta_sufficient[regime] = state$transition_kernels$range_beta_sufficient[regime] -1/sqrt(iter_start + iter +100)
        if(!is.nan(current_U-proposed_U))
        {
          if (log(runif(1)) < current_U-proposed_U)
          {
            # print("sufficient!")
            # print(regime)
            state$transition_kernels$range_beta_sufficient[regime] = state$transition_kernels$range_beta_sufficient[regime] + 4/sqrt(iter_start + iter +100)
            state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
            state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
            state$params$range_beta[] = new_range_beta
            state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
          }
        }
      }
    }
    if(nrow(state$params$range_beta)>1){
      n_regimes = 1+(ncol(state$params$range_beta)>1)
      updated_comps = list(c(1), c(2,3))
      for (regime in seq(n_regimes)){
        ##########################
        # Range beta (ancillary) #
        ##########################
        q = data$covariates$range$chol_crossprod_X_locs %*% state$params$range_beta[,updated_comps[[regime]]] # whitening wrt covariates of the range
        current_U =
          (
            - Bidart::beta_prior_log_dens(beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = state$params$range_log_scale) # normal prior
            + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
          )
        # MALA whitened
        state$momenta$range_beta_ancillary[,updated_comps[[regime]]] = sqrt(.9) * state$momenta$range_beta_ancillary[,updated_comps[[regime]]] + sqrt(.1)*matrix(rnorm(length(q)), nrow = nrow(q))
        p = state$momenta$range_beta_ancillary[,updated_comps[[regime]]]
        # Make a half step for momentum at the beginning
        p = p - exp(state$transition_kernels$range_beta_ancillary[regime]) *
          as.matrix(
            solve(t(data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                  - Bidart::beta_prior_log_dens_derivative(beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                                           beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                                           beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                                           log_scale = state$params$range_log_scale) # normal prior
                  # normal prior derivative                
                  + Bidart::X_PP_crossprod(
                    X = data$covariates$range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
                    Y = # Jacobian of range field wrt range_beta
                      (
                        # natural gradient of obs likelihood wrt range field
                        Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                      left_vector = as.vector(
                                                        Matrix::solve(
                                                          Matrix::t(state$sparse_chol_and_stuff$sparse_chol), 
                                                          - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                                        ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                                          * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                                        )), 
                                                      right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                      NNarray = vecchia_approx$NNarray  
                        )
                      ))
            ))[,updated_comps[[regime]]]/ 2
        
        
        #######testing the gradient
        ##source("Bidart/R/Useful_stuff.R")
        ##d1 =         - Bidart::beta_prior_log_dens_derivative(
        ##  beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
        ##  beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
        ##  beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
        ##  log_scale = state$params$range_log_scale) # normal prior
        ##beta_ = state$params$range_beta 
        ##beta_[2, 2] = beta_[2, 2] + 1/10000
        ##d2 = 10000*(
        ##  -
        ##  Bidart::beta_prior_log_dens(beta = beta_, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
        ##                        beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
        ##                        beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
        ##                        log_scale = state$params$range_log_scale) # normal prior
        ##  + 
        ##  Bidart::beta_prior_log_dens(beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
        ##                        beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
        ##                        beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
        ##                        log_scale = state$params$range_log_scale) # normal prior
        ##)
        ##d1/d2
        #### i = 4
        #### j = 2
        #### 
        #### range_beta_ = state$params$range_beta
        #### range_beta_[i,j] = state$params$range_beta[i,j] + .0001
        #### 
        #### compressed_sparse_chol_and_grad_ = 
        ####   Bidart::compute_sparse_chol(
        ####     num_threads = num_threads, 
        ####     anisotropic = hierarchical_model$anisotropic, 
        ####     sphere = hierarchical_model$sphere, 
        ####     range_beta = range_beta_, NNarray = vecchia_approx$NNarray, 
        ####     locs = data$locs, 
        ####     range_X = data$covariates$range_X$X_locs, 
        ####     nu = hierarchical_model$nu, 
        ####     PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
        ####     compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
        ####   )
        #### sparse_chol_ = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = compressed_sparse_chol_and_grad_[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        #### field_ = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(sparse_chol_, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
        #### 
        #### 10000*(
        ####   + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  field_[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) - # observation ll
        ####     + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
        #### )
        #### 
        #### 
        #### + Bidart::X_PP_crossprod(
        ####   X = data$covariates$range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
        ####   Y = # Jacobian of range field wrt range_beta
        ####     (
        ####       # natural gradient of obs likelihood wrt range field
        ####       Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
        ####                                     left_vector = as.vector(
        ####                                       Matrix::solve(
        ####                                         Matrix::t(state$sparse_chol_and_stuff$sparse_chol), 
        ####                                         - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
        ####                                                       ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
        ####                                         * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
        ####                                       )), 
        ####                                     right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
        ####                                     NNarray = vecchia_approx$NNarray  
        ####       )
        ####     ))[i,j]
        
        
        # Make a full step for the position
        q = q + exp(state$transition_kernels$range_beta_ancillary[regime]) * p
        new_range_beta = state$params$range_beta
        new_range_beta[,updated_comps[[regime]]] = solve(data$covariates$range_X$chol_crossprod_X_locs, q )
        new_compressed_sparse_chol_and_grad = 
          Bidart::compute_sparse_chol(
            num_threads = num_threads, 
            anisotropic = hierarchical_model$anisotropic, 
            sphere = hierarchical_model$sphere, 
            range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
            locs = data$locs, 
            range_X = data$covariates$range_X$X_locs, 
            nu = hierarchical_model$nu, 
            PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
            compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
          )
        new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
        # Make a half step for momentum at the end.
        p = p - exp(state$transition_kernels$range_beta_ancillary[regime]) *
          as.matrix(
            solve(t(data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                  - Bidart::beta_prior_log_dens_derivative
                  (beta = new_range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                    beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                    beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                    log_scale = state$params$range_log_scale) # normal prior
                  #normal prior derivative                
                  + Bidart::X_PP_crossprod(
                    X = data$covariates$range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
                    Y = # Jacobian of range field wrt range_beta
                      (
                        # natural gradient of obs likelihood wrt range field
                        Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                      left_vector = as.vector(
                                                        Matrix::solve(
                                                          Matrix::t(new_sparse_chol), 
                                                          - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                                        ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                                          * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                                        )), 
                                                      right_vector = new_field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                      NNarray = vecchia_approx$NNarray  
                        )
                      )
                  )
            ))[,updated_comps[[regime]]]/ 2
        # Evaluate potential and kinetic energies at start and end of trajectory
        current_K = sum (state$momenta$range_beta_ancillary[,updated_comps[[regime]]] ^2) / 2
        proposed_U =
          (
            - Bidart::beta_prior_log_dens(beta = new_range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = state$params$range_log_scale) # normal prior
            # normal prior 
            + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
          )
        proposed_K = sum(p^2) / 2
        state$transition_kernels$range_beta_ancillary[regime] = state$transition_kernels$range_beta_ancillary[regime]- 1/sqrt(iter_start + iter +100)
        if(!is.nan(current_U-proposed_U+current_K- proposed_K))
        {
          if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
          {
            state$transition_kernels$range_beta_ancillary[regime] = state$transition_kernels$range_beta_ancillary[regime] + 2/sqrt(iter_start + iter +100)
            #print("tatato!")
            state$momenta$range_beta_ancillary[,updated_comps[[regime]]] = p
            state$params$field = new_field
            state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
            state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
            state$params$range_beta[] = new_range_beta
            state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
          }
        }
        
        ###########################
        # Range beta (sufficient) #
        ###########################
        q = data$covariates$range$chol_crossprod_X_locs %*% state$params$range_beta[,updated_comps[[regime]]] # whitening wrt covariates of the range
        current_U =
          (
            - Bidart::beta_prior_log_dens(beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = state$params$range_log_scale) # normal prior
            # normal prior 
            + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
            - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
          )
        
        # MALA whitened
        state$momenta$range_beta_sufficient[,updated_comps[[regime]]] = sqrt(.9) * state$momenta$range_beta_sufficient[,updated_comps[[regime]]] + sqrt(.1)*matrix(rnorm(length(q)), nrow = nrow(q))
        p = state$momenta$range_beta_sufficient[,updated_comps[[regime]]]
        # Make a half step for momentum at the beginning
        p = p - exp(state$transition_kernels$range_beta_sufficient[regime]) *
          as.matrix(
            solve(t(data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                  - Bidart::beta_prior_log_dens_derivative
                  (beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                    beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                    beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                    log_scale = state$params$range_log_scale) # normal prior
                  # normal prior derivative                
                  + Bidart::X_PP_crossprod(
                    X = data$covariates$range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
                    Y = # Jacobian of range field wrt range_beta
                      (# natural gradient of obs likelihood wrt range field
                        Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                      left_vector = as.vector(state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                                      right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                      NNarray = vecchia_approx$NNarray  
                        )
                        - Bidart::log_determinant_derivatives(sparse_chol_and_grad = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
                      ))))[,updated_comps[[regime]]]/ 2
        #### Checking the gradient
        ##      source("Bidart/R/Useful_stuff.R")
        ##       # recomputing current sparse chol
        ##       state$params$range_beta[] = rnorm(length(state$params$range_beta[]))
        ##       state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = 
        ##         compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
        ##                             range_beta = state$params$range_beta, NNarray = vecchia_approx$NNarray, 
        ##                             locs = data$locs, 
        ##                             range_X = data$covariates$range_X$X_locs, 
        ##                             nu = hierarchical_model$nu, 
        ##                             PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
        ##                             compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
        ##         )
        ##       state$sparse_chol_and_stuff$sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        ##       
        ##       # compute gradient using derivative of sparse chol
        ##       d1 = X_PP_crossprod(
        ##         X = data$covariates$range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
        ##         Y = # Jacobian of range field wrt range_beta
        ##           (
        ##             # natural gradient of obs likelihood wrt range field
        ##             Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
        ##                                           left_vector = as.vector(state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
        ##                                           right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
        ##                                           NNarray = vecchia_approx$NNarray  
        ##             )
        ##             - log_determinant_derivatives(sparse_chol_and_grad = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
        ##             ))
        ##       # compute gradient using finite diff
        ##       new_range_beta = state$params$range_beta
        ##       new_range_beta[10,1] = new_range_beta[10,1] + .0001
        ##       new_compressed_sparse_chol_and_grad = 
        ##         compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
        ##                                     range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
        ##                                     locs = data$locs, 
        ##                                     range_X = data$covariates$range_X$X_locs, 
        ##                                     nu = hierarchical_model$nu, 
        ##                                     PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
        ##                                     compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
        ##         )
        ##       new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        ##       new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
        ##       d2 = 10000 * (
        ##         
        ##         (
        ##           + .5* sum((new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
        ##           - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
        ##         )
        ##         -
        ##           (
        ##             + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
        ##             - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
        ##           )
        ##       )
        ##       d1/d2
        # Make a full step for the position
        q = q + exp(state$transition_kernels$range_beta_sufficient[regime]) * p
        new_range_beta = state$params$range_beta
        new_range_beta[,updated_comps[[regime]]] = solve(data$covariates$range_X$chol_crossprod_X_locs, q)
        new_compressed_sparse_chol_and_grad =
          Bidart::compute_sparse_chol(
            num_threads = num_threads, 
            anisotropic = hierarchical_model$anisotropic, 
            sphere = hierarchical_model$sphere, 
            range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
            locs = data$locs, 
            range_X = data$covariates$range_X$X_locs, 
            nu = hierarchical_model$nu, 
            PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
            compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
          )
        new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        # Make a half step for momentum at the end.
        p = p - exp(state$transition_kernels$range_beta_sufficient[regime]) *
          as.matrix(
            solve(t(data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                  - Bidart::beta_prior_log_dens_derivative
                  (beta = new_range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                    beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                    beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                    log_scale = state$params$range_log_scale) # normal prior
                  # normal prior derivative                
                  + Bidart::X_PP_crossprod(
                    X = data$covariates$range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
                    Y = # Jacobian of range field wrt range_beta
                      (# natural gradient of obs likelihood wrt range field
                        Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                      left_vector = as.vector(new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                                      right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                      NNarray = vecchia_approx$NNarray  
                        )
                        - Bidart::log_determinant_derivatives(sparse_chol_and_grad = new_compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
                      ))))[,updated_comps[[regime]]]/ 2
        # Evaluate potential and kinetic energies at start and end of trajectory
        current_K = sum (state$momenta$range_beta_sufficient[,updated_comps[[regime]]] ^2) / 2
        proposed_U =
          (
            - Bidart::beta_prior_log_dens(beta = new_range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = state$params$range_log_scale) # normal prior
            + .5* sum((new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
            - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
          )
        
        state$transition_kernels$range_beta_sufficient[regime] = state$transition_kernels$range_beta_sufficient[regime] -1/sqrt(iter_start + iter +100)
        proposed_K = sum(p^2) / 2
        if(!is.nan(current_U-proposed_U+current_K- proposed_K))
        {
          if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
          {
            #print("tatata!")
            state$transition_kernels$range_beta_sufficient[regime] = state$transition_kernels$range_beta_sufficient[regime] + 2/sqrt(iter_start + iter +100)
            state$momenta$range_beta_sufficient[,updated_comps[[regime]]] = p
            state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
            state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
            state$params$range_beta[] = new_range_beta
            state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
          }
        }
      }
    }
    #############################
    # Variance of the  range PP #
    #############################
    
    if(hierarchical_model$range_PP){
      # ancillary-sufficient ####
      q = state$params$range_log_scale
      current_U =
        (
          + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
          - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
        )
      
      # MALA whitened
      state$momenta$range_log_scale_sufficient = sqrt(.9) * state$momenta$range_log_scale_sufficient + sqrt(.1)*rnorm(length(q))
      p = state$momenta$range_log_scale_sufficient
      # Make a half step for momentum at the beginning
      d_beta_d_scale = 
        Bidart::derivative_field_wrt_scale(
          state$params$range_beta[-seq(data$covariates$range_X$n_regressors),,drop = F], 
          state$params$range_log_scale
        )
      # derivative of potential wrt range beta
      d_potential_d_beta = as.matrix(
        + Bidart::X_PP_crossprod(
          X = data$covariates$range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
          Y = # Jacobian of range field wrt range_beta
            (# natural gradient of obs likelihood wrt range field
              Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                            left_vector = as.vector(state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                            right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                            NNarray = vecchia_approx$NNarray  
              )
              - Bidart::log_determinant_derivatives(sparse_chol_and_grad = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
            ))[-seq(data$covariates$range_X$n_regressors),,drop = F]
      )
      p = p - exp(state$transition_kernels$range_log_scale_sufficient) *
        (
          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
        )/2
      
      # Make a full step for the position
      q = q + exp(state$transition_kernels$range_log_scale_sufficient) * p
      new_range_beta = state$params$range_beta
      new_range_beta[-seq(data$covariates$range_X$n_regressors),] = new_range_beta[-seq(data$covariates$range_X$n_regressors),] %*% 
        solve(chol(Bidart::expmat(state$params$range_log_scale))) %*% chol(Bidart::expmat(q))
      
      new_compressed_sparse_chol_and_grad =
        Bidart::compute_sparse_chol(
          num_threads = num_threads, 
          anisotropic = hierarchical_model$anisotropic, 
          sphere = hierarchical_model$sphere, 
          range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
          locs = data$locs, 
          range_X = data$covariates$range_X$X_locs, 
          nu = hierarchical_model$nu, 
          PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
          compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
        )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      # Make a half step for momentum at the end.
      d_beta_d_scale = 
        Bidart::derivative_field_wrt_scale(
          new_range_beta[-seq(data$covariates$range_X$n_regressors),,drop = F], 
          q
        )
      # derivative of potential wrt range beta
      d_potential_d_beta = as.matrix(
        + Bidart::X_PP_crossprod(
          X = data$covariates$range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
          Y = # Jacobian of range field wrt range_beta
            (# natural gradient of obs likelihood wrt range field
              Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                            left_vector = as.vector(new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                            right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                            NNarray = vecchia_approx$NNarray  
              )
              - Bidart::log_determinant_derivatives(sparse_chol_and_grad = new_compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
            ))[-seq(data$covariates$range_X$n_regressors),,drop = F]
      )
      p = p - exp(state$transition_kernels$range_log_scale_sufficient) *
        (
          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
        )/2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (state$momenta$range_log_scale_sufficient ^2) / 2
      proposed_U =
        (
          + .5* sum((new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
          - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
        )
      proposed_K = sum(p^2) / 2
      state$transition_kernels$range_log_scale_sufficient  = state$transition_kernels$range_log_scale_sufficient - 1/sqrt(iter_start + iter +100)
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
        {
          state$transition_kernels$range_log_scale_sufficient  = state$transition_kernels$range_log_scale_sufficient + 2/sqrt(iter_start + iter +100)
          new_eigen = eigen(Bidart::expmat(q))$values
          if(
            all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
            all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2])))
          {
            state$momenta$range_log_scale_sufficient = p
            
            state$params$range_beta = new_range_beta
            state$params$range_log_scale[] = q
            
            state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
            state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
            state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
          }
        }
      }
      
      
      # sufficient - sufficient ####
      for(i in seq(10))
      {
        new_range_log_scale = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, .1)
        old_eigen = eigen(Bidart::expmat(state$params$range_log_scale))$values
        new_eigen = eigen(Bidart::expmat(new_range_log_scale))$values
        if(
          all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
          all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2]))&
          (
            + Bidart::beta_prior_log_dens(beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = new_range_log_scale)
            - Bidart::beta_prior_log_dens(beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = state$params$range_log_scale) 
            > log(runif(1))
          )
          
        )
        {
          state$params$range_log_scale = new_range_log_scale
        }
      }
      # ancillary - ancillary ####
      q = state$params$range_log_scale
      current_U =
        (
          + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
        )
      # MALA whitened
      state$momenta$range_log_scale_ancillary = sqrt(.9) * state$momenta$range_log_scale_ancillary + sqrt(.1)*rnorm(length(q))
      p = state$momenta$range_log_scale_ancillary
      # Make a half step for momentum at the beginning
      d_beta_d_scale = 
        Bidart::derivative_field_wrt_scale(
          state$params$range_beta[-seq(data$covariates$range_X$n_regressors),,drop = F], 
          state$params$range_log_scale
        )
      # derivative of potential wrt range beta
      d_potential_d_beta = Bidart::X_PP_crossprod(
        X = data$covariates$range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
        Y = # Jacobian of range field wrt range_beta
          (
            # natural gradient of obs likelihood wrt range field
            Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                          left_vector = as.vector(
                                            Matrix::solve(
                                              Matrix::t(state$sparse_chol_and_stuff$sparse_chol), 
                                              - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                            ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                              * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                            )), 
                                          right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                          NNarray = vecchia_approx$NNarray  
            )
          ))[-seq(data$covariates$range_X$n_regressors),,drop = F]
      p = p - exp(state$transition_kernels$range_log_scale_ancillary) *
        (
          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
        )/2
      
      q = q + exp(state$transition_kernels$range_log_scale_ancillary) * p
      new_range_beta = state$params$range_beta
      new_range_beta[-seq(data$covariates$range_X$n_regressors),] = new_range_beta[-seq(data$covariates$range_X$n_regressors),] %*% 
        solve(chol(Bidart::expmat(state$params$range_log_scale))) %*% chol(Bidart::expmat(q))
      
      new_compressed_sparse_chol_and_grad =
        Bidart::compute_sparse_chol(
          num_threads = num_threads,
          anisotropic = hierarchical_model$anisotropic, 
          sphere = hierarchical_model$sphere, 
          range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
          locs = data$locs, 
          range_X = data$covariates$range_X$X_locs, 
          nu = hierarchical_model$nu, 
          PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP,
          compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
        )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
      # Make a half step for momentum at the end.
      d_beta_d_scale = 
        Bidart::derivative_field_wrt_scale(
          new_range_beta[-seq(data$covariates$range_X$n_regressors),,drop = F], 
          q
        )
      # derivative of potential wrt range beta
      d_potential_d_beta = Bidart::X_PP_crossprod(
        X = data$covariates$range_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$range_PP, locs_idx = vecchia_approx$hctam_scol_1,
        Y = # Jacobian of range field wrt range_beta
          (
            # natural gradient of obs likelihood wrt range field
            Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                          left_vector = as.vector(
                                            Matrix::solve(
                                              Matrix::t(new_sparse_chol), 
                                              - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                            ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                              * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                            )), 
                                          right_vector = new_field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                          NNarray = vecchia_approx$NNarray  
            )
          )
      )[-seq(data$covariates$range_X$n_regressors),,drop = F]
      p = p - exp(state$transition_kernels$range_log_scale_ancillary) *
        (
          apply(d_beta_d_scale * array(rep(d_potential_d_beta, length(p)), dim = c(dim(d_potential_d_beta), length(p))), 3, sum)
        )/2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (state$momenta$range_log_scale_ancillary ^2) / 2
      proposed_U =
        (
          + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
        )
      proposed_K = sum(p^2) / 2
      state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary -1/sqrt(iter_start + iter +100)
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
        {
          new_eigen = eigen(Bidart::expmat(q))$values
          if(
            all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
            all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2])))
          {
            
            state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary + 2/sqrt(iter_start + iter +100)
            state$momenta$range_log_scale_ancillary = p
            
            state$params$range_beta[] = new_range_beta
            state$params$range_log_scale[] = q
            state$params$field = new_field
            
            state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
            state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
            state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
          }
        }
      }
      
      # sufficient - ancillary (equivalent to sufficient-sufficient) ####
      for(i in seq(10))
      {
        new_range_log_scale = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, .1)
        old_eigen = eigen(Bidart::expmat(state$params$range_log_scale))$values
        new_eigen = eigen(Bidart::expmat(new_range_log_scale))$values
        if(
          all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
          all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2]))&
          (
            + Bidart::beta_prior_log_dens(beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = new_range_log_scale)
            - Bidart::beta_prior_log_dens(beta = state$params$range_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$range_PP, 
                                          beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                          beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                          log_scale = state$params$range_log_scale) 
            > log(runif(1))
          )
          
        )
        {
          state$params$range_log_scale = new_range_log_scale
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
    state$sparse_chol_and_stuff$noise = Bidart::variance_field(
      beta = state$params$noise_beta, X = data$covariates$noise_X$X, 
      PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP, 
      locs_idx = NULL
    )
    # VEWY IMPOWTANT don't remove or comment
    squared_residuals = as.vector(state$sparse_chol_and_stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2
    # HMC update
    q = data$covariates$noise_X$chol_crossprod_X %*% state$params$noise_beta
    current_U =
      (
        - Bidart::beta_prior_log_dens(beta = state$params$noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
                                      beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                                      beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
                                      log_scale = state$params$noise_log_scale) # normal prior 
        +.5* sum(log(state$sparse_chol_and_stuff$noise)) # det
        +.5*sum(squared_residuals/state$sparse_chol_and_stuff$noise) # observations
      )
    # HMC whitened
    state$momenta$noise_beta = sqrt(.9) * state$momenta$noise_beta + sqrt(.1)*rnorm(length(state$momenta$noise_beta))
    p = state$momenta$noise_beta
    
    # Make a half step for momentum at the beginning
    p = p - exp(state$transition_kernels$noise_beta_mala) *
      (
        + solve(t(data$covariates$noise_X$chol_crossprod_X), # solving by prior chol because of whitening
                - Bidart::beta_prior_log_dens_derivative(beta = state$params$noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
                                                         beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                                                         beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
                                                         log_scale = state$params$noise_log_scale) # normal prior
                + Bidart::X_PP_crossprod(X = data$covariates$noise_X$X, PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP, 
                                         Y = 
                                           (
                                             + .5 # determinant part of normal likelihood
                                             - (squared_residuals/state$sparse_chol_and_stuff$noise)/2 # exponential part of normal likelihood
                                           ))
        ))/ 2
    # checking gradient with finite differences, to update
    
    #### idx = 200
    #### noise_beta_ = state$params$noise_beta 
    #### noise_beta_[idx] = noise_beta_[idx] + 0.0000001
    #### noise_ = Bidart::variance_field(
    ####   beta = noise_beta_, X = data$covariates$noise_X$X, 
    ####   PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP
    #### )
    #### U_ =
    ####   (
    ####     - Bidart::beta_prior_log_dens(beta = noise_beta_, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
    ####                                   beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
    ####                                   beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
    ####                                   log_scale = state$params$noise_log_scale) # normal prior 
    ####     +.5* sum(log(noise_)) # det
    ####     +.5*sum(squared_residuals/noise_) # observations
    ####   )
    #### print(
    ####   (
    ####             - Bidart::beta_prior_log_dens_derivative(beta = state$params$noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
    ####                                                      beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
    ####                                                      beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
    ####                                                      log_scale = state$params$noise_log_scale) # normal prior
    ####             + Bidart::X_PP_crossprod(X = data$covariates$noise_X$X, PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP, 
    ####                                      Y = 
    ####                                        (
    ####                                          + .5 # determinant part of normal likelihood
    ####                                          - (squared_residuals/state$sparse_chol_and_stuff$noise)/2 # exponential part of normal likelihood
    ####                                        ))
    ####     )[idx]
    #### )
    #### print(10000000*(U_- current_U))
    
    
    # Make a full step for the position
    q = q + exp(state$transition_kernels$noise_beta_mala) * p
    new_noise_beta = solve(data$covariates$noise_X$chol_crossprod_X, q)
    new_noise = Bidart::variance_field(beta = new_noise_beta, PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP, 
                                       X = data$covariates$noise_X$X, locs_idx = NULL)
    # Make a half step for momentum at the end
    p = p - exp(state$transition_kernels$noise_beta_mala) *
      (
        + solve(t(data$covariates$noise_X$chol_crossprod_X), # solving by prior sparse chol because of whitening
                - Bidart::beta_prior_log_dens_derivative(beta = new_noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
                                                         beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                                                         beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
                                                         log_scale = state$params$noise_log_scale) # normal prior  
                + Bidart::X_PP_crossprod(X = data$covariates$noise_X$X, PP = hierarchical_model$PP, use_PP = hierarchical_model$noise_PP, 
                                         (
                                           + .5 # determinant part of normal likelihood
                                           - (squared_residuals/new_noise)/2 # exponential part of normal likelihood
                                         ))
        ))/ 2
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (state$momenta$noise_beta ^2) / 2
    proposed_U = 
      (
        - Bidart::beta_prior_log_dens(beta = new_noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
                                      beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                                      beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
                                      log_scale = state$params$noise_log_scale) # normal prior        
        +.5* sum(log(new_noise)) # det
        +.5*sum(squared_residuals/new_noise) # observations
      )
    proposed_K = sum(p^2) / 2
    
    
    state$transition_kernels$noise_beta_mala = state$transition_kernels$noise_beta_mala- 1/sqrt(iter_start + iter +100)
    if(!is.nan(current_U-proposed_U+current_K- proposed_K))
    {
      if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
      {
        state$transition_kernels$noise_beta_mala = state$transition_kernels$noise_beta_mala + 2/sqrt(iter_start + iter +100)
        state$momenta$noise_beta = p
        state$params$noise_beta[] = new_noise_beta
        state$sparse_chol_and_stuff$noise = new_noise
      }
    }
    
    ###################
    # Noise log scale # 
    ###################
    if(hierarchical_model$noise_PP)
    {
      # ancillary -- sufficient ####
      new_noise_log_scale = state$params$noise_log_scale + rnorm(1, 0, exp(state$transition_kernels$noise_log_scale))
      new_noise_beta = state$params$noise_beta
      new_noise_beta[-seq(data$covariates$noise_X$n_regressors)] = new_noise_beta[-seq(data$covariates$noise_X$n_regressors)] *
        exp((new_noise_log_scale - state$params$noise_log_scale)/2)
      new_noise =Bidart::variance_field(beta = new_noise_beta, PP = hierarchical_model$PP, 
                                        use_PP = hierarchical_model$noise_PP, X = data$covariates$noise_X$X, 
                                        locs_idx = NULL)
      ll_ratio = (
        -.5* sum(log(new_noise)) 
        -.5*sum(squared_residuals/new_noise)
        +.5* sum(log(state$sparse_chol_and_stuff$noise)) 
        +.5*sum(squared_residuals/state$sparse_chol_and_stuff$noise)
      )
      if(!is.nan(ll_ratio))
      {
        if(ll_ratio > log(runif(1)))
        {
          if(
            (new_noise_log_scale > hierarchical_model$noise_log_scale_prior[1])&
            (new_noise_log_scale < hierarchical_model$noise_log_scale_prior[2])
          )
          {
            state$params$noise_log_scale = new_noise_log_scale
            state$params$noise_beta = new_noise_beta 
            state$sparse_chol_and_stuff$noise = new_noise
            state$transition_kernels$noise_log_scale = state$transition_kernels$noise_log_scale + 4/sqrt(iter_start + iter +100)
          }
        }
      }
      state$transition_kernels$noise_log_scale = state$transition_kernels$noise_log_scale - 1/sqrt(iter_start + iter +100)
      
      # sufficient -- sufficient ####
      for(i in seq(10))
      {
        new_noise_log_scale = state$params$noise_log_scale + rnorm(1, 0, .1)
        if(
          (
            + Bidart::beta_prior_log_dens(beta = state$params$noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
                                          beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                                          beta_precision = hierarchical_model$beta_priors$noise_beta_precision, 
                                          log_scale = new_noise_log_scale) - 
            Bidart::beta_prior_log_dens(beta = state$params$noise_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$noise_PP, 
                                        beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                                        beta_precision = hierarchical_model$beta_priors$noise_beta_precision, 
                                        log_scale = state$params$noise_log_scale)
          ) > log(runif(1))
          &(new_noise_log_scale > hierarchical_model$noise_log_scale_prior[1])
          &(new_noise_log_scale < hierarchical_model$noise_log_scale_prior[2])
        )
        {
          state$params$noise_log_scale = new_noise_log_scale
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
      sparse_chol_diag_field = state$sparse_chol_and_stuff$sparse_chol %*% Matrix::Diagonal(x = state$params$field)
      q = data$covariates$scale_X$chol_crossprod_X_locs %*% state$params$scale_beta
      current_U =
        (
          - Bidart::beta_prior_log_dens(beta = state$params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                        beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                        beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                        log_scale = state$params$scale_log_scale) # normal prior 
          +0.5*sum(log(state$sparse_chol_and_stuff$scale))# determinant part
          +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/state$sparse_chol_and_stuff$scale))^2)# covmat product part
        )
      # HMC whitened
      state$momenta$scale_beta_sufficient = sqrt(.9) * state$momenta$scale_beta_sufficient + sqrt(.1)*rnorm(length(state$momenta$scale_beta_sufficient))
      p = state$momenta$scale_beta_sufficient
      # Make a half step for momentum at the beginning
      p = p - exp(state$transition_kernels$scale_beta_sufficient_mala) *
        (
          + solve(t(data$covariates$scale_X$chol_crossprod_X_locs), # solving by t chol because of whitening
                  - Bidart::beta_prior_log_dens_derivative(beta = state$params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                                           beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                                           beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                                           log_scale = state$params$scale_log_scale) # normal prior
                  +  Bidart::X_PP_crossprod(X = data$covariates$scale_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, locs_idx = vecchia_approx$hctam_scol_1,
                                            Y = (
                                              .5  # determinant part 
                                              -.5 * sqrt(1/state$sparse_chol_and_stuff$scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/state$sparse_chol_and_stuff$scale)))# natural derivative
                                            ))
          ))/ 2
      
      # Make a full step for the position
      q = q + exp(state$transition_kernels$scale_beta_sufficient_mala) * p
      new_scale_beta = solve(data$covariates$scale_X$chol_crossprod_X_locs, q)
      new_scale =Bidart::variance_field(beta = new_scale_beta, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, 
                                        X = data$covariates$scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
      
      # Make a half step for momentum at the end.
      p = p - exp(state$transition_kernels$scale_beta_sufficient_mala) *
        (
          + solve(t(data$covariates$scale_X$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                  - Bidart::beta_prior_log_dens_derivative(beta = new_scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                                           beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                                           beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                                           log_scale = state$params$scale_log_scale) # normal prior  
                  + Bidart::X_PP_crossprod(X = data$covariates$scale_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, locs_idx = vecchia_approx$hctam_scol_1,
                                           (
                                             .5  # determinant part 
                                             -.5 * sqrt(1/new_scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/new_scale)))# natural derivative
                                           ))
          ))/ 2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (state$momenta$scale_beta_sufficient ^2) / 2
      proposed_U = 
        (
          - Bidart::beta_prior_log_dens(beta = new_scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                        beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                        beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                        log_scale = state$params$scale_log_scale) # normal prior 
          +0.5*sum(log(new_scale))# determinant part
          +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/new_scale))^2)# covmat product part
        )
      proposed_K = sum(p^2) / 2
      state$transition_kernels$scale_beta_sufficient_mala = state$transition_kernels$scale_beta_sufficient_mala - 1/sqrt(iter_start + iter +100)
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
        {
          state$transition_kernels$scale_beta_sufficient_mala = state$transition_kernels$scale_beta_sufficient_mala + 2/sqrt(iter_start + iter +100)
          state$momenta$scale_beta_sufficient = p
          state$params$scale_beta[] = new_scale_beta
          state$sparse_chol_and_stuff$scale = new_scale
        }
      }
      
      #############################
      # scale beta ancillary MALA #
      #############################
      q = data$covariates$scale_X$chol_crossprod_X_locs %*% state$params$scale_beta
      current_U =
        (
          - Bidart::beta_prior_log_dens(beta = state$params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                        beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                        beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                        log_scale = state$params$scale_log_scale) # normal prior 
          +.5 * sum((data$observed_field - state$sparse_chol_and_stuff$lm_fit - state$params$field[vecchia_approx$locs_match] )^2/state$sparse_chol_and_stuff$noise)
        )
      # MALA whitened
      state$momenta$scale_beta_ancillary = sqrt(.9) * state$momenta$scale_beta_ancillary   + sqrt(.1)*rnorm(length(state$momenta$scale_beta_ancillary))
      p = state$momenta$scale_beta_ancillary
      # Make a. half step for momentum at the beginning
      p = p - exp(state$transition_kernels$scale_beta_ancillary_mala) *
        (
          + solve(t(data$covariates$scale_X$chol_crossprod_X_locs), # solving by prior chol because of whitening
                  - Bidart::beta_prior_log_dens_derivative(beta = state$params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                                           beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                                           beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                                           log_scale = state$params$scale_log_scale) # normal prior 
                  + Bidart::X_PP_crossprod(X = data$covariates$scale_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, locs_idx = vecchia_approx$hctam_scol_1,
                                           (.5 * state$params$field * 
                                              as.vector(vecchia_approx$locs_match_matrix %*% 
                                                          ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                                             state$sparse_chol_and_stuff$noise)
                                              )))
          ))/ 2
      # Make a full step for the position
      q = q + exp(state$transition_kernels$scale_beta_ancillary_mala) * p
      new_scale_beta = solve(data$covariates$scale_X$chol_crossprod_X_locs, q)
      new_scale =Bidart::variance_field(beta = new_scale_beta, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, 
                                        X = data$covariates$scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
      new_field  = state$params$field*sqrt(new_scale)/sqrt(state$sparse_chol_and_stuff$scale)
      # Make a half step for momentum at the end.
      p = p - exp(state$transition_kernels$scale_beta_ancillary_mala) *
        (
          + solve(t(data$covariates$scale_X$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                  - Bidart::beta_prior_log_dens_derivative(beta = new_scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                                           beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                                           beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                                           log_scale = state$params$scale_log_scale) # normal prior  
                  + Bidart::X_PP_crossprod(X = data$covariates$scale_X$X_locs, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, locs_idx = vecchia_approx$hctam_scol_1,
                                           (.5 * new_field * 
                                              as.vector(vecchia_approx$locs_match_matrix %*% 
                                                          ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                                             state$sparse_chol_and_stuff$noise)
                                              )))
          ))/ 2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (state$momenta$scale_beta_ancillary ^2) / 2
      proposed_U = 
        (
          - Bidart::beta_prior_log_dens(beta = new_scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                        beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                        beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                        log_scale = state$params$scale_log_scale) # normal prior  
          +.5 * sum((data$observed_field - state$sparse_chol_and_stuff$lm_fit - new_field[vecchia_approx$locs_match] )^2/state$sparse_chol_and_stuff$noise)
        )
      proposed_K = sum(p^2) / 2
      state$transition_kernels$scale_beta_ancillary_mala = state$transition_kernels$scale_beta_ancillary_mala - 1/sqrt(iter_start + iter +100)
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
        {
          state$transition_kernels$scale_beta_ancillary_mala = state$transition_kernels$scale_beta_ancillary_mala + 2/sqrt(iter_start + iter +100)
          state$momenta$scale_beta_ancillary = p
          state$params$scale_beta[] = new_scale_beta
          state$params$field = new_field
          state$sparse_chol_and_stuff$scale = new_scale
        }
      }
      
      ###################
      # scale log scale #
      ###################
      if(hierarchical_model$scale_PP){
        # ancillary -- sufficient  ####
        # a change in hyperprior scale changes (rescales) the scale, which is then compared with the latent field
        new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, exp(state$transition_kernels$scale_log_scale_sufficient)) 
        new_scale_beta = state$params$scale_beta
        new_scale_beta[-seq(data$covariates$scale_X$n_regressors)] = new_scale_beta[-seq(data$covariates$scale_X$n_regressors)] *
          exp((new_scale_log_scale - state$params$scale_log_scale)/2)
        new_scale = Bidart::variance_field(beta = new_scale_beta, PP = hierarchical_model$PP, 
                                           use_PP = hierarchical_model$scale_PP, X = data$covariates$scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
        state$transition_kernels$scale_log_scale_sufficient = state$transition_kernels$scale_log_scale_sufficient - .25/sqrt(iter_start + iter +100)
        if(
          (
            (
              +.5* sum(log(state$sparse_chol_and_stuff$scale)) 
              -.5* sum(log(new_scale)) # log determinant
              +.5*sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
              -.5*sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(new_scale)))^2) # Gaussian density of the latent field
            )
            > log(runif(1))
          )
        )
        {
          if(
            (new_scale_log_scale > hierarchical_model$scale_log_scale_prior[1])&
            (new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
          )
          {
            state$transition_kernels$scale_log_scale_sufficient = state$transition_kernels$scale_log_scale_sufficient + 1/sqrt(iter_start + iter +100)
            state$params$scale_log_scale = new_scale_log_scale
            state$params$scale_beta = new_scale_beta 
            state$sparse_chol_and_stuff$scale = new_scale
          }
        }
        
        
        # sufficient -- sufficient ####
        for(i in seq(4))
        {
          new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, .1)
          if(
            (
              + Bidart::beta_prior_log_dens(beta = state$params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                            beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                            beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                            log_scale = new_scale_log_scale)     
              - Bidart::beta_prior_log_dens(beta = state$params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                            beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                            beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                            log_scale = state$params$scale_log_scale)     
              > log(runif(1))
            )
            &(new_scale_log_scale > hierarchical_model$scale_log_scale_prior[1])
            &(new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
          )
          {
            state$params$scale_log_scale = new_scale_log_scale
          }
        }
        # ancillary -- ancillary ####
        new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, exp(state$transition_kernels$scale_log_scale_ancillary)) 
        new_scale_beta = state$params$scale_beta
        new_scale_beta[-seq(data$covariates$scale_X$n_regressors)] = new_scale_beta[-seq(data$covariates$scale_X$n_regressors)] *
          exp((new_scale_log_scale - state$params$scale_log_scale)/2)
        new_scale =Bidart::variance_field(new_scale_beta, PP = hierarchical_model$PP, use_PP = hierarchical_model$scale_PP, 
                                          X = data$covariates$scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
        new_field = state$params$field * sqrt(new_scale)/sqrt(state$sparse_chol_and_stuff$scale)
        state$transition_kernels$scale_log_scale_ancillary = state$transition_kernels$scale_log_scale_ancillary - .25/sqrt(iter_start + iter +100)
        if(
          (
            -.5* sum((state$sparse_chol_and_stuff$lm_residuals -          new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) 
            +.5* sum((state$sparse_chol_and_stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) 
            > log(runif(1))
          )
        )
        {
          if(
            (new_scale_log_scale >  hierarchical_model$scale_log_scale_prior[1])
            &(new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
          )
          {  
            state$transition_kernels$scale_log_scale_ancillary = state$transition_kernels$scale_log_scale_ancillary + 1/sqrt(iter_start + iter +100)
            state$params$field = new_field
            state$params$scale_log_scale = new_scale_log_scale
            state$params$scale_beta = new_scale_beta 
            state$sparse_chol_and_stuff$scale = new_scale
          }
        }
        # sufficient -- ancillary ####
        for(i in seq(4))
        {
          new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, .1)
          if(
            (
              + Bidart::beta_prior_log_dens(beta = state$params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                            beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                            beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                            log_scale = new_scale_log_scale)     
              - Bidart::beta_prior_log_dens(beta = state$params$scale_beta, n_PP = hierarchical_model$PP$n_PP*hierarchical_model$scale_PP, 
                                            beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                            beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                            log_scale = state$params$scale_log_scale)     
              > log(runif(1))
            )
            &(new_scale_log_scale > hierarchical_model$scale_log_scale_prior[1])
            &(new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
          )
          {
            state$params$scale_log_scale = new_scale_log_scale
          }
        }
        
      }
      
    }
    

    #######################
    # Storing the samples #
    #######################
    if(thinning * iter == floor(thinning*iter))
    {
      for(name in names(state$params))
      {
        params_records[[name]][,,thinning*iter] = state$params[[name]]
      }
    }
  }
  return(list("state" = state, "params_records" = params_records))
}
