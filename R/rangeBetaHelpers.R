
rangeBetaS <- function(state, hierarchical_model, vecchia_approx,
                       range_X, iter, iter_start, num_threads){
  
  range_reparam_mat <- rangeReparamMat(hierarchical_model)
  # initial density gradient ####
  dens_grad <- rangeDensGradS(
    field_log_var = state$params$field_log_var, range_beta = state$params$range_beta,   
    compressed_chol = state$stuff$compressed_chol, sparse_chol = state$stuff$sparse_chol,   
    field = state$params$field, range_log_scale = state$params$range_log_scale, 
    range_reparam_mat = range_reparam_mat, range_X = range_X, 
    vecchia_approx = vecchia_approx, hierarchical_model = hierarchical_model, 
    num_threads = num_threads
  )
  # conditioning matrix and stepsize ####
  stepsize <- exp(state$ker_var$range_beta_sufficient[1])
  state$stuff$range_beta_conditioning_s <- 
    updateFisher(iter = iter, iter_start = iter_start, dens_grad = dens_grad, 
                 fisher_mat = state$stuff$range_beta_conditioning_s)
  cond_mat <- 
    condMat(range_X = range_X, hierarchical_model = hierarchical_model, 
            iter_start = iter_start, iter = iter, 
            fisher_mat = state$stuff$range_beta_conditioning_s)

    # renewing momenta ####
    state$momenta$range_beta_sufficient = 
    renewMomentum(state$momenta$range_beta_sufficient)
    state$momenta$field_log_var_sufficient_grouped = 
      renewMomentum(state$momenta$field_log_var_sufficient_grouped)
    # proposing new parameters ####
    move_forward <- moveForward(
      range_beta = state$params$range_beta, field_log_var = state$params$field_log_var, 
      dens_grad = dens_grad, 
      range_beta_momentum = state$momenta$range_beta_sufficient, 
      field_log_var_momentum = state$momenta$field_log_var_sufficient_grouped,
      stepsize = stepsize, cond_mat = cond_mat)
    new_range_beta <- state$params$range_beta
    new_field_log_var <- state$params$field_log_var
    new_range_beta[] <- move_forward[-1]
    new_field_log_var[] <- move_forward[1]
    # computing new Vecchia from new parameters ####
    vecchia_(
      start_idx = 1,
      log_range = t(computeLogRange(
          range_beta = new_range_beta, PP = hierarchical_model$range$PP,
          vecchia_approx = vecchia_approx, range_X = range_X)), 
      locs = vecchia_approx$t_locs, NNarray = vecchia_approx$NNarray,
      smoothness = hierarchical_model$matern_smoothness,
      compute_derivative = TRUE, num_threads = num_threads, 
      result = state$stuff$proposed_compressed_chol
    )
    state$stuff$proposed_sparse_chol@x <- 
      state$stuff$proposed_compressed_chol[, , 1][vecchia_approx$sparse_chol_x_reorder]
    # computing gradient at proposed parameters ####
    dens_grad_back <- rangeDensGradS(
      field_log_var = new_field_log_var, range_beta = new_range_beta, # changed
      compressed_chol = state$stuff$proposed_compressed_chol, #compressed
      sparse_chol = state$stuff$proposed_sparse_chol,  # changed
      field = state$params$field, range_log_scale = state$params$range_log_scale, 
      range_reparam_mat = range_reparam_mat, 
      range_X = range_X, vecchia_approx = vecchia_approx, 
      hierarchical_model = hierarchical_model, num_threads = num_threads
    )
    
    # deducing proposed momentum  ####
    new_momentum <- momentumBack(
      current_params = c(state$params$field_log_var, state$params$range_beta),
      proposed_params =  c(new_field_log_var, new_range_beta),
      dens_grad_back = dens_grad_back,
      stepsize = stepsize, cond_mat = cond_mat
    )
    # kinetic energies of momenta  ####
    current_momentum_dens <- momentumDens(c(state$momenta$range_beta_sufficient, 
                                    state$momenta$field_log_var_sufficient))
    proposed_momentum_dens <- momentumDens(new_momentum)
    # target densities ####
    current_dens <- rangeDensS(
      field_log_var = state$params$field_log_var, range_beta = state$params$range_beta, # changed
      compressed_chol = state$stuff$compressed_chol, sparse_chol = state$stuff$sparse_chol, # changed
      range_log_scale = state$params$range_log_scale, field = state$params$field, 
      lm_residuals = state$stuff$lm_residuals, noise_var = state$stuff$noise_var, 
      hierarchical_model = hierarchical_model, vecchia_approx = vecchia_approx)
    proposed_dens <- rangeDensS(
      field_log_var = new_field_log_var, range_beta = new_range_beta, # changed
      compressed_chol = state$stuff$proposed_compressed_chol, # changed
      sparse_chol = state$stuff$proposed_sparse_chol, # changed
      range_log_scale = state$params$range_log_scale, field = state$params$field, 
      lm_residuals = state$stuff$lm_residuals, noise_var = state$stuff$noise_var, 
      hierarchical_model = hierarchical_model, vecchia_approx = vecchia_approx)
    
    
    # if(iter %/% 10 == iter / 10){
    #   par(mfrow = c(2, 1))
    #   testGradientS(dens_to_test = current_dens, range_beta = state$params$range_beta, 
    #                 field_log_var = state$params$field_log_var, compressed_chol = state$stuff$compressed_chol, 
    #                 grad_to_test = dens_grad, sparse_chol = state$stuff$sparse_chol, 
    #                 num_threads = num_threads, 
    #                 range_X, hierarchical_model, vecchia_approx, state = state)
    #   testGradientS(dens_to_test = proposed_dens, range_beta = new_range_beta, 
    #                field_log_var = new_field_log_var, compressed_chol = state$stuff$proposed_compressed_chol, 
    #                grad_to_test = dens_grad_back, sparse_chol = state$stuff$proposed_sparse_chol, 
    #                num_threads = num_threads, 
    #                range_X, hierarchical_model, vecchia_approx, state = state)
    # }
    
    # Metropolis ####
    # always lowering stepsize
    state$ker_var$range_beta_sufficient[1] <- updateKernel(
      iter = iter, iter_start = iter_start,
      kernel_value = state$ker_var$range_beta_sufficient[1], mult = -2
    )
    ratio <- proposed_dens - current_dens + proposed_momentum_dens - current_momentum_dens
    if (!is.nan(ratio )) {
      if (log(runif(1)) < ratio) {
        # increasing stepsize if acceptance
        state$ker_var$range_beta_sufficient[1] <- updateKernel(
          iter = iter, iter_start = iter_start,
          kernel_value = state$ker_var$range_beta_sufficient[1], mult = 3
        )
        # updating momenta
        state$momenta$range_beta_sufficient <- new_momentum[-1]
        state$momenta$field_log_var_sufficient_grouped <- new_momentum[1]
        # updating parameters
        state$params$range_beta[] <- new_range_beta
        state$params$field_log_var[] <- new_field_log_var
        # invert positions of proposed and current matrix to avoid recreating objects.
        stuff_compressed <- match(c("compressed_chol", "proposed_compressed_chol"), names(state$stuff))
        stuff_compressed <- stuff_compressed[!is.na(stuff_compressed)]
        names(state$stuff)[stuff_compressed] <- names(state$stuff)[rev(stuff_compressed)]
        stuff_sparse <- match(c("sparse_chol", "proposed_sparse_chol"), names(state$stuff))
        stuff_sparse <- stuff_sparse[!is.na(stuff_sparse)]
        names(state$stuff)[stuff_sparse] <- names(state$stuff)[rev(stuff_sparse)]
      }
    }
    # negating momentum to induce "slippery slope" behavior
    state$momenta$range_beta_sufficient <- -state$momenta$range_beta_sufficient
    state$momenta$field_log_var_sufficient_grouped <- -state$momenta$field_log_var_sufficient_grouped
    
  return(state)
}


rangeBetaA <- function(state, hierarchical_model, vecchia_approx,
                       range_X, iter, iter_start, num_threads){
  range_reparam_mat <- rangeReparamMat(hierarchical_model)
  # initial density gradient ####
  dens_grad <- rangeDensGradA(
field_log_var = state$params$field_log_var
    , range_beta = state$params$range_beta,   
    compressed_chol = state$stuff$compressed_chol, sparse_chol = state$stuff$sparse_chol,   
    field = state$params$field, range_log_scale = state$params$range_log_scale, 
    range_reparam_mat = range_reparam_mat, range_X = range_X, 
    vecchia_approx = vecchia_approx, hierarchical_model = hierarchical_model, 
    num_threads = num_threads, 
    noise_var = state$stuff$noise_var, lm_residuals = state$stuff$lm_residuals
  )
  # conditioning matrix and stepsize ####
  stepsize <- exp(state$ker_var$range_beta_ancillary[1])
  state$stuff$range_beta_conditioning_a <- 
    updateFisher(iter = iter, iter_start = iter_start, dens_grad = dens_grad, 
                 fisher_mat = state$stuff$range_beta_conditioning_a)
  cond_mat <- 
    condMat(range_X = range_X, hierarchical_model = hierarchical_model, 
            iter_start = iter_start, iter = iter, 
            fisher_mat = state$stuff$range_beta_conditioning_a)

    # renewing momenta ####
    state$momenta$range_beta_ancillary = 
    renewMomentum(state$momenta$range_beta_ancillary)
    state$momenta$field_log_var_ancillary_grouped = 
      renewMomentum(state$momenta$field_log_var_ancillary_grouped)
    # proposing new parameters ####
    move_forward <- moveForward(
      range_beta = state$params$range_beta, field_log_var = state$params$field_log_var, 
      dens_grad = dens_grad, 
      range_beta_momentum = state$momenta$range_beta_ancillary, 
      field_log_var_momentum = state$momenta$field_log_var_ancillary_grouped,
      stepsize = stepsize, cond_mat = cond_mat)
    new_range_beta <- state$params$range_beta
    new_field_log_var <- state$params$field_log_var
    new_range_beta[] <- move_forward[-1]
    new_field_log_var[] <- move_forward[1]
    # computing new Vecchia from new parameters ####
    vecchia_(
      start_idx = 1,
      log_range = t(computeLogRange(
          range_beta = new_range_beta, PP = hierarchical_model$range$PP,
          vecchia_approx = vecchia_approx, range_X = range_X)), 
      locs = vecchia_approx$t_locs, NNarray = vecchia_approx$NNarray,
      smoothness = hierarchical_model$matern_smoothness,
      compute_derivative = TRUE, num_threads = num_threads, 
      result = state$stuff$proposed_compressed_chol
    )
    state$stuff$proposed_sparse_chol@x <- 
      state$stuff$proposed_compressed_chol[, , 1][vecchia_approx$sparse_chol_x_reorder]
    # computing new latent field ####
    new_field <-
      exp(.5 * (new_field_log_var[1,1] - state$params$field_log_var[1,1])) *
      as.vector(Matrix::solve(state$stuff$proposed_sparse_chol, 
                              state$stuff$sparse_chol %*% (state$params$field))) 
    # computing gradient at proposed parameters ####
    dens_grad_back <- rangeDensGradA(
      field_log_var = new_field_log_var, range_beta = new_range_beta, # changed
      compressed_chol = state$stuff$proposed_compressed_chol, # changed
      sparse_chol = state$stuff$proposed_sparse_chol,  # changed
      field = new_field,   # changed
      range_log_scale = state$params$range_log_scale, range_reparam_mat = range_reparam_mat, 
      lm_residuals = state$stuff$lm_residuals, noise_var = state$stuff$noise_var, 
      range_X = range_X, vecchia_approx = vecchia_approx, 
      hierarchical_model = hierarchical_model, num_threads = num_threads
    )
    
    # deducing proposed momentum  ####
    new_momentum <- momentumBack(
      current_params = c(state$params$field_log_var, state$params$range_beta),
      proposed_params =  c(new_field_log_var, new_range_beta),
      dens_grad_back = dens_grad_back,
      stepsize = stepsize, cond_mat = cond_mat
    )
    # kinetic energies of momenta  ####
    current_momentum_dens <- momentumDens(c(state$momenta$range_beta_ancillary, 
                                    state$momenta$field_log_var_ancillary))
    proposed_momentum_dens <- momentumDens(new_momentum)
    # target densities ####
    current_dens <- rangeDensA(
      field_log_var = state$params$field_log_var, range_beta = state$params$range_beta, # changed
      compressed_chol = state$stuff$compressed_chol, sparse_chol = state$stuff$sparse_chol, # changed
      field = state$params$field, 
      range_log_scale = state$params$range_log_scale, 
      lm_residuals = state$stuff$lm_residuals, 
      noise_var = state$stuff$noise_var, 
      hierarchical_model = hierarchical_model, 
      vecchia_approx = vecchia_approx)
    proposed_dens <- rangeDensA(
      field_log_var = new_field_log_var, range_beta = new_range_beta, # changed
      compressed_chol = state$stuff$proposed_compressed_chol, # changed
      sparse_chol = state$stuff$proposed_sparse_chol, # changed
      range_log_scale = state$params$range_log_scale, field = new_field, 
      lm_residuals = state$stuff$lm_residuals, noise_var = state$stuff$noise_var, 
      hierarchical_model = hierarchical_model, vecchia_approx = vecchia_approx)
    
    
     #if(iter %/% 10 == iter / 10){
     #  par(mfrow = c(2, 1))
     #  testGradientA(dens_to_test = current_dens, 
     #                range_beta = state$params$range_beta, 
     #                field_log_var = state$params$field_log_var, 
     #                sparse_chol = state$stuff$sparse_chol, 
     #                compressed_chol = state$stuff$compressed_chol, 
     #                field = state$params$field, 
     #                grad_to_test = dens_grad, 
     #                num_threads = num_threads, 
     #                range_X, hierarchical_model, vecchia_approx, state = state)
     #  testGradientA(dens_to_test = proposed_dens, range_beta = new_range_beta, 
     #               field_log_var = new_field_log_var, compressed_chol = state$stuff$proposed_compressed_chol, 
     #               field = new_field, 
     #               grad_to_test = dens_grad_back, sparse_chol = state$stuff$proposed_sparse_chol, 
     #               num_threads = num_threads, 
     #               range_X, hierarchical_model, vecchia_approx, state = state)
     #}
    
    # Metropolis ####
    # always lowering stepsize
    state$ker_var$range_beta_ancillary[1] <- updateKernel(
      iter = iter, iter_start = iter_start,
      kernel_value = state$ker_var$range_beta_ancillary[1], mult = -2
    )
    ratio <- proposed_dens - current_dens + proposed_momentum_dens - current_momentum_dens
    if (!is.nan(ratio )) {
      if (log(runif(1)) < ratio) {
        # increasing stepsize if acceptance
        state$ker_var$range_beta_ancillary[1] <- updateKernel(
          iter = iter, iter_start = iter_start,
          kernel_value = state$ker_var$range_beta_ancillary[1], mult = 3
        )
        # updating momenta
        state$momenta$range_beta_ancillary <- new_momentum[-1]
        state$momenta$field_log_var_ancillary_grouped <- new_momentum[1]
        # updating parameters
        state$params$range_beta[] <- new_range_beta[]
        state$params$field_log_var[] <- new_field_log_var[]
        state$params$field <- new_field
        # invert positions of proposed and current matrix to avoid recreating objects.
        stuff_compressed <- match(c("compressed_chol", "proposed_compressed_chol"), names(state$stuff))
        stuff_compressed <- stuff_compressed[!is.na(stuff_compressed)]
        names(state$stuff)[stuff_compressed] <- names(state$stuff)[rev(stuff_compressed)]
        stuff_sparse <- match(c("sparse_chol", "proposed_sparse_chol"), names(state$stuff))
        stuff_sparse <- stuff_sparse[!is.na(stuff_sparse)]
        names(state$stuff)[stuff_sparse] <- names(state$stuff)[rev(stuff_sparse)]
        #print("ANCILLARYYYY")
      }
    }
    # negating momentum to induce "slippery slope" behavior
    state$momenta$range_beta_ancillary <- -state$momenta$range_beta_ancillary
    state$momenta$field_log_var_ancillary_grouped <- -state$momenta$field_log_var_ancillary_grouped
    
  return(state)
}


momentumDens <- function(momentum)-.5*sum(momentum^2)

moveForward <- function(stepsize, cond_mat, dens_grad, 
                        range_beta_momentum, field_log_var_momentum,
                        range_beta, field_log_var){
  move_forward <-
    as.vector(stepsize * (cond_mat %*% (crossprod(cond_mat, c(dens_grad))))/2) + 
    sqrt(stepsize) * as.vector(cond_mat %*% c( field_log_var_momentum, range_beta_momentum))
  new_range_beta <- range_beta + move_forward[-1]
  new_field_log_var <- field_log_var[] + move_forward[1]
  c(new_field_log_var[], new_range_beta[])
}

rangeReparamMat <- function(hierarchical_model){
  range_reparam_mat <- matrix(2)
  if (hierarchical_model$anisotropic) {
    range_reparam_mat <- matrix(c(
      2, 2, 0,
      2, -2, 0,
      0, 0, 2 * sqrt(2)
    ), 3)
  }
  range_reparam_mat
}

condMat <- function(
    fisher_mat, 
    iter, iter_start, 
    hierarchical_model, range_X){
  cond_mat <- solve(
    as(
      (1 - min(1, max((500 - (iter + iter_start))/200, .05))) *
        fisher_mat / sqrt(sum(fisher_mat^2)) 
      + Matrix::bdiag(
        1, 
        min(1, max((500 - (iter + iter_start))/200, .05)) *
          (diag(rep(1, 1+ 2*hierarchical_model$anisotropic)) %x% 
             as.matrix(range_X$crossprod_X)) / 
          sqrt(3*sum(as.matrix(range_X$crossprod_X)^2)))
      , 
      "sparseMatrix"
    )
  )
  cond_mat <- cond_mat / sqrt(sum(cond_mat^2))
  cond_mat <- t(chol(cond_mat))
  cond_mat
}

updateFisher <- function(iter, iter_start, dens_grad, fisher_mat){
  if(iter + iter_start >200){
    fisher_mat <- 
      ((iter+iter_start)/(iter+iter_start-1)) * fisher_mat + 
      tcrossprod(c(dens_grad)) / (iter+iter_start)
  }
  fisher_mat
}

momentumBack <- function(
  current_params, proposed_params, 
  dens_grad_back, 
  cond_mat, stepsize
){
  solve(
    sqrt(stepsize) * cond_mat,
    current_params
    - proposed_params
    - as.vector(cond_mat %*% (crossprod(cond_mat, c(dens_grad_back)))) * stepsize / 2
  )
}

rangeDensGradS <- function(
    field_log_var, range_beta, field, 
    range_log_scale,
    compressed_chol, sparse_chol, 
    range_reparam_mat, 
    vecchia_approx, hierarchical_model, range_X, 
    num_threads){
  c(
    # log var
    betaPriorLogDensDerivative(
      beta = as.matrix(field_log_var), n_PP = 0, log_scale = 0,
      beta0_mean = hierarchical_model$scale$beta0_mean,
      beta0_var = hierarchical_model$scale$beta0_sd^2
    )
    + .5 * sum((sparse_chol %*% field)^2) / exp(field_log_var[1, 1]) # observation ll
    - .5 * vecchia_approx$n_locs,
    # range
    +betaPriorLogDensDerivative(
      beta = range_beta,
      n_PP = hierarchical_model$range$PP$n_knots,
      beta0_mean = hierarchical_model$range$beta0_mean,
      beta0_var = hierarchical_model$range$beta0_sd^2,
      log_scale = range_log_scale
    ) # normal prior
    # normal prior derivative
    - xPPCrossprod(
      X = range_X$X_locs, PP = hierarchical_model$range$PP,
      permutate_PP_to_obs = FALSE,
      vecchia_approx = vecchia_approx,
      Y = # Jacobian of range field wrt range_beta
        t( # natural gradient of obs likelihood wrt range field
          derivativeSandwiches_(
            vecchia = compressed_chol, # derivative of the (unscaled) NNGP factor
            left_vector = as.vector(sparse_chol %*% (field / exp(.5 * field_log_var[1, 1]))), # left vector = whitened latent field
            right_vector = field / exp(.5 * field_log_var[1, 1]), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
            NNarray = vecchia_approx$NNarray,
            sauce_determinant_chef = TRUE,
            num_threads = num_threads
          )
        )
    ) %*% range_reparam_mat
  )
}

rangeDensS <- function(
    field_log_var, range_beta, field, range_log_scale, 
    compressed_chol, sparse_chol, 
    lm_residuals, noise_var, 
    hierarchical_model, vecchia_approx){
  (
    # field log var prior
    betaPriorLogDens(
      beta = as.matrix(field_log_var), n_PP = 0, log_scale = 0,
      beta0_mean = hierarchical_model$scale$beta0_mean,
      beta0_var = hierarchical_model$scale$beta0_sd^2
    )
    # range
    + betaPriorLogDens(
      beta = range_beta,
      n_PP = hierarchical_model$range$PP$n_knots,
      beta0_mean = hierarchical_model$range$beta0_mean,
      beta0_var = hierarchical_model$range$beta0_sd^2,
      log_scale = range_log_scale
    )
    # data ll
    - .5 * sum((sparse_chol %*% (field / exp(.5 * field_log_var[1, 1])))^2)
    + sum(log(Matrix::diag(sparse_chol)))
    - vecchia_approx$n_locs * (.5 * field_log_var[1, 1])
  )

}

rangeDensGradA <- function(
    field_log_var, range_beta, field, compressed_chol, sparse_chol, 
    range_log_scale, lm_residuals, noise_var,
    vecchia_approx, hierarchical_model, range_X, 
    range_reparam_mat, num_threads){
  
  c(
    # log var
    betaPriorLogDensDerivative(
      beta = as.matrix(field_log_var), n_PP = 0, log_scale = 0,
      beta0_mean = hierarchical_model$scale$beta0_mean,
      beta0_var = hierarchical_model$scale$beta0_sd^2) + 
      .5 * sum(field[vecchia_approx$locs_match] * 
                 (lm_residuals - field[vecchia_approx$locs_match]) / noise_var),
    # field range
    +betaPriorLogDensDerivative(
      beta = range_beta, n_PP = hierarchical_model$range$PP$n_knots,
      beta0_mean = hierarchical_model$range$beta0_mean,
      beta0_var = hierarchical_model$range$beta0_sd^2,
      range_log_scale
    ) # normal prior
    - xPPCrossprod(
      X = range_X$X_locs, vecchia_approx = vecchia_approx,
      permutate_PP_to_obs = FALSE,
      PP = hierarchical_model$range$PP,
      Y = # Jacobian of range field wrt range_beta
        t(
          # natural gradient of obs likelihood wrt range field
          derivativeSandwiches_(
            vecchia = compressed_chol,
            left_vector = as.vector(
              Matrix::solve(
                Matrix::t(sparse_chol),
                -as.vector(vecchia_approx$locs_match_matrix %*% # gradient of  Gaussian observations ll wrt latent field
                             ((field[vecchia_approx$locs_match] - lm_residuals) / noise_var))
                * exp(.5 * field_log_var[1, 1]) 
              )
            ),
            right_vector = field / exp(.5 * field_log_var[1, 1]), 
            NNarray = vecchia_approx$NNarray,
            sauce_determinant_chef = FALSE,
            num_threads = num_threads
          )
        ) %*% range_reparam_mat
    )
  )
}

rangeDensA <- function(
    field_log_var, range_beta, field, compressed_chol, sparse_chol, 
    range_log_scale, lm_residuals, noise_var,
    vecchia_approx, hierarchical_model, range_X, 
    range_reparam_mat){
  (
    # field log var prior
   +betaPriorLogDens(
     beta = as.matrix(field_log_var), n_PP = 0, log_scale = 0,
     beta0_mean = hierarchical_model$scale$beta0_mean,
     beta0_var = hierarchical_model$scale$beta0_sd^2
   )
   # field range prior
   + betaPriorLogDens(
     beta = range_beta,
     n_PP = hierarchical_model$range$PP$n_knots,
     beta0_mean = hierarchical_model$range$beta0_mean,
     beta0_var = hierarchical_model$range$beta0_sd^2,
     range_log_scale
   )
    # data ll
    - .5 * sum((lm_residuals - field[vecchia_approx$locs_match])^2 / noise_var) # observation ll
  )
  
}


testGradientS <- function(
    range_beta, field_log_var, 
    compressed_chol, sparse_chol, 
    grad_to_test, dens_to_test, num_threads,
    range_X, hierarchical_model, vecchia_approx, 
    state
){
  empirical_grad = c()
  sparse_chol_grad_test = sparse_chol
  compressed_chol_grad_test = compressed_chol
  for(i in seq(length(range_beta)+1)){
    range_beta_grad_test <- range_beta 
    field_log_var_grad_test <- field_log_var 
    if(i==1) field_log_var_grad_test = field_log_var_grad_test + 1e-6
    if(i>1) range_beta_grad_test[i-1] = range_beta_grad_test[i-1] + 1e-6
    
    # re-computing Vecchia approx
    log_range <- computeLogRange(
      range_beta = range_beta_grad_test,
      PP = hierarchical_model$range$PP,
      vecchia_approx = vecchia_approx,
      range_X = range_X
    )
    vecchia_(start_idx = 1,
             log_range = t(log_range), locs = vecchia_approx$t_locs,
             NNarray = vecchia_approx$NNarray,
             smoothness = hierarchical_model$matern_smoothness,
             compute_derivative = TRUE, num_threads = num_threads, result = compressed_chol_grad_test
    )
    sparse_chol_grad_test@x <- compressed_chol_grad_test[,,1][vecchia_approx$sparse_chol_x_reorder]
    
    test_dens <- rangeDensS(
      field_log_var = field_log_var_grad_test, range_beta = range_beta_grad_test, # changed
      compressed_chol = compressed_chol_grad_test, # changed
      sparse_chol = sparse_chol_grad_test, # changed
      range_log_scale = state$params$range_log_scale, field = state$params$field, 
      lm_residuals = state$stuff$lm_residuals, noise_var = state$stuff$noise_var, 
      hierarchical_model = hierarchical_model, vecchia_approx = vecchia_approx)
    empirical_grad = c(empirical_grad,  (test_dens - dens_to_test)*1e6)
  }
  plot(grad_to_test, empirical_grad,  main  = "sufficient")
  abline(a=0, b=1)
  abline(h=0)
  abline(v=0)
}
  
  

testGradientA <- function(
    range_beta, field_log_var, 
    compressed_chol, sparse_chol, 
    field,
    grad_to_test, dens_to_test, num_threads,
    range_X, hierarchical_model, vecchia_approx, 
    state
){
  empirical_grad = c()
  sparse_chol_grad_test = sparse_chol
  compressed_chol_grad_test = compressed_chol
  for(i in seq(length(range_beta)+1)){
    range_beta_grad_test <- range_beta 
    field_log_var_grad_test <- field_log_var 
    if(i==1) field_log_var_grad_test[] = field_log_var_grad_test[] + 1e-6
    if(i>1) range_beta_grad_test[i-1] = range_beta_grad_test[i-1] + 1e-6
    
    # re-computing Vecchia approx
    log_range <- computeLogRange(
      range_beta = range_beta_grad_test,
      PP = hierarchical_model$range$PP,
      vecchia_approx = vecchia_approx,
      range_X = range_X
    )
    vecchia_(start_idx = 1,
             log_range = t(log_range), locs = vecchia_approx$t_locs,
             NNarray = vecchia_approx$NNarray,
             smoothness = hierarchical_model$matern_smoothness,
             compute_derivative = TRUE, num_threads = num_threads, 
             result = compressed_chol_grad_test
    )
    sparse_chol_grad_test@x <- compressed_chol_grad_test[,,1][vecchia_approx$sparse_chol_x_reorder]
    field_grad_test <-
      exp(.5*(field_log_var_grad_test[1,1] - field_log_var[1,1])) *
      as.vector(Matrix::solve(sparse_chol_grad_test, sparse_chol %*% field))
    test_dens <- rangeDensA(
      field_log_var = field_log_var_grad_test, 
      range_beta = range_beta_grad_test, # changed
      compressed_chol = compressed_chol_grad_test, # changed
      field = field_grad_test, sparse_chol = sparse_chol_grad_test, # changed
      range_log_scale = state$params$range_log_scale, 
      lm_residuals = state$stuff$lm_residuals, noise_var = state$stuff$noise_var, 
      hierarchical_model = hierarchical_model, vecchia_approx = vecchia_approx, 
      range_reparam_mat = range_reparam_mat, range_X = range_X)
    empirical_grad = c(empirical_grad,  (test_dens - dens_to_test)*1e6)
  }
  plot(grad_to_test, empirical_grad,  main  = "ancillary")
  abline(a=0, b=1)
  abline(h=0)
  abline(v=0)
}
  
  