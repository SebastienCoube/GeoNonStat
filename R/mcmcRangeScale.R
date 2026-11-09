
# Samples the log-scale parameter of the range PP coeffs, using 
# ancillary parametrization of the PP coeffs 
# and sufficient parametrization of the latent field 
# (nested interweaving)
rangeLogScaleS <- function(state, hierarchical_model, vecchia_approx,
                           range_X, iter, iter_start, num_threads){
  
  range_reparam_mat <- rangeReparamMat(hierarchical_model)
  # initial density gradient ####
  grad_PP <- rangePPLogLikGradS(
    field_log_var = state$params$field_log_var,
    compressed_chol = state$stuff$compressed_chol, sparse_chol = state$stuff$sparse_chol,
    field = state$params$field, 
    range_reparam_mat = range_reparam_mat, 
    vecchia_approx = vecchia_approx, hierarchical_model = hierarchical_model,
    num_threads = num_threads)
  dens_grad <- logScaleDensGrad(
    PP_coeff = getPPCoeff(state$params$range_beta, hierarchical_model$range$PP$n_knots),
    grad_PP_coeff = grad_PP)
  # Step size and conditioning matrix ####
  stepsize <- exp(state$ker_var$range_log_scale_sufficient[1])
  cond_mat <- diag(rep(1, 1+hierarchical_model$anisotropic))
  # always lowering stepsize
  state$ker_var$range_log_scale_sufficient[1] <- updateKernel(
    iter = iter, iter_start = iter_start,
    kernel_value = state$ker_var$range_log_scale_sufficient[1], mult = -2
  )
  # renewing momenta ####
  state$momenta$range_log_scale_sufficient <-
    renewMomentum(state$momenta$range_log_scale_sufficient)
  # proposing new parameters ####
  new_range_log_scale <- moveForward(
    current_params= state$params$range_log_scale,
    dens_grad = dens_grad,
    momentum =state$momenta$range_log_scale_sufficient,
    stepsize = stepsize, cond_mat = cond_mat)
  
  # checking that uniform prior is respected ####
  if(all(inBounds(hierarchical_model$range$log_scale_bounds, new_range_log_scale))){
    # Updating PP coefficients ####
    new_range_beta <- movePPCoeff(
      beta = state$params$range_beta, 
      log_scale = state$params$range_log_scale, 
      new_log_scale = new_range_log_scale, 
      n_knots = hierarchical_model$range$PP$n_knots)
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
    grad_PP_back <- rangePPLogLikGradS(
      field_log_var = state$params$field_log_var,
      compressed_chol = state$stuff$proposed_compressed_chol, 
      sparse_chol = state$stuff$proposed_sparse_chol,
      field = state$params$field, 
      range_reparam_mat = range_reparam_mat, 
      vecchia_approx = vecchia_approx, hierarchical_model = hierarchical_model,
      num_threads = num_threads
    )
    dens_grad_back <- logScaleDensGrad(
      PP_coeff = getPPCoeff(new_range_beta, hierarchical_model$range$PP$n_knots),
      grad_PP_coeff = grad_PP_back)
    # deducing proposed momentum  ####
    new_momentum <- momentumBack(
      current_params = state$params$range_log_scale,
      proposed_params =  new_range_log_scale,
      dens_grad_back = dens_grad_back,
      stepsize = stepsize, cond_mat = cond_mat
    )
    
    # kinetic energies of momenta  ####
    current_momentum_dens <- momentumDens(state$momenta$range_log_scale_sufficient)
    proposed_momentum_dens <- momentumDens(new_momentum)
    # target densities ####
    current_dens <- rangeLogLikS(
      field_log_var = state$params$field_log_var, 
      sparse_chol = state$stuff$sparse_chol, 
      field = state$params$field)
    proposed_dens <- rangeLogLikS(
      field_log_var = state$params$field_log_var, 
      sparse_chol = state$stuff$proposed_sparse_chol, 
      field = state$params$field)
    
    #if(iter %/% 10 == iter / 10){
    #  par(mfrow = c(3, 2))
    #  plot(
    #    state$params$range_log_scale, 
    #    moveForward(stepsize, cond_mat, dens_grad_back, 
    #                new_range_log_scale, momentum = new_momentum),
    #    xlab = "current range log scale", ylab = "moving back proposed range log scale"
    #  )
    #  abline(a=0, b=1)
    #  plot(
    #    state$params$range_beta, 
    #    movePPCoeff(n_knots = hierarchical_model$range$PP$n_knots, 
    #                beta = new_range_beta, 
    #                log_scale = new_range_log_scale, 
    #                new_log_scale = moveForward(stepsize, cond_mat, dens_grad_back, 
    #                                                  new_range_log_scale, momentum = new_momentum)),
    #    xlab = "current range beta", ylab = "moving back proposed range beta"
    #  )
    #  abline(a=0, b=1)
    #  testrangePPLogLikGradS(
    #    range_beta = state$params$range_beta, 
    #    grad_to_test = grad_PP, 
    #    dens_to_test = current_dens, 
    #    num_threads = num_threads, 
    #    hierarchical_model = hierarchical_model, range_X = range_X,
    #    vecchia_approx = vecchia_approx, state = state)
    #  testrangePPLogLikGradS(
    #    range_beta = new_range_beta, 
    #    grad_to_test = grad_PP_back, 
    #    dens_to_test = proposed_dens, 
    #    num_threads = num_threads, 
    #    hierarchical_model = hierarchical_model, range_X = range_X,
    #    vecchia_approx = vecchia_approx, state = state)
    #  testGradientLogScaleS(
    #    range_log_scale = state$params$range_log_scale, 
    #    range_beta = state$params$range_beta, 
    #    grad_to_test = dens_grad, 
    #    dens_to_test = current_dens, 
    #    num_threads = num_threads, 
    #    hierarchical_model = hierarchical_model, range_X = range_X,
    #    vecchia_approx = vecchia_approx, state = state)
    #  testGradientLogScaleS(
    #    range_log_scale = new_range_log_scale, range_beta = new_range_beta, 
    #    grad_to_test = dens_grad_back, dens_to_test = proposed_dens, 
    #    num_threads = num_threads, hierarchical_model = hierarchical_model, range_X = range_X,
    #    vecchia_approx = vecchia_approx, state = state)
    #}
    
    # Metropolis ####
    ratio <- proposed_dens - current_dens + proposed_momentum_dens - current_momentum_dens
    if (!is.nan(ratio)) {
      if (log(runif(1)) < ratio) {
        #print("Log Scale Sufficient!!!")
        # increasing stepsize if acceptance
        state$ker_var$range_log_scale_sufficient[1] <- updateKernel(
          iter = iter, iter_start = iter_start,
          kernel_value = state$ker_var$range_log_scale_sufficient[1], mult = 3
        )
        # updating momenta
        state$momenta$range_log_scale_sufficient <- new_momentum
        # updating parameters
        state$params$range_log_scale[] <- new_range_log_scale
        state$params$range_beta[] <- new_range_beta
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
    state$momenta$range_log_scale_sufficient <- -state$momenta$range_log_scale_sufficient
  }
  return(state)
}

# Samples the log-scale parameter of the range PP coeffs, using 
# ancillary parametrization of the PP coeffs 
# and ancillary parametrization of the latent field 
# (nested interweaving)
rangeLogScaleA <- function(state, hierarchical_model, vecchia_approx,
                       range_X, iter, iter_start, num_threads){

  range_reparam_mat <- rangeReparamMat(hierarchical_model)
  # initial density gradient ####
  grad_PP <- rangePPLogLikGradA(
    field = state$params$field, lm_residuals = state$stuff$lm_residuals, 
    noise_var = state$stuff$noise_var,
    compressed_chol = state$stuff$compressed_chol, sparse_chol = state$stuff$sparse_chol,
    range_reparam_mat = range_reparam_mat, field_log_var = state$params$field_log_var,
    vecchia_approx = vecchia_approx, hierarchical_model = hierarchical_model,
    num_threads = num_threads)
  dens_grad <- logScaleDensGrad(
    PP_coeff = getPPCoeff(state$params$range_beta, hierarchical_model$range$PP$n_knots),
    grad_PP_coeff = grad_PP)
  # Step size and conditioning matrix ####
  stepsize <- exp(state$ker_var$range_log_scale_ancillary[1])
  cond_mat <- diag(rep(1, 1+hierarchical_model$anisotropic))
  # always lowering stepsize
  state$ker_var$range_log_scale_ancillary[1] <- updateKernel(
    iter = iter, iter_start = iter_start,
    kernel_value = state$ker_var$range_log_scale_ancillary[1], mult = -2
  )
  # renewing momenta ####
  state$momenta$range_log_scale_ancillary <-
    renewMomentum(state$momenta$range_log_scale_ancillary)
  # proposing new parameters ####
  new_range_log_scale <- moveForward(
    current_params= state$params$range_log_scale,
    dens_grad = dens_grad,
    momentum =state$momenta$range_log_scale_ancillary,
    stepsize = stepsize, cond_mat = cond_mat)

  # checking that uniform prior is respected ####
  if(all(inBounds(hierarchical_model$range$log_scale_bounds, new_range_log_scale))){
  # Updating PP coefficients ####
    new_range_beta <- movePPCoeff(
      beta = state$params$range_beta,
      log_scale = state$params$range_log_scale,
      new_log_scale = new_range_log_scale,
      n_knots = hierarchical_model$range$PP$n_knots)
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
    # computing new field ####
    new_field <- 
      as.matrix(Matrix::solve(state$stuff$proposed_sparse_chol,
                              state$stuff$sparse_chol %*% (state$params$field)))
    # computing gradient at proposed parameters ####
    grad_PP_back <-rangePPLogLikGradA(
      field = new_field, lm_residuals = state$stuff$lm_residuals, 
      noise_var = state$stuff$noise_var,
      compressed_chol = state$stuff$proposed_compressed_chol, sparse_chol = state$stuff$proposed_sparse_chol,
      range_reparam_mat = range_reparam_mat, field_log_var = state$params$field_log_var,
      vecchia_approx = vecchia_approx, hierarchical_model = hierarchical_model,
      num_threads = num_threads)
    dens_grad_back <- logScaleDensGrad(
      PP_coeff = getPPCoeff(new_range_beta, hierarchical_model$range$PP$n_knots),
      grad_PP_coeff = grad_PP_back)
    # deducing proposed momentum  ####
    new_momentum <- momentumBack(
      current_params = state$params$range_log_scale,
      proposed_params =  new_range_log_scale,
      dens_grad_back = dens_grad_back,
      stepsize = stepsize, cond_mat = cond_mat
    )
    # kinetic energies of momenta  ####
    current_momentum_dens <- momentumDens(state$momenta$range_log_scale_ancillary)
    proposed_momentum_dens <- momentumDens(new_momentum)
    # target densities ####
    current_dens <- rangeLogLikA(
      field = state$params$field,
      lm_residuals = state$stuff$lm_residuals, noise_var = state$stuff$noise_var, 
      vecchia_approx = vecchia_approx
    )
    proposed_dens <- rangeLogLikA(
      field = new_field,
      lm_residuals = state$stuff$lm_residuals, noise_var = state$stuff$noise_var, 
      vecchia_approx = vecchia_approx
    )

    # if(iter %/% 10 == iter / 10){
    #   par(mfrow = c(3, 2))
    #   plot(
    #     state$params$range_log_scale,
    #     moveForward(stepsize, cond_mat, dens_grad_back,
    #                 new_range_log_scale, momentum = new_momentum),
    #     xlab = "current range log scale", ylab = "moving back proposed range log scale"
    #   )
    #   abline(a=0, b=1)
    #   plot(
    #     state$params$range_beta,
    #     movePPCoeff(n_knots = hierarchical_model$range$PP$n_knots,
    #                 beta = new_range_beta,
    #                 log_scale = new_range_log_scale,
    #                 new_log_scale = moveForward(stepsize, cond_mat, dens_grad_back,
    #                                                   new_range_log_scale, momentum = new_momentum)),
    #     xlab = "current range beta", ylab = "moving back proposed range beta"
    #   )
    #   abline(a=0, b=1)
    #   testrangePPLogLikGradA(
    #     range_beta = state$params$range_beta, field = state$params$field, 
    #     sparse_chol = state$stuff$sparse_chol,
    #     grad_to_test = grad_PP,
    #     dens_to_test = current_dens,
    #     num_threads = num_threads,
    #     hierarchical_model = hierarchical_model, range_X = range_X,
    #     vecchia_approx = vecchia_approx, state = state)
    #   testrangePPLogLikGradA(
    #     range_beta = new_range_beta, field = new_field, 
    #     sparse_chol = state$stuff$proposed_sparse_chol,
    #     grad_to_test = grad_PP_back,
    #     dens_to_test = proposed_dens,
    #     num_threads = num_threads,
    #     hierarchical_model = hierarchical_model, range_X = range_X,
    #     vecchia_approx = vecchia_approx, state = state)
    #   testGradientLogScaleA(
    #     field = state$params$field, 
    #     range_log_scale = state$params$range_log_scale, 
    #     range_beta = state$params$range_beta, 
    #     grad_to_test = dens_grad, 
    #     dens_to_test = current_dens, 
    #     num_threads = num_threads, 
    #     hierarchical_model = hierarchical_model, range_X = range_X,
    #     vecchia_approx = vecchia_approx, state = state)
    #   testGradientLogScaleA(
    #     field = state$params$field, 
    #     range_log_scale = new_range_log_scale, range_beta = new_range_beta, 
    #     grad_to_test = dens_grad_back, dens_to_test = proposed_dens, 
    #     num_threads = num_threads, hierarchical_model = hierarchical_model, range_X = range_X,
    #     vecchia_approx = vecchia_approx, state = state)
    #   
    # }

     # Metropolis ####
     ratio <- proposed_dens - current_dens + proposed_momentum_dens - current_momentum_dens
     if (!is.nan(ratio)) {
       if (log(runif(1)) < ratio) {
         # increasing stepsize if acceptance
         state$ker_var$range_log_scale_ancillary[1] <- updateKernel(
           iter = iter, iter_start = iter_start,
           kernel_value = state$ker_var$range_log_scale_ancillary[1], mult = 3
         )
         # updating momenta
         state$momenta$range_log_scale_ancillary <- new_momentum
         # updating parameters
         state$params$range_log_scale[] <- new_range_log_scale
         state$params$range_beta[] <- new_range_beta
         state$params$field <- new_field
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
     state$momenta$range_log_scale_ancillary <- -state$momenta$range_log_scale_ancillary
  }
  return(state)
}

# Gradient of the field log likelihood with respect to the PP coefficients, +
# with sufficient parametrization of the latent field
rangePPLogLikGradS <- function(
    field_log_var, field, 
    compressed_chol, sparse_chol,
    range_reparam_mat,
    vecchia_approx, hierarchical_model,
    num_threads){
  # normal prior derivative, without range X because only PP
  - xPPCrossprod(
    X = NULL, PP = hierarchical_model$range$PP,
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
}


# Gradient of the field log likelihood with respect to the PP coefficients, +
# with ancillary parametrization of the latent field
rangePPLogLikGradA <- function(
    field, lm_residuals, noise_var, 
    compressed_chol, sparse_chol, 
    range_reparam_mat, field_log_var,
    vecchia_approx, hierarchical_model,
    num_threads){
  # normal prior derivative, without range X because only PP
  - xPPCrossprod(
    X = NULL, vecchia_approx = vecchia_approx,
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
}

