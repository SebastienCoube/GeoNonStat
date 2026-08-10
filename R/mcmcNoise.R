

noiseBetaDensGrad <- function(noise_beta, noise_var, 
                          squared_residuals, noise_log_scale, 
                          hm_noise, noise_X, vecchia_approx){
  return(
    betaPriorLogDensDerivative(
      beta = noise_beta, n_PP = hm_noise$PP$n_knots,
      beta0_mean = hm_noise$beta0_mean,
      beta0_var = hm_noise$beta0_sd^2,
      log_scale = noise_log_scale
    ) # normal prior
    - xPPCrossprod(
      X = noise_X$X, hm_noise$PP, vecchia_approx = vecchia_approx, permutate_PP_to_obs = TRUE,
      Y =
        (
          +.5 # determinant part of normal likelihood
          - (squared_residuals / noise_var) / 2 # exponential part of normal likelihood
        )
    )
    )
}

noiseBetaDens <- function(noise_beta, noise_var,
                      hm_noise, noise_log_scale, squared_residuals){
  return(
    betaPriorLogDens(
      beta = noise_beta,
      n_PP = hm_noise$PP$n_knots,
      beta0_mean = hm_noise$beta0_mean,
      beta0_var = hm_noise$beta0_sd^2,
      log_scale = noise_log_scale
    ) # normal prior
    + noiseLogLik(noise_var, squared_residuals)
  )
}

noiseLogLik <- function(noise_var, squared_residuals){
  return(
    - .5 * sum(log(noise_var)) # det
    - .5 * sum(squared_residuals / noise_var) # observations
  )
}

noisePPLogLikGrad <- function(noise_var, 
                             squared_residuals, 
                             hm_noise, vecchia_approx){
  return(
    - xPPCrossprod(
      X = NULL, hm_noise$PP, vecchia_approx = vecchia_approx, permutate_PP_to_obs = TRUE,
      Y =
        (
          +.5 # determinant part of normal likelihood
          - (squared_residuals / noise_var) / 2 # exponential part of normal likelihood
        )
    )
  )
}

testGradientNoise <- function(
    noise_beta, noise_var, 
    grad_to_test, dens_to_test, squared_residuals,
    noise_X, hm_noise, vecchia_approx, noise_log_scale, 
    main
){
  empirical_grad = c()
  for(i in seq(length(noise_beta))){
    noise_beta_grad_test <- noise_beta 
    noise_beta_grad_test[i] = noise_beta_grad_test[i] + 1e-6
    
    # re-computing Vecchia approx
    noise_var_test <- exp(xPPMultRight(
      X = noise_X$X, PP = hm_noise$PP,
      vecchia_approx = vecchia_approx, Y = noise_beta_grad_test,
      permutate_PP_to_obs = TRUE
    ))[]
    
    test_dens <- noiseBetaDens(
      noise_beta = noise_beta_grad_test, noise_var = noise_var_test, 
      squared_residuals = squared_residuals, noise_log_scale = noise_log_scale,
      hm_noise = hm_noise
    )
    empirical_grad = c(empirical_grad,  (test_dens - dens_to_test)*1e6)
  }
  plot(grad_to_test, empirical_grad,  main  =  main)
  abline(a=0, b=1)
  abline(h=0)
  abline(v=0)
}




noiseBeta <- function(state, hm_noise, noise_X, vecchia_approx, iter, iter_start) {
  stepsize <- exp(state$ker_var$noise_beta)
  # residuals 
  squared_residuals <- as.matrix(state$stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2
  # initial gradient 
  dens_grad <- noiseBetaDensGrad(
    noise_var = state$stuff$noise_var, noise_beta = state$params$noise_beta, 
    noise_X = noise_X, 
    vecchia_approx = vecchia_approx, 
    hm_noise = hm_noise, 
    squared_residuals = squared_residuals, 
    noise_log_scale = state$params$noise_log_scale)
  # conditioning matrix update
  cond_mat <- noise_X$t_chol_solve_crossprod_X
  grad_record_for_Fisher <- list()
  # pre-allocation
  new_noise_var <- state$stuff$noise_var
  
  for (mala_rep in seq(2)) {
    # for conditioning matrix update
    grad_record_for_Fisher[[mala_rep]] <- as.vector(dens_grad)
    # updating momentum
    state$momenta$noise_beta <- renewMomentum(state$momenta$noise_beta, .9)
    # proposing parameter ####
    new_noise_beta <- moveForward(
      position = state$params$noise_beta, momentum = state$momenta$noise_beta, 
      dens_grad = dens_grad, cond_mat = cond_mat, stepsize = stepsize)
    new_noise_var[] <- exp(xPPMultRight(
      X = noise_X$X, PP = hm_noise$PP,
      vecchia_approx = vecchia_approx, Y = new_noise_beta,
      permutate_PP_to_obs = TRUE
    ))[]
    # caculating new gradient ####
    new_dens_grad <- noiseBetaDensGrad(
      noise_var = new_noise_var, noise_beta = new_noise_beta, 
      noise_X = noise_X, vecchia_approx = vecchia_approx, 
      hm_noise = hm_noise, 
      noise_log_scale = state$params$noise_log_scale, squared_residuals = squared_residuals)
    
    # deducing proposed momentum  ####
    new_momentum <- momentumBack(
      current_params = state$params$noise_beta,
      proposed_params =  new_noise_beta,
      dens_grad_back = new_dens_grad,
      stepsize = stepsize, cond_mat = cond_mat
    )
    # kinetic energies of momenta  ####
    current_momentum_dens <- momentumDens(state$momenta$noise_beta)
    proposed_momentum_dens <- momentumDens(new_momentum)
    # target densities ####
    current_dens <- noiseBetaDens(
      noise_beta = state$params$noise_beta, noise_var = state$stuff$noise_var, 
      squared_residuals = squared_residuals, 
      noise_log_scale = state$params$noise_log_scale,
      hm_noise = hm_noise)
    proposed_dens <- noiseBetaDens(
      noise_beta = new_noise_beta, noise_var = new_noise_var, 
      squared_residuals = squared_residuals, noise_log_scale = state$params$noise_log_scale,
      hm_noise = hm_noise)
    # always lowering kernel var ####
    state$ker_var$noise_beta_mala <- updateKernel(
      iter = iter,
      iter_start = iter_start,
      kernel_value = state$ker_var$noise_beta_mala,
      mult = -.6
    )
    
    ## Testing gradient ####
    #if(iter %/%10 == iter /10){
    #  par(mfrow = c(2,1))
    #  testGradientNoise(
    #    noise_beta = state$params$noise_beta, noise_var  =state$stuff$noise_var, 
    #    grad_to_test = dens_grad, dens_to_test = current_dens,
    #    noise_X = noise_X, hm_noise=  hm_noise, vecchia_approx = vecchia_approx, 
    #    noise_log_scale = state$params$noise_log_scale, squared_residuals = squared_residuals, 
    #    main = paste("Forward, iter =", iter, "MALA = ", mala_rep)
    #  )
    #  testGradientNoise(
    #    noise_beta = new_noise_beta, noise_var  = new_noise_var, 
    #    grad_to_test = new_dens_grad, dens_to_test = proposed_dens,
    #    noise_X = noise_X, hm_noise=  hm_noise, vecchia_approx = vecchia_approx, 
    #    noise_log_scale = state$params$noise_log_scale, squared_residuals = squared_residuals,
    #    main = paste("Back, iter =", iter, "MALA = ", mala_rep)
    #  )
    #}
    
    # Metropolis ####
    ratio <- proposed_dens - current_dens + proposed_momentum_dens - current_momentum_dens
    if (!is.nan(ratio)) {
      if (log(runif(1)) < ratio) {
        state$ker_var$noise_beta_mala <- updateKernel(
          iter = iter,
          iter_start = iter_start,
          kernel_value = state$ker_var$noise_beta_mala,
          mult = 1
        )
        dens_grad[] <- new_dens_grad[]
        state$momenta$noise_beta[] <- new_momentum[]
        state$params$noise_beta[] <- new_noise_beta
        state$stuff$noise_var[] <- new_noise_var[]
      }
    }
    # negating momentum ####
    state$momenta$noise_beta[] <- -state$momenta$noise_beta[]
  }
  
  
  return(state)
}



#' Samples the log-scale parameter of the noise variance PP coeffs, using 
#' ancillary parametrization of the PP coeffs 
#' (interweaving)
noiseLogScale <- function(state, hierarchical_model, vecchia_approx,
                           noise_X, iter, iter_start){
  # residuals 
  squared_residuals <- as.matrix(state$stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2
  # initial density gradient ####
  grad_PP <- noisePPLogLikGrad(noise_var = state$stuff$noise_var,
    squared_residuals = squared_residuals, hm_noise = hierarchical_model$noise, vecchia_approx = vecchia_approx)
  dens_grad <- logScaleDensGrad(
    PP_coeff = getPPCoeff(state$params$noise_beta, hierarchical_model$noise$PP$n_knots),
    grad_PP_coeff = grad_PP)
  # Step size and conditioning matrix ####
  stepsize <- exp(state$ker_var$noise_log_scale)
  cond_mat <- matrix(1)
  # always lowering stepsize
  state$ker_var$noise_log_scale <- updateKernel(
    iter = iter, iter_start = iter_start,
    kernel_value = state$ker_var$noise_log_scale, mult = -2
  )
  # renewing momenta ####
  state$momenta$noise_log_scale <-
    renewMomentum(state$momenta$noise_log_scale)
  # proposing new parameters ####
  new_noise_log_scale <- moveForward(
    position= state$params$noise_log_scale,
    dens_grad = dens_grad,
    momentum =state$momenta$noise_log_scale,
    stepsize = stepsize, cond_mat = cond_mat)
  
  # checking that uniform prior is respected ####
  if(all(inBounds(hierarchical_model$noise$log_scale_bounds, new_noise_log_scale))){
    # Updating PP coefficients ####
    new_noise_beta <- movePPCoeff(
      beta = state$params$noise_beta, 
      log_scale = state$params$noise_log_scale, 
      new_log_scale = new_noise_log_scale, 
      n_knots = hierarchical_model$noise$PP$n_knots)
    # computing new noise var from new parameters ####
    new_noise_var <- as.vector(exp(xPPMultRight(
      X = noise_X$X, PP = hierarchical_model$noise$PP,
      vecchia_approx = vecchia_approx, Y = new_noise_beta,
      permutate_PP_to_obs = TRUE
    )))
    # computing gradient at proposed parameters ####
    grad_PP_back <- noisePPLogLikGrad(
      noise_var = new_noise_var,
      squared_residuals = squared_residuals, hm_noise = hierarchical_model$noise, vecchia_approx = vecchia_approx)
    dens_grad_back <- logScaleDensGrad(
      PP_coeff = getPPCoeff(new_noise_beta, hierarchical_model$no$PP$n_knots),
      grad_PP_coeff = grad_PP_back)
    # deducing proposed momentum  ####
    new_momentum <- momentumBack(
      current_params = state$params$noise_log_scale,
      proposed_params =  new_noise_log_scale,
      dens_grad_back = dens_grad_back,
      stepsize = stepsize, cond_mat = cond_mat
    )
    
    # kinetic energies of momenta  ####
    current_momentum_dens <- momentumDens(state$momenta$noise_log_scale)
    proposed_momentum_dens <- momentumDens(new_momentum)
    # target densities ####
    current_dens <- noiseLogLik(
      squared_residuals = squared_residuals, noise_var = state$stuff$noise_var)
    proposed_dens <- noiseLogLik(
      squared_residuals = squared_residuals, noise_var = new_noise_var)
    # 
    # if(iter %/% 10 == iter / 10){
    #   par(mfrow = c(3, 2))
    #   plot(
    #     state$params$noise_log_scale, 
    #     moveForward(stepsize, cond_mat, dens_grad_back, 
    #                 new_noise_log_scale, momentum = new_momentum),
    #     xlab = "current noise log scale", ylab = "moving back proposed noise log scale"
    #   )
    #   abline(a=0, b=1)
    #   plot(
    #     state$params$noise_beta, 
    #     movePPCoeff(n_knots = hierarchical_model$noise$PP$n_knots, 
    #                 beta = new_noise_beta, 
    #                 log_scale = new_noise_log_scale, 
    #                 new_log_scale = moveForward(stepsize, cond_mat, dens_grad_back, 
    #                                                   new_noise_log_scale, momentum = new_momentum)),
    #     xlab = "current noise beta", ylab = "moving back proposed noise beta"
    #   )
    #   abline(a=0, b=1)
    #   testNoisePPLogLikGrad(
    #     noise_beta = state$params$noise_beta, 
    #     grad_to_test = grad_PP, 
    #     dens_to_test = current_dens, 
    #     hierarchical_model = hierarchical_model, noise_X = noise_X,
    #     vecchia_approx = vecchia_approx, state = state)
    #   testNoisePPLogLikGrad(
    #     noise_beta = new_noise_beta, 
    #     grad_to_test = grad_PP_back, 
    #     dens_to_test = proposed_dens, 
    #     hierarchical_model = hierarchical_model, noise_X = noise_X,
    #     vecchia_approx = vecchia_approx, state = state)
    # }
    
     # Metropolis ####
     ratio <- proposed_dens - current_dens + proposed_momentum_dens - current_momentum_dens
     if (!is.nan(ratio)) {
       if (log(runif(1)) < ratio) {
         # increasing stepsize if acceptance
         state$ker_var$noise_log_scale <- updateKernel(
           iter = iter, iter_start = iter_start,
           kernel_value = state$ker_var$noise_log_scale, mult = 3
         )
         # updating momenta
         state$momenta$noise_log_scale <- new_momentum
         # updating parameters
         state$params$noise_log_scale[] <- new_noise_log_scale
         state$params$noise_beta[] <- new_noise_beta
         # updating stuff
         state$stuff$noise_var <- new_noise_var
       }
     }
     # negating momentum to induce "slippery slope" behavior
     state$momenta$noise_log_scale <- -state$momenta$noise_log_scale
  }
  return(state)
}

# Testing function for the gradient of PP
testNoisePPLogLikGrad <- function(
    noise_beta, field, 
    grad_to_test ,
    dens_to_test ,
    hierarchical_model, noise_X,
    vecchia_approx, state){
  squared_residuals <- as.matrix(state$stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2
  empirical_grad = c()
    for(i in seq(nrow(noise_beta) - hierarchical_model$noise$PP$n_knots+1, nrow(noise_beta))){
      # moving noise beta 
      noise_beta_test <- noise_beta
      noise_beta_test[i] <- noise_beta_test[i]+1e-6
      # re-computing noise_var
      noise_var_grad_test <- exp(xPPMultRight(
        X = noise_X$X, PP = hierarchical_model$noise$PP,
        vecchia_approx = vecchia_approx, Y = noise_beta_test,
        permutate_PP_to_obs = TRUE
      ))
      # re-computing density
      test_dens <- noiseLogLik(noise_var_grad_test, squared_residuals)
      # empitical gradient
      empirical_grad = c(empirical_grad,  (test_dens - dens_to_test)*1e6)
    }
  plot(grad_to_test, empirical_grad,  main  = "PP sufficient gradient test")
  abline(a=0, b=1)
  abline(h=0)
  abline(v=0)
}
