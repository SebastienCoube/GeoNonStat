

densGradNoise <- function(noise_beta, noise_var, 
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

densNoise <- function(noise_beta, noise_var,
                      hm_noise, noise_log_scale, squared_residuals){
  return(
    betaPriorLogDens(
      beta = noise_beta,
      n_PP = hm_noise$PP$n_knots,
      beta0_mean = hm_noise$beta0_mean,
      beta0_var = hm_noise$beta0_sd^2,
      log_scale = noise_log_scale
    ) # normal prior
    - .5 * sum(log(noise_var)) # det
    - .5 * sum(squared_residuals / noise_var) # observations
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
    
    test_dens <- densNoise(
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




updateNoiseBetaFisher <- function(state, hm_noise, noise_X, vecchia_approx, iter, iter_start) {
  stepsize <- exp(state$ker_var$noise_beta)
  # residuals 
  squared_residuals <- as.matrix(state$stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2
  # initial gradient 
  dens_grad <- densGradNoise(
    noise_var = state$stuff$noise_var, noise_beta = state$params$noise_beta, 
    noise_X = noise_X, 
    vecchia_approx = vecchia_approx, 
    hm_noise = hm_noise, 
    squared_residuals = squared_residuals, 
    noise_log_scale = state$params$noise_log_scale)
  # conditioning matrix update
  cond_mat <- condMat(
    empirical_fisher = state$stuff$noise_beta_empirical_fisher, 
    prior_fisher = noise_X$crossprod_X, iter, iter_start)
  grad_record_for_Fisher <- list()
  # pre-allocation
  new_noise_var <- state$stuff$noise_var
  
  for (mala_rep in seq(20)) {
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
    new_dens_grad <- densGradNoise(
      noise_var = new_noise_var, noise_beta = new_noise_beta, 
      noise_X = noise_X, vecchia_approx = vecchia_approx, hm_noise = hm_noise, 
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
    current_dens <- densNoise(
      noise_beta = state$params$noise_beta, noise_var = state$stuff$noise_var, 
      squared_residuals = squared_residuals, 
      noise_log_scale = state$params$noise_log_scale,
      hm_noise = hm_noise)
    proposed_dens <- densNoise(
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
  
  # updating conditioning matrix ####
  if (iter + iter_start > 50) {
    state$stuff$noise_beta_empirical_fisher[] <-
      updateFisher(empirical_fisher = state$stuff$noise_beta_empirical_fisher, 
                   dens_grad = do.call(cbind, grad_record_for_Fisher), 
                   iter = iter, iter_start = iter_start)
      
  }
  
  
  
  return(list("state" = state, "squared_residuals" = squared_residuals))
}
