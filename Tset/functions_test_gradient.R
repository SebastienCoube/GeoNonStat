# Testing function for the sufficient gradient of range log scale
testGradientLogScaleS <- function(
    range_log_scale,
    range_beta,
    grad_to_test, dens_to_test, num_threads,
    range_X, hierarchical_model, vecchia_approx,
    state){
  empirical_grad = c()
  sparse_chol_grad_test = state$stuff$sparse_chol
  compressed_chol_grad_test = state$stuff$compressed_chol
  for(i in seq(length(range_log_scale))){
    # moving range log scale
    range_log_scale_grad_test <- range_log_scale 
    range_log_scale_grad_test[i]  <- range_log_scale_grad_test[i] + 1e-6
    # moving range beta 
    range_beta_test <- movePPCoeff(
      beta = range_beta, log_scale = range_log_scale, 
      new_log_scale = range_log_scale_grad_test, 
      n_knots = hierarchical_model$range$PP$n_knots)
    # re-computing Vecchia approx
    log_range <- computeLogRange(
      range_beta = range_beta_test,
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
    test_dens <-rangeLogLikS(
      field = state$params$field, 
      sparse_chol = sparse_chol_grad_test, 
      field_log_var = state$params$field_log_var
    )
    empirical_grad = c(empirical_grad,  (test_dens - dens_to_test)*1e6)
  }
  plot(grad_to_test, empirical_grad,  main  = "Range log scale sufficient gradient test")
  abline(a=0, b=1)
  abline(h=0)
  abline(v=0)
}

# Testing function for the ancillary gradient of range log scale
testGradientLogScaleA <- function(
    field,
    range_log_scale,
    range_beta,
    grad_to_test, dens_to_test, num_threads,
    range_X, hierarchical_model, vecchia_approx,
    state){
  empirical_grad = c()
  sparse_chol_grad_test = state$stuff$sparse_chol
  compressed_chol_grad_test = state$stuff$compressed_chol
  for(i in seq(length(range_log_scale))){
    # moving range log scale
    range_log_scale_grad_test <- range_log_scale 
    range_log_scale_grad_test[i]  <- range_log_scale_grad_test[i] + 1e-6
    # moving range beta 
    range_beta_test <- movePPCoeff(
      beta = range_beta, log_scale = range_log_scale, 
      new_log_scale = range_log_scale_grad_test, 
      n_knots = hierarchical_model$range$PP$n_knots)
    # re-computing Vecchia approx
    log_range <- computeLogRange(
      range_beta = range_beta_test,
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
    field_grad_test <- Matrix::solve(sparse_chol_grad_test, sparse_chol %*% field)
    test_dens <- rangeLogLikA(
      field = field_grad_test, lm_residuals = state$stuff$lm_residuals, 
      vecchia_approx = vecchia_approx, noise_var = state$stuff$noise_var
    )
    empirical_grad = c(empirical_grad,  (test_dens - dens_to_test)*1e6)
  }
  plot(grad_to_test, empirical_grad,  main  = "Range log scale sufficient gradient test")
  abline(a=0, b=1)
  abline(h=0)
  abline(v=0)
}

# Testing function for the sufficient gradient of PP
testrangePPLogLikGradS <- function(
    range_beta,
    grad_to_test, dens_to_test, num_threads,
    range_X, hierarchical_model, vecchia_approx,
    state){
  empirical_grad = c()
  sparse_chol_grad_test = state$stuff$sparse_chol
  compressed_chol_grad_test = state$stuff$compressed_chol
  for(j in seq_len(ncol(range_beta))){
    for(i in seq(nrow(range_beta) - hierarchical_model$range$PP$n_knots+1, nrow(range_beta))){
      # moving range beta 
      range_beta_test <- range_beta
      range_beta_test[i,j] <- range_beta_test[i,j]+1e-6
      # re-computing Vecchia approx
      log_range <- computeLogRange(
        range_beta = range_beta_test,
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
      test_dens <- rangeLogLikS(
        field = state$params$field, 
        sparse_chol = sparse_chol_grad_test, 
        field_log_var = state$params$field_log_var
      )
      empirical_grad = c(empirical_grad,  (test_dens - dens_to_test)*1e6)
    }}
  plot(grad_to_test, empirical_grad,  main  = "PP sufficient gradient test")
  abline(a=0, b=1)
  abline(h=0)
  abline(v=0)
}

# Testing function for the ancillary gradient of PP
testrangePPLogLikGradA <- function(
    range_beta, field, 
    sparse_chol,
    grad_to_test ,
    dens_to_test ,
    num_threads  ,
    hierarchical_model, range_X,
    vecchia_approx, state){
  empirical_grad = c()
  sparse_chol_grad_test = state$stuff$sparse_chol
  compressed_chol_grad_test = state$stuff$compressed_chol
  for(j in seq_len(ncol(range_beta))){
    for(i in seq(nrow(range_beta) - hierarchical_model$range$PP$n_knots+1, nrow(range_beta))){
      # moving range beta 
      range_beta_test <- range_beta
      range_beta_test[i,j] <- range_beta_test[i,j]+1e-6
      # re-computing Vecchia approx
      log_range <- computeLogRange(
        range_beta = range_beta_test,
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
      field_grad_test <- Matrix::solve(sparse_chol_grad_test, sparse_chol %*% field)
      test_dens <- rangeLogLikA(
        field = field_grad_test, lm_residuals = state$stuff$lm_residuals, 
        vecchia_approx = vecchia_approx, noise_var = state$stuff$noise_var
      )
      empirical_grad = c(empirical_grad,  (test_dens - dens_to_test)*1e6)
    }}
  plot(grad_to_test, empirical_grad,  main  = "PP sufficient gradient test")
  abline(a=0, b=1)
  abline(h=0)
  abline(v=0)
}



# Testing function for the ancillary gradient of range log scale
testGradientLogScaleA <- function(
    field,
    range_log_scale,
    range_beta,
    grad_to_test, dens_to_test, num_threads,
    range_X, hierarchical_model, vecchia_approx,
    state){
  empirical_grad = c()
  sparse_chol_grad_test = state$stuff$sparse_chol
  compressed_chol_grad_test = state$stuff$compressed_chol
  for(i in seq(length(range_log_scale))){
    # moving range log scale
    range_log_scale_grad_test <- range_log_scale 
    range_log_scale_grad_test[i]  <- range_log_scale_grad_test[i] + 1e-6
    # moving range beta 
    range_beta_test <- movePPCoeff(
      beta = range_beta, log_scale = range_log_scale, 
      new_log_scale = range_log_scale_grad_test, 
      n_knots = hierarchical_model$range$PP$n_knots)
    # re-computing Vecchia approx
    log_range <- computeLogRange(
      range_beta = range_beta_test,
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
    field_grad_test <- Matrix::solve(sparse_chol_grad_test, sparse_chol %*% field)
    test_dens <- rangeLogLikA(
      field = field_grad_test, vecchia_approx = vecchia_approx, 
      lm_residuals = state$stuff$lm_residuals, noise_var = state$stuff$noise_var
    )
    empirical_grad = c(empirical_grad,  (test_dens - dens_to_test)*1e6)
  }
  plot(grad_to_test, empirical_grad,  main  = "Range log scale sufficient gradient test")
  abline(a=0, b=1)
  abline(h=0)
  abline(v=0)
}
