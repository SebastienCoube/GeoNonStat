plot(
  round(vecchia(
    log_range = matrix(log(matern_range), 1, vecchia_approx$n_locs), 
    locs = vecchia_approx$t_locs, 
    NNarray = vecchia_approx$NNarray, 
    num_threads = 5, smoothness = 1.5, compute_derivative = F
  )[,1,], 4), 
  round(t(GpGp::vecchia_Linv(
    covparms = c(1, (matern_range), 1.5, .000001), 
    locs = vecchia_approx$locs, 
    NNarray = t(vecchia_approx$NNarray), 
    covfun_name = "matern_isotropic"
  )), 4)
)