observed_locs = cbind(runif(400000), 1)
X = cbind(observed_locs[,1], rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)))

NNarray = GpGp::find_ordered_nn(observed_locs, 12)
Linv = GpGp::vecchia_Linv(c(1, .0005, .0001), "matern15_isotropic", locs = observed_locs, NNarray)
w = GpGp::fast_Gp_sim_Linv(Linv, NNarray)
reg_coeffs = rnorm(ncol(X)+1)

observed_field = c(cbind(1, X) %*% reg_coeffs + w + rnorm(nrow(observed_locs)))

plot(observed_locs[,1], observed_field)
points(observed_locs[,1], w, col = 2, cex = .5, pch = 16)


vecchia_approx = createVecchia(observed_locs)
PP = createPP(vecchia_approx)

mygns = GeoNonStat(
  vecchia_approx = vecchia_approx, 
  observed_field = c(observed_field), X = as.data.frame(X), 
  matern_smoothness = 1.5, anisotropic = F, 
  n_chains = 3, 
  noise_X = NULL, range_X = NULL, scale_X = NULL, noise_PP = NULL, range_PP = NULL, scale_PP = NULL, seed = 1
)

samples1 = MCMC(
  covariates = mygns$covariates, 
  observed_field = mygns$observed_field, 
  hierarchical_model =  mygns$hierarchical_model, 
  vecchia_approx = mygns$vecchia_approx, 
  n_iterations_update = 100, 
  state = mygns$states$chain_1,
  num_threads = 5, 
  iter_start = 1, 
  seed = 1
    )
samples2 = MCMC(
  covariates = mygns$covariates, 
  observed_field = mygns$observed_field, 
  hierarchical_model =  mygns$hierarchical_model, 
  vecchia_approx = mygns$vecchia_approx, 
  n_iterations_update = 100, 
  state = mygns$states$chain_2,
  num_threads = 5, 
  iter_start = 1, 
  seed = 1
    )
samples3 = MCMC(
  covariates = mygns$covariates, 
  observed_field = mygns$observed_field, 
  hierarchical_model =  mygns$hierarchical_model, 
  vecchia_approx = mygns$vecchia_approx, 
  n_iterations_update = 100, 
  state = mygns$states$chain_3,
  num_threads = 5, 
  iter_start = 1, 
  seed = 1
    )

samples3$params_records[[1]]$scale_beta

mygns.run(100)

