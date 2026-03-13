# devtools::load_all("GeoNonStat")
devtools::load_all()
set.seed(2)
nlocs = 30000
observed_locs = cbind(runif(nlocs), runif(nlocs))
vecchia_approx = createVecchia(observed_locs, m = 6)

# fixed effects
X = as.data.frame(cbind(observed_locs[,1], rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs))))

#X_range = as.data.frame(observed_locs)

PP_range = createPP(vecchia_approx, knots = 16, matern_range = .3)
PP_noise = createPP(vecchia_approx, knots = 30, matern_range = .25)

gns_simulator = createGnsSimulator(
  vecchia_approx = vecchia_approx,
  X = X, 
  noise_X = X, noise_PP = PP_noise, 
  range_X = NULL,range_PP = PP_range, 
  anisotropic = T
)

# case with nice non-stationarity and little noise. 
coeff_list = createGnsSimulatorParameters(
  gns_simulator, 
  range_PP_log_var =c(-1,-1), 
  noise_intercept = 1, 
  noise_PP_log_var = 1, 
  field_log_var = 2)


# case with no non-stationarity and little noise. 
# case with nice non-stationarity and crazy noise. 
# case with no non-stationarity and crazy noise. 

## coeff_list$range_X_coeff[1] = -5

fake_data = simulateGnsData(
  gns_simulator = gns_simulator, 
  gns_params = coeff_list)

gnsDemo = GeoNonStat(
  vecchia_approx = vecchia_approx, 
  observed_field = fake_data$observed_field, X = X, 
  matern_smoothness = 1.5, anisotropic = T, 
  n_chains = 3, 
  noise_X = NULL, range_X = NULL,
  noise_PP = PP_noise,
  range_PP = PP_range
)

future::plan(future::multisession, workers=3)
set.seed(123)
processedGnsDemo = runParallelGeoNonStatMcmc(
  object = gnsDemo, n_chains_in_parallel = 1,
  n_threads_per_chain = 6, n_iterations = 15, seed = 1)

set.seed(123)
mise_a_jour_mcmc_parallele2 = runParallelGeoNonStatMcmc(
  object = gnsDemo, n_chains_in_parallel = 1,
  n_threads_per_chain = 6, n_iterations = 15, seed = 1)

identical(mise_a_jour_mcmc_parallele1, mise_a_jour_mcmc_parallele2)

set.seed(123)
mise_a_jour_mcmc_parallele3 = runParallelGeoNonStatMcmc(
  object = gnsDemo, n_chains_in_parallel = 3,
  n_threads_per_chain = 6, n_iterations = 15, seed = 1)

set.seed(123)
mise_a_jour_mcmc_parallele4 = runParallelGeoNonStatMcmc(
  object = gnsDemo, n_chains_in_parallel = 3,
  n_threads_per_chain = 6, n_iterations = 15, seed = 1)

identical(mise_a_jour_mcmc_parallele3, mise_a_jour_mcmc_parallele4)



# parallel run
# profvis::profvis({
t <- Sys.time()
mise_a_jour_mcmc_parallele = runParallelGeoNonStatMcmc(
  object = gnsDemo, n_chains_in_parallel = 1,
  n_threads_per_chain = 6, n_iterations = 15, seed = 1)
print(Sys.time() - t)
# })
