library(GeoNonStat)
set.seed(2)
nobs = 200000
nlocs = 200000
observed_locs = cbind(runif(nlocs), runif(nlocs))[c(seq(nlocs), sample(seq(nlocs), nobs - nlocs, T)),]

vecchia_approx = createVecchia(observed_locs, m = 6, round_locs = 0)

# fixed effects
X = as.data.frame(cbind(vecchia_approx$locs[vecchia_approx$locs_match,1], vecchia_approx$observed_locs[,1] + rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs))))

#X_range = as.data.frame(observed_locs)

PP_range = createPP(vecchia_approx, knots = 10, matern_range = .5)
PP_noise = createPP(vecchia_approx, knots = 100, matern_range = .15)


gns_simulator = createGnsSimulator(
  vecchia_approx = vecchia_approx,
  X = X, 
  noise_X = X, noise_PP = PP_noise
)


field_log_var=  0
# case with nice non-stationarity and little noise. 
range_log_scale = c(0, 0)
noise_intercept = 2
noise_PP_log_var = -.5
# case with no non-stationarity and little noise. 
# case with nice non-stationarity and crazy noise. 
# case with no non-stationarity and crazy noise. 
coeff_list = createGnsSimulatorParameters(
  gns_simulator, range_PP_log_var =range_log_scale, 
  noise_intercept = noise_intercept, noise_PP_log_var = noise_PP_log_var, 
  field_log_var = field_log_var)
coeff_list$range_X_coeff[1,1] = -3
coeff_list$noise_X_coeff[-1] = .1*rnorm(4)
set.seed(2)
fake_data = simulateGnsData(gns_simulator = gns_simulator, gns_params = coeff_list)
plotPointillistPainting(vecchia_approx$locs, fake_data$latent_field, cex= 1, pch= 15)



vecchia_approx_0.01 = createVecchia(observed_locs, m = 6, round_locs = 0.01)
PP_noise_0.01 = createPP(vecchia_approx_0.01, knots = 100, matern_range = .15)
X = as.data.frame(cbind(
  vecchia_approx_0.01$locs[vecchia_approx_0.01$locs_match,1], 
  vecchia_approx_0.01$observed_locs[,1] + rnorm(nrow(observed_locs)), 
  rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs))))
geo_non_stat_0.01 = GeoNonStat(
  vecchia_approx = vecchia_approx_0.01, 
  observed_field = fake_data$observed_field, X = X, 
  matern_smoothness = 1.5,
  n_chains = 3, 
  noise_X = X, range_X = NULL,
  noise_PP = PP_noise_0.01
)

future::plan(strategy = "multisession", workers = 3)
geo_non_stat_0.01 = GeoNonStatMcmc(
  object = geo_non_stat_0.01, n_chains_in_parallel = 3, 
  n_threads_per_chain = 7, n_iterations = 1000, seed = 1)
tracePlots(geo_non_stat_0.01, keep = "range_beta", burn_in = .2)
tracePlots(geo_non_stat_0.01, keep = "field_log_var", burn_in = .2)


vecchia_approx_0.005 = createVecchia(observed_locs, m = 6, round_locs = 0.005)
PP_noise_0.005 = createPP(vecchia_approx_0.005, knots = 100, matern_range = .15)
X = as.data.frame(cbind(vecchia_approx_0.005$locs[vecchia_approx_0.005$locs_match,1], vecchia_approx_0.005$observed_locs[,1] + rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs))))
geo_non_stat_0.005 = GeoNonStat(
  vecchia_approx = vecchia_approx_0.005, 
  observed_field = fake_data$observed_field, X = X, 
  matern_smoothness = 1.5,
  n_chains = 3, 
  noise_X = NULL, range_X = NULL,
  noise_PP = PP_noise_0.005
)
future::plan(strategy = "multisession", workers = 3)
geo_non_stat_0.005 = GeoNonStatMcmc(
  object = geo_non_stat_0.005, n_chains_in_parallel = 3, 
  n_threads_per_chain = 7, n_iterations = 2000, seed = 1)
tracePlots(geo_non_stat_0.005, keep = "range_beta", burn_in = .2)
tracePlots(geo_non_stat_0.005, keep = "field_log_var", burn_in = .2)



vecchia_approx_0.001 = createVecchia(observed_locs, m = 6, round_locs = 0.001)
PP_noise_0.001 = createPP(vecchia_approx_0.001, knots = 100, matern_range = .15)
X = as.data.frame(cbind(vecchia_approx_0.001$locs[vecchia_approx_0.001$locs_match,1], vecchia_approx_0.001$observed_locs[,1] + rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs))))
geo_non_stat_0.001 = GeoNonStat(
  vecchia_approx = vecchia_approx_0.001, 
  observed_field = fake_data$observed_field, X = X, 
  matern_smoothness = 1.5,
  n_chains = 3, 
  noise_X = NULL, range_X = NULL,
  noise_PP = PP_noise_0.001
)
future::plan(strategy = "multisession", workers = 3)
geo_non_stat_0.001 = GeoNonStatMcmc(
  object = geo_non_stat_0.001, n_chains_in_parallel = 3, 
  n_threads_per_chain = 7, n_iterations = 3000, seed = 1)
tracePlots(geo_non_stat_0.001, keep = "range_beta", burn_in = .2)
tracePlots(geo_non_stat_0.001, keep = "field_log_var", burn_in = .2)
