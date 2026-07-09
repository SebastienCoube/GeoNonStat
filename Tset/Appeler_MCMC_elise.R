library(GeoNonStat)
set.seed(2)
nobs = 50000
nlocs = 2500
observed_locs = cbind(runif(nlocs), runif(nlocs))[c(seq(nlocs), sample(seq(nlocs), nobs - nlocs, T)),]

vecchia_approx = createVecchia(observed_locs, m = 6, round_locs = .005)

# fixed effects
X = as.data.frame(cbind(vecchia_approx$locs[vecchia_approx$locs_match,1], vecchia_approx$observed_locs[,1] + rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs))))

#X_range = as.data.frame(observed_locs)

PP_range = createPP(vecchia_approx, knots = 10, matern_range = .5)
PP_noise = createPP(vecchia_approx, knots = 100, matern_range = .15)


gns_simulator = createGnsSimulator(
  vecchia_approx = vecchia_approx,
  X = X, 
  noise_X = X, noise_PP = PP_noise, 
  range_X = NULL,range_PP = PP_range, 
  anisotropic = T
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
coeff_list$range_X_coeff[1,1] = -5
coeff_list$noise_X_coeff[-1] = .1*rnorm(4)
set.seed(2)
fake_data = simulateGnsData(gns_simulator = gns_simulator, gns_params = coeff_list)
plotPointillistPainting(vecchia_approx$locs, fake_data$latent_field, cex= .5, pch= 15)

plotPointillistPainting(vecchia_approx$observed_locs, fake_data$observed_field, cex= 1, pch= 15)

plotPointillistPainting(vecchia_approx$observed_locs, fake_data$noise, cex= 1, pch= 15)

geo_non_stat = GeoNonStat(
  vecchia_approx = vecchia_approx, 
  observed_field = fake_data$observed_field, X = X, 
  matern_smoothness = 1.5, anisotropic = T, 
  n_chains = 3, 
  noise_X = NULL, range_X = NULL,
  noise_PP = PP_noise,
  range_PP = PP_range
)

# parallel run
future::plan(strategy = "multisession", workers = 3)
geo_non_stat = GeoNonStatMcmc(
  object = geo_non_stat, n_chains_in_parallel = 3, 
  n_threads_per_chain = 2, n_iterations = 200, seed = 1)

# plotting and diagnostics
tracePlots(geo_non_stat, keep = "range_beta", burn_in = .0)
tracePlots(geo_non_stat, keep = "field_log_var", burn_in = .0)
tracePlots(geo_non_stat, keep = "field_log_var", burn_in = .2)
tracePlots(geo_non_stat, keep = "noise_beta", burn_in = 0)
tracePlots(geo_non_stat, keep = "range_log_scale", burn_in = .0)
MCMC_diags = mcmcDiags(geo_non_stat, burn_in = .2)
MCMC_diags$diags["field_log_var. ",] 

# prediction
new_locs = as.matrix(expand.grid(seq(0, 1, .01), seq(0, 1, .01))) 
new_locs = rbind(new_locs, new_locs[rev(seq_len(nrow(new_locs))),])
extended_vecchia_approx = extendVecchia(vecchia_approx, new_locs)

# prediction of 1 log - noise variance realization
extended_PP_noise = extendPP(
  PP = PP_noise, 
  extended_vecchia_approx = extended_vecchia_approx)
noise_beta = geo_non_stat$states$chain_2$params$noise_beta

predicted_noise = predict1Noise(
  noise_beta, 
  extended_PP_noise = extended_PP_noise, 
  extended_vecchia_approx, 
  new_noise_X = matrix(1, nrow(new_locs)))
par(mfrow = c(2,3))
plotPointillistPainting(vecchia_approx$observed_locs, fake_data$log_noise_var_field, main = "true log noise var", pch = 16, cex=  1)
plotPointillistPainting(vecchia_approx$observed_locs, log(geo_non_stat$states$chain_2$stuff$noise_var), main = "sampled log noise var", pch = 16, cex=  1)
plotPointillistPainting(extended_vecchia_approx$new_observed_locs, predicted_noise, main = "predicted log noise variance")
# prediction of 1 realization of the latent field + range field

range_beta = geo_non_stat$states$chain_2$params$range_beta
extended_PP_range = extendPP(
  PP = PP_range, 
  extended_vecchia_approx = extended_vecchia_approx)
complete_range_X_locs = list(X_locs = matrix(1, extended_vecchia_approx$n_locs))
field = geo_non_stat$states$chain_2$params$field

predicted_field = predict1Field(
  range_beta = range_beta, field = field, 
  extended_PP_range = extended_PP_range, 
  extended_vecchia_approx, 
  complete_range_X_locs, geo_non_stat$hierarchical_model, 
  num_threads = 8)

plotPointillistPainting(vecchia_approx$locs, fake_data$latent_field, main = "true latent field", pch = 16, cex=  1)
plotPointillistPainting(vecchia_approx$locs, field, main = "sampled latent field", pch = 16, cex=  1)
plotPointillistPainting(
  extended_vecchia_approx$new_locs, 
  predicted_field$field, 
  main = "predicted latent field")

plotPointillistPainting(vecchia_approx$locs, fake_data$log_range_field[,1], main = "true log range", pch = 16, cex=  1)
plotPointillistPainting(extended_vecchia_approx$new_locs, predicted_field$log_range[,1], main = "predicted log range")
plotPointillistPainting(vecchia_approx$locs, fake_data$log_range_field[,2], main = "true log aniso 1", pch = 16, cex=  1)
plotPointillistPainting(extended_vecchia_approx$new_locs, predicted_field$log_range[,2], main = "predicted log aniso 1")
plotPointillistPainting(vecchia_approx$locs, fake_data$log_range_field[,3], main = "true log aniso 2", pch = 16, cex=  1)
plotPointillistPainting(extended_vecchia_approx$new_locs, predicted_field$log_range[,3], main = "predicted log aniso 2")


predict.GeoNonStat(
  object = geo_non_stat, 
  new_locs = new_locs, 
  new_X = NULL, 
  new_noise_X = matrix(1, nrow(new_locs)), 
  new_range_X = complete_range_X_locs, 
  num_threads=2)
  
