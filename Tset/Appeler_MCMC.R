library(GeoNonStat)
set.seed(2)
nobs = 10000
nlocs = 4000
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
  noise_X = X, noise_PP = PP_noise, 
  range_X = NULL,range_PP = PP_range, 
  anisotropic = T
)


field_log_var=  0
# case with nice non-stationarity and little noise. 
range_log_scale = c(0, 0)
noise_intercept = -2
noise_PP_log_var = -.5
# case with no non-stationarity and little noise. 
# case with nice non-stationarity and crazy noise. 
# case with no non-stationarity and crazy noise. 
coeff_list = createGnsSimulatorParameters(
  gns_simulator, range_PP_log_var =range_log_scale, 
  noise_intercept = noise_intercept, noise_PP_log_var = noise_PP_log_var, 
  field_log_var = field_log_var)
coeff_list$range_X_coeff[1,1] = -4.5
coeff_list$noise_X_coeff[-1] = .1*rnorm(4)
set.seed(2)
fake_data = simulateGnsData(gns_simulator = gns_simulator, gns_params = coeff_list)
plotPointillistPainting(vecchia_approx$locs, fake_data$latent_field, cex= .5, pch= 15)

plotPointillistPainting(vecchia_approx$observed_locs, fake_data$observed_field, cex= 1, pch= 15)

plotPointillistPainting(vecchia_approx$observed_locs, fake_data$log_noise_var_field, cex= 1, pch= 15)

geo_non_stat = GeoNonStat(
  vecchia_approx = vecchia_approx, 
  observed_field = fake_data$observed_field, X = X, 
  matern_smoothness = 1.5, anisotropic = T, 
  n_chains = 3, 
  noise_X = NULL, range_X = NULL,
  noise_PP = PP_noise,
  range_PP = PP_range
)


summary(geo_non_stat)

# ## # #Run MCMC chain for 40 iterations
#  list2env(geo_non_stat, environment())
#  state = geo_non_stat$states$chain_1
#  #state$params$range_log_scale[1] = 0
#  #state$params$range_log_scale[2] = -3
#  num_threads = 10
#  n_iterations = 300
#  iter_start = 1
#  iter = 1
#  seed=  1
#  range_X = covariates$range_X
#  
#  samples = runOneChainMcmc(
#   covariates = geo_non_stat$covariates, observed_field = geo_non_stat$observed_field, 
#   hierarchical_model = geo_non_stat$hierarchical_model, vecchia_approx = geo_non_stat$vecchia_approx, 
#   state = state, n_iterations = 59, num_threads = 8, iter_start = 1, seed = 1
#  )
#  
#  state = samples$state
#  samples = runOneChainMcmc(
#   covariates = geo_non_stat$covariates, observed_field = geo_non_stat$observed_field, 
#   hierarchical_model = geo_non_stat$hierarchical_model, vecchia_approx = geo_non_stat$vecchia_approx, 
#   state = state, n_iterations = 200, num_threads = 8, iter_start = 60, seed = 1
#  )
#  
#  state = samples$state
#  params = samples$state$params
#  params_records = samples$params_records
#  iter = length(samples$params_records)
#
# 
# plotPointillistPainting(vecchia_approx$locs, params$field, main = "sampled", pch = 16, cex=  1)
# plotPointillistPainting(vecchia_approx$locs, fake_data$latent_field, main = "true", pch = 16, cex=  1)
# plotPointillistPainting(vecchia_approx$observed_locs, log(samples$state$stuff$noise_var), main = "sampled log noise var", pch = 16, cex=  1)
# plotPointillistPainting(vecchia_approx$observed_locs, fake_data$log_noise_var_field, main = "true log noise var", pch = 16, cex=  1)
#  
# plot(params$field, fake_data$latent_field)
# abline(a=0, b=1)
#  
# par(mar = rep(2, 4))
#  par(mfrow = c(3,3))
#  for(i in seq(nrow(params$range_beta))){
#    for(j in seq(ncol(params$range_beta))){
#      plot(c(rbind(coeff_list$range_X_coeff, coeff_list$range_PP_coeff)[i,j], sapply(params_records, function(x)x$range_beta[i,j])[seq(iter)]), ylab = "", main  = paste(row.names(params$range_beta)[i], c("range", "aniso 1", "aniso 2")[j]))
#      abline(h = rbind(coeff_list$range_X_coeff, coeff_list$range_PP_coeff)[i,j])
#    }}
#  
#  
#  range_log_scale_samples = rbind(range_log_scale, t(sapply(params_records, function(x)x$range_log_scale))[seq(iter-1),])
#  plot(range_log_scale_samples[,1], main = "range log scale")
#  abline(h = range_log_scale[1])
#  plot(range_log_scale_samples[,2], main = "aniso log scale")
#  abline(h = range_log_scale[2])
#  
#  noise_log_scale_samples = c(noise_PP_log_var, sapply(params_records, function(x)x$noise_log_scale)[seq(iter/10, iter-1)])
#  plot(noise_log_scale_samples, main = "noise log scale")
#  abline(h = noise_PP_log_var)
#  
#  plot(c(coeff_list$noise_X_coeff[1], sapply(params_records, function(x)x$noise_beta[1])[seq(iter/10, iter-1)]), main = "noise intercept")
#  abline(h = coeff_list$noise_X_coeff[1])
#  for(i in c(1,2,20, 50, 100, 250, 500, 800)){
#  plot(c(coeff_list$noise_PP_coeff[i], sapply(params_records, function(x)x$noise_beta[i+1])[seq(iter-1)]), main = paste("noise PP", i))
#  abline(h = coeff_list$noise_PP_coeff[i])
#  abline(h = 0)
#  }
#  
#  
#  field_log_var_samples = c(field_log_var, sapply(params_records, function(x)x$field_log_var)[seq(iter/10, iter-1)])
#  plot(field_log_var_samples, main = "field log var")
#  abline(h = field_log_var)

# parallel run
future::plan(strategy = "multisession", workers = 3)
geo_non_stat = GeoNonStatMcmc(
  object = geo_non_stat, n_chains_in_parallel = 3, 
  n_threads_per_chain = 7, n_iterations = 600, seed = 1)

print("DOOOONE")

# plotting and diagnostics
tracePlots(geo_non_stat, keep = "beta", burn_in = .3)
tracePlots(geo_non_stat, keep = "range_beta", burn_in = .3)
tracePlots(geo_non_stat, keep = "field_log_var", burn_in = .3)
tracePlots(geo_non_stat, keep = "noise_beta", burn_in = .3)
tracePlots(geo_non_stat, keep = "range_log_scale", burn_in = .3)
tracePlots(geo_non_stat, keep = "noise_log_scale", burn_in = .3)
MCMC_diags = printMcmcDiags(geo_non_stat, burn_in = .3)
MCMC_diags$diags["field_log_var. ",] 

# prediction
new_locs = as.matrix(expand.grid(seq(0, 1, .01), seq(0, 1, .01))) 
new_locs = rbind(new_locs, new_locs[rev(seq_len(nrow(new_locs))),])
extended_vecchia_approx = extendVecchia(vecchia_approx, new_locs)

# prediction of 1 log - noise variance realization
extended_PP_noise = extendPP(
  PP = PP_noise, 
  extended_vecchia_approx = extended_vecchia_approx)
noise_beta = geo_non_stat$states$chain_1$params$noise_beta

predicted_noise = predict1Noise(
  noise_beta, 
  extended_PP_noise = extended_PP_noise, 
  extended_vecchia_approx, 
  new_noise_X = matrix(1, nrow(new_locs)))

par(mfrow = c(2,3))


plotPointillistPainting(
  vecchia_approx$observed_locs, fake_data$log_noise_var_field, 
  main = "true log noise var", pch = 16, cex=  1)

plotPointillistPainting(
  vecchia_approx$observed_locs, 
  log(geo_non_stat$states$chain_1$stuff$noise_var), 
  main = "sampled log noise var", pch = 16, cex=  1)

plotPointillistPainting(
  extended_vecchia_approx$new_locs, predicted_noise, 
  main = "predicted log noise variance")


# prediction of 1 realization of the latent field + range field


range_beta = geo_non_stat$states$chain_1$params$range_beta
extended_PP_range = extendPP(
  PP = PP_range, 
  extended_vecchia_approx = extended_vecchia_approx)
complete_range_X_locs = list(X_locs = matrix(1, extended_vecchia_approx$n_locs))
field = geo_non_stat$states$chain_1$params$field

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

