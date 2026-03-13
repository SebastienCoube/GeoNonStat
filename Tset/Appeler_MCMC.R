library(GeoNonStat)
set.seed(2)
nlocs = 20000
observed_locs = cbind(runif(nlocs), runif(nlocs))
vecchia_approx = createVecchia(observed_locs, m = 6)

# fixed effects
X = as.data.frame(cbind(observed_locs[,1], rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs))))

#X_range = as.data.frame(observed_locs)

PP_range = createPP(vecchia_approx, knots = 16, matern_range = .3)
PP_noise = createPP(vecchia_approx, knots = 1000, matern_range = .05)


gns_simulator = createGnsSimulator(
  vecchia_approx = vecchia_approx,
  X = X, 
  noise_X = X, noise_PP = PP_noise, 
  range_X = NULL,range_PP = PP_range, 
  anisotropic = T
)


field_log_var=  0
# case with nice non-stationarity and little noise. 
range_log_scale = c(-1, -1)
noise_intercept = 1
noise_PP_log_var = 0
# case with no non-stationarity and little noise. 
# case with nice non-stationarity and crazy noise. 
# case with no non-stationarity and crazy noise. 
coeff_list = createGnsSimulatorParameters(gns_simulator, range_PP_log_var =range_log_scale, noise_intercept = noise_intercept, noise_PP_log_var = noise_PP_log_var, field_log_var = field_log_var)
coeff_list$range_X_coeff[1] = -4.5


fake_data = simulateGnsData(gns_simulator = gns_simulator, gns_params = coeff_list)
dev.off()


plotPointillistPainting(vecchia_approx$locs, fake_data$latent_field, cex= 1, pch= 15)

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


# #Run MCMC chain for 40 iterations
list2env(geo_non_stat, environment())
state = geo_non_stat$states$chain_1
#state$params$range_log_scale[1] = 0
#state$params$range_log_scale[2] = -3
num_threads = 10
n_iterations = 40
iter_start = 1
iter = 1
seed=  1
range_X = covariates$range_X

samples = runGeoNonStatMcmc(
 covariates = geo_non_stat$covariates, observed_field = geo_non_stat$observed_field, 
 hierarchical_model = geo_non_stat$hierarchical_model, vecchia_approx = geo_non_stat$vecchia_approx, 
 state = state, n_iterations = 400, num_threads = 8, iter_start = 1, seed = 1
)

params = samples$state$params
params_records = samples$params_records
iter = length(samples$params_records)


plotPointillistPainting(vecchia_approx$locs, params$field, main = "sampled", pch = 16, cex=  1)
plotPointillistPainting(vecchia_approx$locs, fake_data$latent_field, main = "true", pch = 16, cex=  1)
plotPointillistPainting(vecchia_approx$observed_locs, log(samples$state$stuff$noise_var), main = "sampled log noise var", pch = 16, cex=  1)
plotPointillistPainting(vecchia_approx$observed_locs, fake_data$log_noise_var_field, main = "true log noise var", pch = 16, cex=  1)

 par(mar = rep(2, 4))
 par(mfrow = c(3,3))
 for(i in seq(nrow(params$range_beta))){
   for(j in seq(ncol(params$range_beta))){
     plot(c(rbind(coeff_list$range_X_coeff, coeff_list$range_PP_coeff)[i,j], sapply(params_records, function(x)x$range_beta[i,j])[seq(iter)]), ylab = "", main  = paste(row.names(params$range_beta)[i], c("range", "aniso 1", "aniso 2")[j]))
     abline(h = rbind(coeff_list$range_X_coeff, coeff_list$range_PP_coeff)[i,j])
   }}
 
 range_log_scale_samples = rbind(range_log_scale, t(sapply(params_records, function(x)x$range_log_scale))[seq(iter-1),])
 plot(range_log_scale_samples[,1], main = "range log scale")
 abline(h = range_log_scale[1])
 plot(range_log_scale_samples[,2], main = "aniso log scale")
 abline(h = range_log_scale[2])
 
 noise_log_scale_samples = c(noise_PP_log_var, sapply(params_records, function(x)x$noise_log_scale)[seq(iter-1)])
 plot(noise_log_scale_samples, main = "noise log scale")
 abline(h = noise_PP_log_var)
 
 plot(c(noise_PP_log_var, sapply(params_records, function(x)x$noise_beta[1])[seq(iter-1)]), main = "noise intercept")
 
 
 field_log_var_samples = c(field_log_var, sapply(params_records, function(x)x$field_log_var)[seq(iter-1)])
 plot(field_log_var_samples, main = "field log var")
 abline(h = field_log_var)
 
 
# # # Current i.e. last  state of the MCMC chain
# # samples$state
# # # e.g. current utils 
# # samples$state$stuff
# # # e.g. current model parameters 
# # samples$state$params
# # # History of the samples
# # samples$params_records[[1]]$beta
# 
# # on veut arriver à ça... mais en version pas º(00)º : 
# 
# # parallel run
# t1 = Sys.time()
# mise_a_jour_mcmc_parallele = runParallelGeoNonStatMcmc(
#   object = geo_non_stat, n_chains_in_parallel = 3, 
#   n_threads_per_chain = 10, n_iterations = 30, seed = 1)
# print(Sys.time()-t1)
# # update of the state and the records
# for(chain_idx in seq(length(mygns$states))){
#   mygns$states[[chain_idx]] = mise_a_jour_mcmc_parallele[[chain_idx]][[1]]
#   mygns$records[[chain_idx]] = c(mygns$records[[chain_idx]], mise_a_jour_mcmc_parallele[[chain_idx]][[2]])
# }
# 
# 
# t1 = Sys.time()
# mise_a_jour_mcmc_parallele = runParallelGeoNonStatMcmc(
#   object = geo_non_stat, n_chains_in_parallel = 3, 
#   n_threads_per_chain = 3, n_iterations = 30, seed = 1)
# print(Sys.time()-t1)
# # update of the state and the records
# for(chain_idx in seq(length(mygns$states))){
#   mygns$states[[chain_idx]] = mise_a_jour_mcmc_parallele[[chain_idx]][[1]]
#   mygns$records[[chain_idx]] = c(mygns$records[[chain_idx]], mise_a_jour_mcmc_parallele[[chain_idx]][[2]])
# }
# 
# 
# t1 = Sys.time()
# mise_a_jour_mcmc_parallele = runParallelGeoNonStatMcmc(
#   object = geo_non_stat, n_chains_in_parallel = 3, 
#   n_threads_per_chain = 5, n_iterations = 30, seed = 1)
# print(Sys.time()-t1)
# # update of the state and the records
# for(chain_idx in seq(length(mygns$states))){
#   mygns$states[[chain_idx]] = mise_a_jour_mcmc_parallele[[chain_idx]][[1]]
#   mygns$records[[chain_idx]] = c(mygns$records[[chain_idx]], mise_a_jour_mcmc_parallele[[chain_idx]][[2]])
# }
# 
# 
# records = geo_non_stat$records
# 
# PrintMcmcDiags(geo_non_stat,burn_in = .2)
# TracePlots(geo_non_stat, .2)
# #TracePlots(geo_non_stat, .3, who = c("range_beta", "range_log_scale"))
# 
# estimation = estimate(geo_non_stat, .2)
# 
# elppd_train = elppdTrain(geo_non_stat, .2)
# 
# 
# new_locs = rbind(observed_locs[1,], as.matrix(expand.grid(seq(-.1, 1.1, .01), seq(-.1, 1.1, .01))))
# 