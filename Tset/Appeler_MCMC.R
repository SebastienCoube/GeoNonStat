
set.seed(2)
nlocs = 5000
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


field_log_var=  2
# case with nice non-stationarity and little noise. 
range_log_scale = c(-1,-1)
noise_intercept = 1
noise_PP_log_var = -1
# case with no non-stationarity and little noise. 
# case with nice non-stationarity and crazy noise. 
# case with no non-stationarity and crazy noise. 
coeff_list = createGnsSimulatorParameters(gns_simulator, range_PP_log_var =range_log_scale, noise_intercept = noise_intercept, noise_PP_log_var = noise_PP_log_var, field_log_var = field_log_var)
coeff_list$range_X_coeff[1] = -3.5


fake_data = simulateGnsData(gns_simulator = gns_simulator, gns_params = coeff_list)
dev.off()


plot_pointillist_painting(vecchia_approx$locs, fake_data$latent_field, cex= 0.5, pch= 15)

plot_pointillist_painting(vecchia_approx$observed_locs, fake_data$observed_field)
plot_pointillist_painting(vecchia_approx$observed_locs, fake_data$noise)
plot_pointillist_painting(vecchia_approx$observed_locs, fake_data$log_noise_var_field)


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
# samples = MessyMessyMcmc(
#   covariates = mygns$covariates, observed_field = mygns$observed_field, 
#   hierarchical_model = mygns$hierarchical_model, vecchia_approx = mygns$vecchia_approx, 
#   state = mygns$states$chain_1, n_iterations_update = 40, num_threads = 10, iter_start = 1, seed = 1
# )
# 
# # Current i.e. last  state of the MCMC chain
# samples$state
# # e.g. current utils 
# samples$state$stuff
# # e.g. current model parameters 
# samples$state$params
# # History of the samples
# samples$params_records[[1]]$beta

# on veut arriver à ça... mais en version pas º(00)º : 

# parallel run
t1 = Sys.time()
mise_a_jour_mcmc_parallele = run_parallel_version_goret(
  object = geo_non_stat, n_chains_in_parallel = 3, 
  n_threads_per_chain = 10, n_iterations = 30, seed = 1)
print(Sys.time()-t1)
# update of the state and the records
for(chain_idx in seq(length(mygns$states))){
  mygns$states[[chain_idx]] = mise_a_jour_mcmc_parallele[[chain_idx]][[1]]
  mygns$records[[chain_idx]] = c(mygns$records[[chain_idx]], mise_a_jour_mcmc_parallele[[chain_idx]][[2]])
}


t1 = Sys.time()
mise_a_jour_mcmc_parallele = run_parallel_version_goret(
  object = geo_non_stat, n_chains_in_parallel = 3, 
  n_threads_per_chain = 3, n_iterations = 30, seed = 1)
print(Sys.time()-t1)
# update of the state and the records
for(chain_idx in seq(length(mygns$states))){
  mygns$states[[chain_idx]] = mise_a_jour_mcmc_parallele[[chain_idx]][[1]]
  mygns$records[[chain_idx]] = c(mygns$records[[chain_idx]], mise_a_jour_mcmc_parallele[[chain_idx]][[2]])
}


t1 = Sys.time()
mise_a_jour_mcmc_parallele = run_parallel_version_goret(
  object = geo_non_stat, n_chains_in_parallel = 3, 
  n_threads_per_chain = 5, n_iterations = 30, seed = 1)
print(Sys.time()-t1)
# update of the state and the records
for(chain_idx in seq(length(mygns$states))){
  mygns$states[[chain_idx]] = mise_a_jour_mcmc_parallele[[chain_idx]][[1]]
  mygns$records[[chain_idx]] = c(mygns$records[[chain_idx]], mise_a_jour_mcmc_parallele[[chain_idx]][[2]])
}


records = geo_non_stat$records

PrintMcmcDiags(geo_non_stat,burn_in = .2)
TracePlots(geo_non_stat, .2)
#TracePlots(geo_non_stat, .3, who = c("range_beta", "range_log_scale"))

estimation = Estimate(geo_non_stat, .2)

elppd_train = ElppdTrain(geo_non_stat, .2)


new_locs = rbind(observed_locs[1,], as.matrix(expand.grid(seq(-.1, 1.1, .01), seq(-.1, 1.1, .01))))
