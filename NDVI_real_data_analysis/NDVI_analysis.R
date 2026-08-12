load(file ="NDVI_real_data_analysis/data_cleaned_small_expanded.RData")

train_test_split <- createSplitData(
  data = list(y = data_cleaned_small$NDVI), 
  locs = cbind(data_cleaned_small$scaled_x, data_cleaned_small$scaled_y), 
  n_clust = 1000, 
  prop_test = c(.15, .02), 
  round_locs = .004)

plotPointillistPainting(
  train_test_split$train$locs, 
  train_test_split$train$y)

# creating model
vecchia_approx <- createVecchia(train_test_split$train$locs, m = 6)
range_PP <- createPP(vecchia_approx, .16)
noise_PP <- createPP(vecchia_approx, .08)
geo_non_stat_full_monty <- GeoNonStat(
  vecchia_approx = vecchia_approx, observed_field = train_test_split$train$y,
  noise_PP = noise_PP, range_PP = range_PP, n_chains = 3, anisotropic = T)


 run_time = Sys.time()
 future::plan(strategy = "multisession", workers = 3)
 geo_non_stat_full_monty <- GeoNonStatMcmc(
   geo_non_stat_full_monty, n_iterations = 600, 
   n_chains_in_parallel = 3, n_threads_per_chain = 8)
 run_time = run_time - Sys.time()
 
 
tracePlots(geo_non_stat_full_monty, keep = "beta", burn_in = 0) 
tracePlots(geo_non_stat_full_monty, keep = "range_beta", burn_in = 0) 
tracePlots(geo_non_stat_full_monty, keep = "range_log_scale", burn_in = 0) 
tracePlots(geo_non_stat_full_monty, keep = "field_log_var", burn_in = 0) 
tracePlots(geo_non_stat_full_monty, keep = "noise_log_scale", burn_in = 0) 
mcmcDiags(geo_non_stat_full_monty, .3)

# 
# geo_non_stat_full_monty$states$chain_1$ker_var$range_beta_sufficient
# geo_non_stat_full_monty$states$chain_1$ker_var$range_beta_ancillary
# 
# geo_non_stat_full_monty$states$chain_1$ker_var
# geo_non_stat_full_monty$states$chain_2$ker_var
# geo_non_stat_full_monty$states$chain_3$ker_var
# 
# solve(geo_non_stat_full_monty$states$chain_1$stuff$range_beta_empirical_fisher_s)
# 
# 
# samples = runOneChainMcmc(
#  covariates = geo_non_stat_full_monty$covariates, observed_field = geo_non_stat_full_monty$observed_field, 
#  hierarchical_model = geo_non_stat_full_monty$hierarchical_model, vecchia_approx = geo_non_stat_full_monty$vecchia_approx, 
#  state = geo_non_stat_full_monty$states$chain_1, n_iterations = 300, num_threads = 5, iter_start = 70
# )
# 
# 
