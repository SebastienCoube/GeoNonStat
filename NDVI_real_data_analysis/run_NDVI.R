library(GeoNonStat)

set.seed(1)
number_NNGP_neighbors <- 10
n_MCMC_chains <- 3
n_threads <- 8

out_dir <- "NDVI_real_data_analysis"

# Loading the data set ####
load(file.path(out_dir, "data_cleaned_small_expanded.RData"))

# Splitting between train and test ####
train_test_split <- GeoNonStat::splitData(
  data = list(observed_field = data_cleaned_small$NDVI),
  locs = cbind(data_cleaned_small$scaled_x, data_cleaned_small$scaled_y),
  n_clust = 1000, round_locs = .004, prop_test = c(.2, .1)
)
saveRDS(train_test_split, file.path(out_dir, "train_test_split.RDS"))

plotPointillistPainting(train_test_split$train$locs,
                        train_test_split$train$observed_field)

# Model structures ####
vecchia_approx <- createVecchia(train_test_split$train$locs, m = number_NNGP_neighbors)
PP_large <- createPP(vecchia_approx = vecchia_approx, matern_range = .16)
PP_small <- createPP(vecchia_approx = vecchia_approx, matern_range = .08)

# Run design ####
run_design <- matrix(F, 4, 3)
colnames(run_design) <- c("heterosk noise", "nonstat range", "aniso range")
run_design[-1, "heterosk noise"] <- T
run_design[-c(1, 2), "nonstat range"] <- T
run_design[-c(1, 2, 3), "aniso range"] <- T
print(run_design)

# Run, predict and score according to the design ####
scores <- NULL

for (i in seq(3, 4)) {
  noise_PP <- NULL
  noise_log_scale_bounds = NULL
  if (run_design[i, "heterosk noise"]) {
    noise_PP <- PP_small
    noise_log_scale_bounds <- c(0, 5)
  }
  range_PP <- NULL
  if (run_design[i, "nonstat range"]) range_PP <- PP_large
  aniso <- F
  if (run_design[i, "aniso range"]) aniso <- T
  
  # initialize
  geo_non_stat <- GeoNonStat(
    vecchia_approx = vecchia_approx, data_list = train_test_split$train,
    noise_PP = noise_PP, range_PP = range_PP, anisotropic = aniso,
    n_chains = n_MCMC_chains, matern_smoothness = 1.5, 
    noise_log_scale_bounds = noise_log_scale_bounds
  )
  
  # run MCMC (in a separate R session)
  model_save_name <- file.path(out_dir, paste0("run", i, ".RDS"))
  diags_save_name <- file.path(out_dir, paste0("model", i, ".pdf"))
  future::plan(strategy = "multisession", workers = length(geo_non_stat$states))
  geo_non_stat <- GeoNonStat::automaticMcmc(
    object = geo_non_stat,
    n_threads_per_chain = n_threads,
    satisfying_ESS = 100, burn_in = .3,
    satisfying_Gelman_Rubin = 1.1, iter_per_step = 100,
    save_checkpoint_model = model_save_name,
    save_checkpoint_plots = diags_save_name
  )
  future::plan(strategy = "sequential")
  # 
  # 
  # # predict on test sets
  # prediction <- predict.GeoNonStat(
  #   object = geo_non_stat,
  #   new_data_lists = train_test_split[c("test_close", "test_far")],
  #   num_threads = n_threads, burn_in = .3
  # )
  # saveRDS(prediction, file.path(out_dir, paste0("prediction", i)))
  # 
  # # estimation (train scores)
  # estimation <- estimate(geo_non_stat)
  # saveRDS(estimation, file.path(out_dir, paste0("estimation", i)))
  # 
  # scores <- rbind(scores, c(estimation$scores$aggregated,
  #                           prediction$test_close$scores$aggregated,
  #                           prediction$test_far$scores$aggregated))
  # colnames(scores) <- paste(rep(c("train", "close", "far"), each = 4),
  #                           colnames(prediction$test_close$scores$per_obs))
  # print(scores)
  # rm(geo_non_stat, prediction, estimation)
  # gc()
}

saveRDS(scores, file.path(out_dir, "scores.RDS"))