library(GeoNonStat)


Rendre noise X NULL en cas de full stationary
# loading the data set and performing basic visualization ####
synthetic <- readRDS("Tset/synthetic_data_set.RDS")
plotPointillistPainting(synthetic$vecchia_approx$locs, synthetic$hidden_fields$latent_field)
plotPointillistPainting(synthetic$observed$locs, synthetic$observed$observed_field)

# Splitting between train and test ####
train_test_split <- GeoNonStat::splitData(
  data = synthetic$observed, locs = synthetic$observed$locs, n_clust = 1000, round_locs = .00, prop_test = c(.2,.1))

# Model structures ####
vecchia_approx <- createVecchia(train_test_split$train$locs, m = 10)
PP_long_range <- createPP(vecchia_approx = vecchia_approx, matern_range = .25)
PP_short_range <- createPP(vecchia_approx = vecchia_approx, matern_range = .1)

# Run ####
# design
run_design <- matrix(F, 4, 3)
colnames(run_design) <- c("heterosk noise", "nonstat range", "aniso range")
run_design[-1,"heterosk noise"] <- T
run_design[-c(1,2),"nonstat range"] <- T
run_design[-c(1,2,3),"aniso range"] <- T
print(run_design)

# running according to design
res <- list()
for(i in seq(nrow(run_design))){
  noise_PP = NULL
  if(run_design[i,"heterosk noise"])noise_PP = PP_short_range
  range_PP = NULL
  if(run_design[i,"nonstat range"])range_PP = PP_long_range
  aniso = F
  if(run_design[i,"aniso range"])aniso = T
  
  # initialize
  geo_non_stat <- GeoNonStat(
    vecchia_approx = vecchia_approx, 
    observed_field = train_test_split$train$observed_field, 
    X = train_test_split$train$X, 
    noise_X = train_test_split$train$noise_X, noise_PP = noise_PP, 
    range_PP = range_PP, anisotropic = aniso, 
    n_chains = 3, matern_smoothness = 1.5
  )
  # run
  future::plan(strategy = "multisession", workers = 3)
  geo_non_stat <- automaticMcmc(
    object = geo_non_stat, 
    n_threads_per_chain =  5, satisfying_ESS = 100, 
    satisfying_Gelman_Rubin = 1.10, burn_in = .3, 
    iter_per_step = 200, verbose = T)
  GeoNonStat::tracePlots(geo_non_stat, .0)
  # Criteria
  scores <- list()
  scores$train <- trainScores(geo_non_stat, burn_in = .3)
  scores$close <- testScores(
    geo_non_stat, burn_in = .3, num_threads = 8, 
    test_noise_X = train_test_split$test_close$noise_X, test_range_X = NULL, 
    test_X = train_test_split$test_close$X, 
    test_observed_field = train_test_split$test_close$observed_field, 
    test_locs = train_test_split$test_close$locs)
  scores$far <- testScores(
    geo_non_stat, burn_in = .3, num_threads = 8, 
    test_noise_X = train_test_split$test_far$noise_X, test_range_X = NULL, 
    test_X = train_test_split$test_far$X, 
    test_observed_field = train_test_split$test_far$observed_field, 
    test_locs = train_test_split$test_far$locs)
  saveRDS(object = geo_non_stat, file = paste("geo_non_stat", i,".RDS", sep=""))
  res[[i]] <- list(scores = scores)
  gc()
}

# Model comparison #### 

criteria <- rbind(
c(res[[1]]$scores$train$scores, res[[1]]$scores$far$scores, res[[1]]$scores$close$scores),
c(res[[2]]$scores$train$scores, res[[2]]$scores$far$scores, res[[2]]$scores$close$scores),
c(res[[3]]$scores$train$scores, res[[3]]$scores$far$scores, res[[3]]$scores$close$scores),
c(res[[4]]$scores$train$scores, res[[4]]$scores$far$scores, res[[4]]$scores$close$scores)
)
colnames(criteria) <- c(
  paste("train", names(res[[4]]$scores$train$scores)),
  paste("far", names(res[[4]]$scores$train$scores)),
  paste("close", names(res[[4]]$scores$train$scores))
  )
criteria <- cbind(run_design, criteria)
