library(GeoNonStat)
set.seed(1)
nobs = 750
nlocs = 700
observed_locs = cbind(runif(nlocs), runif(nlocs))[c(seq(nlocs), sample(seq(nlocs), nobs - nlocs, T)),]

vecchia_approx = createVecchia(observed_locs, m = 10)

# fixed effects
X = as.data.frame(cbind(vecchia_approx$observed_locs[,1] + rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs))))

PP = createPP(vecchia_approx, matern_range = .25, plot = F)

gns_simulator = createGnsSimulator(
  vecchia_approx = vecchia_approx,
  X = X, 
  noise_X = X, noise_PP = PP, 
  range_X = NULL,range_PP = PP, 
  anisotropic = T, matern_smoothness = 1.5
)


set.seed(1)
# case with nice non-stationarity and little noise. 
range_log_scale = c(-1,-1)
field_log_var=  1
noise_intercept = -2
noise_PP_log_var = 2
# case with no non-stationarity and little noise. 
# case with nice non-stationarity and crazy noise. 
# case with no non-stationarity and crazy noise. 
coeff_list = createGnsSimulatorParameters(
  gns_simulator, range_PP_log_var =range_log_scale, 
  noise_intercept = noise_intercept, noise_PP_log_var = noise_PP_log_var, 
  field_log_var = field_log_var)
coeff_list$X_coeff[] <- 0
coeff_list$range_X_coeff[1,1] = -2
fake_data = simulateGnsData(gns_simulator = gns_simulator, gns_params = coeff_list)

plotPointillistPainting(vecchia_approx$locs, fake_data$hidden_fields$latent_field, main = "latent", cex = 8)
plotPointillistPainting(vecchia_approx$observed_locs, fake_data$observed$observed_field, main = "observed", cex = 8)
plotPointillistPainting(vecchia_approx$observed_locs, fake_data$hidden_fields$log_noise_var_field, main = "noise var", cex = 8)


gnsDemo = GeoNonStat(
  vecchia_approx = vecchia_approx,
  data_list = fake_data$observed[c("X", "noise_X", "observed_field")], 
  matern_smoothness = 1.5, anisotropic = ,
  n_chains = 3,
  noise_PP = PP,
  range_PP = PP
)
save(gnsDemo, file="data/gnsDemo.rda", compress="xz")
processedGnsDemo = GeoNonStatMcmc(
  object = gnsDemo,
  n_threads_per_chain = 5,
  n_iterations = 22
)
tracePlots(processedGnsDemo)
save(processedGnsDemo, file="data/processedGnsDemo.rda", compress="xz")
