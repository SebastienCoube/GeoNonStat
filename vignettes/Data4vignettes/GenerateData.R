# Small working example

library(GeoNonStat)
set.seed(2)
nobs = 10000
nlocs = 10000
observed_locs = cbind(runif(nlocs), runif(nlocs))[c(seq(nlocs), sample(seq(nlocs), nobs - nlocs, T)),]

vecchia_approx = createVecchia(observed_locs, m = 6)

# fixed effects
X = as.data.frame(cbind(vecchia_approx$locs[vecchia_approx$locs_match,1], vecchia_approx$observed_locs[,1] + rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs))))

PP_range = createPP(vecchia_approx, matern_range = .5, plot = F)
plot(PP_range)
PP_noise = createPP(vecchia_approx, matern_range = .2)
plot(PP_noise)

gns_simulator = createGnsSimulator(
  vecchia_approx = vecchia_approx,
  X = X, 
  noise_X = X, noise_PP = PP_noise, 
  range_X = NULL,range_PP = PP_range, 
  anisotropic = T, matern_smoothness = 1.5
)

field_log_var=  0
range_log_scale = c(0, 0)
noise_intercept = 0
noise_PP_log_var = -.5

coeff_list = createGnsSimulatorParameters(
  gns_simulator, range_PP_log_var =range_log_scale, 
  noise_intercept = noise_intercept, noise_PP_log_var = noise_PP_log_var, 
  field_log_var = field_log_var)

coeff_list$range_X_coeff[] = 0
coeff_list$range_X_coeff[1,1] = -4.5
coeff_list$range_PP_coeff[,-1] = 0
coeff_list$noise_X_coeff[-1] = .1*rnorm(4)
set.seed(2)

fake_data = simulateGnsData(gns_simulator = gns_simulator, gns_params = coeff_list)

plotPointillistPainting(vecchia_approx$locs, fake_data$hidden_fields$latent_field, cex= .7, pch= 15)
plotPointillistPainting(vecchia_approx$observed_locs, fake_data$observed$observed_field, cex= 1, pch= 15)
plotPointillistPainting(vecchia_approx$observed_locs, fake_data$hidden_fields$log_noise_var_field, cex= 1, pch= 15)
fake_data$observed$range_X <- NULL
saveRDS(object = fake_data$observed, "vignettes/Data4vignettes/smallExample.RDS")
saveRDS(object = fake_data, "vignettes/Data4vignettes/smallExampleFull.RDS")




