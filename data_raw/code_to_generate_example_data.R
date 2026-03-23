set.seed(2)
nlocs = 3000
observed_locs = cbind(runif(nlocs), runif(nlocs))
vecchia_approx = createVecchia(observed_locs, m = 6)

# fixed effects
X = data.frame(
  V1 = observed_locs[,1],
  V2 = rnorm(nrow(observed_locs)),
  V3 = rnorm(nrow(observed_locs)),
  V4 = rnorm(nrow(observed_locs))
)

PP_range = createPP(vecchia_approx,
                    knots = 16,
                    matern_range = .3,
                    plot=FALSE)
PP_noise = createPP(vecchia_approx,
                    knots = 30,
                    matern_range = .25,
                    plot=FALSE)

gns_simulator = createGnsSimulator(
  vecchia_approx = vecchia_approx,
  X = X,
  noise_X = X, noise_PP = PP_noise,
  range_X = NULL,range_PP = PP_range,
  anisotropic = TRUE
)

coeff_list = createGnsSimulatorParameters(
  gns_simulator,
  range_PP_log_var =c(-1,-1),
  noise_intercept = 1,
  noise_PP_log_var = -1,
  field_log_var = 2,
  verbose = FALSE)
coeff_list$range_X_coeff[1] = -5

fake_data = simulateGnsData(
  gns_simulator = gns_simulator,
  gns_params = coeff_list)

gnsDemo = GeoNonStat(
  vecchia_approx = vecchia_approx,
  observed_field = fake_data$observed_field, X = X,
  matern_smoothness = 1.5, anisotropic = T,
  n_chains = 3,
  noise_X = NULL, range_X = NULL,
  noise_PP = PP_noise,
  range_PP = PP_range
)
save(gnsDemo, file="data/gnsDemo.rda", compress="xz")




set.seed(123)
processedGnsDemo = GeoNonStatMcmc(
  object = gnsDemo,
  n_chains_in_parallel = 1,
  n_threads_per_chain = 1,
  n_iterations = 15,
  seed = 1
)
save(processedGnsDemo, file="data/processedGnsDemo.rda", compress="xz")
