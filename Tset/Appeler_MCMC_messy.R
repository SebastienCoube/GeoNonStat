
# observed_locs = cbind(runif(40000), 1)
# 
# X = cbind(observed_locs[,1], rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)))
# 
# NNarray = GpGp::find_ordered_nn(observed_locs, 12)
# Linv = GpGp::vecchia_Linv(c(1, .0005, .0001), "matern15_isotropic", locs = observed_locs, NNarray)
# w = GpGp::fast_Gp_sim_Linv(Linv, NNarray)
# reg_coeffs = rnorm(ncol(X)+1)
# 
# observed_field = c(cbind(1, X) %*% reg_coeffs + w + rnorm(nrow(observed_locs)))
# 
# plot(observed_locs[,1], observed_field)
# points(observed_locs[,1], w, col = 2, cex = .5, pch = 16)
# 
# vecchia_approx = createVecchia(observed_locs)
# PP = createPP(vecchia_approx)
# 
# mygns = GeoNonStat(
#   vecchia_approx = vecchia_approx, 
#   observed_field = c(observed_field), X = as.data.frame(X), 
#   matern_smoothness = 1.5, anisotropic = F, 
#   n_chains = 3, 
#   noise_X = NULL, range_X = NULL, scale_X = NULL, noise_PP = NULL, 
#   range_PP = NULL, scale_PP = NULL, seed = 1
# )
# 
# names(mygns)
# 
# 
# list2env(mygns$states$chain_1, environment())
# list2env(mygns, environment())
# 
# n_iterations_update = 200
# num_threads = 5
# iter_start = 1
# iter = 1
# seed=123 
# 
# # # test only beta suff
# # 
# # params$field = w[vecchia_approx$hctam_scol_1]
# # stuff$noise_var[] = 1
# # 
# # source("R/Gibbs.R")
# # 
# # plot( params_records$beta[1,,])
# # abline(h = reg_coeffs[1])
# # plot(params_records$beta[2,,])
# # abline(h = reg_coeffs[2])
# # plot(params_records$beta[3,,])
# # abline(h = reg_coeffs[3])
# # plot(params_records$beta[4,,])
# # abline(h = reg_coeffs[4])
# # plot(params_records$beta[5,,])
# # abline(h = reg_coeffs[5])
# 
# 
# # # test field only
# # stuff$noise_var[] = 1
# # params$range_beta = reg_coeffs
# # source("R/Gibbs.R")
# # plot(observed_locs[,1], w, col = 1)
# # points(vecchia_approx$locs[,1], params$field, col = 2, pch = 16, cex = .5)
# 
# # test field and centered covariates
# stuff$noise_var[] = 1
# source("R/Gibbs.R")
# plot(observed_locs[,1], w, col = 1)
# points(vecchia_approx$locs[,1], params$field, col = 2, pch = 16, cex = .5)
# 
# 
# plot( params_records$beta[1,,])
# abline(h = reg_coeffs[1])
# plot(params_records$beta[2,,])
# abline(h = reg_coeffs[2])
# plot(params_records$beta[3,,])
# abline(h = reg_coeffs[3])
# plot(params_records$beta[4,,])
# abline(h = reg_coeffs[4])
# plot(params_records$beta[5,,])
# abline(h = reg_coeffs[5])



###############
# Anisotropic #
###############

set.seed(2)
observed_locs = cbind(runif(10000), runif(10000))[sample(seq(10000),size = 10000, replace = F),]
vecchia_approx = createVecchia(observed_locs, m = 6)

# fixed effects
X = as.data.frame(cbind(observed_locs[,1], rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs))))

#X_range = as.data.frame(observed_locs)
  
PP_range = createPP(vecchia_approx, knots = 16, matern_range = .3)
PP_noise = createPP(vecchia_approx, knots = 500, matern_range = .06)


GNSSimulator = GeoNonStatSimulator(
  vecchia_approx = vecchia_approx,
  X = X, 
  noise_X = X, noise_PP = PP_noise, 
  range_X = NULL,range_PP = PP_range, 
  anisotropic = T
)

coeff_list = create_coeff_list(GNSSimulator)
coeff_list$noise_PP_coeff[] = 0*rnorm(length(coeff_list$noise_PP_coeff))
coeff_list$noise_X_coeff[1] = -4


range_log_scale = c(-1,-1)
coeff_list$range_PP_coeff[,1] = exp(range_log_scale[1])*rnorm(nrow(coeff_list$range_PP_coeff))
coeff_list$range_PP_coeff[,-1] = exp(range_log_scale[2])*rnorm(2*nrow(coeff_list$range_PP_coeff))
coeff_list$range_X_coeff[1,1] = -4

fake_data = simulate(
  GNSSimulator = GNSSimulator, 
  coeff_list = coeff_list)
dev.off()


plot_pointillist_painting(vecchia_approx$locs, fake_data$latent_field, cex= 1)

plot_pointillist_painting(vecchia_approx$observed_locs, fake_data$observed_field)
plot_pointillist_painting(vecchia_approx$observed_locs, fake_data$noise)
plot_pointillist_painting(vecchia_approx$observed_locs, fake_data$log_noise_var_field)


mygns = GeoNonStat(
  vecchia_approx = vecchia_approx, 
  observed_field = fake_data$observed_field, X = X, 
  matern_smoothness = 1.5, anisotropic = T, 
  n_chains = 3, 
  noise_X = NULL, range_X = NULL,
  noise_PP = PP_noise,
  range_PP = PP_range
)

plot_pointillist_painting(mygns$vecchia_approx$locs, mygns$states$chain_1$params$field)
plot_pointillist_painting(mygns$vecchia_approx$observed_locs, mygns$observed_field)


n_iterations_update = 2000
num_threads = 12
iter_start = 1
iter = 1
seed=123 
list2env(mygns, environment())
state = mygns$states$chain_1
list2env(state, envir = environment())
remove(state)
#cl = parallel::makeCluster(num_threads)

params_records = list()

params$field_log_var = 0
params$range_log_scale = range_log_scale
stuff$noise_var[] = exp(fake_data$log_noise_var_field)

ker_var$range_beta_ancillary


source("Tset/MCMC_messy.R")

dev.off()
plot_pointillist_painting(vecchia_approx$locs, fake_data$latent_field, cex= .5)
plot_ellipses(vecchia_approx$locs[seq(100),], log_range = covariates$range_X$X_locs[seq(100),]%*%params$range_beta, add = T, shrink =  sqrt(8*1.5))

plot_pointillist_painting(vecchia_approx$observed_locs, X_range[,5], cex= .5)
