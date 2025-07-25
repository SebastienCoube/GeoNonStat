

observed_locs = cbind(runif(40000), 1)

X = cbind(observed_locs[,1], rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)))

NNarray = GpGp::find_ordered_nn(observed_locs, 12)
Linv = GpGp::vecchia_Linv(c(1, .0005, .0001), "matern15_isotropic", locs = observed_locs, NNarray)
w = GpGp::fast_Gp_sim_Linv(Linv, NNarray)
reg_coeffs = rnorm(ncol(X)+1)

observed_field = c(cbind(1, X) %*% reg_coeffs + w + rnorm(nrow(observed_locs)))

plot(observed_locs[,1], observed_field)
points(observed_locs[,1], w, col = 2, cex = .5, pch = 16)


vecchia_approx = createVecchia(observed_locs)
PP = createPP(vecchia_approx)

mygns = GeoNonStat(
  vecchia_approx = vecchia_approx, 
  observed_field = c(observed_field), X = as.data.frame(X), 
  matern_smoothness = 1.5, anisotropic = F, 
  n_chains = 3, 
  noise_X = NULL, range_X = NULL, scale_X = NULL, noise_PP = NULL, 
  range_PP = NULL, scale_PP = NULL, seed = 1
)

names(mygns)


list2env(mygns$states$chain_1, environment())
list2env(mygns, environment())

n_iterations_update = 200
num_threads = 5
iter_start = 1
iter = 1
seed=123 

# # test only beta suff
# 
# params$field = w[vecchia_approx$hctam_scol_1]
# stuff$noise_var[] = 1
# 
# source("R/Gibbs.R")
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


# # test field only
# stuff$noise_var[] = 1
# params$range_beta = reg_coeffs
# source("R/Gibbs.R")
# plot(observed_locs[,1], w, col = 1)
# points(vecchia_approx$locs[,1], params$field, col = 2, pch = 16, cex = .5)

# test field and centered covariates
stuff$noise_var[] = 1
source("R/Gibbs.R")
plot(observed_locs[,1], w, col = 1)
points(vecchia_approx$locs[,1], params$field, col = 2, pch = 16, cex = .5)


plot( params_records$beta[1,,])
abline(h = reg_coeffs[1])
plot(params_records$beta[2,,])
abline(h = reg_coeffs[2])
plot(params_records$beta[3,,])
abline(h = reg_coeffs[3])
plot(params_records$beta[4,,])
abline(h = reg_coeffs[4])
plot(params_records$beta[5,,])
abline(h = reg_coeffs[5])








###############
# Anisotropic #
###############

observed_locs = cbind(runif(40000), runif(40000))
vecchia_approx = createVecchia(observed_locs)

# fixed effects
X = cbind(observed_locs[,1], rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)), rnorm(nrow(observed_locs)))
X = process_covariates(X = as.data.frame(X), vecchia_approx = vecchia_approx, PP = PP)
reg_coeffs = rnorm(ncol(X$X))

range_X = observed_locs
PP = createPP(vecchia_approx)
range_X = process_covariates(X = as.data.frame(range_X), vecchia_approx = vecchia_approx, PP = PP)
range_beta = matrix(0, range_X$n_regressors + PP$n_knots, 3)
range_beta[] = .5*rnorm(length(range_beta))
range_beta[1,1] = -4

compressed_sparse_chol = compute_sparse_chol(
  range_beta = range_beta, vecchia_approx = vecchia_approx, 
  range_X = range_X, PP = PP, matern_smoothness = 1.5, 
  compute_derivative = F, num_threads = 6)
decompressed_chol = decompress_chol(vecchia_approx, compressed_sparse_chol)

w = as.vector(Matrix::solve(decompressed_chol, rnorm(nrow(vecchia_approx$locs))))
plot_pointillist_painting(vecchia_approx$locs, w, cex = .5)

noise_var = 1



mygns = GeoNonStat(
  vecchia_approx = vecchia_approx, 
  observed_field = c(observed_field), X = as.data.frame(X), 
  matern_smoothness = 1.5, anisotropic = F, 
  n_chains = 3, 
  noise_X = NULL, range_X = NULL, scale_X = NULL, noise_PP = NULL, 
  range_PP = NULL, scale_PP = NULL, seed = 1
)

names(mygns)


list2env(mygns$states$chain_1, environment())
list2env(mygns, environment())

n_iterations_update = 200
num_threads = 5
iter_start = 1
iter = 1
seed=123 
