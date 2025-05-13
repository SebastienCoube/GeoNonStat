

set.seed(100)
# spatial locations 
locs = cbind(runif(10000), runif(10000))

# Predictive Process
PP = GeoNonStat::createPP(
  observed_locs = locs[seq(10000),], # spatial sites
  matern_range = .1,
  n_PP = 50, # number of knots
  m = 15 # number of NNGP parents
)

# range parameters
range_beta = rbind(
  c(-4, 0, 0), # intercept
  matrix(.5*rnorm(150), 50)
) %*% diag(c(1, 2,.5))
NNarray_aniso = GpGp::find_ordered_nn(locs[seq(10000),], 10)

# GP coefficients
chol_precision = 
  GeoNonStat::compute_sparse_chol(
    range_beta = range_beta, # range parameters
    NNarray = NNarray_aniso, # Vecchia approx DAG
    locs = locs[seq(10000),], # spatial coordinates
    range_X = matrix(1, 10000), # covariates for the range (just an intercept)
    PP = PP, use_PP = T, # predictive process
    nu = 1.5, # smoothness
    anisotropic = T # anisotropy
  )[[1]]
# putting coefficients in precision Cholesly

chol_precision = Matrix::sparseMatrix(
  x = chol_precision[!is.na(NNarray_aniso)], 
  i = row(NNarray_aniso)[!is.na(NNarray_aniso)], 
  j = (NNarray_aniso)[!is.na(NNarray_aniso)], 
  triangular = T
)
# sampling the anisotropic process
seed_vector = rnorm(10000)
aniso_latent_field = as.vector(Matrix::solve(chol_precision, seed_vector)) 
aniso_observed_field = aniso_latent_field + .8*rnorm(10000)

GeoNonStat::plot_pointillist_painting(locs, aniso_latent_field)
  
MCMC_NNGP = GeoNonStat::initialize(
  observed_locs = locs[seq(10000),], observed_field = aniso_observed_field, 
  nu = 1.5, n_chains = 5,
  range_PP = T, PP = PP, # use PP for range
  anisotropic = F # Covariance will be anisotropic
)

t1 = Sys.time()
MCMC_samples = GeoNonStat::update(
  MCMC_NNGP, # model MCMC_NNGP$state
  n_iterations = 1, seed = 1, 
  num_threads = 1 # practical settings
)
Sys.time()-t1




MCMC_NNGP = Bidart::mcmc_nngp_initialize_nonstationary(
  observed_locs = locs[seq(10000),], observed_field = aniso_observed_field, 
  nu = 1.5, n_chains = 5,
  range_PP = T, PP = PP, # use PP for range
  anisotropic = F # Covariance will be anisotropic
)

t1 = Sys.time()
MCMC_samples = Bidart::mcmc_nngp_update_Gaussian(
  data = MCMC_NNGP$data, hierarchical_model = MCMC_NNGP$hierarchical_model, 
  vecchia_approx = MCMC_NNGP$vecchia_approx, state = MCMC_NNGP$states$chain_1, 
  n_iterations_update = 1, thinning = 1, iter_start = 1, num_threads = 1
)
Sys.time()-t1


