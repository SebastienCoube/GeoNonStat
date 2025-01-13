#' @export
get_array_summary = function(samples)
{
  if(any(is.na(samples)))message("there are NAs in the samples")
  if(is.null(samples))return(NULL)
  if(dim(samples)[2]>1 |  length(dim(samples))==2){
      out = apply(X = samples, MARGIN = c(1, 2), FUN = function(x)c(mean(x, na.rm = T), quantile(x, c(0.025, 0.5, 0.975), na.rm = T), sd(x, na.rm = T)))
  }
  if(dim(samples)[2]==1){
    cl = parallel::makeCluster(parallel::detectCores()-1)
    out = parallel::parApply(cl = cl, X = samples, MARGIN = c(1), FUN = function(x)c(mean(x, na.rm = T), quantile(x, c(0.025, 0.5, 0.975), na.rm = T), sd(x, na.rm = T)))
    parallel::stopCluster(cl)
    out = array(out, dim = c(nrow(out), ncol(out), 1))
  }
  dimnames(out)[[1]] = c("mean", "q0.025", "median", "q0.975", "sd")
  out
}

#' Predicts the latent field at unobserved locations
#' @param mcmc_nngp_list A mcmc_nngp_list object generated using mcmc_nngp_initialize and ran using mcmc_nngp_run
#' @param predicted_locs A set of predicted locations disjunct from the locations of mcmc_nngp_list
#' @param X_range_pred Covariates for the range observed at predicted_locs
#' @param X_scale_pred Covariates for the marginal variance observed at predicted_locs
#' @param burn_in MCMC burn-in
#' @param num_threads_per_chain Number of OMP threads for the predictions
#' @param  lib.loc Location of libraries if installed in manually specified folder
#' @export
predict_latent_field = function(mcmc_nngp_list, predicted_locs, X_range_pred = NULL, X_scale_pred = NULL, burn_in = .5, num_threads_per_chain = 1, lib.loc = NULL, parallel = T)
{
  # Sanity checks
  if(length(intersect(split(predicted_locs, row(predicted_locs)), split(mcmc_nngp_list$data$locs, row(mcmc_nngp_list$data$locs)) )) > 0) stop("the predicted locations must contain none of the original observed locations, use estimate instead")
  if(!is.null(X_range_pred))if(nrow(X_range_pred)!=nrow(predicted_locs))stop("X_range_pred must have the same Number of rows as predicted_locs")
  if(!is.null(X_scale_pred))if(nrow(X_scale_pred)!=nrow(predicted_locs))stop("X_scale_pred must have the same Number of rows as predicted_locs")
  if(!is.null(X_scale_pred))if(!is.null(X_range_pred))if(nrow(X_scale_pred)!=nrow(X_range_pred))stop("X_scale_pred and X_range_pred must have the same Number of rows")
  n_cores = min(parallel::detectCores()-1, length(mcmc_nngp_list$records))
  # Dealing with duplicated data
  duplicated_idx = duplicated(predicted_locs)
  order_idx = GpGp::order_maxmin(predicted_locs[!duplicated_idx,])
  X_range_pred_unique   = X_range_pred[!duplicated_idx,][order_idx,]
  X_scale_pred_unique   = X_scale_pred[!duplicated_idx,][order_idx,]
  predicted_locs_unique = predicted_locs[!duplicated_idx,][order_idx,]
  # surjection from duplicated locs to unique locs
  locs_match = match(split(predicted_locs, row(predicted_locs)), split(predicted_locs_unique, row(predicted_locs_unique)))
  # injection from unique locs to duplicated locs
  hctam_scol_1 = sapply(split(seq(nrow(predicted_locs)), locs_match), function(x)x[1])
  # adding intercept to prediction regressors
  X_range_pred_unique = cbind(rep(1, nrow(predicted_locs_unique)), X_range_pred_unique)
  X_scale_pred_unique = cbind(rep(1, nrow(predicted_locs_unique)), X_scale_pred_unique)
  # binding regressors of observed and predicted locations 
  X_range = rbind(as.matrix(mcmc_nngp_list$data$covariates$range_X$X_locs), as.matrix(X_range_pred_unique))
  X_scale = rbind(as.matrix(mcmc_nngp_list$data$covariates$scale_X$X_locs), as.matrix(X_scale_pred_unique))
  # constructing NNarray of observed and predicted locations
  locs = rbind(mcmc_nngp_list$data$locs, predicted_locs_unique)
  NNarray = rbind(mcmc_nngp_list$vecchia_approx$NNarray, GpGp::find_ordered_nn(locs, m = ncol(mcmc_nngp_list$vecchia_approx$NNarray)-1)[-seq(mcmc_nngp_list$vecchia_approx$n_locs),])
  NNarray_non_NA = !is.na(NNarray)
  sparse_chol_row_idx = row(NNarray)[NNarray_non_NA]
  sparse_chol_column_idx = NNarray[NNarray_non_NA]
  # constructing NNarray of observed and predicted locations for parameter fields
  if(mcmc_nngp_list$hierarchical_model$scale_PP | mcmc_nngp_list$hierarchical_model$range_PP) 
  {
    PP = mcmc_nngp_list$hierarchical_model$PP
    PP_ = mcmc_nngp_list$hierarchical_model$PP
    # appending predicted locations to NNarray
    NNarray_p = GpGp::find_ordered_nn(rbind(PP$knots, PP$unique_reordered_locs, predicted_locs_unique), m = ncol(PP_$NNarray)-1)
    NNarray_non_NA_p = !is.na(NNarray_p)
    sparse_chol_row_idx_p = row(NNarray_p)[NNarray_non_NA_p]
    sparse_chol_column_idx_p = NNarray_p[NNarray_non_NA_p]
    PP_$unique_reordered_locs = rbind(PP$unique_reordered_locs, predicted_locs_unique)
    PP_$sparse_chol =  Matrix::sparseMatrix(
      i = sparse_chol_row_idx_p, 
      j = sparse_chol_column_idx_p, 
      x = GpGp::vecchia_Linv(
        covparms = c(1, PP$matern_range, .0001),
        covfun_name = "matern15_isotropic",
        NNarray = NNarray_p, 
        locs = rbind(PP_$knots, PP_$unique_reordered_locs))[NNarray_non_NA_p],
      triangular = T
    )
    PP_$idx = c(PP$idx, nrow(PP$unique_reordered_locs)+locs_match)
    #### testing PP
    ## seed_vector =  rnorm(PP_$n_PP + nrow(PP_$unique_reordered_locs))
    ## Bidart::plot_pointillist_painting(PP_$unique_reordered_locs[PP_$idx,], Bidart::X_PP_mult_right(PP = PP_, use_PP = T, Y = seed_vector[seq(PP$n_PP)]), cex = .3)
    ## Bidart::plot_pointillist_painting(PP$unique_reordered_locs[PP$idx,], Bidart::X_PP_mult_right(PP = PP, use_PP = T, Y = seed_vector[seq(PP$n_PP)]), cex = .3)
  }
  
  # burn in
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
  # useful stuff
  n_samples = length(kept_iterations)
  hierarchical_model = mcmc_nngp_list$hierarchical_model
  vecchia_approx = mcmc_nngp_list$vecchia_approx
  # samples function
  get_samples = function(chain)
  {
    res = list()
    # creating arrays
    # fields
    res$field = array(0, dim = c(nrow(predicted_locs_unique), 1, n_samples))
    res$log_scale = array(0, dim = c(nrow(predicted_locs_unique), 1, n_samples))
    res$log_range = array(0, dim = c(nrow(predicted_locs_unique), dim(chain$range_beta)[2], n_samples))
    # looping over saved observations
    i_predict = 20
    for(i_predict in seq(n_samples))
    {
      # computing NNGP factor
      sparse_chol = 
        Matrix::sparseMatrix(
          i = sparse_chol_row_idx, 
          j = sparse_chol_column_idx, 
          x = 
            Bidart::compute_sparse_chol(
              range_beta = matrix(chain$range_beta[,,i_start + i_predict], nrow =dim(chain$range_beta)[1]), 
              NNarray = NNarray, locs = locs, 
              anisotropic = hierarchical_model$anisotropic, 
              sphere = hierarchical_model$sphere,
              PP = PP_, use_PP = hierarchical_model$range_PP, 
              range_X = X_range, 
              compute_derivative = F, 
              nu = hierarchical_model$nu, 
              locs_idx = c(vecchia_approx$hctam_scol_1, vecchia_approx$n_obs+hctam_scol_1), 
              num_threads = num_threads_per_chain
            )[[1]][NNarray_non_NA], 
          triangular = T
        )
      res$log_range[,,i_predict] = 
        (
          Bidart::X_PP_mult_right(
            X = X_range, 
            PP = PP_, use_PP = hierarchical_model$range_PP, 
            locs_idx = c(vecchia_approx$hctam_scol_1, vecchia_approx$n_obs+hctam_scol_1), 
            Y = matrix(chain$range_beta[,,i_start + i_predict], nrow =dim(chain$range_beta)[1])
          )[-seq(vecchia_approx$n_locs),,drop = F]
        )
      scale_field = Bidart::variance_field(
        beta = chain$scale_beta[,,i_start + i_predict], PP = PP_, use_PP = hierarchical_model$scale_PP,  
        locs_idx =  c(vecchia_approx$hctam_scol_1, vecchia_approx$n_obs+hctam_scol_1),
        X = X_scale
      )
      res$log_scale[,,i_predict] = 
        (
          log(scale_field[-seq(vecchia_approx$n_locs)])
        )
      
      # predicting latent field #######
      res$field[,,i_predict] = 
        (
          sqrt(scale_field) * 
            as.vector(Matrix::solve(
              sparse_chol, 
              c(
                as.vector(sparse_chol[1:vecchia_approx$n_locs, 1:vecchia_approx$n_locs] %*% (chain$field[,,i_start + i_predict]/sqrt(scale_field[1:vecchia_approx$n_locs]))), 
                rnorm(nrow(predicted_locs_unique))
              )
            ))
        )[-seq(vecchia_approx$n_locs)]
    }
    return(res)
  }
  # parallelization
  if(parallel){
    cl = parallel::makeCluster(n_cores, outfile = NULL)
    parallel::clusterExport(cl = cl, varlist = setdiff(ls(envir = environment()),c("cl", "mcmc_nngp_list")), envir = environment())
    if(!is.null(lib.loc)){
      parallel::clusterEvalQ(cl = cl, expr = library(GpGp, lib.loc =   lib.loc))
      parallel::clusterEvalQ(cl = cl, expr = library(Bidart, lib.loc = lib.loc))
      parallel::clusterEvalQ(cl = cl, expr = library(expm, lib.loc =   lib.loc))
      parallel::clusterEvalQ(cl = cl, expr = library(Matrix, lib.loc =   lib.loc))
    }
    predicted_samples = 
      parallel::parLapply(
                          cl,
                          X = mcmc_nngp_list$records, 
                          fun = get_samples
      )
  }
  if(!parallel){
    predicted_samples = 
      lapply(
        X = mcmc_nngp_list$records,
        FUN = get_samples
      )
  }
  summaries = list()
  predicted_samples_ = list()
  for(name in names(predicted_samples[[1]]))
  {
    predicted_samples_[[name]] = Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(predicted_samples, function(y)y[[name]]))
  }
  summaries = parallel::parLapply(cl = cl, X = predicted_samples_, fun = Bidart::get_array_summary)
  for(name in names(predicted_samples[[1]]))
  {
    predicted_samples_[[name]] = predicted_samples_[[name]][locs_match,,,drop=F]
    summaries[[name]] = summaries[[name]][,locs_match,,drop=F]
  }
  parallel::stopCluster(cl)
  return(list("predicted_locs_unique" = predicted_locs_unique,"predicted_samples" = predicted_samples_, "summaries" = summaries))
}

#' Predicts the fixed effects, whatever the locations
#' @param mcmc_nngp_list A mcmc_nngp_list object generated using mcmc_nngp_initialize and ran using mcmc_nngp_run
#' @param X_pred New covariates
#' @param burn_in MCMC burn-in
#' @param  lib.loc Location of libraries if installed in manually specified folder
#' @export
predict_fixed_effects = function(mcmc_nngp_list, X_pred = NULL, burn_in = .5,  lib.loc = NULL)
{
  n_cores = min(parallel::detectCores()-1, length(mcmc_nngp_list$records))
  if(!is.null(X_pred)){if(ncol(X_pred)!=ncol(mcmc_nngp_list$data$covariates$X$arg))stop("The Number of provided covariates does not match")}
  # adding intercept to prediction regressors
  if(!is.null(X_pred))X_pred = cbind(rep(1, nrow(X_pred)), X_pred)
  if(is.null(X_pred))X_pred = matrix(1, 1, 1)
  colnames(X_pred)[1] = "(Intercept)"
  # converting to matrix
  X_pred  =as.matrix(X_pred)
  # burn in
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
  n_samples = length(kept_iterations)
  cl = parallel::makeCluster(n_cores, outfile = NULL)
  hierarchical_model = mcmc_nngp_list$hierarchical_model
  vecchia_approx = mcmc_nngp_list$vecchia_approx
  parallel::clusterExport(cl = cl, varlist = setdiff(ls(envir = environment()),c("cl", "mcmc_nngp_list")), envir = environment())
  
  if(!is.null(lib.loc)){
    parallel::clusterEvalQ(cl = cl, expr = library(GpGp, lib.loc =   lib.loc))
    parallel::clusterEvalQ(cl = cl, expr = library(Bidart, lib.loc = lib.loc))
    parallel::clusterEvalQ(cl = cl, expr = library(expm, lib.loc =   lib.loc))
    parallel::clusterEvalQ(cl = cl, expr = library(Matrix, lib.loc =   lib.loc))
  }
  predicted_samples = parallel::parLapply(
    cl = cl,
    X = mcmc_nngp_list$records, 
    fun = function(chain)
    {
      res = array(0, dim = c(nrow(X_pred), 1, n_samples))
      # looping over saved observations
      for(i_predict in seq(n_samples))
      {
        res[,,i_predict] = X_pred %*% chain$beta[,,i_start + i_predict]
      }
      return(res)
    })
  parallel::stopCluster(cl)
  summaries = list()
  predicted_samples_ = list()
  library(abind())
  predicted_samples_ = do.call("abind", predicted_samples)
  summaries = Bidart::get_array_summary(predicted_samples_)
  return(list("predicted_samples" = predicted_samples_, "summaries" = summaries))
}

#' Predicts the noise variance at unobserved locations
#' @param mcmc_nngp_list A mcmc_nngp_list object generated using mcmc_nngp_initialize and ran using mcmc_nngp_run
#' @param predicted_locs A set of predicted locations, who can be the same as the locations in mcmc_nngp_list
#' @param X_noise_pred Covariates for the noise observed at predicted_locs
#' @param burn_in MCMC burn-in
#' @export
predict_noise = function(mcmc_nngp_list, X_noise_pred = NULL, burn_in = .5, predicted_locs = NULL)
{
  if(mcmc_nngp_list$hierarchical_model$noise_PP & is.null(predicted_locs))stop("no predicted locs were provided, yet there is a GP prior for the noise")
  if(!is.null(X_noise_pred)){if(identical(X_noise_pred, "No covariates were provided")){X_noise_pred = NULL}}
  if(!is.null(X_noise_pred))if(nrow(predicted_locs)!= nrow(X_noise_pred))stop("X_noise_pred and predicted_locs should have the same Number of rows")
  if(is.null(X_noise_pred)&(!identical(mcmc_nngp_list$data$covariates$noise_X$arg,"No covariates were provided"))) stop("No covariates were provided for prediction, while covariates were provided for fit")
  if(!is.null(X_noise_pred)&(identical(mcmc_nngp_list$data$covariates$noise_X$arg, "No covariates were provided"))) stop("Covariates were provided for prediction, while no covariates were provided for fit")
  if(!is.null(X_noise_pred)) if(ncol(X_noise_pred)!=ncol(mcmc_nngp_list$data$covariates$noise_X$arg))stop("The Number of provided covariates does not match")
  # adding intercept to prediction regressors
  if(!is.null(X_noise_pred))X_noise_pred = cbind(rep(1, nrow(X_noise_pred)), X_noise_pred)
  if(is.null(X_noise_pred))X_noise_pred = matrix(1, nrow(predicted_locs), 1)
  colnames(X_noise_pred) = colnames(mcmc_nngp_list$data$covariates$noise_X$X)

  # GP if needed
  if(mcmc_nngp_list$hierarchical_model$noise_PP) 
  {
    # Dealing with duplicated locations
    duplicated_idx = duplicated(predicted_locs)
    order_idx = GpGp::order_maxmin(predicted_locs[!duplicated_idx,])
    predicted_locs_unique = predicted_locs[!duplicated_idx,][order_idx,]
    # surjection from duplicated locs to unique locs
    locs_match = match(split(predicted_locs, row(predicted_locs)), split(predicted_locs_unique, row(predicted_locs_unique)))
    # injection from unique locs to duplicated locs
    hctam_scol_1 = sapply(split(seq(nrow(predicted_locs)), locs_match), function(x)x[1])
    # constructing NNarray of observed and predicted locations
    ### locs = rbind(mcmc_nngp_list$data$locs, predicted_locs_unique)
    ### NNarray = rbind(mcmc_nngp_list$vecchia_approx$NNarray, GpGp::find_ordered_nn(locs, m = ncol(mcmc_nngp_list$vecchia_approx$NNarray)-1)[-seq(mcmc_nngp_list$vecchia_approx$n_locs),])
    ### NNarray_non_NA = !is.na(NNarray)
    ### sparse_chol_row_idx = row(NNarray)[NNarray_non_NA]
    ### sparse_chol_column_idx = NNarray[NNarray_non_NA]
    
    PP = mcmc_nngp_list$hierarchical_model$PP
    PP_ = mcmc_nngp_list$hierarchical_model$PP
    # appending predicted locations to NNarray
    NNarray_p = GpGp::find_ordered_nn(rbind(PP$knots, PP$unique_reordered_locs, predicted_locs_unique), m = ncol(PP_$NNarray)-1)
    NNarray_non_NA_p = !is.na(NNarray_p)
    sparse_chol_row_idx_p = row(NNarray_p)[NNarray_non_NA_p]
    sparse_chol_column_idx_p = NNarray_p[NNarray_non_NA_p]
    PP_$unique_reordered_locs = rbind(PP$unique_reordered_locs, predicted_locs_unique)
    PP_$sparse_chol =  Matrix::sparseMatrix(
      i = sparse_chol_row_idx_p, 
      j = sparse_chol_column_idx_p, 
      x = GpGp::vecchia_Linv(
        covparms = c(1, PP$matern_range, .0001),
        covfun_name = "matern15_isotropic",
        NNarray = NNarray_p, 
        locs = rbind(PP_$knots, PP_$unique_reordered_locs))[NNarray_non_NA_p],
      triangular = T
    )
    PP_$idx = c(PP$idx, nrow(PP$unique_reordered_locs)+locs_match)
  }
  # burn in
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
  # parallelization
  n_samples = length(kept_iterations)
  predicted_samples = lapply(
    X = mcmc_nngp_list$records, 
    FUN = 
      function(chain)
      {
        res = array(0, dim = c(nrow(X_noise_pred), 1, n_samples))
        # looping over saved observations
        for(i_predict in seq(n_samples))
        {
          res[,,i_predict] = Bidart::X_PP_mult_right(
            X = as.matrix(rbind(mcmc_nngp_list$data$covariates$noise_X$X, X_noise_pred)), 
            PP = PP_, 
            use_PP = mcmc_nngp_list$hierarchical_model$noise_PP, 
            locs_idx = seq(nrow(rbind(mcmc_nngp_list$data$covariates$noise_X$X, X_noise_pred))), 
            Y = chain$noise_beta[,,i_start + i_predict]
          )[-seq(mcmc_nngp_list$vecchia_approx$n_obs),]
        }
        return(res)
      })
  summaries = list()
  predicted_samples_ = list()
  predicted_samples_ = Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = predicted_samples)
  summaries = Bidart::get_array_summary(predicted_samples_)
  return(list("predicted_samples" = predicted_samples_, "summaries" = summaries, "predicted_locs" = predicted_locs))
}

#' Estimates the parameters from mcmc_nngp_list
#' @param mcmc_nngp_list A mcmc_nngp_list object generated using mcmc_nngp_initialize and ran using mcmc_nngp_run
#' @param burn_in MCMC burn-in
#' @param get_samples return Samples of the latent field 
#' @param lib.loc Location of libraries if installed in manually specified folder
#' @export
estimate_parameters = function(mcmc_nngp_list, burn_in = .5, get_samples = F, lib.loc = NULL)
{
  # burn in
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
  summaries = list()
  samples = list()
  n_cores = min(parallel::detectCores()-1, length(mcmc_nngp_list$records))
  cl = parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("kept_iterations"), envir = environment())
  # matched latent field 
  x = Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = parallel::parLapply(cl, mcmc_nngp_list$records, function(y)array(y$field[,,kept_iterations], dim = c(dim(y$field)[1], 1 , length(kept_iterations)))))
  summaries$field_at_observed_locs = Bidart::get_array_summary(x)[,mcmc_nngp_list$vecchia_approx$locs_match,,drop=F]
  if(get_samples)samples$field_at_observed_locs = x[mcmc_nngp_list$vecchia_approx$locs_match,,,drop=F]
  for(name in names(mcmc_nngp_list$records[[1]]))
  {
    parallel::clusterExport(cl, "name", envir = environment())
    x = Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = parallel::parLapply(cl, mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations)))))
    if(get_samples)samples[[name]] = x
    summaries[[name]] = Bidart::get_array_summary(x)
  }
  parallel::stopCluster(cl)
  list("summaries" = summaries, "samples" = samples)
}

#' Deviance Information Criterion
#' @param mcmc_nngp_list A mcmc_nngp_list object generated using mcmc_nngp_initialize and ran using mcmc_nngp_run
#' @param burn_in MCMC burn-in
#' @export
DIC = function(mcmc_nngp_list, burn_in = .5)
{
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  mean_field = rep(0, mcmc_nngp_list$vecchia_approx$n_locs)
  mean_noise = rep(0, mcmc_nngp_list$vecchia_approx$n_obs)
  mean_beta  = rep(0, ncol(mcmc_nngp_list$data$covariates$X$X))
  ll = c()
  for(i in kept_iterations)
  {
    for(j in seq(length(mcmc_nngp_list$states)))
    {
      noise = Bidart::variance_field(
        beta = mcmc_nngp_list$records[[j]]$noise_beta[,,i], 
        PP = mcmc_nngp_list$hierarchical_model$PP, use_PP = mcmc_nngp_list$hierarchical_model$noise_PP, 
        X = mcmc_nngp_list$data$covariates$noise_X$X, 
        locs_idx = NULL
        )
      field = mcmc_nngp_list$records[[j]]$field[,,i]
      beta = mcmc_nngp_list$records[[j]]$beta[,,i]
      ll = c(ll, 
             sum(dnorm(mcmc_nngp_list$data$observed_field - field[mcmc_nngp_list$vecchia_approx$locs_match]-mcmc_nngp_list$data$covariates$X$X%*%beta, 0, sqrt(noise), log = T))
             )
      mean_field = mean_field + field 
      mean_noise = mean_noise + noise 
      mean_beta  = mean_beta  + beta  
    }
  }
  mean_field = mean_field / (length(kept_iterations)*length(mcmc_nngp_list$states)) 
  mean_noise = mean_noise / (length(kept_iterations)*length(mcmc_nngp_list$states)) 
  mean_beta  = mean_beta  / (length(kept_iterations)*length(mcmc_nngp_list$states))
  ll_ = sum(dnorm(mcmc_nngp_list$data$observed_field - mean_field[mcmc_nngp_list$vecchia_approx$locs_match]-mcmc_nngp_list$data$covariates$X$X%*%mean_beta, 0, sqrt(mean_noise), log = T))
  -2 * (2 * mean(ll) - ll_)
}

#' Log-density of observations with respect to the model.
#' Gives a training score when the observations are used in the model.
#' Gives a validation score when the observations are not used in the model.
#' @param mcmc_nngp_list A mcmc_nngp_list object generated using mcmc_nngp_initialize and ran using mcmc_nngp_run
#' @param burn_in MCMC burn-in
#' @param latent_field_samples Samples of the latent field, either given by predict_latent_field or estimate_parameters
#' @param log_noise_samples Samples of the log noise variance, given by predict_noise
#' @param fixed_effects_samples Samples of the fixed effects, given by predict_fixed_effects
#' @export
log_score_Gaussian = function(observed_field, latent_field_samples, log_noise_samples, fixed_effects_samples)
{
  nsamples =  dim(latent_field_samples)[3]
  samples = dnorm(
    observed_field, 
    latent_field_samples + fixed_effects_samples, 
    exp(.5 * log_noise_samples)
  )
  log_score_per_obs = log(apply(samples, 1, mean))
  log_score_total = sum(log_score_per_obs)
  list(
    "per_obs" =   log_score_per_obs,
    "total" = log_score_total
  )
}

### #' Log-density of observations with respect to the model.
### #' Gives a training score when the observations are used in the model.
### #' Gives a validation score when the observations are not used in the model.
### #' @param mcmc_nngp_list A mcmc_nngp_list object generated using mcmc_nngp_initialize and ran using mcmc_nngp_run
### #' @param burn_in MCMC burn-in
### #' @param latent_field_samples Samples of the latent field, either given by predict_latent_field or estimate_parameters
### #' @param log_noise_samples Samples of the log noise variance, given by predict_noise
### #' @param fixed_effects_samples Samples of the fixed effects, given by predict_fixed_effects
### #' @export
### scoring = function(observed_field, latent_field_samples, log_noise_samples, fixed_effects_samples)
### {
###   nsamples =  dim(latent_field_samples)[3]
###   
###   response_pred = array(0, c(dim(latent_field_samples), ceiling(1000 / dim(latent_field_samples)[3])))
###   for(iii in seq(dim(response_pred)[4])){
###     response_pred[,,,iii] = latent_field_samples + fixed_effects_samples +  rnorm(length(latent_field_samples)) * exp(.5*log_noise_samples)
###   }
###   
###   cl = parallel::makeCluster(parallel::detectCores())
###   densities = parallel::parApply(cl, response_pred, 1, function(x)approxfun(density(x)))
###   densities_eval = log(mapply(function(z, dens)dens(z), z = observed_field, dens = densities))
###   remove(densities)
###   
###   list(
###   "per_obs" =   apply(samples, 1, mean), 
###   "total" = 
###     sum(
###       samples
###     )/nsamples
###   )
### }

