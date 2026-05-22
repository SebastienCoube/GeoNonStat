
extendVecchia <- function(vecchia_approx, new_locs){
   new_vecchia_approx <- list()
   new_vecchia_approx$new_observed_locs <- new_locs
   new_vecchia_approx$observed_locs <- rbind(vecchia_approx$observed_locs, new_locs)
   #new_vecchia_approx$observed_locs[seq(vecchia_approx$n_obs),] - vecchia_approx$observed_locs
   #new_vecchia_approx$observed_locs[-seq(vecchia_approx$n_obs),] - new_locs
   new_vecchia_approx$previous_n_obs <- vecchia_approx$n_obs
   new_vecchia_approx$previous_n_locs <-  vecchia_approx$n_locs 
   new_vecchia_approx$n_obs <- nrow(new_vecchia_approx$observed_locs)
   new_vecchia_approx$locs <- unique(rbind(vecchia_approx$locs, new_locs))
   #new_vecchia_approx$locs[seq(vecchia_approx$n_locs),] - vecchia_approx$locs
   new_vecchia_approx$new_locs <- new_vecchia_approx$locs[-seq(vecchia_approx$n_locs),]
   new_vecchia_approx$n_locs <- nrow(new_vecchia_approx$locs)
   new_vecchia_approx$t_locs <- t(new_vecchia_approx$locs)
   new_vecchia_approx$locs_match <- match(split(new_vecchia_approx$observed_locs, row(new_vecchia_approx$observed_locs)), split(new_vecchia_approx$locs, row(new_vecchia_approx$locs)))
   new_vecchia_approx$locs_match_matrix <- Matrix::sparseMatrix(i = new_vecchia_approx$locs_match, j = seq_len(new_vecchia_approx$n_obs), x = 1)
   # doing reversed operation : for a given unrepeated location, tell which observations correspond
   new_vecchia_approx$hctam_scol <- split(seq_len(new_vecchia_approx$n_obs), new_vecchia_approx$locs_match)
   # extracting NNarray =  nearest neighbours for Vecchia approximation
   unobserved_new_locs = new_vecchia_approx$locs[-seq(vecchia_approx$n_locs),]
   new_vecchia_approx$NNarray <- cbind(
     vecchia_approx$NNarray,  
     t(cbind(
       seq(vecchia_approx$n_locs + 1, new_vecchia_approx$n_locs),
       FNN::get.knnx(
         data = vecchia_approx$locs, query = unobserved_new_locs, 
         k = nrow(vecchia_approx$NNarray)-1)$nn.index
     ))
   )
   new_vecchia_approx$NNNoNA <- !is.na(new_vecchia_approx$NNarray)
   
   sparse_mat <- Matrix::sparseMatrix(
     x = seq_len(sum(new_vecchia_approx$NNNoNA)),
     i = col(new_vecchia_approx$NNarray )[new_vecchia_approx$NNNoNA],
     j = new_vecchia_approx$NNarray [new_vecchia_approx$NNNoNA],
     triangular = TRUE
   )
   new_vecchia_approx$sparse_chol_x_reorder <- seq_along(new_vecchia_approx$NNarray)[new_vecchia_approx$NNNoNA][match(sparse_mat@x, seq_len(sum(new_vecchia_approx$NNNoNA)))]
   return(new_vecchia_approx) 
 }

#extended_PP = extendPP(PP, extended_vecchia_approx)
extendPP <- function(PP, extended_vecchia_approx){
  if(is.null(PP))return(NULL)
  createPP(
    vecchia_approx = extended_vecchia_approx, 
    matern_range = PP$matern_range, knots = PP$knots, reorder_knots = F,
    plot = F)
}


predict1Noise <- function(noise_beta, extended_PP_noise, extended_vecchia_approx, new_noise_X){
  PP_effect <- xPPMultRight(X = NULL, PP = extended_PP_noise, 
                           vecchia_approx = extended_vecchia_approx, 
                           permutate_PP_to_obs = T, 
                           Y = noise_beta[-seq_len(ncol(new_noise_X))])
  
  PP_effect = PP_effect[-seq_len(extended_vecchia_approx$previous_n_obs)]
  X_effect <- new_noise_X %*% noise_beta[seq_len(ncol(new_noise_X))]
  res <- PP_effect + X_effect
  return(res)
}

predictNoise <- function(geo_non_stat, extended_vecchia_approx, new_noise_X, burn_in){
  extended_PP_noise <- extendPP(geo_non_stat$hierarchical_model$noise$PP)
  new_noise_X <- cbind(rep(1, extended_vecchia_approx$n_obs - 
                             extended_vecchia_approx$previous_n_obs), 
                       new_noise_X)
  # apply to MCMC samples of noise_beta from geo_non_stat
}

predict1Field <- function(range_beta, 
                          field, 
                          extended_PP_range, 
                          extended_vecchia_approx, 
                          complete_range_X_locs, 
                          hierarchical_model, 
                          num_threads){
  compressed_chol <- 
    array(0, c(nrow(extended_vecchia_approx$NNarray), extended_vecchia_approx$n_locs, 1))
  log_range <- computeLogRange(
    range_beta = range_beta,
    PP = extended_PP_range,
    vecchia_approx = extended_vecchia_approx,
    range_X = complete_range_X_locs
  )
  vecchia_(
    log_range = t(log_range), locs = extended_vecchia_approx$t_locs,
    NNarray = extended_vecchia_approx$NNarray,
    smoothness = hierarchical_model$matern_smoothness,
    compute_derivative = FALSE, num_threads = num_threads, result = compressed_chol
  )
  pred_field_sample <- 
    (
      -GpGp::Linv_mult(
        Linv = t(compressed_chol[,,1]), 
        z = c(field, rep(0, nrow(extended_vecchia_approx$new_locs))), 
        NNarray = t(extended_vecchia_approx$NNarray))[-seq(extended_vecchia_approx$previous_n_locs)] 
      + rnorm(extended_vecchia_approx$n_locs - extended_vecchia_approx$previous_n_locs)
    ) / compressed_chol[1,-seq(extended_vecchia_approx$previous_n_locs),1]
  res = list(
    "log_range" = log_range[-seq(extended_vecchia_approx$previous_n_locs),], 
    "field" = pred_field_sample
  )
  
  return(res)
}


#' Predict a 'GeoNonStat' object on new locations
#'
#' @param object an object of class \code{GeoNonStat}
#' @param new_locs new locations
#' @param new_X new X 
#' @param new_noise_X new noise X
#' @param new_range_X new range X
#' @param num_threads number of threads. Integer.
#' @param ... additional arguments (unused)
#' @rdname GeoNonStat
#' @export
#' @method predict GeoNonStat
predict.GeoNonStat <-  function(object, new_locs, new_X, new_noise_X, new_range_X, num_threads, ...) {
  # Checks 
  # new_X au même format que X
  # new_noise_X au même format que noise_X
  # new_range_X au même format que range_X
  
  extended_vecchia_approx <- extendVecchia(object$vecchia_approx, new_locs)
  
  extended_PP_noise <- extendPP(object$hierarchical_model$noise$PP, extended_vecchia_approx)
  extended_PP_range <- extendPP(object$hierarchical_model$range$PP, extended_vecchia_approx) 
  
  estimated_beta<- aggregateRecords(object, keep="noise_beta", keep_separate_chains = TRUE)
  # estimated_field_range<- aggregateRecords(object, keep=c("range_beta", "field"), 
                                           # keep_separate_chains = TRUE)
  
  predicted_noises <- lapply(
    object$states, function(x) {
      apply(x[["noise_beta"]], 1, function(y) {
        predict1Noise(y, extended_PP_noise, extended_vecchia_approx, new_noise_X)
      })
    }
  )
  predicted_field <- lapply(
    object$states, function(x) {
        predict1Field(x$params$range_beta,
                      x$params$field, 
                      extended_PP_range, 
                      extended_vecchia_approx, 
                      new_range_X, 
                      object$hierarchical_model,
                      num_threads)
    }
  )
  
  # Attention à l'ordre des locs pour représenter après. Il faudra sûrement les sortir.
  
  # predicted_field = predict1Field(
  #   range_beta = range_beta, field = field, 
  #   extended_PP_range = extended_PP_range, 
  #   extended_vecchia_approx, 
  #   complete_range_X_locs, geo_non_stat$hierarchical_model, 
  #   num_threads = 8)
  
  # predictedY <- lapply( 
  #                        )
  
  # print(summary des predictions)
  # return( predicted values)
  # A voir si on retourne le summary. 
  # Sebastien est pour. 
  
  return(list("predicted_noises" = predicted_noises, "predicted_field" = predicted_field))
}
