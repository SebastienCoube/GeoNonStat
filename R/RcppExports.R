# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

nonstat_vecchia_Linv <- function(log_range, covfun_name, sphere, locs, NNarray, num_threads, compute_derivative) {
    .Call('_GeoNonStat_nonstat_vecchia_Linv', PACKAGE = 'GeoNonStat', log_range, covfun_name, sphere, locs, NNarray, num_threads, compute_derivative)
}

derivative_sandwich <- function(derivative, left_vector, right_vector, NNarray) {
    .Call('_GeoNonStat_derivative_sandwich', PACKAGE = 'GeoNonStat', derivative, left_vector, right_vector, NNarray)
}

log_determinant_derivative <- function(derivative, compressed_sparse_chol, NNarray) {
    .Call('_GeoNonStat_log_determinant_derivative', PACKAGE = 'GeoNonStat', derivative, compressed_sparse_chol, NNarray)
}

nonstat_covmat <- function(log_range, covfun_name, locs) {
    .Call('_GeoNonStat_nonstat_covmat', PACKAGE = 'GeoNonStat', log_range, covfun_name, locs)
}

