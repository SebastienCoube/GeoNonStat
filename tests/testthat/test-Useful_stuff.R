set.seed(123)
vecchia_approx = createVecchia(cbind(runif(1000), runif(1000)), 10, ncores = 1)

set.seed(123)
PP = createPP(vecchia_approx, plot=FALSE)
X = matrix(rnorm(500), 1000)
Y = matrix(rnorm(30*nrow(X)), nrow(X))


test_that("derivative_sandwiches produce expected output", {
  # TODO
})

test_that("log_determinant_derivatives produce expected output", {
  # TODO
})

test_that("expmat produce expected output", {
  coords <- c(1,2,3,4,5,6)
  expect_error(
    res <- expmat(coords), 
    NA
  )
  expect_true(is(res, "array"))
  expect_identical(dim(res), c(3L,3L))
  
  expected_mat <- matrix(
    c(47859.7700011562, 55641.9937959158, 62462.1209161336, 
      55641.9937959158, 64689.8331149192, 72618.8933675331, 
      62462.1209161336, 72618.8933675331, 81519.8962286822),
    nrow=3)
  expect_equal(res, expected_mat)
  expect_error(
    res <- expmat(coords, eps=0), 
    NA
  )
  diag(expected_mat) <- diag(expected_mat) - 0.0001
  expect_equal(res, expected_mat)
})

test_that("symmat produce expected output", {
  # 1x1 matrix
  coords <- c(1)
  expect_error(
    res <- symmat(coords), 
    NA
  )
  expect_true(is(res, "matrix"))
  expect_identical(dim(res), c(1L,1L))
  expect_equal(res, matrix(c(1)))
  
  # 2x2 matrix
  coords <- c(1,2,3)
  expect_error(
    res <- symmat(coords), 
    NA
  )
  expect_true(is(res, "matrix"))
  expect_identical(dim(res), c(2L,2L))
  expect_equal(res, 
               matrix(c(1,3,3,2),
                      nrow=2))
  
  # 3x3 matrix
  coords <- c(1,2,3,4,5,6)
  expect_error(
    res <- symmat(coords), 
    NA
  )
  expect_true(is(res, "matrix"))
  expect_identical(dim(res), c(3L,3L))
  expect_equal(res, 
               matrix(c(1,4,5,4,2,6,5,6,3),
                      nrow=3))
  
  # Vector incompatible
  coords <- c(1,2,3,4,5)
  expect_error(
    res <- symmat(coords), 
    "length of coords incompatible with a symetric matrix"
  )
})


# 
# test_that("variance_field produce expected output", {
#   # Here we test only structure since results are directly extracted from X_PP_mult_right
#   set.seed(123)
#   locs = cbind(runif(100), runif(100))
#   n_PP = 50
#   PP = createPP(locs, c(1, .1, 1.5, 0), knots = n_PP, m = 15)
#   X = matrix(rnorm(10*nrow(PP$unique_reordered_locs)), ncol = 10)
# 
#   # When using PP
#   expect_error(
#     res <- variance_field(beta = rnorm(n_PP), 
#                           PP = PP, 
#                           use_PP = TRUE, 
#                           X = X),
#     NA
#   )
#   
#   expect_true(is(res, "numeric"))
#   expect_length(res, 100)
#   
#   # When not using PP
#   expect_error(
#     res <- variance_field(beta = rnorm(n_PP), 
#                           PP = PP, 
#                           use_PP = FALSE, 
#                           X = X),
#     NA
#   )
#   expect_true(is(res, "numeric"))
#   expect_length(res, 100)
# })

# test_that("compute_sparse_chol produce expected output", {
#   set.seed(123)
#   locs = cbind(seq(100)/10, 0)
#   NNarray = GpGp::find_ordered_nn(locs, 10)
#   
#   expect_error(
#     res <- compute_sparse_chol(
#       range_beta = matrix(.5/sqrt(2),1,1), 
#       NNarray = NNarray, 
#       locs = locs,
#       nu = 2
#     ),
#     "nu must be equal to 0.5 or 1.5"
#   )
#   
#   set.seed(123)
#   expect_error(
#     res <- compute_sparse_chol(
#               range_beta = matrix(.5/sqrt(2),1,1), 
#               NNarray = NNarray, 
#               locs = locs,
#               use_PP = F, 
#               num_threads = 1, 
#               anisotropic = F,
#               range_X = matrix(1, nrow(locs), 1), 
#               nu = 1.5
#             ),
#     NA
#   )
#   expect_type(res, "list")
#   expect_length(res, 2)
#   
#   expect_true(is(res[[1]], "array"))
#   expect_identical(dim(res[[1]]), c(100L, 11L))
#   expect_equal(mean(res[[1]]), 0.01397719, tolerance = 1e-7)
#   expect_equal(res[[1]][2,1:3], c(13.8705735, -13.8310216, 0.00000))
#   
#   expect_type(res[[2]], "list")
#   expect_true(is(res[[2]][[1]], "array")) 
#   expect_identical(dim(res[[2]][[1]]), c(100L, 11L, 11L))
#   expect_equal(mean(res[[2]][[1]]), -0.00044410, tolerance = 1e-6)
#   expect_equal(res[[2]][[1]][2,1:3,1], c(6.11473435, 6.11473435, 0.000000))
# })


test_that("beta_prior_log_dens produce expected output", {
  # TODO
})

test_that("beta_prior_log_dens_derivative produce expected output", {
  # TODO
})


test_that("derivative_chol_expmat produce expected output", {
  expect_error( 
    res <- derivative_chol_expmat(c(1,2,3, 3,2,4)),
    NA
  )
  expect_true(is(res, "array"))
  expect_identical(dim(res), c(3L,3L,6L))
  expect_equal(
    apply(res, 1, mean),
    c(8.5554350147, -0.0108160777, 0.0003808579)
  )
  expect_equal(
    apply(res, 2, mean),
    c(2.343741568, 3.071973047, 3.129285180)
  )
})

test_that("derivative_field_wrt_scale produce expected output", {
  field <- array(c(1,2,3), dim=c(1,3))
  coords <- c(1,2,3, 3,2,4)
  expect_error( 
    res <- derivative_field_wrt_scale(field, coords),
    NA
  )
  expect_true(is(res, "array"))
  expect_identical(dim(res), c(1L,3L,6L))
  expect_equal(
    apply(res, 1, mean),
    c(0.312839647)
  )
  expect_equal(
    apply(res, 2, mean),
    c(0.261804197, 0.335998026, 0.340716718)
  )
})



test_that("X_PP_crossprod Simple crossprod if no PP, vecchia unused", {
  expect_error(
    res1 <- X_PP_crossprod(X = X, PP = NULL, Y = Y, vecchia_approx = vecchia_approx),
    NA
  )
  expect_true(is(res1, "matrix"))
  expect_identical(dim(res1), c(1L, 30L))
  expect_identical(res1,crossprod(X, Y))
  
  # Simple corssprod, vecchia unused. 
  expect_error(
    res2 <- X_PP_crossprod(X = X, PP = NULL, Y = Y, vecchia_approx = NULL),
    NA
  )
  expect_identical(res1,res2)
})

test_that("X_PP_crossprod with PP", {
  expect_error(
    res2 <- X_PP_crossprod(X = X, PP = PP, Y = Y, vecchia_approx = vecchia_approx),
    NA
  )
  expect_true(is(res2, "matrix"))
  expect_identical(dim(res2), c(101L, 30L))
  
  # First line is crossprod
  expect_identical(res2[1,],crossprod(X, Y)[1,])
  expect_equal(colMeans(res2)[1:5], 
               c(0.1999806, -0.8748530, 0.9684133, -0.7444476, 0.1350033),
               tolerance = 1e-5)
})

test_that("X_PP_crossprod with PP and permute obs", {
  expect_error(
    res3 <- X_PP_crossprod(X = X, PP = PP, Y = Y, vecchia_approx = vecchia_approx, permutate_PP_to_obs = TRUE),
    NA
  )
  expect_true(is(res3, "matrix"))
  expect_identical(dim(res3), c(101L, 30L))
  
  # First line is crossprod
  expect_identical(res3[1,],crossprod(X, Y)[1,])
  expect_equal(colMeans(res3)[1:5], 
               c(0.02263887, -1.04065836, 0.84462225, -0.70220738, 0.36781311),
               tolerance = 1e-5)
})

# res1 <- X_PP_mult_right(X = X, Y = covariate_coefficients, vecchia_approx = vecchia_approx)
test_that("X_PP_mult_right with PP and permute obs", {
  set.seed(123)
  expect_error(
    resmr <- X_PP_mult_right(X = X, PP = PP, Y = Y, vecchia_approx = vecchia_approx, permutate_PP_to_obs = FALSE),
    "Y should have"
  )
  # No PP
  expect_error(
    resmr <- X_PP_mult_right(X = X, PP = NULL, Y = matrix(c(0.5, 1, 1.1, 0.2), nrow = 1), vecchia_approx = vecchia_approx, permutate_PP_to_obs = FALSE),
    NA
  )
  expect_true(is(resmr, "matrix"))
  expect_identical(dim(resmr), c(1000L, 4L))
  expect_equal(mean(resmr), 0.01392224, tolerance = 1e-5)
  
  # No X
  expect_error(
    resmr <- X_PP_mult_right(X = NULL, PP = PP, Y = matrix(rnorm(200), nrow = 100), vecchia_approx = vecchia_approx, permutate_PP_to_obs = FALSE),
    NA
  )
  expect_true(is(resmr, "Matrix"))
  expect_identical(dim(resmr), c(1000L, 2L))
  mean(as.matrix(resmr))
  expect_equal(mean(as.matrix(resmr)), -0.1908332, tolerance = 1e-5)
  
  # X and PP
  expect_error(
    resmr <- X_PP_mult_right(X = X, PP = PP, Y =  matrix(rnorm(202), nrow = 101), vecchia_approx = vecchia_approx, permutate_PP_to_obs = FALSE),
    NA
  )
  expect_true(is(resmr, "Matrix"))
  expect_identical(dim(resmr), c(1000L, 2L))
  expect_equal(mean(as.matrix(resmr)), 0.01075514, tolerance = 1e-5)
  
})


