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

test_that("variance_field produce expected output", {
  set.seed(123)
  locs = cbind(runif(100), runif(100))
  n_PP = 50
  PP = get_PP(locs, c(1, .1, 1.5, 0), n_PP = n_PP, m = 15)
  X = matrix(rnorm(10*nrow(PP$unique_reordered_locs)), ncol = 10)

  # When using PP
  expect_error(
    res <- variance_field(beta = rnorm(n_PP), 
                          PP = PP, 
                          use_PP = TRUE, 
                          X = X),
    NA
  )
  
  expect_type(res, "double")
  expect_length(res, 100)
  expect_equal(mean(res), 50.4898)
  expect_equal(res[1], 64.84969)
  
  # When not using PP
  expect_error(
    res <- variance_field(beta = rnorm(n_PP), 
                          PP = PP, 
                          use_PP = FALSE, 
                          X = X),
    NA
  )
  expect_type(res, "double")
  expect_length(res, 100)
  expect_equal(mean(res), 67.298166)
  expect_equal(res[1], 6.136415)
})

test_that("compute_sparse_chol produce expected output", {
  set.seed(123)
  locs = cbind(seq(100)/10, 0)
  NNarray = GpGp::find_ordered_nn(locs, 10)
  
  expect_error(
    res <- compute_sparse_chol(
      range_beta = matrix(.5/sqrt(2),1,1), 
      NNarray = NNarray, 
      locs = locs,
      nu = 2
    ),
    "nu must be equal to 0.5 or 1.5"
  )
  
  expect_error(
    res <- compute_sparse_chol(
              range_beta = matrix(.5/sqrt(2),1,1), 
              NNarray = NNarray, 
              locs = locs,
              use_PP = F, 
              num_threads = 1, 
              anisotropic = F,
              range_X = matrix(1, nrow(locs), 1), 
              nu = 1.5
            ),
    NA
  )
  expect_type(res, "list")
  expect_length(res, 2)
  
  expect_true(is(res[[1]], "array"))
  expect_identical(dim(res[[1]]), c(100L, 11L))
  expect_equal(mean(res[[1]]), 0.01397719)
  expect_equal(res[[1]][2,1:3], c(13.8705735, -13.8310216, 0.00000))
  
  expect_type(res[[2]], "list")
  expect_type(res[[2]][[1]], "array") 
  expect_identical(dim(res[[2]][[1]]), c(100L, 11L, 11L))
  expect_equal(mean(res[[2]][[1]]), -0.00044410, tolerance = 1e-6)
  expect_equal(res[[2]][[1]][2,1:3,1], c(6.11473435, 6.11473435, 0.000000))
})

test_that("get_PP produce expected output", {
  obs_locs <- matrix(rnorm(10), ncol=2)
  expect_error(
    res <- get_PP(obs_locs, matern_range=c(1, 1.1, 1.5, 0), n_PP=4),
    NA
  )
  expect_named(res, c('knots', 'unique_reordered_locs', 'idx', 
                      'lonlat', 'm', 'matern_range', 
                      'sparse_chol', 'NNarray', 'n_PP'))
  
})