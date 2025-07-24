set.seed(123)
vecchia_approx = createVecchia(cbind(runif(1000), runif(1000)), 10, ncores = 1)

test_that("createPP works (automatic)", {
  set.seed(123)
  expect_message(
    pepito <- createPP(vecchia_approx, plot=TRUE), 
    "number of knots set to 100")
  expect_s3_class(pepito, "PP")
  expect_named(pepito, c("knots", "matern_range", "sparse_chol","n_knots", "vecchia_locs"))
  
  # Format
  expect_identical(sapply(pepito, function(x) class(x)[1]),
                   c(
                     "knots" = "matrix",
                     "matern_range" = "numeric",
                     "sparse_chol" = "dtCMatrix",
                     "n_knots" = "integer",
                     "vecchia_locs" = "matrix")
                   )
  expect_identical(dim(pepito$knots), c(100L, 2L))
  expect_identical(dim(pepito$sparse_chol), c(1100L, 1100L))
  expect_identical(dim(pepito$vecchia_locs), c(1000L, 2L))
  
  # values
  expect_identical(rownames(pepito$knots)[1:5], c("74", "50", "14", "12", "4"))
  expect_equal(colMeans(pepito$knots), c(0.478447, 0.457153), tolerance = 1e-5)
  expect_equal(sum(pepito$sparse_chol[1:10,1:10]), 4.2186, tolerance = 1e-5)
  expect_equal(pepito$matern_range, 0.287012, tolerance = 1e-5)
  expect_identical(pepito$n_knots, 100L)
  expect_identical(pepito$vecchia_locs, vecchia_approx$locs)
})


test_that("createPP works matern_range only", {
  set.seed(123)
  expect_message(
    pepito <- createPP(vecchia_approx, matern_range = .1, plot=FALSE), 
    "number of knots set to 100")
  expect_s3_class(pepito, "PP")
  expect_named(pepito, c("knots", "matern_range", "sparse_chol","n_knots", "vecchia_locs"))
  
  # values
  expect_identical(rownames(pepito$knots)[1:5], c("74", "50", "14", "12", "4"))
  expect_equal(colMeans(pepito$knots), c(0.478447, 0.457153), tolerance = 1e-5)
  expect_equal(sum(pepito$sparse_chol[1:10,1:10]), 8.00045, tolerance = 1e-5)
  expect_equal(pepito$matern_range, .1, tolerance = 1e-5)
})

test_that("createPP works matern_range and nb knots", {
  set.seed(123)
  expect_message(
    pepito <- createPP(vecchia_approx, knots = 50, matern_range = .1, plot=FALSE), 
    "knot placement done by default using k-means")
  expect_s3_class(pepito, "PP")
  # Format 
  expect_identical(dim(pepito$knots), c(50L, 2L))
  expect_identical(dim(pepito$sparse_chol), c(1050L, 1050L))
  
  # values
  expect_identical(rownames(pepito$knots)[1:5], c("34", "33", "8",  "18", "6"))
  expect_equal(colMeans(pepito$knots), c(0.496751, 0.469081), tolerance = 1e-5)
  expect_equal(sum(pepito$sparse_chol[1:10,1:10]), 7.91042, tolerance = 1e-5)
  expect_equal(pepito$matern_range, .1, tolerance = 1e-5)
})

test_that("createPP works knots as matrix", {
  set.seed(123)
  myknots <- as.matrix(expand.grid(seq(-.05, 1.05, .05), seq(-.05, 1.05, .05)))
  expect_message(
    pepito <- createPP(vecchia_approx, knots = myknots, matern_range = .1), 
    "This is great")
  expect_s3_class(pepito, "PP")
  # Format 
  expect_identical(dim(pepito$knots), dim(myknots))
  expect_identical(dim(pepito$sparse_chol), c(1529L, 1529L))
  
  # values
  expect_identical(rownames(pepito$knots)[1:5], c("383", "283", "412", "249", "300"))
  expect_equal(colMeans(pepito$knots), colMeans(myknots))
  expect_equal(sum(pepito$sparse_chol[1:10,1:10]), 6.9976, tolerance = 1e-5)
  expect_equal(pepito$matern_range, .1, tolerance = 1e-5)
})

test_that("var_loss_percentage.PP gives expected results", {
  set.seed(123)
  suppressMessages({
    pepito0 <- createPP(vecchia_approx, plot=FALSE)
    pepito1 <- createPP(vecchia_approx, plot=FALSE, matern_range = 0.1)
    pepito2 <- createPP(vecchia_approx, plot=FALSE, matern_range = 0.01)
  })
  
  # Var loss ok
  expect_message(
    res <- var_loss_percentage.PP(pepito0), 
    "This is great")
  expect_true(is(res, "numeric"))
  expect_equal(mean(res), 0.53340, tolerance = 1e-5)
  
  # Var loss almost ok
  expect_message(
    res <- var_loss_percentage.PP(pepito1), 
    "This is fairly good")
  expect_equal(mean(res), 5.61732, tolerance = 1e-5)
  
  # Var loss bad
  expect_message(
    res <- var_loss_percentage.PP(pepito2), 
    "his is quite a bit of loss")
  expect_equal(mean(res), 85.522, tolerance = 1e-5)
})

test_that("summary.PP gives expected results", {
  set.seed(123)
  pepito <- createPP(vecchia_approx, plot=FALSE)
  expect_output(summary(pepito),
                 "Object of class 'PP' with 100 knots, based on 1000 locations, matérn range = 0.2870118")
})

set.seed(123)
PP = createPP(vecchia_approx, plot=FALSE)
X = matrix(rnorm(500), 1000)
Y = matrix(rnorm(30*nrow(X)), nrow(X))

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
