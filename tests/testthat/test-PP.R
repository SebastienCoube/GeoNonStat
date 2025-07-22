set.seed(123)
vecchia_approx = createVecchia(observed_locs  = cbind(runif(1000), runif(1000)), 10)

test_that("createPP works (automatic)", {
  set.seed(123)
  expect_message(
    pepito <- createPP(vecchia_approx), 
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
  expect_identical(dim(pepito$vecchia_locs), c(1000, 2))
  
  # values
  expect_identical(rownames(pepito$knots)[1:5], c("74", "50", "12", "96", "4"))
  expect_equal(colMeans(pepito$knots), c(0.4948819, 0.4584636), tolerance = 1e-5)
  expect_equal(sum(pepito$sparse_chol[1:10,1:10]), 3.752212, tolerance = 1e-5)
  expect_equal(pepito$matern_range, 0.2581213, tolerance = 1e-5)
  expect_identical(pepito$n_knots, 100L)
  expect_identical(pepito$vecchia_locs, vecchia_approx$locs)
})


test_that("createPP works matern_range only", {
  set.seed(123)
  expect_message(
    pepito <- createPP(vecchia_approx, matern_range = .1), 
    "number of knots set to 100")
  expect_s3_class(pepito, "PP")
  expect_named(pepito, c("knots", "matern_range", "sparse_chol","n_knots", "vecchia_locs"))
  
  # values
  expect_identical(rownames(pepito$knots)[1:5], c("74", "50", "12", "96", "4"))
  expect_equal(colMeans(pepito$knots), c(0.4948819, 0.4584636), tolerance = 1e-5)
  expect_equal(sum(pepito$sparse_chol[1:10,1:10]), 7.787548, tolerance = 1e-5)
  expect_equal(pepito$matern_range, .1, tolerance = 1e-5)
})

test_that("createPP works matern_range and nb knots", {
  set.seed(123)
  expect_message(
    pepito <- createPP(vecchia_approx, knots = 50, matern_range = .1), 
    "knot placement done by default using k-means")
  expect_s3_class(pepito, "PP")
  # Format 
  expect_identical(dim(pepito$knots), c(50L, 2L))
  expect_identical(dim(pepito$sparse_chol), c(1050L, 1050L))
  
  # values
  expect_identical(rownames(pepito$knots)[1:5], c("34", "33", "8",  "18", "6"))
  expect_equal(colMeans(pepito$knots), c(0.5000674, 0.4690878), tolerance = 1e-5)
  expect_equal(sum(pepito$sparse_chol[1:10,1:10]), 7.86248, tolerance = 1e-5)
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
  expect_equal(mean(res), 0.5867833, tolerance = 1e-5)
  
  # Var loss almost ok
  expect_message(
    res <- var_loss_percentage.PP(pepito1), 
    "This is fairly good")
  expect_equal(mean(res), 4.618685, tolerance = 1e-5)
  
  # Var loss bad
  expect_message(
    res <- var_loss_percentage.PP(pepito2), 
    "his is quite a bit of loss")
  expect_equal(mean(res), 83.91222, tolerance = 1e-5)
})

test_that("summary.PP gives expected results", {
  set.seed(123)
  pepito <- createPP(vecchia_approx, plot=FALSE)
  expect_output(summary(pepito),
                 "Object of class 'PP' with 100 knots, based on 1000 locations, matérn range = 0.2870118")
})
