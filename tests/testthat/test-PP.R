#' observed_locs = cbind(runif(1000), runif(1000))
#' observed_locs = observed_locs[ceiling(nrow(observed_locs)*runif(3000)),]
#' # automatic 
#' pepito = PP(observed_locs)
#' # choosing manually Matérn range, too small wrt number of knots
#' pepito = PP(observed_locs, matern_range = .1)
#' # choosing manually Matérn range, way too small wrt number of knots
#' pepito = PP(observed_locs, matern_range = .01)
#' # choosing manually number of knots in order to adjust to Matérn range
#' pepito = PP(observed_locs, knots = 200, matern_range = .1)
#' # choosing manually number of knots, but picking too few for default Matérn range
#' pepito = PP(observed_locs, knots = 10)
#' # choosing manually Matérn range in order to adjust to the number of knots
#' pepito = PP(observed_locs, knots = 10, matern_range = .5)
#' # inputing an user-specified grid of knots
#' pepito = PP(observed_locs, knots = as.matrix(expand.grid(seq(-.05, 1.05, .1), seq(-.05, 1.05, .1))))
#' # inputing an user-specified grid of knots in order to adjust to small Matérn range
#' pepito = PP(observed_locs, knots = as.matrix(expand.grid(seq(-.05, 1.05, .1), seq(-.05, 1.05, .1))), matern_range = .1)

set.seed(123)
observed_locs = cbind(runif(1000), runif(1000))
observed_locs = observed_locs[ceiling(nrow(observed_locs)*runif(3000)),]

test_that("PP works (automatic", {
  set.seed(123)
  expect_message(
    pepito <- PP(observed_locs), 
    "number of knots set by default to 100")
  expect_true(is(pepito, "PP"))
  expect_named(pepito, c('knots', 'unique_reordered_locs', 'idx', 
                         'm', 'matern_range', 'sparse_chol', 'NNarray', 'n_knots'))
  
  # Format
  expect_identical(sapply(pepito, function(x) class(x)[1]),
                   c(
                     "knots" = "matrix",
                     "unique_reordered_locs" = "matrix",
                     "idx" = "integer",
                     "m" = "numeric",
                     "matern_range" = "numeric",
                     "sparse_chol" = "dtCMatrix",
                     "NNarray" = "matrix",
                     "n_knots" = "integer")
                   )
  expect_identical(dim(pepito$knots), c(100L, 2L))
  expect_identical(dim(pepito$unique_reordered_locs), c(951L, 2L))
  expect_identical(length(pepito$idx), 3000L)
  expect_identical(dim(pepito$NNarray), c(1051L, 11L))
  
  # values
  expect_identical(rownames(pepito$knots)[1:5], c("10", "95", "85", "57", "42"))
  expect_equal(mean(pepito$knots[,1]), 0.5299115, tolerance = 1e-5)
  expect_equal(mean(pepito$knots[,2]), 0.4845792, tolerance = 1e-5)
  
  expect_equal(mean(pepito$unique_reordered_locs[,1]), 0.4964153, tolerance = 1e-5)
  expect_equal(mean(pepito$unique_reordered_locs[,2]), 0.4976046, tolerance = 1e-5)
  
  expect_equal(sum(pepito$sparse_chol[1:10,1:10]), 3.47036, tolerance = 1e-5)
  expect_equal(mean(pepito$NNarray[,2], na.rm = TRUE), 232.0133, tolerance = 1e-5)
  expect_identical(pepito$idx[4], 550L)
  expect_identical(pepito$m, 10)
  expect_equal(pepito$matern_range, 0.268488, tolerance = 1e-5)
  expect_identical(pepito$n_knots, 100L)
})
