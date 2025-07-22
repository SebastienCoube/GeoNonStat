set.seed(123)
obs_locs <- matrix(rnorm(200), ncol=2)
vecchia_approx = createVecchia(obs_locs, ncores=1)
pepito = createPP(vecchia_approx, plot=FALSE)

test_that("plot_knots.PP doesn't produce errors", {
  expect_error(
    res <- plot_knots.PP(pepito),
    NA
  )
  mar_var = var_loss_percentage.PP(pepito)
  expect_error(
    res <- plot_knots.PP(pepito, mar_var_loss = mar_var),
    NA
  )
  expect_error(
    res <- plot_knots.PP(pepito, mar_var_loss = mar_var, show_knots = FALSE),
    NA
  )
})

test_that("plot.PP doesn't produce errors", {
  expect_message(
    res <- plot.PP(pepito, mar_var_loss = FALSE),
    NA
  )
  expect_message(
    res <- plot.PP(pepito, mar_var_loss = TRUE),
    "0.2% of marginal variance"
  )
  expect_message(
    res <- plot.PP(pepito, mar_var_loss = TRUE, separate=TRUE),
    "0.2% of marginal variance"
  )
  mar_var = var_loss_percentage.PP(pepito)
  expect_error(
    res <- plot_knots.PP(pepito, mar_var_loss = mar_var),
    NA
  )
  expect_error(
    res <- plot_knots.PP(pepito, mar_var_loss = mar_var, show_knots = FALSE),
    NA
  )
})
# 
# test_that("compare_PP_NNGP produce expected output", {
#   obs_locs <- matrix(rnorm(20), ncol=2)
#   ex_PP <- PP(obs_locs, matern_range=c(1, 1.1, 1.5, 0), knots=4)
#   op <- par("mfrow")
#   expect_error(
#     res <- compare_PP_NNGP(ex_PP),
#     NA
#   )
#   # Test that par is correctly reset
#   testthat::expect_identical(op, par("mfrow"))
# })
# 
# test_that("plot_ellipses don't produce errors", {
#   locs <- matrix(rnorm(20), ncol=2)
#   log_range <- matrix(rnorm(30), ncol=3)
#   if(!is.null(dev.list())) dev.off()
#   expect_error(
#     plot_ellipses(locs, log_range, shrink=0.1, add=TRUE),
#     "plot.new has not been called yet"
#   )
#   expect_error(
#     plot_ellipses(locs, log_range, shrink=0.1),
#     NA
#   )
# })
# 
# test_that("get_colors produce expected output", {
#   expect_error(
#     res<- get_colors(c(1,2,3)), 
#           NA
#   )
#   expect_identical(res, c("#FF0000", "#FF9B00", "#FFFF9E"))
# })
# 
# test_that("plot_pointillist_painting produce expected errors", {
#   locs=c(2,3,6)
#   field=c(1,2,3)
#   if(!is.null(dev.list())) dev.off()
#   expect_error(
#     plot_pointillist_painting(locs, field, add=TRUE),
#     "plot.new has not been called yet"
#   )
#   expect_error(
#     plot_pointillist_painting(locs, field),
#     NA
#   )
# })
# 
# test_that("pointillist_colorscale produce expected output", {
#   op <- par("mfrow")
#   expect_error(
#     pointillist_colorscale(c(1,2,3)),
#     NA
#   )
#   # Test that par is correctly reset
#   testthat::expect_identical(op, par("mfrow"))
# })
# 
# test_that("array_cumsum produce expected output", {
#   myarray <- array(abs(rnorm(8)), dim = c(2, 2, 2))
#   expect_error(
#     res <- array_cumsum(myarray),
#     NA
#   )
#   expect_true(class(res)=="array")
#   expect_identical(dim(res), dim(myarray))
#   expect_equivalent(res[,,1],
#                     array(c(0.33867004,0.05762468 ,0.195550,1.094133), dim=c(2,2)))
#   expect_equivalent(res[,,2],
#                     array(c(0.8683650, 0.2959837, 0.2992697, 1.5791371), dim=c(2,2)))
#   
#   myarray2 <- array(abs(rnorm(100)), dim = c(10, 10, 2))
#   expect_error(
#     res2 <- array_cumsum(myarray2),
#     NA
#   )
#   expect_true(class(res2)=="array")
#   expect_identical(dim(res2), dim(myarray2))
#   
# })
# 
# 
# test_that("array_multiply_3 produce expected output", {
#   x <- array(1:12, dim = c(2, 2, 3))
#   M <- Matrix::Matrix(matrix(1:9, nrow = 3, ncol = 3), sparse = TRUE)
#   expect_error(
#     res <- array_multiply_3(x, M),
#     NA
#   )
#   expect_identical(dim(res), dim(x))
#   expect_identical(res[,,1], array(c(38,44,50,56), dim=c(2,2)))
#   expect_identical(res[,,2], array(c(83,98,113,128), dim=c(2,2)))
#   expect_identical(res[,,3], array(c(128,152,176,200), dim=c(2,2)))
# })
# 
# test_that("array_multiply_2 produce expected output", {
#   x <- array(1:24, dim = c(2, 2, 4))
#   M <- matrix(1:4, nrow = 2, ncol = 2)
#   expect_error(
#     res <- array_multiply_2(x, M),
#     NA
#   )
#   # expect_identical(dim(res), dim(x))
# })
# 
# test_that("array_multiply_1 produce expected output", {
#   M <- matrix(1:4, nrow = 2, ncol = 2)
#   x <- array(1:8, dim = c(2, 2, 2)) 
#   expect_error(
#     res <- array_multiply_1(x, M),
#     NA
#   )
#   expect_identical(dim(res), c(2L,2L,2L))
#   expect_identical(res[,,1], array(c(7, 10, 15, 22), dim=c(2,2)))
#   expect_identical(res[,,2], array(c(23, 34, 31, 46), dim=c(2,2)))
# })
# 
# 
# test_that("grb_diags_field produce expected output format", {
#   myarrays <- list(
#     array(rnorm(100000)+10, dim = c(10, 10, 1000)),
#     array(rnorm(100000)+10, dim = c(10, 10, 1000)),
#     array(rnorm(100000)+10, dim = c(10, 10, 1000))
#   )
#   expect_error(
#     res <- grb_diags_field(myarrays, iterations = seq(1000)*10),
#     NA
#   )
#   expect_type(res, "list")
#   expect_named(res, c("iterations", "mean", "var", "PSRF", "PSRF_quantiles"))
#   
#   expect_true(is(res$iterations, "numeric"))
#   expect_true(is(res$mean, "list"))
#   expect_true(is(res$var, "list"))
#   expect_true(is(res$PSRF, "array"))
#   expect_true(is(res$PSRF_quantiles, "array"))
#   
#   expect_length(res$iterations, 500)
#   expect_length(res$mean, 3)
#   expect_length(res$var, 3)
#   expect_length(res$PSRF, 50000)
#   expect_length(res$PSRF_quantiles, 2000)
#   
#   target_dim <- c(10,10,500)
#   expect_identical(lapply(res$mean), rep(list(target_dim)))
#   expect_identical(lapply(res$var), rep(list(target_dim)))
#   expect_identical(dim(res$PSRF), target_dim)
#   expect_identical(dim(res$PSRF), target_dim)
#   expect_identical(dim(res$PSRF_quantiles), c(4,500))
#   
#   # Test with small data to test for output values
# })
# 
# 
# test_that("grb_diags_field produce expected output values", {
#   myarrays <- list(
#     array(rnorm(100)+10, dim = c(10, 10, 10)),
#     array(rnorm(100)+10, dim = c(10, 10, 10)),
#     array(rnorm(100)+10, dim = c(10, 10, 10))
#   )
#   set.seed(123)
#   expect_error(
#     res <- grb_diags_field(myarrays, iterations = seq(10)*10),
#     NA
#   )
#   expect_type(res, "list")
#   expect_named(res, c("iterations", "mean", "var", "PSRF", "PSRF_quantiles"))
#   expect_identical(mean(res$iterations), 80)
#   expect_equal(sapply(res$mean, mean), c(9.835054, 10.146240, 10.129054), tolerance = 1e-5)
#   expect_equal(sapply(res$var, mean), c(2.298606e-15, 1.792936e-15, 1.315688e-15), tolerance = 1e-10)
#   res$PSRF[is.infinite(res$PSRF)] <- 0
#   expect_equal(sum(res$PSRF)/1e29, 2.70271258397, tolerance = 1e-10)
#   res$PSRF_quantiles[is.infinite(res$PSRF_quantiles)] <- 0
#   expect_equal(sum(res$PSRF_quantiles)/1e14, 2.8483681325, tolerance = 1e-10)
# })
# 
# test_that("plot_PSRF produce no error", {
#   myarrays = list(
#     array(rnorm(600), dim = c(4, 2, 100)),
#     array(rnorm(600), dim = c(4, 2, 100)),
#     array(rnorm(600), dim = c(4, 2, 100))
#   )
#   res = grb_diags_field(record_arrays = myarrays, 
#                          iterations = seq(100), 
#                          burn_in = .5, 
#                          starting_proportion = .5)
#   expect_error(
#     plot_PSRF(res, varname = "Example"),
#     NA
#   )
#   res$PSRF[1,1,1] <- 1
#   expect_error(
#     plot_PSRF(res),
#     NA
#   )
# })
# 
# test_that("plot_log_scale produce no error", {
#   myarrays = list(
#       array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(6, 1, 1000)),
#       array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(6, 1, 1000)),
#       array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(6, 1, 1000))
#   )
#   expect_error(
#     res <- plot_log_scale(log_scale_arrays = myarrays, 
#                         iterations = seq(100), 
#                         starting_proportion = .5,
#                         varname = "toto"),
#     NA
#   )
#   expect_type(res, "list")
#   myarrays = list(
#     array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(1, 6, 1000)),
#     array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(1, 6, 1000)),
#     array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(1, 6, 1000))
#   )
#   expect_error(
#     res <- plot_log_scale(log_scale_arrays = myarrays, 
#                           iterations = seq(100), 
#                           starting_proportion = .5,
#                           varname = "toto"),
#     NA
#   )
#   expect_null(res)
# })
# 
# test_that("plot_beta produce no error", {
#   beta_arrays = list(
#     array(rnorm(400), dim = c(2, 2, 100)),
#     array(rnorm(400), dim = c(2, 2, 100)),
#     array(rnorm(400), dim = c(2, 2, 100))
#   )
#   
#   expect_error(
#     res <- plot_beta(beta_arrays, 
#                      iterations=seq(100), 
#                      starting_proportion = .5, 
#                      varname = "Example", 
#                      var_names = c(1, 2)),
#     NA
#   )
# })
# 
# test_that("diagnostic_plots produce no error", {
#   # TODO
# })
# 
# test_that("ESS produce expected_output", {
#   # TODO
# })
