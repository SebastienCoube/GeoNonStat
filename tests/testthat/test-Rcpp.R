test_that("vecchia works", {
  set.seed(100000)
  n <- 400L
  # spatial locations
  locs <- cbind(runif(n), runif(n))
  log_range <- matrix(0, n, 1)
  log_range[, 1] <- -10 + 3 * locs[, 1]
  NNarray <- GpGp::find_ordered_nn(locs, 10)
  NNarray[is.na(NNarray)] <- 0

  # isotropic, with smoothness = 1.5
  set.seed(1)
  expect_error(
    res <- vecchia(
      log_range = t(log_range),
      locs = t(locs),
      NNarray = t(NNarray),
      num_threads = 1,
      compute_derivative = F,
      smoothness = 1.5
    ),
    NA
  )
  expect_true(is(res, "array"))
  expect_identical(dim(res), c(11L, 1L, n))
  print(mean(res))
  expect_equal(mean(res), 0.054481, tolerance = 1e-5)


  # nominally anisotropic, with smoothness = 1.5
  set.seed(1)
  log_range <- cbind(log_range, log_range, 0)
  expect_error(
    res <- vecchia(
      log_range = t(log_range),
      locs = t(locs),
      NNarray = t(NNarray),
      num_threads = 1,
      compute_derivative = F,
      smoothness = 1.5
    ),
    NA
  )
  expect_true(is(res, "array"))
  expect_identical(dim(res), c(11L, 1L, n))
  print(mean(res))
  expect_equal(mean(res), 0.054496, tolerance = 1e-5)

  # actually anisotropic, with smoothness = 1.5
  set.seed(1)
  log_range[, 2] <- -10 + 3 * locs[, 2]
  log_range[, 3] <- locs[, 2] - 2 * locs[, 1]
  tatato <- vecchia(
    log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = F, smoothness = 1.5
  )
  # GeoNonStat::plot_pointillist_painting(locs, GpGp::fast_Gp_sim_Linv (t(tatato[,1,]), NNarray))

  # isotropic, with smoothness = 0.5
  set.seed(1)
  log_range <- matrix(0, n, 1)
  log_range[, 1] <- -10 + 3 * locs[, 1]
  tatato <- vecchia(
    log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = F, smoothness = 0.5
  )
  # GeoNonStat::plot_pointillist_painting(locs, GpGp::fast_Gp_sim_Linv (t(tatato[,1,]), NNarray))

  # nominally anisotropic, with smoothness = 0.5
  set.seed(1)
  log_range <- cbind(log_range, log_range, 0)
  tatato <- vecchia(
    log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = F, smoothness = 0.5
  )
  # GeoNonStat::plot_pointillist_painting(locs, GpGp::fast_Gp_sim_Linv (t(tatato[,1,]), NNarray))

  # actually anisotropic, with smoothness = 0.5
  set.seed(1)
  log_range[, 2] <- -10 + 3 * locs[, 2]
  log_range[, 3] <- locs[, 2] - 2 * locs[, 1]
  tatato <- vecchia(
    log_range = t(log_range), locs = t(locs), NNarray = t(NNarray), num_threads = 1, compute_derivative = F, smoothness = 0.5
  )
  # GeoNonStat::plot_pointillist_painting(locs, GpGp::fast_Gp_sim_Linv (t(tatato[,1,]), NNarray))
})
