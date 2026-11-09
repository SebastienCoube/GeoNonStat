test_that("predict return expected error if no MCMC", {
  locs <- cbind(rnorm(200), rnorm(200))
  X <- data.frame(matrix(rnorm(200*4), ncol=4))
  colnames(X) <- paste0("V", 1:4)
  testdata <- list("locs" = locs, "X" = X, noise_X=X)  
  expect_error(
    testpredictions <- predict(gnsDemo, list(testdata)),
    "MCMC doesn't seem to have"
  )
})


test_that("predict return expected error if no MCMC", {
  locs <- cbind(rnorm(200), rnorm(200))
  X <- data.frame(matrix(rnorm(200*4), ncol=4))
  colnames(X) <- paste0("V", 1:4)
  testdata <- list("locs" = locs, "X" = X, noise_X=X)  
  expect_error(
    testpredictions <- predict(processedGnsDemo, list(testdata)),
    NA
  )
})
