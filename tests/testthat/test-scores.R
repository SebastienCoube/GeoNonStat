test_that("splitData produce expected results", {
  set.seed(123)
  # create repeted locations
  n_sites <- 1000
  n_rep = 2
  sites <- cbind(x = runif(n_sites), y = runif(n_sites))
  locs  <- sites[sample(seq_len(n_sites), n_sites, replace = TRUE), , drop = FALSE]
  rownames(locs) <- paste0("l", seq_len(n_sites))
  nlocs     <- nrow(locs)
  datalist <- list(
    titi = rnorm(nlocs),
    toto = matrix(rnorm(nlocs * 2), ncol = 2)
  )
  expect_error(
    out <- splitData(locs = locs, datalist = datalist, plot=FALSE, n_clust_far = 50, prop_test_far = 0.1),
    NA
  )
  expect_named(out, c("train", "test_close", "test_far", "discarded"),
               ignore.order = TRUE)
  expect_true(!is.null(attr(out, "info")))
  for (b in out) {
    expect_true("locs" %in% names(b))
    expect_true(all(names(datalist) %in% names(b)))
    n <- nrow(b$locs)
    expect_equal(NROW(b$titi), n)
    expect_equal(NROW(b$toto), n)
  }
  
  total <- sum(vapply(out, function(b) nrow(b$locs), numeric(1)))
  expect_equal(total, nlocs)
  # No duplicated locations  
  expect_equal(sum(vapply(out, function(b) length(unique(rownames(b$locs))), numeric(1))), n_sites)
  expect_equal(length(unique(unlist(sapply(out, function(b) rownames(b$locs))))), n_sites)
})


test_that("splitData produce expected results if not enough spatial clusters", {
  set.seed(123)
  # create repeted locations
  n_sites <- 1000
  n_rep = 2
  sites <- cbind(x = runif(n_sites), y = runif(n_sites))
  locs  <- sites[sample(seq_len(n_sites), n_sites, replace = TRUE), , drop = FALSE]
  nlocs     <- nrow(locs)
  datalist <- list(
    titi = rnorm(nlocs),
    toto = matrix(rnorm(nlocs * 2), ncol = 2)
  )
  expect_warning(
    out <- splitData(locs = locs, datalist = datalist, plot=FALSE, n_clust_far = 5, prop_test_far = 0.1),
    "Not enough observations to produce far test sample"
  )
 
  expect_named(out, c("train", "test_close"),
               ignore.order = TRUE)
  expect_true(!is.null(attr(out, "info")))
  for (b in out) {
    expect_true("locs" %in% names(b))
    expect_true(all(names(datalist) %in% names(b)))
    n <- nrow(b$locs)
    expect_equal(NROW(b$titi), n)
    expect_equal(NROW(b$toto), n)
  }
  
  total <- sum(vapply(out, function(b) nrow(b$locs), numeric(1)))
  expect_equal(total, nlocs)
})

