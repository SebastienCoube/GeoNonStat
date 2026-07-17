# inBounds tests
test_that("inBounds works with open and closed intervals", {
  expect_true(inBounds(c(0, 1), 0.5))
  expect_false(inBounds(c(0, 1), 0))
  expect_true(inBounds(c(0, 1), 0, includeBounds = TRUE))
  expect_false(inBounds(c(0, 1), 1, includeBounds = FALSE))
  expect_true(inBounds(c(0, 1), 1, includeBounds = c(TRUE, TRUE)))
  expect_equal(inBounds(c(0, 1), c(-1, 0, 0.5, 1, 2), includeBounds = TRUE),
               c(FALSE, TRUE, TRUE, TRUE, FALSE))
})

test_that("inBounds rejects invalid arguments", {
  expect_error(inBounds(c(0, 1, 2), 0.5), "bounds should be a numeric of length 2")
  expect_error(inBounds(c(0, 1), 0.5, includeBounds = c(TRUE, TRUE, TRUE)), "includeBounds should be a logical of length 1 or 2")
})

# subsetData tests
test_that("subsetData subsets data.frames, vectors and handles NULL", {
  data <- list(
    df = data.frame(a = 1:4, b = 5:8),
    vec = letters[1:4],
    empty = NULL
  )
  out <- subsetData(data, c(2, 4))
  expect_equal(out$df, data.frame(a = c(2, 4), b = c(6, 8)), ignore_attr=TRUE)
  expect_equal(out$vec, c("b", "d"))
  expect_null(out$empty)
})

# symmat and expmat tests
test_that("symmat builds symmetric matrices and rejects invalid lengths", {
  coords <- c(1, 2, 3, 4, 5, 6)
  mat <- symmat(coords)
  expect_equal(mat, matrix(c(1, 4, 5, 4, 2, 6, 5, 6, 3), 3, 3))
  expect_error(symmat(1:2), "length of coords incompatible with a symetric matrix")
})

test_that("expmat returns an exponential matrix with diagonal jitter", {
  coords <- c(1, 2, 3, 4, 5, 6)
  res <- expmat(coords, eps = 1e-6)
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(3, 3))
  expect_true(all(diag(res) > 0))
})

# transposeList tests
test_that("transposeList transposes a named list of sublists", {
  input <- list(
    A = list(x = 1, y = 2),
    B = list(x = 3, y = 4)
  )
  out <- transposeList(input)
  expect_equal(names(out), c("x", "y"))
  expect_equal(out$x, list("A"=1, "B"=3))
  expect_equal(out$y, list("A"=2, "B"=4))
})

# naiveGreedyColoring tests

test_that("naiveGreedyColoring colors a simple symmetric sparse matrix", {
  M <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 3),
    j = c(1, 2, 1, 2, 3),
    x = c(1, 1, 1, 1, 1),
    dims = c(3, 3)
  )
  cols <- naiveGreedyColoring(M)
  expect_equal(length(cols), 3)
  expect_equal(cols, c(1, 2, 1))
})

# xPPMultRight and xPPCrossprod minimal structure tests

test_that("xPPMultRight rejects both X and PP NULL", {
  vecchia_approx <- list(n_locs = 1, locs_match = 1)
  expect_error(xPPMultRight(
    X = NULL, 
    PP = NULL, 
    vecchia_approx = vecchia_approx, 
    Y = matrix(1, 1, 1)), "X and PP can't be both NULL")
})

test_that("xPPCrossprod rejects mismatched row counts", {
  expect_error(xPPCrossprod(
    X = matrix(1, 2), 
    Y = matrix(1, 1), 
    vecchia_approx = NULL), "X and Y should have the same number of rows")
})
