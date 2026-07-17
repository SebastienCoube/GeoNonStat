library(testthat)

if (!requireNamespace("GeoNonStat", quietly = TRUE)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required to run tests from source")
  }
  pkgload::load_all(".")
}

test_check("GeoNonStat")
