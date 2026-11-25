#' @title
#' GeoNonStat demo data
#'
#' @description
#' A S3 GeoNonStat object with nice non-stationarity and little noise.
#' generated for demonstration and examples purpose.
#' @format a GeoNonStat object with 3000 observed locations and field.
#' No MCMC has run on this object.
#'
#' @details
#' The data have been generated with the code in 
#' `data_raw/code_to_generate_example_data.R`
#'
#'
#' @return
#' A \code{\link{GeoNonStat}}.
#' @examples
#' gnsDemo
"gnsDemo"


#' @title
#' GeoNonStat demo data
#'
#' @description
#' A S3 GeoNonStat object with nice non-stationarity and little noise.
#' generated for demonstration and examples purpose. MCMC has already run on
#' this object.
#' @format a GeoNonStat object with 3000 observed locations and field.
#' 15 iterations of MCMC have run on this object.
#'
#' @details
#' The data have been generated with the code in 
#' `data_raw/code_to_generate_example_data.R`
#'
#' @return
#' A \code{\link{GeoNonStat}} with processed MCMC.
#' @examples
#' processedGnsDemo
"processedGnsDemo"
