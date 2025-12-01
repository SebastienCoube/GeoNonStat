#' Title
#'
#' @param nlocs numeric, number of locations
#' @param nvars numeric, number of explicatives variables
#' @param nppr numeric, number of PP in PP range
#' @param nppn numeric, number of PP in PP noise
#'
#' @returns
#' @export
#'
#' @examples
initStateParams <- function(stateParams){
  lapply(stateParams, function(x){
    if(is.matrix(x)){
      return(matrix(NA_real_, nrow = nrow(x), ncol = ncol(x)))
    } else if(is.numeric(x)) {
      return(rep(NA_real_, length(x)))
    }
  })
}