#' bin
#'
#' @param x a vector of numbers
#'
#' Function for calculating presence-absence from a vector. It does not take into account NAs.
#' @export

bin <- function( x ){
  sum( any( x > 0 ) )
}
