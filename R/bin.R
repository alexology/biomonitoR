#' bin
#'
#' @description
#' Function for calculating presence-absence from a vector. It does not take into account NAs.
#'
#' @param x A vector of numbers.
#'
#' @export

bin <- function(x) {
  sum(any(x > 0))
}
