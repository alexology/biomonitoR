#' Check a vector for presence-absence
#'
#' @description
#' Function for calculating presence-absence from a vector. It does not take into account `NA`.
#'
#' @param x A vector of numbers.
#'
#' @export
#'
#' @examples
#'
#' # There is at least a presence
#' bin(c(0, 1, 10))
#'
#' # There are not presence
#' bin(c(0, 0))

bin <- function(x) {
  sum(any(x > 0))
}
