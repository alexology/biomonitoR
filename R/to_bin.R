#' Abundance to presence-absence
#'
#' @description
#' This function transforms abundance to presence-absence.
#'
#' @param x A `data.frame`
#'
#' @details
#' `to_bin` will transform values greter than 0 to 1.
#'
#' @examples
#'
#' data(macro_ex)
#'
#' to_bin(macro_ex)
#'
#' @export

to_bin <- function(x) {
  x_num <- unlist(lapply(x, is.numeric))
  x_bin <- x[, x_num]
  x_bin[x_bin > 0] <- 1
  x.df <- data.frame(x[, !x_num, drop = FALSE], x_bin, check.names = FALSE)
  class(x.df[, !x_num]) <- class(x[, !x_num])
  x.df
}
