#' add_bias_to_traits
#'
#' @description
#' Add bias to traits to avoid duplicates.
#'
#' @param x Result of `average_traits()` or a `data.frame` with a similar structure.
#' @param fuzzy Are you working with fuzzy traits? Default to `TRUE`.
#' @param col_blocks A vector that contains the number of modalities for each trait
#' @param SD The amount of bias.
#'
#' @details This function works by adding a small bias to traits.
#' The bias is added with the [stats::rnorm()] function.
#' Fuzzy data are prepared with the function [ade4::prep.fuzzy()] of the `ade4` package.
#'
#' @export
#' @importFrom stats rnorm
#' @importFrom ade4 prep.fuzzy
#'
#' @examples
#'
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' data_ts <- assign_traits(data_agr)
#' data_ts_av <- average_traits(data_ts)
#'
#' data_ts_av[1:5, 1:5]
#'
#' # set col_blocks
#' col_blocks <- c(8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8)
#'
#' add_bias_to_traits(data_ts_av, fuzzy = TRUE, col_blocks = col_blocks)[1:5, 1:5]
#' add_bias_to_traits(data_ts_av, fuzzy = TRUE, col_blocks = col_blocks, SD = 0.01)[1:5, 1:5]
add_bias_to_traits <- function(x, fuzzy = TRUE, col_blocks = NULL, SD = 0.001) {

  # simulate modalities of duplicated traits
  # at first a small bias is addedd to the original traits with rnorm
  # traits with 0 values are preserved by multypling biased values by ( x > 0 )
  # to be sure to have positive values we take the abolute values and data
  # are forced to sum to 1 by dividing for sum( abs( simDat ) )

  if (fuzzy & is.null(col_blocks)) stop("Please set col_blocks")
  if (!fuzzy & !is.null(col_blocks)) stop("col_blocks will be ignored because fuzzy = FALSE")

  simTraits <- function(x, ...) {
    simDat <- x + rnorm(length(x), ...)
    simDat <- simDat * (x > 0)
    abs(simDat) / sum(abs(simDat))
  }


  i <- unlist(lapply(x, is.numeric)) & !colnames(x) %in% "Taxonomic_distance"
  x[i] <- x[i] + abs(rnorm(prod(dim(x[i])), mean = 0, sd = SD))



  if (fuzzy) {
    x[i] <- as.data.frame(prep.fuzzy(x[i], col.blocks = col_blocks))
  }

  x
}
