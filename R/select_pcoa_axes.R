#' @title Select optimal number of pcoa axes
#'
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("experimental") }
#'
#' Select the number of the axis for the calculation of functional richness.
#'
#' @param x A distance matrix.
#' @param method Quality measure of the functional space. Can be `cor`, `legendre` and `maire`.
#' @param tresh the treshold for selecting the optimal number of axes. By default is 0.7 for `cor` and `legendre` and
#' 0.01 for `maire`.
#' @param nbdim The maximum number of dimension when the option is set to `maire`. Default to 15.
#'
#' @details This function provides some information to select the number of axes for the calculation of
#' functional richness. As implemented in `biomonitoR`, functional richness requires a pcoa to be computed from
#' the trait distances. For calculating functional richness a subset of the number of axes originating
#' from the pcoa is kept. The choiche of the number of axis can affect the calculation and some
#' information about how much the subset represents the original distance matrix is desirable. To this purpose, 3 measures
#' are provide. The `cor` option represents the correlation between the distance of the points in the
#' reduced space and the original distance matrix. The distance of points in the reduced space is calculated
#' with the `dist()` function from the `stats` package on the coordinates of the selected number of axis.
#' A p-value based on a mantel statistic is also provided.
#' The `legendre` option calculates the R^2 like ratio as described by Legendre and Legendre (2008)
#' and implemented in the `FD` package.
#' The `maire` option calculate the quality of the functional space according to Maire et al. (2015).
#' The function selectPcoaAxis also provides information about the euclidean property of the trait
#' distance. The lack of this property can lead to biases and some solutions are proposed here and
#' described in Legendre and Legendre (2008). Please take a look also to the help of the `FD`
#' function `dbFD` for further information.
#'
#' @keywords aggregate_taxa
#' @export
#' @importFrom dplyr bind_rows
#' @seealso [aggregate_taxa]
#' @examples
#'
#' library(ade4)
#'
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' data_ts <- assign_traits(data_agr)
#'
#' # averaging
#' data_ts_av <- average_traits(data_ts)
#'
#' col_blocks <- c(8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8)
#'
#' rownames(data_ts_av) <- data_ts_av$Taxa
#' traits_prep <- prep.fuzzy(data_ts_av[, -1], col.blocks = col_blocks)
#'
#' traits_dist <- ktab.list.df(list(traits_prep))
#' traits_dist <- dist.ktab(traits_dist, type = "F")
#'
#' select_pcoa_axes(traits_dist, method = "cor", tresh = 0.7)
#' select_pcoa_axes(traits_dist, method = "legendre", tresh = 0.7)
#' select_pcoa_axes(traits_dist, method = "maire", tresh = 0.01)
#' @references Maire, E., Grenouillet, G., Brosse, S., & Villeger, S. (2015).
#'   How many dimensions are needed to accurately assess functional diversity?
#'   A pragmatic approach for assessing the quality of functional spaces. Global
#'   Ecology and Biogeography, 24(6), 728-740.
#' @references Legendre, P. and L. Legendre (1998) Numerical Ecology.
#' 2nd English edition. Amsterdam: Elsevier.





select_pcoa_axes <- function(x, method = "legendre", tresh = NULL, nbdim = 15) {
  if (identical(method, "cor") & is.null(tresh)) {
    tresh <- 0.7
  }
  if (identical(method, "legendre") & is.null(tresh)) {
    tresh <- 0.7
  }
  if (identical(method, "maire") & is.null(tresh)) {
    tresh <- 0.01
  }

  if (suppressWarnings(is.euclid(x)) & identical(method, "cor")) {
    res <- pcoaQuality(x, "none", method = "cor", tresh = tresh)
  }

  if (suppressWarnings(is.euclid(x)) & identical(method, "legendre")) {
    res <- pcoaQuality(x, "none", method = "legendre", tresh = tresh)
  }

  if (suppressWarnings(is.euclid(x)) & identical(method, "maire")) {
    res <- pcoaQuality(x, "none", method = "maire", tresh = tresh)
  }

  if (!suppressWarnings(is.euclid(x)) & identical(method, "cor")) {
    res1 <- pcoaQuality(x, "none", method = "cor", tresh = tresh)
    res2 <- pcoaQuality(x, "cailliez", method = "cor", tresh = tresh)
    res3 <- pcoaQuality(x, "lingoes", method = "cor", tresh = tresh)
    res4 <- pcoaQuality(x, "sqrt", method = "cor", tresh = tresh)
    res5 <- pcoaQuality(x, "quasi", method = "cor", tresh = tresh)
    res <- suppressWarnings(bind_rows(res1, res2, res3, res4, res5))
  }

  if (!suppressWarnings(is.euclid(x)) & identical(method, "legendre")) {
    res1 <- pcoaQuality(x, "none", method = "legendre", tresh = tresh)
    res2 <- pcoaQuality(x, "cailliez", method = "legendre", tresh = tresh)
    res3 <- pcoaQuality(x, "lingoes", method = "legendre", tresh = tresh)
    res4 <- pcoaQuality(x, "sqrt", method = "legendre", tresh = tresh)
    res5 <- pcoaQuality(x, "quasi", method = "legendre", tresh = tresh)
    res <- suppressWarnings(bind_rows(res1, res2, res3, res4, res5))
  }


  if (!suppressWarnings(is.euclid(x)) & identical(method, "maire")) {
    res1 <- pcoaQuality(x, "none", method = "maire", nbdim = nbdim, tresh = tresh)
    res2 <- pcoaQuality(x, "cailliez", method = "maire", nbdim = nbdim, tresh = tresh)
    res3 <- pcoaQuality(x, "lingoes", method = "maire", nbdim = nbdim, tresh = tresh)
    res4 <- pcoaQuality(x, "sqrt", method = "maire", nbdim = nbdim, tresh = tresh)
    res5 <- pcoaQuality(x, "quasi", method = "maire", nbdim = nbdim, tresh = tresh)
    res <- suppressWarnings(bind_rows(res1, res2, res3, res4, res5))
  }


  res
}
