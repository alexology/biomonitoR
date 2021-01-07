#' selectPcoaAxes
#'
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("deprecated") }
#'
#' Select the number of the axis for the calculation of functional richness.
#'
#' @param x a distance matrix.
#' @param method the measure of the quality of the functional space. Can be `cor`, `legendre` and `maire`.
#' @param tresh the treshold for selecting the axis. By default is 0.7 for `cor` and `legendre` and
#' 0.01 for `maire`.
#' @param nbdim the maximum number of dimension when the option is set to `maire`. Default to 15.
#' @details This function provide some information to select the number of axis for the calculation of
#' functional richness. As implemented in biomonitoR, functional richness require a pcoa to be computed from
#' the trait distances. For calculating functional richness a subset of the number of axis originating
#' from the pcoa is kept. The choiche of the number of axis can affect the calculation and some
#' information about how much the subset represents the original distance matrix. To this purpose, 3 measures
#' are provide. The `cor` option represents the correlation between the distance of the points in the
#' reduced space and distance matrix. The distance of points in the reduced space is calculated
#' with the `dist` function on the coordinates of the selected number of axis.
#' A p-value based on a mantel statistic is also provided.
#' The `legendre` option calculates the R^2 like ratio as described by Legendre and Legendre (2008)
#' and implemented in the FD package.
#' The `maire` option calculate the quality of the functional space according to Maire et al. (2015).
#' The function selectPcoaAxis provide also information about the euclidean property of the trait
#' distance. The lack of this property can lead to biases and some solutions are proposed here and
#' described in Legendre and Legendre (2008). Please take a look also to the help of the FD
#' function dbFD for further information.
#'
#' @keywords aggregatoR
#' @export
#' @importFrom dplyr bind_rows
#' @seealso \code{\link{aggregatoR}}
#'
#' @references Maire, E., Grenouillet, G., Brosse, S., & Villeger, S. (2015).
#'   How many dimensions are needed to accurately assess functional diversity?
#'   A pragmatic approach for assessing the quality of functional spaces. Global
#'   Ecology and Biogeography, 24(6), 728-740.
#' @references Legendre, P. and L. Legendre (1998) Numerical Ecology.
#' 2nd English edition. Amsterdam: Elsevier.





selectPcoaAxes <- function(x, method = "legendre", tresh = NULL, nbdim = 15) {
  .Deprecated("select_pcoa_axes", package = "biomonitoR")

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
