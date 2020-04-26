#' eptd
#'
#' This function calculates the log10(Sel_EPTD + 1) metric, where EPTD stands for Ephemeroptera, Plecoptera, Trichoptera and Diptera.
#' @param x results of function aggregatoR
#' @keywords ept
#' @details log10(Sel_EPTD + 1) the base-10 logarithm of the abundance of the selected EPTD families plus 1. Selected EPTD families are Heptageniidae, Ephemeridae, Leptophlebiidae, Brachycentridae, Goeridae, Polycentropodidae, Limnephilidae, Odontoceridae, Dolichopodidae, Stratiomyidae, Dixidae, Empididae, Athericidae and Nemouridae.
#' This metric is part of the italian STAR_ICMi index, where it is supposed to be relate to habitat integrity.
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' eptd(data.agR)


eptd <- function (x){

  # check if the object d is of class "biomonitoR"

  # check if the object x is of class "biomonitoR"
  classCheck(x, group = "mi")

  eptd_fam <- c("Heptageniidae", "Ephemeridae", "Leptophlebiidae", "Brachycentridae", "Goeridae", "Polycentropodidae", "Limnephilidae", "Odontoceridae", "Dolichopodidae", "Stratiomyidae", "Dixidae", "Empididae", "Athericidae", "Nemouridae")
  x_fam <- x[["Family"]]
  x_eptd <- x_fam[which(x_fam$Family %in% eptd_fam), , drop=F]
  temp <- log10(apply(x_eptd[ , -1 , drop = FALSE ], 2 , sum) + 1 )

  return( temp )

}
