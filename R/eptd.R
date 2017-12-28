#' eptd
#'
#' This function calculates the log10(Sel_EPTD + 1) metric, where EPTD stands for Ephemeroptera, Plecoptera Oligochaeta and Diptera.
#' @param x results of function aggregatoR
#' @keywords ept
#' @details log10(Sel_EPTD + 1) the base-10 logarithm of the abundance of the selected EPTD families plus 1.
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' eptd(data.agR)


eptd <- function (x){

  # check if the object d is of class "biomonitoR"

  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }

  eptd_fam <- c("Heptageniidae", "Ephemeridae", "Leptophlebiidae", "Brachycentridae", "Goeridae", "Polycentropodidae", "Limnephilidae", "Odontoceridae", "Dolichopodidae", "Stratiomyidae", "Dixidae", "Empididae", "Athericidae", "Nemouridae")
  x_fam <- x[["Family"]]
  x_eptd <- x_fam[which(x_fam$Family %in% eptd_fam), , drop=F]
  temp <- log10(apply(x_eptd[ , -1], 2 , sum)+1)

  return( temp )

}
