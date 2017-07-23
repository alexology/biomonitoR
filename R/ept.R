#' ept
#'
#' This function calculates EPT richnes
#' @param x results of function aggregatoR
#' @keywords ept
#' @details
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' ept(data.bio)
#' ept()

ept <- function (x){
  x_ept <- x[["Order"]]
  ept_taxa <- x_ept[which(x_ept$Order == "Plecoptera" |
              x_ept$Order == "Ephemeroptera" | x_ept$Order == "Trichoptera"),]
  if(nrow(ept_taxa)==0){
    ept_taxa[1,-1] <- rep(0, ncol(ept_taxa)-1)
    return(ept_taxa[,-1])
  } else {
    temp <- colSums(ept_taxa[,-1,drop=F] > 1)
    return(temp)
  }
}
