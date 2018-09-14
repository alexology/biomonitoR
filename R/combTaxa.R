#' combTaxa
#'
#' This function returns all the combinations of Taxa at the desired taxonomic resolution.
#' @param x results of function aggregatoR.
#' @param ntaxa number of Taxa to choose.
#' @param taxLev taxonimc level on which the calculation has to be made.
#' @keywords aggregatoR
#' @importFrom utils combn
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' combTaxa(data.agR, taxLev = "Family")



combTaxa <- function(x, ntaxa = 2, taxLev = "Family"){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  df <- x[[taxLev]]

  if("unassigned" %in% df[,1]){
    z <- which(df[,1]=="unassigned")
    df <- df[-z,] # remove unassigned row from the species count
  }

  cbn <- combn(nrow(df), ntaxa)
  df.sub <- lapply(seq(ncol(cbn)), function(x) df[cbn[,x],])

  df.agg <- lapply(df.sub, FUN = agg_fun)

  return(do.call(rbind, df.agg))
}
