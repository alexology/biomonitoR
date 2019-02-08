#' combTaxa
#'
#' This function returns all the combinations of Taxa at the desired taxonomic resolution.
#' @param x results of function aggregatoR.
#' @param ntaxa number of Taxa to choose.
#' @param taxLev taxonimc level on which the calculation has to be made.
#' @details This function is intended to help the user to identify the best subset of taxa that correlate well with environmental variables.
#' Metrics based on only certain taxa (EPT, 1-GOLD, etc) are currently used in biomonitoring and their are usually based on autoecological knowledge
#' of the selected taxa. However, the relationship between an environmental variable and a subset of taxa could exist and not detected due to the actual limitations
#' about the autoecological knowledge of a group of organsism. With more than 4 taxa the calculations should become infeasible, expecially when the number
#' of taxa in the user dataset is high. For instance the number of combinations taken 20 at time is nearly 20000.
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
