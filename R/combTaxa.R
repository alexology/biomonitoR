#' combTaxa
#'
#' This function returns all the combinations of n taxa at the desired taxonomic resolution.
#' @param x result of the function aggregatoR.
#' @param ntaxa number of taxa.
#' @param taxLev taxonomic level on which the calculation has to be made.
#' @details This function is intended to help the user to identify the best subset of taxa that correlates with environmental variables.
#' Metrics based on a specific group of taxa (EPT, 1-GOLD, etc) are currently used in biomonitoring. They rely on available information on the responses
#' of the selected taxa to target strssors. The relationship between an environmental variable and a subset of taxa could exist and not detected due to the scarce
#' availaibility of autoecological data. With more than 4 taxa the calculations should become infeasible, expecially when the number
#' of taxa in the user dataset is high.
#' @keywords aggregatoR
#' @importFrom utils combn
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' combTaxa(data.agR, taxLev = "Family")



combTaxa <- function( x, ntaxa = 2 , taxLev = "Taxa" ){

  # check if the object x is of class "biomonitoR"
  classCheck( x )

  # get the data.frame at the desired taxonomic level
  df <- x[[ taxLev ]]

  # remove unassigned row from the species count if present
  if( "unassigned" %in% df[ , 1 ] ){
    z <- which( df[ , 1 ] == "unassigned" )
    df <- df[ -z , ] # remove unassigned row from the species count
  }

  # list all the combination of x taxa taken n at time
  cbn <- combn( nrow( df ) , ntaxa )
  df.sub <- lapply( seq( ncol( cbn ) ) , function( x ) df[ cbn[ , x ] , ] )

  # sum the abundances of the n-th combination of df.sub and change the label
  df.agg <- lapply( df.sub , FUN = agg_fun )

  # create a data.frame to store the results
  return( do.call( rbind, df.agg ) )
}
