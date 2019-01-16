#' @name abundance
#' @title abundance
#'
#' @description This function allow the calculation of abundances with the possibility to exlude unassigned taxon.
#' @aliases  richness abu
#' @param x results of function aggregatoR
#' @param taxLev taxonomic level on which the calculation has to be made.
#' @param unassigned Does unassigned taxa need to be taken into account? If yes set unassigned to TRUE, otherwise FALSE.
#' @keywords abundance
#' @details If unassigned equal to TRUE calculated abundances are equal among taxonomic levels.
#' @export
#' @export richness
#' @export abu
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' abundance(data.agR, taxLev = "Family" )
#' abundance(data.agR, taxLev = "Family" , unassigned = TRUE)

abundance <- function( x , taxLev = "Family" , unassigned = TRUE ){
  classCheck(x)

  tax <- x[[ taxLev ]]

  if( unassigned == TRUE ){
    if( "unassigned" %in% tax[,1] ){
      z <- which( tax$Taxon == "unassigned" )
      tax <- tax[ -z , ] # remove unassigned row from the species count
    }
  }

  res <- apply( tax[ , -1 , drop = F ], 2, FUN = sum )
  return( res )

}
