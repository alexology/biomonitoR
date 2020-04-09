#' @name abundance
#' @title abundance
#'
#' @description This function allows the calculation of abundances with the possibility to exclude unassigned taxa at the desired taxonomic level.
#' @param x result of the function aggregatoR
#' @param taxLev taxonomic level on which the calculation has to be made.
#' @param unassigned Does unassigned taxa need to be taken into account? If yes set unassigned to TRUE, otherwise FALSE.
#' @keywords abundance
#' @details If unassigned is set to TRUE abundances are equal among taxonomic levels.
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' abundance(data.agR, taxLev = "Family" )
#' abundance(data.agR, taxLev = "Family" , unassigned = TRUE)

abundance <- function( x , taxLev = "Family" , unassigned = FALSE ){

  # check if the object d is of class "biomonitoR"
  classCheck( x )

  # get the desired taxomic level
  tax <- x[[ taxLev ]]

  if( unassigned == FALSE ){
    if( "unassigned" %in% tax[ , 1 ] ){
      z <- which( tax[ , 1 ] == "unassigned" ) # find the row corresponding to the unassigned
      tax <- tax[ -z , ] # remove unassigned row from the species count
    }
  }

  # calculate abundance from the desired taxonomic level. The -1 is to remove the taxa column
  res <- apply( tax[ , -1 , drop = F ], 2 , FUN = sum )
  return( res )

}
