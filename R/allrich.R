#' @name allrich
#' @title Function for calculating richness at different taxonomic levels
#'
#' @description This function allows the calculation of richness at different taxonomic levels and of the total abundance.
#' @aliases  richness
#' @param x results of function aggregatoR
#' @param taxLev taxonomic level on which the calculation has to be made.
#' @keywords allrich, richness
#' @export
#' @export richness
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' richness(data.agR, taxLev = "Family" )


allrich <- function( x ){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  n.tree <- names( x[[ "Tree" ]][ , sapply(x[[ "Tree" ]], is.factor) ] )
  res <- data.frame( matrix( NA , ncol = length( n.tree ),
                                nrow = sum( ! names( x[[ "Tree" ]]) %in% n.tree ) ) )
  for( i in 1 : length( n.tree ) ){
    res[ , i ] <- richness( x , taxLev = n.tree[ i ])
  }

  colnames( res ) <- n.tree
  rownames( res ) <- colnames( x[[ "Tree" ]] )[ ! colnames( x[[ "Tree" ]] ) %in% n.tree ]
  return( res )

}

