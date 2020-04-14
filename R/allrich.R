#' @name allrich
#' @title Function for calculating richness at different taxonomic levels
#'
#' @description This function allows the calculation of richness at different taxonomic levels and of the total abundance.
#' @aliases  richness
#' @param x result of the function aggregatoR
#' @param taxLev taxonomic level on which the calculation has to be made.
#' @details The function `richness` returns a warning meassage if unassigned taxa are detected at the desired taxonomic level. Unassigned taxa can be present
#' because some taxa are entered at higher taxonomic levels than those at which the calculation is performed.
#' For instance this problem appears when the calculation is done at genus level and some taxa has been entered at family level.
#' Unassigned taxa can be present also because of missing levels in the reference database at the desired taxonomic level.
#' This problem is present, for example, for Mollusca in the default database for macroinvertebrates.
#' Mollusca orders are not present because there is no consensus among taxonomist (see the [www.molluscabase.org](https://www.molluscabase.org/index.php) website).
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
  classCheck( x )

  # extract the names of all the taxonomic levels on which perform the calculation
  n.tree <- names( x[[ "Tree" ]][ , sapply( x[[ "Tree" ]] , is.factor ) ] )

  # initialize the data.frame to store the results
  res <- data.frame( matrix( NA , ncol = length( n.tree ), nrow = sum( ! names( x[[ "Tree" ]]) %in% n.tree ) ) )

  # a for loop to calculate richness at different taxonomic levels
  # while filling the res data.frame
  for( i in 1 : length( n.tree ) ){
    res[ , i ] <- richness( x , taxLev = n.tree[ i ] )
  }

  # assign column and row names and return the results
  colnames( res ) <- n.tree
  rownames( res ) <- colnames( x[[ "Tree" ]] )[ ! colnames( x[[ "Tree" ]] ) %in% n.tree ]
  res

}

