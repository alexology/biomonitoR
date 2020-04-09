#' abuTax
#'
#' This function calculates the absolute or relative abundance of a taxon or of a set taxa.
#' @param x result of the function aggregatoR.
#' @param taxa a taxon or a vector of taxa.
#' @param rel if TRUE calculates relative abundance. Default to FALSE.
#' @details This function does not check for ambiguous. For instance if the vector of taxa contains both the genus *Baetis* and the
#' family Baetidae, the abundance of the genus *Baetis* will be double counted. Check the function \code{\link{ambiguousSolver}} for a solution
#' to this problem.
#' @keywords aggregatoR
#' @export
#' @seealso \code{\link{aggregatoR}} \code{\link{ambiguousSolver}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' abuTax(data.agR, taxa = "Ephemeroptera")
#' abuTax(data.agR, taxa = c("Setodes", "Orthocladiinae"), rel = TRUE)

abuTax <- function ( x, taxa = NULL, rel = FALSE )
{

  # check if the object d is of class "biomonitoR"
  classCheck( x )

  # check if the taxa argument is empty or contains NULL strings
  if ( is.null( taxa ) == T || ( any( taxa == "" ) & length( taxa ) == 1 ) ) {
    stop( "Please provide a taxon name" )
  }

  # transform the taxa argument to character
  taxa <- as.character( taxa )

  # get the list of all the taxa present in the x object including the information stored
  # in the taxonomic tree
  df <- x[["Tree"]][, 1:10 ]
  df.vec <- unique( as.character( unlist( df ) ) )
  df.vec <- df.vec[ ! df.vec  %in% ""  ]

  # get the name of the x object
  x.name <- deparse( substitute( x ) )

  # find any taxa not present in the dataset provided by the user. If no taxa belongs to the list of all the taxa present
  # in the x object the function will be stopped, otherwise it will provide the list of missing taxa
  taxa.sub <- taxa[ ! taxa %in% df.vec ]
  if( length( taxa.sub ) > 0 ){
    print( paste( "The following taxon were not find in the ", x.name ," database and has been excluded: ", taxa.sub , sep = "" ) )
    taxa <- taxa[ taxa %in% df.vec ]
    if( length( taxa ) == 0 ){
      stop( "None of the taxa provided were found in the ", x.name ," database" )
    }
  }

  # list the row and column numbers at which the desired taxa are stored
  taxind <- data.frame( row = numeric(), col = numeric() )
  for ( i in 1:length( taxa ) ) {
    temp <- which( df == taxa[i] , arr.ind = T )
    taxind <- rbind( temp, taxind )
  }
  taxcol <- unique( taxind[ , "col" ] )
  taxgroup <- names( df )[ taxcol ]
  taxgroup <- unique( taxgroup )

  # subset the data.frame of the taxonomic levels of the taxa provided by the user from the x object and rbind them
  # alternatively we could merge the data.frame of all the taxonomic levels, probably faster for
  # small dataset

  taxsub <- x[ taxgroup ]
  for ( i in 1:length( taxsub ) ){
    colnames( taxsub[[i]] )[1] <- "Taxon"
  }
  taxsub <- do.call( "rbind" , taxsub )
  rownames( taxsub ) <- NULL

  # calculate the abundance of the selected taxon or taxa. As always -1 remove the first column containing the taxa list
  abucum <- apply( taxsub[ which( taxsub$Taxon %in% taxa ), -1 , drop = FALSE ] , 2 , sum )
  if (rel == TRUE) {
    abucum <- abucum / abu( x )
  }

  abuperc.v <- as.vector( t( abucum ) )
  names( abuperc.v ) <- names( abucum )
  abuperc.v
}
