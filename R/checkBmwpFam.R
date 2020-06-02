#' @importFrom stats aggregate
checkBmwpFam <- function( DF , famNames , stNames ){
  # initialize the Taxon object with which the aggregation will be done
  # this is needed to avoid warning from R CMD
  Taxon <- NULL

  # transform taxa names to characted
  DF[ , "Taxon" ] <- as.character( DF[ , "Taxon" ] )

  # subset the taxa that need to be changed
  famCheck <- DF[ which( DF[ , "Taxon" ] %in% famNames[ , "Taxon" ] ) , "Taxon" ]

  # if taxa that need to be changed are present chenge them, otherwise return the orginal DF
  if( length( famCheck ) > 0 ){
    for( i in 1 : length( famCheck ) ){
      taxName <- famCheck[ i ]
      DF$Taxon[ DF$Taxon == taxName ] <- as.character( subset( famNames, Taxon == taxName )[ , "Correct_Taxon" ] )
    }

    DF <- aggregate( . ~ Taxon , DF , sum )
    names( DF ) <- c( "Taxon" , stNames )
    return( DF )
  }
  else{ DF <- DF
  return( DF )
  }
}
