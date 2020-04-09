abu <- function( x ){

  # check if the object d is of class "biomonitoR"
  classCheck( x )

  # calculate abundance from the highest taxonomic level. At this level
  # there should not be loss of information. The -1 is to remove the
  # taxa column
  nord <- apply( x[[ "Phylum" ]][ , -1 , drop = F ], 2, FUN = sum )
  return( nord )
}

