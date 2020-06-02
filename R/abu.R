abu <- function( x ){

  # calculate abundance from the highest taxonomic level. At this level
  # there should not be loss of information. The -1 is to remove the
  # taxa column. Two option because of the presence-absence issue

  res <- apply( x[[ "Phylum" ]][ , -1 , drop = F ] , 2 , FUN = sum )

  res

}

