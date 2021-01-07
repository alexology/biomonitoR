#' @describeIn allrich Taxonomic richness

richness <- function( x , tax_lev = "Taxa" ){

  # check if the object x is of class "biomonitoR"
  classCheck( x )

  # get the data.frame at the desired taxonomic level
  res <- x[[ tax_lev ]]

  # remove unassigned row from the species count if present
  if( "unassigned" %in% res[ , 1 ] ){
    z <- which( res[ , tax_lev ] == "unassigned" )
    res <- res[ -z , ] # remove unassigned row from the species count
    # a warning telling that richness calculation can be baised because unassigned taxa were found
    # unassigned taxa are present because of their entry as higher taxonomic levels
    # or because some levels are missing at the desired taxonomic level (e.g. Orders for Mollusca)
    warning(paste( "Unassigned taxa at " , tax_lev , " level, richness calculation could be biased", sep = "" ) )
  }

  # calculate richness
  nres <- apply( res[ , -1 , drop = F ] , 2 , FUN = function( x ){ length( x[ x > 0 ] ) } )
  nres

  # # Taxa cannot have empty rows so richness value is returned
  # if( tax_lev == "Taxa" ) ( return( nres ) )
  #
  # ## the following bunch of code is to check if the taxonomic level at which the calculation have been made
  # ## contains empty rows (e.g. because some data was entered at family level when calculating genus richness)
  #
  # # extract the taxonomic tree discarding columns corresponding to samples (numeric columns, this is the reason why we keep only non-numeric columns)
  # x.tree <- x[[ "Tree" ]][ , sapply( x[[ "Tree" ]] , is.factor ) ]
  #
  # # remove the column called taxa
  # x.tree <- x.tree[ , ! colnames( x.tree ) %in% "Taxa" ]
  #
  # # count and identify non null columns
  # null.col <- apply( x.tree , 2 , function( x ) sum( x == "" ) )
  # null.col <- rev( null.col == nrow( x.tree ) )
  #
  # # find the first non null taxonomic level
  # fnn.tl <- which( null.col == FALSE )[ 1 ]
  #
  #
  # if( tax_lev == names( fnn.tl ) ){
  #   return( nres )
  # } else {
  #   keep.col <- which( colnames( x.tree ) == tax_lev ) : which( colnames( x.tree ) == names( fnn.tl ) )
  #   keep.row <- apply( x.tree[ , keep.col , drop = F ], 2 , function( x )( which( x  != "" ) )  )
  #   def.length <- length( unique( unlist( keep.row ) ) )
  #   if( def.length <= sum( x.tree[ , tax_lev ] != "" ) ){
  #     return( nres )
  #   } else {
  #     warning(paste( "Unassigned taxa at " , tax_lev , " level, richness calculation could be biased", sep = "" ) )
  #     return( nres )
  #   }
  # }
}
