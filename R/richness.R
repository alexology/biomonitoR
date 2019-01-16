#' @describeIn allrich Taxonomic richness

richness <- function( x , taxLev = "Family" ){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  res <- x[[ taxLev ]]
  if( "unassigned" %in% res[ , 1 ] ){
    z <- which( res[ , taxLev ] == "unassigned" )
    res <- res[-z,] # remove unassigned row from the species count
  }
  nres <- apply(res[,-1, drop = F], 2, FUN = function( x ){ length( x[ x > 0 ] ) } )

  if( taxLev == "Taxa" ) ( return( nres ) )

  x.tree <- x[[ "Tree" ]][ , sapply(x[[ "Tree" ]], is.factor) ]

  # remove the column called taxa
  x.tree <- x.tree[ , ! colnames( x.tree ) %in% "Taxa" ]

  # identify non null columns
  null.col <- apply( x.tree , 2 , function(x) sum( x == "" ) )
  null.col <- rev( null.col == nrow( x.tree ) )

  # find the first non null taxonomic level
  fnn.tl <- which( null.col == FALSE )[1]


  if( taxLev == names( fnn.tl ) ){
    return( nres )
  } else {
    keep.col <- which( colnames( x.tree ) == taxLev ) : which( colnames( x.tree ) == names( fnn.tl ) )
    keep.row <- apply( x.tree[ , keep.col , drop = F], 2 , function( x )( which( x  != "" ) )  )
    def.length <- length( unique( unlist( keep.row ) ) )
    if( def.length <= sum( x.tree[ , taxLev ] != "") ){
      return( nres )
    } else {
      warning(paste( "missing taxonomic levels in ", taxLev , ", richness calculation could be biased", sep = "") )
      return( nres )
    }
  }
}
