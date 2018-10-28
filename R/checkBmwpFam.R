#' @importFrom stats aggregate
checkBmwpFam <- function(df, famNames, stNames){
  Taxon <- NULL
  df[ , "Taxon" ] <- as.character( df[ , "Taxon" ] )
  famCheck <- df[ which( df[ , "Taxon" ] %in% famNames[ , "Taxon" ] ), "Taxon" ]
  if( length( famCheck ) > 0 ){

    for( i in 1 : length( famCheck ) ){
      taxName <- famCheck[ i ]
      df$Taxon[ df$Taxon == taxName ] <- as.character( subset( famNames, Taxon == taxName)[,2])
    }

    df <- aggregate( . ~ Taxon , df , sum )
    names( df ) <- c( "Taxon" , stNames)
    return(df)
  }
  else{df <- df
  return(df)
  }
}
