capWords <- function( x ) {
  # see ?tolower
  cap <- paste( toupper( substring( x , 1 , 1 ) ), tolower( substring( x , 2 ) ),
        sep = "" , collapse = " " )
}


trim <- function( x )( as.factor( trimws( as.character( x ) ) ) )

traitS <- function( x , y , z , w ){
  DFtaxa <- as.character( x[  11 ] )
  trait.final <- data.frame( Taxa_db = character() , Traits_real = character() , z[-c(1:nrow(z)) , ] , stringsAsFactors = FALSE)


  taxa_tree <- y[ y[  , "Taxa" ] == DFtaxa , 1:10 ]
  taxa_tree <- rev( taxa_tree[ taxa_tree != "" ] )
  for( i in 1:length( taxa_tree) ){
    if( any( z$Taxa == taxa_tree[ i ] ) ){
      temp <- z[ z$Taxa == taxa_tree[ i ] , ]
      name.temp <- w[ z$Taxa == taxa_tree[ i ] ]
      nr.temp <- nrow( temp )
      nr.trait.final <- nrow(trait.final)
      trait.final[ (nr.trait.final + 1):(nr.trait.final + nr.temp )  , -c( 1 : 2 ) ] <- rbind( trait.final[ , - 1 ] , temp )
      trait.final[ (nr.trait.final + 1):(nr.trait.final + nr.temp )  , 1 ] <- rep( DFtaxa, nr.temp )
      trait.final[ (nr.trait.final + 1):(nr.trait.final + nr.temp )  , 2 ] <- name.temp
      if( nrow(temp) > 0 ){ break }
    } else{ next }
  }
  return( trait.final )
}
