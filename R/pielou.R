#' @describeIn allindices Pielou's evenness index

pielou <- function( x , base = exp( 1 ) , tax_lev = "Taxa" ){

  # check if the object x is of class "biomonitoR"
  classCheck( x )

  # get the data.frame at the desired taxonomic level
  DF <-  x[[ tax_lev ]]

  if( inherits( x , "bin" ) ){
    DF <- to_bin( DF )
  }

  # remove unassigned row from the species count if present
  if( "unassigned" %in% DF[ , 1 ] ){
    z <- which( DF[ , 1 ] == "unassigned" )
    DF <- DF[ -z , ]
  }

  # apply the function for calculating the Pielou's evenness index
  rich <- apply( DF[ , -1 , drop = F ] , 2 , FUN = function( x ){ sum( x > 0 ) } )
  Pie <- apply( DF[ , -1 , drop = F ] , 2 , FUN = function( x ){ Pi( x , index = "Shannon" , base = base ) } ) / log( rich, base = base )
  Pie
}
