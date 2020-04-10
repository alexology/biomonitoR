#' @describeIn allindices Pielou's evenness index

pielou <- function( x , base = exp( 1 ) , taxLev = "Taxa" ){

  # check if the object x is of class "biomonitoR"
  classCheck( x )

  # get the data.frame at the desired taxonomic level
  df <-  x[[ taxLev ]]

  # remove unassigned row from the species count if present
  if( "unassigned" %in% df[ , 1 ] ){
    z <- which( df[ , 1 ] == "unassigned" )
    df <- df[ -z , ]
  }

  # apply the function for calculating the Pielou's evenness index
  rich <- apply( df[ , -1 , drop = F ] , 2 , FUN = function( x ){ length( x[ x > 0 ] ) } )
  Pie <- apply( df[ , -1 , drop = F ] , 2 , FUN = function( x ){ Pi( x , index = "Shannon" , base = base ) } ) / log( rich, base = base )
  Pie
}
