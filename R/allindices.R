#' @describeIn shannon calculate all indices at once

allindices <- function(x, taxLev = "Family", base = 2){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  st.names <- names( x[["Taxa"]][ , -1 , drop = FALSE ] )

  f1 <- function( x ) ( shannon( x , taxLev = taxLev, base = base ) )
  f2 <- function( x ) ( berpar( x , taxLev = taxLev) )
  f3 <- function( x ) ( brill( x , taxLev = taxLev) )
  f4 <- function( x ) ( invberpar( x , taxLev = taxLev) )
  f5 <- function( x ) ( invsimpson( x , taxLev = taxLev) )
  f6 <- function( x ) ( margalef( x , taxLev = taxLev) )
  f7 <- function( x ) ( mcintosh( x , taxLev = taxLev) )
  f8 <- function( x ) ( menhinick( x , taxLev = taxLev) )
  f9 <- function( x ) ( odum( x , taxLev = taxLev, base = base ) )
  f10 <- function( x ) ( pielou( x , taxLev = taxLev, base = base ) )
  f11 <- function( x ) ( simpson( x , taxLev = taxLev) )
  f12 <- function( x ) ( esimpson(x , taxLev = taxLev) )

  funs <- list( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12)
  indices <- lapply(funs, function(f) f(data.agR))
  # from list to data.frame
  res <- data.frame( t( matrix( unlist( indices ) , ncol = length( st.names ), byrow = T ) ) )
  rownames( res ) <- st.names
  colnames( res ) <- c( "shannon", "berpar", "brill", "invberpar", "invsimpson",
                        "margalef", "mcintosh", "menhinick", "odum", "pielou",
                        "simpson", "esimpson")
  res

}
