#' @name allindices
#' @title Functions for calculating diversity, eveness and dominance indices
#'
#' @description Functions for calculating shannon, simpson, margalef, menhinick, pielou and other indices.
#' @aliases  shannon simpson margalef menhinick pielou berpar brill esimpson invberpar invsimpson mcintosh allindices
#' @param x results of function aggregatoR
#' @param base the base of the logarithm
#' @param taxLev taxonimc level on which the calculation has to be made.
#' @keywords shannon, simpson, margalef, menhinick, pielou
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' allindices(data.agR)
#' shannon(data.agR)
#'
#' # natural logarithm
#' shannon(data.agR, base=exp(1))
#' @export simpson
#' @export esimpson
#' @export invsimpson
#' @export margalef
#' @export menhinick
#' @export pielou
#' @export berpar
#' @export invberpar
#' @export brill
#' @export mcintosh
#' @export shannon

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
  f9 <- function( x ) ( pielou( x , taxLev = taxLev, base = base ) )
  f10 <- function( x ) ( simpson( x , taxLev = taxLev) )
  f11 <- function( x ) ( esimpson(x , taxLev = taxLev) )

  funs <- list( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11)
  indices <- lapply(funs, function(f) f(x))
  # from list to data.frame
  res <- data.frame( t( matrix( unlist( indices ) , ncol = length( st.names ), byrow = T ) ) )
  rownames( res ) <- st.names
  colnames( res ) <- c( "shannon", "berpar", "brill", "invberpar", "invsimpson",
                        "margalef", "mcintosh", "menhinick", "pielou",
                        "simpson", "esimpson")
  res

}
