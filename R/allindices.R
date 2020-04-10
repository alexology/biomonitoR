#' @name allindices
#' @title Functions for calculating diversity, eveness and dominance indices.
#'
#' @description Functions for calculating shannon, simpson, margalef, menhinick, pielou and other indices.
#' @aliases  shannon simpson margalef menhinick pielou berpar brill esimpson invberpar invsimpson mcintosh fisher allindices
#' @param x result of the function aggregatoR
#' @param base the base of the logarithm
#' @param taxLev taxonomic level on which the calculation has to be made.
#' @details Shannon index:
#' \deqn{H'=-\sum_{i=1}^{S} p_{i}\log(p_{i})}
#' Pielou index:
#' \deqn{J'=\frac{H'}{\log(S)}}
#' Margalef diversity:
#' \deqn{D_{mg}=\frac{(S - 1)}{\log(N)}}
#' Menhinick diversity:
#' \deqn{D_{mn}=\frac{S}{\sqrt(N)}}
#' Brillouin index:
#' \deqn{HB=\frac{\log(N) -\sum(\log(n_{i})) }{N}}
#' Simpson's index is calculated as D. biomonitoR returns 1 - D and 1/D when `index` is set to `Simpson` or `Invsimpson`.
#' \deqn{D=\sum_{i=1}^{S} p_{i}^{2}}
#' Simpson's evenness:
#' \deqn{D_{1/D}=\frac{1/D}{\sqrt(S)}}
#' Berger-Parker index:
#' \deqn{d=\frac{N_{max}}{N}}
#' The inverse of this index is also provided when `index` is set to `Invberpar`. McIntosh's diversity:
#' \deqn{D=\frac{N-U}{N-\sqrt(N)}}
#' where
#' \deqn{U=\sqrt(\sum_{i=1}^{S} n_{i}^{2})}
#' Fisher alpha is calculated as follow:
#' \deqn{\alpha=\frac{N(1-x)}{x}}
#' where x is estimated from the iterative solution of:
#' \deqn{\frac{S}{N}=-\frac{N(1-x)}{x} \log(1-x)}
#'
#'
#' p_i is the proportion of individuals found in the i-th species, n_i is the number of individuals found in the i-th species, S the species richness,
#'  N the number of individuals and N_max the number of individuals in the most abundant species. All the indices are calculated according to Magurran (2004).
#' @keywords shannon, simpson, margalef, menhinick, pielou
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @references Magurran, A. E. (2004). Measuring biological diversity. Blackwell Science ltd.
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' allindices(data.agR)
#' shannon(data.agR)
#'
#' # base 2
#' shannon(data.agR, base = 2)
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
#' @export fisher

allindices <- function( x, taxLev = "Taxa", base = exp( 1 ) ){

  # check if the object x is of class "biomonitoR"
  classCheck( x )

  st.names <- names( x[[ "Taxa" ]][ , -1 , drop = FALSE ] )

  f1 <- function( x ) ( shannon( x , taxLev = taxLev , base = base ) )
  f2 <- function( x ) ( berpar( x , taxLev = taxLev ) )
  f3 <- function( x ) ( brill( x , taxLev = taxLev ) )
  f4 <- function( x ) ( invberpar( x , taxLev = taxLev ) )
  f5 <- function( x ) ( invsimpson( x , taxLev = taxLev ) )
  f6 <- function( x ) ( margalef( x , taxLev = taxLev ) )
  f7 <- function( x ) ( mcintosh( x , taxLev = taxLev ) )
  f8 <- function( x ) ( menhinick( x , taxLev = taxLev ) )
  f9 <- function( x ) ( pielou( x , taxLev = taxLev, base = base ) )
  f10 <- function( x ) ( simpson( x , taxLev = taxLev ) )
  f11 <- function( x ) ( esimpson( x , taxLev = taxLev ) )
  f12 <- function( x ) ( fisher( x , taxLev = taxLev ) )

  funs <- list( f1 , f2 , f3 , f4 , f5 , f6 , f7 , f8 , f9 , f10 , f11 , f12 )
  indices <- lapply( funs , function( f ) f( x ) )
  # from list to data.frame
  res <- data.frame( t( matrix( unlist( indices ) , ncol = length( st.names ), byrow = T ) ) )
  rownames( res ) <- st.names
  colnames( res ) <- c( "shannon", "berpar", "brill", "invberpar", "invsimpson",
                        "margalef", "mcintosh", "menhinick", "pielou",
                        "simpson", "esimpson", "fisher")
  res

}
