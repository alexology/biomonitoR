#' bioco
#'
#' bioco function provides site-specific biocontamination index (SBCI), abundance contamination index (ACI) and richness contamination index (RCI) at family rank according to those proposed by Arbaciauskas et al. (2008).
#'
#' @param x result of function aggregatoR.
#' @param alien a vector containing the alien taxa name. Only taxa up to family level will be considered.
#' @param dfref custom reference database if other than default is used. Should be the same of that used in asBiomonitor.
#' @details Biocontamination of sampling sites was assessed using a site-specific biocontamination index (SBCI) derived from two metrics: abundance contamination index (ACI) and richness contamination index (RCI) at family/ordinal rank. These indices were calculated as:
#' \deqn{ACI = Na/Nt}
#' where Na and Nt are numbers of specimens of alien taxa and total specimens in a sample, respectively, and
#' \deqn{RCI = Ta/Tt}
#' where Ta is the total number of alien families/orders, and Tt is the total number of identified families/orders. With values of ACI and RCI, the site-specific biocontamination index (SBCI) can then be derived from an easy double entry matrix (Arbaciauskas et al. 2008). Five classes of biocontamination ranging from 0 ("no" contamination) to 4 ("severe" contamination) are defined.
#' A global alien species dataset cannot be provided because alien species may vary among different biogeographical regions and countries. For this reason users need to define and properly upload their own alien reference database based on the studied area they want to study.
#' Applications and examples of these indices are available in Cuk et al. (2019) for Croatian rivers and MacNeil et al. 2010 (The Isle of Man).
#' @keywords aggregatoR
#' @references Arbaciauskas, K., Semenchenko, V., Grabowski, M., Leuven, R.S.E.W., Paunovic, M., Son, M.O., Csanyi, B., Gumuliauskaite, S., Konopacka, A., Nehring, S., Van der Velde, G., Vezhnovetz, V., Panov, V.E. (2008). Assessment of biocontamination of benthic macroinvertebrate communities in European inland waterways. Aquatic Invasions 3 (2): 211-230.
#' @references Cuk, R., Milisa, M., Atanackovic, A., Dekic, S., Blazekovic, L., & Zganec, K. (2019). Biocontamination of benthic macroinvertebrate assemblages in Croatian major rivers and effects on ecological quality assessment. Knowledge & Management of Aquatic Ecosystems, (420), 11.
#' @references MacNeil, C., Briffa, M., Leuven, R.S., Gell, F.R. and Selman, R. (2010). An appraisal of a biocontamination assessment method for freshwater macroinvertebrate assemblages; a practical way to measure a significant biological pressure? Hydrobiologia 638:151-159
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data( macro_ex )
#' data.bio <- asBiomonitor( macro_ex )
#' data.agR <- aggregatoR( data.bio )
#' # Toy example:
#' alien <- c( "Laccobius" , "Setodes bulgaricus" , "Caenidae" )
#' bioco( data.agR , alien = alien )

bioco <- function( x , alien = NULL , dfref = NULL ){

  #  check if the object x is of class "biomonitoR"
  classCheck( x )

  # format names with the first letter as capital letter and the others lowercase

  alien <- sapply( alien , capWords , USE.NAMES = FALSE)

  # if dfref = NULL check if x is of class biomonitor mi, mf or fi, otherwise dfref is needed

  if( is.null( dfref )){
    if( ! any( class( x ) %in% c( "mi" , "mf" , "fi") ) ) stop( ( "Custom reference database you used in aggregatoRis needed" ) )
    if( any( class( x ) %in% "mi"  ) ) ( ref <- mi_ref)
    if( any( class( x ) %in% "mf"  ) ) ( ref <- mf_ref)
    if( any( class( x ) %in% "fi"  ) ) ( stop("Fish reference database database not implemented yet") )
    if( ! any( class( x ) %in% c( "mi" , "mf" , "fi" ) ) ) ( stop("No default reference dataset available for class = custom") )
    } else { ref = dfref }

  # check
  if( is.null( alien ) ) ( stop( "Please provide a vector containing the names of alien taxa" ) )

  alien <- trimws( alien )
  st.names <- names( x[[ 1 ]][ -1 ] )
  DF <- ref[ , 1:10 ]

  # check if alien contains taxa don't present in the reference database
  taxa.vec <- as.character( unique( unlist( DF ) ) )
  if( any( ! alien %in% taxa.vec ) ){
    absent <- alien[ ! alien %in% taxa.vec ]
    absent <- paste( absent , collapse = ", ")
    mes <- paste( "The following taxa are absent from the reference database:" , absent , sep = " " )
    stop( mes )
  }

  # Position of taxon in the df data.frame
  taxind <- data.frame( row = numeric() , col = numeric( ) )
  for(i in 1:length( alien ) ){
    temp <- which( DF == alien[ i ], arr.ind = T )
    taxind <- rbind( temp , taxind )
  }

  getAlienAll <- c()
  for( i in 1:nrow( taxind ) ){
    a <- taxind[ i , 1 ]
    b <- taxind[ i , 2 ]:10
    temp <- as.character( unlist( DF[ a , b ] )  )
    getAlienAll <- c( getAlienAll , temp )
  }

  getAlienAll <- unique( getAlienAll )
  getAlienAll <- getAlienAll[ getAlienAll != "" ]

  x.taxa <- x[[ "Tree" ]]
  x.taxa <- x.taxa[ x.taxa[ , "Taxa" ] %in% getAlienAll , ]
  x.taxa <- x.taxa[ , colnames( x.taxa ) %in% c( "Family" , st.names ) ]

  abu.alien <- apply( x.taxa[ , -1 , drop = FALSE ] , 2 , sum )
  tax.alien <- apply( x.taxa[ , -1 , drop = FALSE ] , 2 , function( x ) sum( x > 0 ) )

  aci <- round( abu.alien / abundance( x , taxLev = "Taxa" ) , 2 )
  rci <- round( tax.alien / richness( x , taxLev = "Family" ) , 2 )
  cl.lim <- c( 1 , 0.5 , 0.2 , 0.1 , 0.01 , 0)
  cl.lab <- c( 0:4 )
  cl.abu <- cut( aci , cl.lim  , cl.lab , right = TRUE , include.lowest = T )
  cl.tax <- cut( rci , cl.lim  , cl.lab , right = TRUE , include.lowest = T )
  cl <- data.frame( as.numeric( as.character( cl.abu ) ) , as.numeric( as.character( cl.tax ) ) )
  data.frame( aci = aci , rci = rci , sbci = apply(  cl , 1 , max  ) )
}


