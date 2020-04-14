#' aggregatoR
#'
#' This function prepares data for further calculations.
#' @param x results of function asBiomonitoR
#' @param FUN the function to be applied for aggregating to higher taxonomic levels. Default to `sum`.
#' @keywords aggregatoR
#' @importFrom stats aggregate
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#'
#' # example for macrophytes
#' data(oglio)
#'
#' oglio.asb <- asBiomonitor(oglio, group = "mf")
#' oglio.agg <- aggregatoR(oglio.asb)
#' richness(oglio.agg , taxLev = "Species")
#' richness(oglio.agg , taxLev = "Genus")
#' richness(oglio.agg , taxLev = "Family")

aggregatoR <- function ( x , FUN = sum )
{

  # check if the object x is of class "biomonitoR"
  classCheck( x )

  # The following part is to aggregate the dataset provided by the user at different taxonomic levels
  tx <- c( "Phylum" , "Class" , "Subclass" , "Order" , "Family" , "Subfamily" , "Tribus" , "Genus" , "Species" , "Subspecies" , "Taxa" )
  stz <- x[ ! ( names( x ) %in% tx ) ]
  check.pa <- any( as.data.frame( stz )  != 1 )
  stz_n <- names( stz )
  phy.agg <- aggregate( stz , by = list( x$Phylum ) , FUN )
  if( ! check.pa ) ( phy.agg <- data.frame( phy.agg[ , 1  , drop = FALSE ] , ( phy.agg[ , -1 , drop = FALSE ] > 0 ) * 1 ) )
  levels( phy.agg$Group.1 )[ levels( phy.agg$Group.1 ) == "" ] <- "unassigned"
  names( phy.agg ) <- c( "Phylum" , stz_n )
  cla.agg <- aggregate( stz , by = list( x$Class ) , FUN )
  if( ! check.pa ) ( cla.agg <- data.frame( cla.agg[ , 1 , drop = FALSE], ( cla.agg[ , -1, drop = FALSE ] > 0 ) * 1 ) )
  levels( cla.agg$Group.1 )[ levels( cla.agg$Group.1 ) == "" ] <- "unassigned"
  names( cla.agg ) <- c( "Class" , stz_n )
  scla.agg <- aggregate( stz , by = list( x$Subclass ), FUN )
  if( ! check.pa ) ( scla.agg <- data.frame( scla.agg[ , 1 , drop = FALSE ] , ( scla.agg[ , -1, drop = FALSE ] > 0 ) * 1 ) )
  levels( scla.agg$Group.1 )[ levels( scla.agg$Group.1 ) == "" ] <- "unassigned"
  names( scla.agg ) <- c( "Subclass" , stz_n )
  ord.agg <- aggregate( stz , by = list( x$Order ), FUN )
  if( ! check.pa ) ( ord.agg <- data.frame( ord.agg[ , 1 , drop = FALSE ] , ( ord.agg[ , -1 , drop = FALSE ] > 0 ) * 1 ) )
  levels( ord.agg$Group.1 )[ levels( ord.agg$Group.1 ) == "" ] <- "unassigned"
  names( ord.agg ) <- c( "Order" , stz_n )
  fam.agg <- aggregate( stz , by = list( x$Family ) , FUN )
  if( ! check.pa ) ( fam.agg <- data.frame( fam.agg[ , 1 , drop = FALSE ], ( fam.agg[ , -1, drop = FALSE ] > 0 ) * 1 ) )
  levels( fam.agg$Group.1 )[ levels( fam.agg$Group.1 ) == "" ] <- "unassigned"
  names( fam.agg ) <- c( "Family" , stz_n )
  sfam.agg <- aggregate( stz , by = list( x$Subfamily ) , FUN )
  if( ! check.pa ) ( sfam.agg <- data.frame( sfam.agg[ , 1 , drop = FALSE ], ( sfam.agg[ , -1, drop = FALSE ] > 0 ) * 1 ) )
  levels( sfam.agg$Group.1 )[ levels( sfam.agg$Group.1 ) == "" ] <- "unassigned"
  names( sfam.agg ) <- c( "Subfamily" , stz_n )
  tri.agg <- aggregate( stz , by = list( x$Tribus ) , FUN )
  if( ! check.pa ) ( tri.agg <- data.frame( tri.agg[ , 1 , drop = FALSE ], ( tri.agg[ , -1, drop = FALSE ] > 0 ) * 1 ) )
  levels(tri.agg$Group.1)[ levels( tri.agg$Group.1 ) == "" ] <- "unassigned"
  names( tri.agg ) <- c( "Tribus" , stz_n )
  gen.agg <- aggregate( stz , by = list( x$Genus ), FUN )
  if( ! check.pa ) ( gen.agg <- data.frame( gen.agg[ , 1 , drop = FALSE ], ( gen.agg[ , -1, drop = FALSE ] > 0 ) * 1 ) )
  levels( gen.agg$Group.1 )[ levels( gen.agg$Group.1 ) == "" ] <- "unassigned"
  names( gen.agg ) <- c( "Genus" , stz_n )
  spe.agg <- aggregate( stz , by = list( x$Species ), FUN )
  if( ! check.pa ) ( spe.agg <- data.frame( spe.agg[ , 1 , drop = FALSE ], ( spe.agg[ , -1, drop = FALSE ] > 0 ) * 1 ) )
  levels(spe.agg$Group.1)[ levels( spe.agg$Group.1 ) == "" ] <- "unassigned"
  names( spe.agg ) <- c( "Species" , stz_n )
  sspe.agg <- aggregate(stz , by = list( x$Subspecies ) , FUN )
  if( ! check.pa ) ( sspe.agg <- data.frame( sspe.agg[ , 1 , drop = FALSE ], ( sspe.agg[ , -1, drop = FALSE ] > 0 ) * 1 ) )
  levels( sspe.agg$Group.1 )[ levels( sspe.agg$Group.1 ) == "" ] <- "unassigned"
  names( sspe.agg ) <- c( "Subspecies" , stz_n )
  tax.agg <- aggregate( stz, by = list( x$Taxa ), FUN )
  if( ! check.pa ) ( tax.agg <- data.frame( tax.agg[ , 1 , drop = FALSE ], ( tax.agg[ , -1, drop = FALSE ] > 0 ) * 1 ) )
  levels( tax.agg$Group.1 )[ levels( tax.agg$Group.1 ) == "" ] <- "unassigned"
  names( tax.agg ) <- c( "Taxon" , stz_n )
  tree.agg <- aggregate( stz , by = list( x$Phylum , x$Class , x$Subclass , x$Order , x$Family , x$Subfamily ,
                                       x$Tribus , x$Genus , x$Species , x$Subspecies , x$Taxa ) , FUN )
  names(tree.agg) <- c( tx , stz_n )
  temp <- list( phy.agg , cla.agg , scla.agg , ord.agg , fam.agg , sfam.agg , tri.agg , gen.agg , spe.agg , sspe.agg , tax.agg , tree.agg )
  names( temp ) <- c( "Phylum" , "Class" , "Subclass" , "Order" , "Family" , "Subfamily" , "Tribus" , "Genus" , "Species" , "Subspecies" , "Taxa" , "Tree" )
  class( temp ) <- class( x )
  temp
}
