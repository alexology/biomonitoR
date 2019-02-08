#' totInfo
#'
#' Function to calculate the overall richness at different taxonomic level. 
#' @param x results of function asBiomonitoR
#' @param taxalist if true returns the list of taxa for each taxonomic level
#' @keywords totInfo
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' totInfo( data.agR )


totInfo <- function( x , taxalist = FALSE ){
  
  # check if the object x is of class "biomonitoR"
  classCheck(x)
  
  x <- x[[ "Tree" ]]
  tx <- c("Phylum", "Class", "Subclass", "Order", "Family", "Subfamily", "Tribus", "Genus", "Species", "Subspecies", "Taxa")
  stz1 <- x[ !( names( x ) %in% tx ) ]
  stz <- apply( stz1 , 1 , sum )
  stz[stz>0] <- 1
  stz <- as.logical( stz )
  treeo <- x[ stz , names( x ) %in% tx ]
  tlist <- apply(treeo , 2 , unique )
  if( taxalist ) { tlist }
  else { if( any( stz1  > 1 ) ){
    c( unlist( lapply( tlist, function( x ) ( sum( x != "" ) ) ) ) , Abundance = sum( stz1 ) )
  } else{ unlist( lapply( tlist, function( x ) ( sum( x != "" ) ) ) ) }
     }
}

