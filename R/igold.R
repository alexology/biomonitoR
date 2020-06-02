#' igold
#'
#' This function calculates the 1 - GOLD metric, where GOLD stands for Gastropoda, Oligochaeta and Diptera. This metric should decrease with increasing organic pollution (Pinto et al., 2004).
#' @param x results of function aggregatoR
#' @param traceB if set to `TRUE` a list as specified below will be returned.
#' @keywords ept
#' @details The metric 1 - GOLD is calculated as 1 minus the relative abundance of Gastropoda, Oligochaeta and Diptera. If a custom database is provided (see \code{\link{aggregatoR}}) please be sure that Gastropoda are at Class, Oligochaeta at Sublclass and Diptera at Order level, otherwise the gold calculation will be meaningless.
#' If this is the case please see [abuTax].
#'
#' @return If `traceB` is set to `TRUE` a list with the following elements will be returned:
#' \itemize{
#'  \item `results` Results of the `igold` index.
#'  \item `taxa_df` The data.frame used for the calculation containing the abundance of the GOLD taxa.
#' }
#'
#' @importFrom stats aggregate
#' @references Pinto, P., Rosado, J., Morais, M., & Antunes, I. (2004). Assessment methodology for southern siliceous basins in Portugal. In Integrated Assessment of Running Waters in Europe (pp. 191-214). Springer, Dordrecht.
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' igold(data.agR)


igold <- function( x , traceB = FALSE ){

  # check if the object x is of class "biomonitoR"
  classCheck( x )

  # basic operation on th asBiomonitor object. Extract Tree and samples from the Tree element of the asBiomonitor object
  x_gold <- x[[ "Tree" ]]
  tx <- c( "Phylum" , "Class" , "Subclass" , "Order" , "Family" , "Subfamily" , "Tribus" , "Genus" , "Species" , "Subspecies" , "Taxa" )
  # get sample names
  st_names <- names( x_gold[ ! ( names( x_gold ) %in% tx ) ] )
  gold_taxa <- x_gold[ x_gold$Class == "Gastropoda" | x_gold$Subclass == "Oligochaeta" | x_gold$Order == "Diptera" , , drop=F]

  gold_temp <- gold_taxa[ , st_names , drop = FALSE ]

  if( inherits( x , "bin" ) ){
    gold_temp[ gold_temp > 0 ] <- 1
  }

  temp <- 1 - apply( gold_temp[ , , drop = FALSE ], 2 , sum ) / abundance( x , "Taxa" , unassigned = TRUE )
  if( traceB == FALSE ){
    temp
  } else {
    list( results = temp , taxa_df = gold_taxa[ , c( "Taxa" , st_names ) , drop = FALSE ] )
  }

}
