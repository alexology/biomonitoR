#' convertovegan
#'
#' An utility function to export data to other formats.
#' @param x results of aggregatoR function
#' @param taxLev taxonomic level of interest. Possible choices are Phylum, Class, Subclass, Order, Family, Subfamily, Tribus, Genus, Species, Subspecies, Taxa
#' @keywords convertovegan
#' @details `convertovegan` converts data to a vegan style format, with sites on rows and taxa on columns. `convertovegan` converts data to
#' the format of biotic package. For Oligochaeta, if both family and above family level taxonomic information are concurrently present,
#' information will be provided at Oligochaeta level.
#' @export
#' @export convertobiotic
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.cv <- convertovegan(data.agR, taxLev = "Family")


convertovegan <- function(x, taxLev = "Family"){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  # vector of possible taxonomic levels
  taxa.vec <- c("Phylum",	"Class",	"Subclass",	"Order",	"Family",	"Subfamily",	"Tribus",
                "Genus",	"Species",	"Subspecies",	"Taxa")

  if(!taxLev %in% taxa.vec){
    stop("Please provide a valid taxLev")
  }

  x <- x[[ taxLev ]]
  st.names <- names( x )[ -1 ]
  taxa.names <- x[ , 1 ]

  if( "unassigned" %in% taxa.names ){
    # position of unassigned
    un.numb <- which( "unassigned" %in% taxa.names )
    x <- x[ -un.numb , ]
    x.t <- as.data.frame( t (x [, st.names ] ) )
    names( x.t ) <- taxa.names[ -un.numb ]
    return( x.t )
  }

  else{
    x.t <- as.data.frame( t (x [, st.names ] ) )
    names( x.t ) <- taxa.names
    return( x.t )
  }
}
