#' convertovegan
#'
#' An utility function to export data to other formats.
#' @param x result of the aggregatoR function
#' @param taxLev taxonomic level of interest. Possible choices are Phylum, Class, Subclass, Order, Family, Subfamily, Tribus, Genus, Species, Subspecies, Taxa
#' @keywords convertovegan
#' @details `convertovegan` converts data to the `vegan` package style, with sites on rows and taxa on columns. `convertobiotic` converts data to
#' the format of `biotic` package. This function will extract information at the family level, with an exception for Oligochaeta.
#' In `convertobiotic` if both family and above family level taxonomic information of Oligochaeta are concurrently present,
#' information will be provided at Oligochaeta level.
#' `biomonitoR` automatically identify the taxonomic level at which Oligochaeta are stored (subclass is the default in `biomonitoR`).
#' Be careful, for all the other taxa entered at taxonomic level higher than family will be discarded (e.g. Trombidiformes).
#' @export
#' @export convertobiotic
#' @seealso \code{\link{asBiomonitor}}

convertovegan <- function(x, taxLev = "Family"){

  .Deprecated("convert_to_vegan", package = "biomonitoR")

  # check if the object x is of class "biomonitoR"
  classCheck( x )

  # useful for transforming data to 0-1 later
  if( inherits( x , "bin" ) ){
    BIN <- TRUE
  } else { BIN <- FALSE }

  # vector of possible taxonomic levels
  taxa.vec <- c( "Phylum" ,	"Class" ,	"Subclass" ,	"Order" ,	"Family" ,	"Subfamily" ,	"Tribus" ,
                "Genus" ,	"Species" ,	"Subspecies" ,	"Taxa" )

  # taxLev must be one of the tax.vec levels
  if( ! taxLev %in% taxa.vec ){
    stop( "Please provide a valid taxLev" )
  }

  # get the data.frame at the specified taxonomic level
  x <- x[[ taxLev ]]
  st.names <- names( x )[ -1 ]
  taxa.names <- x[ , 1 ]

  # transform the data.frame from abundance to presence-absence if needed
  if( BIN ){
    x <- to_bin( x )
  }

  # if unassigned are present remove them, otherwise return the x after having changed column names
  if( "unassigned" %in% taxa.names ){
    # position of unassigned
    un.numb <- which( "unassigned" %in% taxa.names )
    x <- x[ -un.numb , ]
    x.t <- as.data.frame( t ( x[ , st.names ] ) )
    names( x.t ) <- taxa.names[ -un.numb ]
    x.t
  }

  else{
    x.t <- as.data.frame( t ( x[ , st.names ] ) )
    names( x.t ) <- taxa.names
    x.t
  }
}
