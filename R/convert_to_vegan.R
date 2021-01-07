#' convert_to_vegan
#'
#' @description
#' Utility functions to export data to other formats.
#'
#' @param x Result of thefunction `aggregate_taxa()`.
#' @param tax_lev taxonomic level of interest. Possible choices are Phylum, Class, Subclass, Order, Family, Subfamily, Tribus, Genus, Species, Subspecies, Taxa
#'
#' @keywords convert_to_vegan
#' @details `convert_to_vegan()` converts data to the `vegan` package style, with sites on rows and taxa on columns. `convert_to_biotic()` converts data to
#' the format of `biotic` package. This function will extract information at the family level, with an exception for Oligochaeta.
#' In `convert_to_biotic()` if both family and above family level taxonomic information of Oligochaeta are concurrently present,
#' information will be provided for Oligochaeta only.
#' `biomonitoR` automatically identify the taxonomic level at which Oligochaeta are stored (subclass is the default in `biomonitoR`).
#' Be careful, for all the other taxa entered at taxonomic level higher than family will be discarded (e.g. Trombidiformes).
#' @export
#' @export convert_to_biotic
#' @seealso \code{\link{aggregate_taxa}}
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' data_cv <- convert_to_vegan(data_agr, tax_lev = "Family")


convert_to_vegan <- function(x, tax_lev = "Family"){

  # check if the object x is of class "biomonitoR"
  classCheck( x )

  # useful for transforming data to 0-1 later
  if( inherits( x , "bin" ) ){
    BIN <- TRUE
  } else { BIN <- FALSE }

  # vector of possible taxonomic levels
  taxa.vec <- c( "Phylum" ,	"Class" ,	"Subclass" ,	"Order" ,	"Family" ,	"Subfamily" ,	"Tribus" ,
                 "Genus" ,	"Species" ,	"Subspecies" ,	"Taxa" )

  # tax_lev must be one of the tax.vec levels
  if( ! tax_lev %in% taxa.vec ){
    stop( "Please provide a valid tax_lev" )
  }

  # get the data.frame at the specified taxonomic level
  x <- x[[ tax_lev ]]
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
