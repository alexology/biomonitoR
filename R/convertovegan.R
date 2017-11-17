#' convertovegan
#'
#' An utility function to export data to vegan like format
#' @param x results of aggregatoR function
#' @param taxLev taxonomic level of interest. Possible choices are Phylum, Class, Subclass, Order, Family, Subfamily, Tribus, Genus, Species, Subspecies, Taxa
#' @keywords convertovegan
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.cv <- convertovegan(data.agR, taxLev = "Family")


convertovegan <- function(x, taxLev = "Family"){
  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }

  # vector of possible taxonomic levels
  taxa.vec <- c("Phylum",	"Class",	"Subclass",	"Order",	"Family",	"Subfamily",	"Tribus",
                "Genus",	"Species",	"Subspecies",	"Taxa")

  if(!taxLev %in% taxa.vec){
    stop("Please provide a valid taxLev")
  }

  x <- x[[ taxLev ]]
  st.names <- names(x)[ -1 ]
  taxa.names <- x[ , 1 ]

  if("unassigned" %in% taxa.names){
    # position of unassigned
    un.numb <- which("unassigned" %in% taxa.names)
    x <- x[-un.numb,]
    x.t <- as.data.frame( t (x [, st.names ] ) )
    names(x.t) <- taxa.names[-un.numb]
    return(x.t)
  }

  else{
    x.t <- as.data.frame( t (x [, st.names ] ) )
    names(x.t) <- taxa.names
    return(x.t)
  }
}
