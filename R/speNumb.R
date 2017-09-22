#' speNumb
#'
#' Functions for calculating species, genus, family and order Richness and abundance
#' @aliases  speNumb  genNumb famNumb ordNumb abu
#' @param x results of function aggregatoR
#' @keywords speNumb, genNumb, famNumb, ordNumb
#' @details By now only species, genus and family richness calculation are reliable. This is because order assignment for order in the reference database is not completely covered. Unassigned taxon are exluded from the calculations.
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' genNumb(data.agR)
#' @export taxNumb
#' @export genNumb
#' @export famNumb
#' @export ordNumb
#' @export claNumb


speNumb <- function(x){

  # check if the object d is of class "biomonitoR"

  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }


  spe <- x[["Species"]]
  if("unassigned" %in% spe[,1]){
    z <- which(spe$Species=="unassigned")
    spe <- spe[-z,] # remove unassigned row from the species count
  }
  nspe <- apply(spe[,-1], 2, FUN=function(x){length(x[x>0])})
  return(nspe)
}








