#' combTaxa
#'
#' This function calculates the absolute or relative abundance of a Taxon or of a set Taxa. 
#' @param x results of function aggregatoR.
#' @param taxa a Taxon or a vector of taxa.
#' @param rel if TRUE calculates relative abundance. default =F.
#' @keywords aggregatoR
#' @details
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' combTaxa(data.agr, taxLev = "Family")



combTaxa <- function(x, ntaxa = 2, taxLev = "Family", rel = F){
  
  # check if the object d is of class "biomonitoR"
  
  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }
  
  df <- x[[taxLev]]
  
  if("unassigned" %in% df[,1]){
    z <- which(df[,1]=="unassigned")
    df <- df[-z,] # remove unassigned row from the species count
  }
  
  cbn <- combn(nrow(df), ntaca)
  df.sub <- lapply(seq(ncol(cbn)), function(x) df[cbn[,x],])
  
  if(rel = T){
    df.agg <- lapply(df.sub, FUN = agg_fun)
  } else {
    df.agg <- lapply(df.sub, FUN = agg_fun)
  }
    
  return(do.call(rbind, df.agg))
}