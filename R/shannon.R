#' shannon
#'
#' Functions for calculating shannon, simpson, margalef and menhinick indexes.
#' @aliases  shannon simpson margalef menhinick
#' @param x results of function aggregatoR
#' @param base the base of the logarithm
#' @param taxaLev taxonimc level on which the calculation has to be made.
#' @keywords shannon, simpson, margalef, menhinick
#' @details
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' shannon(data.agR)
#' @export simpson
#' @export margalef
#' @export margalef
#' @export menhinick


shannon <- function(x, base=2, taxLev="Family"){
  
  # check if the object d is of class "biomonitoR"
  
    if (class(z) != "biomonitoR") {
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      return("Object x is not an object of class biomonitoR")
    }
  
  if(class(x)!="biomonitoR"){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }
  df <-  x[[taxLev]]
  if("unassigned" %in% df[,1]){
    z <- which(df[,1]=="unassigned")
    df<- df[-z,] # remove unassigned row from the species count
  }
  if(ncol(df)==2){
    sha <- Pi(df[,-1], index = "Shannon", base = base)
  }
  else{
    sha <- apply(df[,-1], 2, FUN = function(x){Pi( x, index = "Shannon", base = base )})
  }
  return( round(sha, 3) )
}

