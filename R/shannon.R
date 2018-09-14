#' shannon
#'
#' Functions for calculating shannon, simpson, margalef, menhinick, pielou and other indices.
#' @aliases  shannon simpson margalef menhinick pielou berpar brill esimpson invberpar invsimpson mcintosh odum
#' @param x results of function aggregatoR
#' @param base the base of the logarithm
#' @param taxLev taxonimc level on which the calculation has to be made.
#' @keywords shannon, simpson, margalef, menhinick, pielou
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' shannon(data.agR)
#'
#' # natural logarithm
#' shannon(data.agR, base=exp(1))
#' @export simpson
#' @export esimpson
#' @export invsimpson
#' @export margalef
#' @export menhinick
#' @export pielou
#' @export berpar
#' @export invberpar
#' @export brill
#' @export mcintosh
#' @export odum

shannon <- function(x, base=2, taxLev="Family"){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

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
  return( sha )
}

