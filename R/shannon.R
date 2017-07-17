#' shannon
#'
#' Functions for calculating shannon, simpson, margalef and menhinick indexes.
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


shannon <- function(x, base=2, taxLev="Family"){
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
    sha <- Pi(df[,-1], index = "Shannon", base = 2)
  }
  else{
    sha <- apply(df[,-1], 2, FUN = function(x){Pi( x, index = "Shannon", base = 2 )})
  }
  return(sha)
}

simpson <- function(x, taxLev = "Family"){
  if(class(x) != "biomonitoR"){
    stop("Object x is not a data.frame of class biomonitoR")
  }
  df <-  x[[taxLev]]
  if("unassigned" %in% df[,1]){
    z <- which(df[,1]=="unassigned")
    df<- df[-z,] # remove unassigned row from the species count
  }
  if(ncol(df)==2){
    sim <- Pi(df[,-1], index = "Simpson")
  }
  else{
    sim <- apply(df[,-1], 2, FUN = function(x){Pi( x, index = "Simpson")})
  }
  return(sim)
}


margalef <- function(x, taxLev = "Family"){
  if(class(x) != "biomonitoR"){
    stop("Object x is not a data.frame of class biomonitoR")
  }
  df <-  x[[taxLev]]
  if("unassigned" %in% df[,1]){
    z <- which(df[,1]=="unassigned")
    df<- df[-z,] # remove unassigned row from the species count
  }
  if(ncol(df)==2){
    marg <- Pi(df[,-1], index = "Margalef")
  }
  else{
    marg <- apply(df[,-1], 2, FUN = function(x){Pi( x, index = "Margalef")})
  }
  return(marg)
}


menhinick <- function(x, taxLev = "Family"){
  if(class(x) != "biomonitoR"){
    stop("Object x is not a data.frame of class biomonitoR")
  }
  df <-  x[[taxLev]]
  if("unassigned" %in% df[,1]){
    z <- which(df[,1]=="unassigned")
    df<- df[-z,] # remove unassigned row from the species count
  }
  if(ncol(df)==2){
    men <- Pi(df[,-1], index = "Menhinick")
  }
  else{
    men <- apply(df[,-1], 2, FUN = function(x){Pi( x, index = "Menhinick")})
  }
  return(men)
}

