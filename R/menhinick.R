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
