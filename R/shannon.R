#' @describeIn allindices Shannon index

shannon <- function(x, base=2, taxLev="Family"){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  df <-  x[[taxLev]]
  if("unassigned" %in% df[,1]){
    z <- which(df[,1]=="unassigned")
    df<- df[-z,] # remove unassigned row from the species count
  }

  sha <- apply(df[ , -1 , drop = FALSE], 2, FUN = function(x){Pi( x, index = "Shannon", base = base )})
  return( sha )
}

