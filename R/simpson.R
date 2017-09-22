simpson <- function(x, taxLev = "Family"){
  
  # check if the object d is of class "biomonitoR"
  
    if (class(z) != "biomonitoR") {
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      return("Object x is not an object of class biomonitoR")
    }
    
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
  return( round(sim,3))
}
