pielou <- function(x, base=2, taxLev="Family"){

  # check if the object d is of class "biomonitoR"

  if (class(x) != "biomonitoR") {
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

  rich <- apply(df[,-1, drop = F], 2, FUN=function(x){length(x[x>0])})

  if(ncol(df)==2){
    pie <- Pi(df[,-1], index = "Shannon", base = base) / log(rich, base = base)
  }
  else{
    pie <- apply(df[,-1], 2, FUN = function(x){Pi( x, index = "Shannon", base = base )}) / log(rich, base = base)
  }
  return( pie )
}
