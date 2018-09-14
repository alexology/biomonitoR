pielou <- function(x, base = 2, taxLev="Family"){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

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

  pie <- apply(df[ ,-1, drop = F ], 2, FUN = function(x){Pi( x, index = "Shannon", base = base )}) / log(rich, base = base)
  
  return( pie )
}
