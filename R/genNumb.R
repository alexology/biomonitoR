genNumb <- function(x){
  
  # check if the object d is of class "biomonitoR"
  
  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }
  
  
  gen <- x[["Genus"]]
  if("unassigned" %in% gen[,1]){
    z <- which(gen$Genus=="unassigned")
    gen <- gen[-z,] # remove unassigned row from the species count
  }
  ngen <- apply(gen[,-1], 2, FUN=function(x){length(x[x>0])})
  return(ngen)
}
