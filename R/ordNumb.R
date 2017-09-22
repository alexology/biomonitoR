ordNumb <- function(x){
  
  # check if the object d is of class "biomonitoR"
  
  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }
  
  ord <- x[["Order"]]
  if("unassigned" %in% ord[,1]){
    z <- which(ord$Order=="unassigned")
    ord <- ord[-z,] # remove unassigned row from the species count
  }
  nord <- apply(ord[,-1], 2, FUN=function(x){length(x[x>0])})
  return(nord)
}
