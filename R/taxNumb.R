taxNumb <- function(x){

  # check if the object d is of class "biomonitoR"

  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }


  tax <- x[["Taxa"]]
  if("unassigned" %in% tax[,1]){
    z <- which(tax$Taxon=="unassigned")
    tax <- tax[-z,] # remove unassigned row from the species count
  }
  tax <- apply(tax[,-1], 2, FUN=function(x){length(x[x>0])})
  return(tax)
}
