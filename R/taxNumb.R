taxNumb <- function(x){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  tax <- x[["Taxa"]]
  if("unassigned" %in% tax[,1]){
    z <- which(tax$Taxon=="unassigned")
    tax <- tax[-z,] # remove unassigned row from the species count
  }
  tax <- apply(tax[,-1, drop = F], 2, FUN=function(x){length(x[x>0])})
  return(tax)
}
