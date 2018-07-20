genNumb <- function(x){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  gen <- x[["Genus"]]
  if("unassigned" %in% gen[,1]){
    z <- which(gen$Genus=="unassigned")
    gen <- gen[-z,] # remove unassigned row from the species count
  }
  ngen <- apply(gen[,-1, drop = F], 2, FUN=function(x){length(x[x>0])})
  return(ngen)
}
