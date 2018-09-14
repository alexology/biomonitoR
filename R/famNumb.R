famNumb <- function(x){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  fam <- x[["Family"]]
  if("unassigned" %in% fam[,1]){
    z <- which(fam$Family=="unassigned")
    fam <- fam[-z,] # remove unassigned row from the species count
  }
  nfam <- apply(fam[,-1, drop = F], 2, FUN=function(x){length(x[x>0])})
  return(nfam)
}
