famNumb <- function(x){

  # check if the object d is of class "biomonitoR"

  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }


  fam <- x[["Family"]]
  if("unassigned" %in% fam[,1]){
    z <- which(fam$Family=="unassigned")
    fam <- fam[-z,] # remove unassigned row from the species count
  }
  nfam <- apply(fam[,-1, drop = F], 2, FUN=function(x){length(x[x>0])})
  return(nfam)
}
