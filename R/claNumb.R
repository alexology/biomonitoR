claNumb <- function(x){

  # check if the object d is of class "biomonitoR"

  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }

  cla <- x[["Class"]]
  if("unassigned" %in% cla[,1]){
    z <- which(cla$Class=="unassigned")
    cla <- cla[-z,] # remove unassigned row from the species count
  }
  ncla <- apply(cla[,-1, drop = F], 2, FUN=function(x){length(x[x>0])})
  return(ncla)
}
