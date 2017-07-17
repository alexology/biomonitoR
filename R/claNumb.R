claNumb <- function(x){
  cla <- x[["Class"]]
  if("unassigned" %in% cla[,1]){
    z <- which(cla$Class=="unassigned")
    cla <- cla[-z,] # remove unassigned row from the species count
  }
  ncla <- apply(cla[,-1], 2, FUN=function(x){length(x[x>0])})
  return(ncla)
}
