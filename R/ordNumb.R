ordNumb <- function(x){
  ord <- x[["Order"]]
  if("unassigned" %in% ord[,1]){
    z <- which(ord$Order=="unassigned")
    ord <- ord[-z,] # remove unassigned row from the species count
  }
  nord <- apply(ord[,-1], 2, FUN=function(x){length(x[x>0])})
  return(nord)
}
