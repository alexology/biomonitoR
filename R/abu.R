abu <- function(x){
  nord <- apply(x[["Order"]][,-1], 2, FUN=sum)
  return(nord)
}
