asBiomonitor <- function(x){
  if(checkNames(x)==TRUE){
    temp <- merge(ref,x, by="Taxa")
  }
  else {
    return("Wrong taxa name are present: use rename function to correct the names")
  }
  return(temp)
}
