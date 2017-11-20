# help function for aggTaxa

agg_fun <- function(x, rel = F){
  lab <- paste(x[,1], collapse = "_")
  temp <- apply(x[,-1], 2, FUN = sum)
  df <- data.frame(Taxa = lab, t(temp))
  return(df)
}