rename <- function(x, write.table=F, custom = F){
  dfName <- deparse(substitute(x))
  if(custom == F){
    tempNames <- suggestNames(x)
  } else {
    tempNames <- suggestNames(x, custom = T)
  }
  
  wrong <- tempNames$wrongNames
  correct <- tempNames$correctNames
  result <- x$Taxa
  n <- length(wrong)
  if (length(wrong)!=length(correct)) {
    stop("pattern and replacement do not have the same length.")
  }
  for(i in 1:n){
    result <- gsub(wrong[i],correct[i],result)
  }
  x$Taxa <- result
  if(write.table==T){
    dfNameChange <- paste0(dfName, "_mod",".txt")
    write.table(x, file=dfNameChange, quote=F, append=F)
  }
  x$Taxa <- sub("_", " ", x$Taxa)
  return(x)
}
