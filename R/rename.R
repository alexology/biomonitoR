rename <- function(x, write.table = FALSE, customx = FALSE, group = "none"){

  dfName <- deparse(substitute(x))
  if(customx == F){
    tempNames <- suggestNames(x, group = group)
  } else {
    tempNames <- suggestNames(x, custom = T, group = group)
  }

  if(length(tempNames) == 0){
    message("All names are correct")
    return(x)
  }
  else{
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
  if(write.table == TRUE){
    dfNameChange <- paste0(dfName, "_mod",".txt")
    write.table(x, file=dfNameChange, quote=FALSE, append=FALSE)
  }
  x$Taxa <- sub("_", " ", x$Taxa)
  return(x)
  }
}
