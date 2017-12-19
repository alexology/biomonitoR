#' quickRename
#'
#' This function allow the user to change the wron taxa names.
#' @param x a data.frame as specified in details
#' @keywords asBiomonitor
#' @details data.frame must have a column called "Taxa" where put species, genus or family names. See data(macro_ex) for an example dataset.\cr
#' The function rename will suggest correct name and allow the user to insert a name (Enter taxon name).
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' macro_ex.mod <- quickRename(macro_ex)



quickRename <- function(x, write.table=F){
  dfName <- deparse(substitute(x))
  tempNames <- suggestUserNames(x)
  if(length(tempNames) == 0){
    message("All names are correct")
    return(x)
  }
  else{
  wrong <- tempNames$wrongNames
  correct <- tempNames$correctNames
  result <- as.character(x$Taxa)
  n <- length(wrong)
  if (length(wrong)!=length(correct)) {
    stop("pattern and replacement do not have the same length.")
  }
  for(i in 1:n){
    result[which(result == wrong[i])] <- correct[i]
  }
  x$Taxa <- result
  if(write.table==T){
    dfNameChange <- paste0(dfName, "_mod",".txt")
    write.table(x, file=dfNameChange, quote=F, append=F)
  }
  x$Taxa <- factor(sub("_", " ", x$Taxa))
  return(x)
  }
}
