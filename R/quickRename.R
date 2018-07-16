#' quickRename
#'
#' This function allow the user to change the wron taxa names.
#' @param x a data.frame as specified in details
#' @param group otic group of interest. Possible values are "mi" for macroinvertebrates and "mf" for macrophytes.
#' @param write.table if T quickRename will save a csv file with changes provided by the user
#' @keywords asBiomonitor
#' @details data.frame must have a column called "Taxa" where put species, genus or family names. See data(macro_ex) for an example dataset.\cr
#' The function rename will suggest correct name and allow the user to insert a name (Enter taxon name).
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' macro_ex.mod <- quickRename(macro_ex, group = "mi")



quickRename <- function(x, group = NULL, write.table=FALSE){

  if(is.null(group) == TRUE){
    stop("Please provide a valid group")
  }

  dfName <- deparse(substitute(x))
  tempNames <- suggestUserNames(x, groups = group)
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
  if(write.table == TRUE){
    dfNameChange <- paste0(dfName, "_mod",".txt")
    write.table(x, file=dfNameChange, quote = FALSE, append = FALSE)
  }
  x$Taxa <- factor(sub("_", " ", x$Taxa))
  return(x)
  }
}
