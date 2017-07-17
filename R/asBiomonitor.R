#' asBiomonitor
#'
#' This function prepares data for further calculations.
#' @param x a data.frame as specified in details
#' @keywords asBiomonitor
#' @details data.frame must have a column called "Taxa" where put species, genus or family names. See data(macro_ex) for an example dataset.\cr
#' asBiomonitor check the correctness of taxa names in the data.frame provided by the user. If names are correct the function will process the data.frame to a biomonitor object, else it suggest to rename the wrong names with the function \code{\link{rename}}.
#' @export
#' @seealso \code{\link{rename}}
#' @examples
#' data(macro_ex)
#' asBiomonitor(macro_ex)
#' # Select the correct name "Ametropus"
#' macro_ex.mod <- rename(macro_ex)
#' asBiomonitoring(macro_ex.mod)



asBiomonitor <- function(x){
  if(checkNames(x)==TRUE){
    x <- aggregate(. ~ Taxa, x, FUN=sum)
    userTaxa <- x$Taxa
    userTaxaCap <- sapply(userTaxa,capWords,USE.NAMES=F)
    x$Taxa <- userTaxaCap
    temp <- merge(ref, x, by="Taxa", all=F)
  }
  else {
    return("Wrong taxa name are present: use rename function to correct the names")
  }
  class(temp) <- "biomonitoR"
  return(temp)
}
