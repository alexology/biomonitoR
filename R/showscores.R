#' showscores
#'
#' This function print the scores used for the calculation of several indices in biomonitoR.
#' @param x name of the scores to print. Allowed names are "aspt_b", "aspt_i", "aspt_h", "life_scores_spe" and "life_scores_fam"
#' @param writecsv if TRUE the scores are saved in the working directory.
#' @keywords asBiomonitor
#' @export
#' @seealso \code{\link{aspt}}, \code{\link{bmwp}}, \code{\link{life}}
#' @examples
#' showscores("aspt_i", writecsv = F)

showscores <- function(x, writecsv = F){
  # list of allowed scores
  score.list <- c("aspt_b", "aspt_i", "aspt_h", "life_scores_spe", "life_scores_fam")
  if(!x %in% score.list){
    stop("provide a valid name")
  }
  
  else {
    score.obj <- get(x)
    if(writecsv == F){
      return(score.obj)
    } else {
      write.csv(score.obj, paste(x, ".csv", sep =""))
      return(score.obj)
    }
  }
}