#' showscores
#'
#' This function print the scores used for the calculation of several indices in biomonitoR.
#' @param x name of the scores to print. Allowed names are "aspt_b", "aspt_i", "aspt_h", "life_scores_fam", "whpt_scores_fam", "psi_scores_fam" and "psi_scores_fam"
#' @param writecsv if TRUE the scores are saved in the working directory.
#' @keywords asBiomonitor
#' @importFrom utils write.csv
#' @export
#' @seealso \code{\link{aspt}}, \code{\link{bmwp}}, \code{\link{life}}, \code{\link{whpt}}
#' @examples
#' showscores("aspt_i", writecsv = FALSE)

showscores <- function(x, writecsv = FALSE){
  # list of allowed scores
  score.list <- c("aspt_b", "aspt_i", "aspt_h", "life_scores_fam", "whpt_scores_fam", "psi_scores_fam", "epsi_scores_fam")
  if(!x %in% score.list){
    stop("provide a valid name")
  }

  else {
    score.obj <- get(x)
    if(writecsv == FALSE){
      return(score.obj)
    } else {
      write.csv(score.obj, paste(x, ".csv", sep =""))
      return(score.obj)
    }
  }
}
