#' aggregatoR
#'
#' Functions for calculating BMWP and ASPT
#' @param x results of function aggregatoR
#' @param method a,b or i. See details.
#' @keywords aggregatoR
#' @details
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' bmwp(data.agR)



bmwpx <- function( d , method = "a") {
  numb <- c(which(names(d)=="Tree"), which(names(d)=="Taxa")) # position of the Tree element in the list to remove
  x <- d[-numb]
  # y is the reference data.set for bmwp calculation
  st.names <- names(x[[1]][-1]) # names of sampled sites
    if(method == "a") (y <- aspt_h)
  if(method == "b") {y <- aspt_b
  z <- bfam_acc}
  if(method == "i") {y <- aspt_i
  z <- ifam_acc}
  if(method == "b" || method == "i") (x <- checkBmwpFam(df=x, famNames=z, stNames=st.names))

  for(i in 1:length(x)){
    colnames(x[[i]])[1] <- "Taxon"
  }
  
  df <- do.call( "rbind" , x )
  rownames( df ) <- NULL
  df <- data.frame( df[ , 1 , drop =F ], (df[ , -1 ] > 0 ) * 1 )
  tot.mer <- merge( y , df )
  
  # check if merge results provided valid data.frame
  if( nrow(tot.mer) == 0 ){
    opt <- options( show.error.messages = T )
    on.exit( options( opt ) )
    return("No valid taxon provided")
  }
  else {
    tot.st <- which(names(tot.mer)%in%st.names)
    tot.bmwp <- apply(tot.mer$Value*tot.mer[ , tot.st, drop=F], 2, sum)
  }
  return( tot.bmwp )
}
