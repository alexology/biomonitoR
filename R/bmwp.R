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



bmwp <- function( x , method = "a") {
  # y is the reference data.set for bmwp calculation
  st.names <- names(x[[1]][-1]) # names of sampled sites
  if(method == "a") (y <- aspt_h)
  if(method == "b") {y <- aspt_b
  z <- bfam_acc}
  if(method == "i") {y <- aspt_i
  z <- ifam_acc}
  if(method == "b" || method == "i") (checkBmwpFam(df=x, famNames=z, stNames=st.names))

  x.bin <- lapply(x, function(x){data.frame( x[,1,drop=F], (x[,-1]>0)*1)})

  # merge families
  fam.mer <- merge( y[["Family"]], x.bin[["Family"]] )
  colnames(fam.mer)[1] <- "Taxa"

  # merge order
  ord.mer <- merge( y[["Order"]], x.bin[["Order"]] )
  colnames(ord.mer)[1] <- "Taxa"

  # merge subclasses
  cla.mer <- merge( y[["Class"]], x.bin[["Class"]] )
  colnames(cla.mer)[1] <- "Taxa"

  # rbind merges
  tot.mer <- rbind(fam.mer, ord.mer, cla.mer)

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
  return(tot.bmwp)
}
