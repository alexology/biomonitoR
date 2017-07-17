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
#' aspt()
#' @export aspt


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
  fam.mer <- merge(y[["Family"]],x.bin[["Family"]])

  # merge order
  ord.mer <- merge(y[["Order"]],x.bin[["Order"]])

  # merge subclasses
  cla.mer <- merge(y[["Class"]],x.bin[["Class"]])

  # check if merge results provided valid data.frame
  if(nrow(fam.mer)==0&&nrow(ord.mer)==0&&nrow(scla.mer)==0){
    opt <- options(show.error.messages=T)
    on.exit(options(opt))
    return("No valid taxon provided")
  }
  else {
    fam.st <- which(names(fam.mer)%in%st.names)
    ord.st <- which(names(ord.mer)%in%st.names)
    cla.st <- which(names(cla.mer)%in%st.names)
    fam.bmwp <- apply(fam.mer$Value*fam.mer[ , fam.st, drop=F], 2, sum)
    ord.bmwp <- apply(ord.mer$Value*ord.mer[ , ord.st, drop=F], 2, sum)
    cla.bmwp <- apply(cla.mer$Value*cla.mer[ , cla.st, drop=F], 2, sum)
    bmwp.res <- fam.bmwp+ord.bmwp+cla.bmwp
  }
  return(bmwp.res)
}
