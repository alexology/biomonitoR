#' aggregatoR
#'
#' This function calculates the Average Score Per Taxon following Armitage et al. (1983), Davy-Bowker et al. (2007) and Alba-Tercedor & Sanchez-Ortega (1988) formulations.
#' @param x results of function aggregatoR
#' @param method the formulation of BMWP needed to calculate ASPT. Possible choises are "a" (Armitage et al. 1983), "b" (Davy-Bowker et al. 2007) and i (Alba-Tercedor & Sanchez Ortega, 1988)  
#' @keywords aggregatoR
#' @details
#' @references ALBA-TERCEDOR,  J.  &  A.  SÁNCHEZ-ORTEGA. 1988.  Un  método  rápido  y  simple  para  evaluar  la calidad biológica de las aguas corrientes basado en el de Hellawell (1978). Limnetica, 4: 51-56.
#' @references Armitage, P. D., Moss, D., Wright, J. F., & Furse, M. T. (1983). The performance of a new biological water quality score system based on macroinvertebrates over a wide range of unpolluted running-water sites. Water research, 17(3), 333-347.
#' @references Davy-Bowker J., Clarke R., Corbin T., Vincent H, Pretty J., Hawczak A., Blackburn J., Murphy J., Jones I., 2008. River Invertebrate Classification Tool. Final report. WFD72C. SNIFFER. 276 pp
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' aspt(data.agR)

aspt <- function( d , method = "a") {
  
  # check if the object d is of class "biomonitoR"
  
    if (class(d) != "biomonitoR") {
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      return("Object x is not an object of class biomonitoR")
    }
  
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
    ntaxa <- colSums(tot.mer[,-c(1:2)] == 1)
    tot.aspt <- apply(tot.mer$Value*tot.mer[ , tot.st, drop=F], 2, sum)/ntaxa
  }
  return( round(tot.aspt, 3) )
}