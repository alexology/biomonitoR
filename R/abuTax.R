#' abuTax
#'
#' This function calculates the absolute or relative abundance of a Taxon or of a set Taxa. 
#' @param x results of function aggregatoR.
#' @param taxa a Taxon or a vector of taxa.
#' @param rel if TRUE calculates relative abundance. default =F.
#' @keywords aggregatoR
#' @details
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' abuTax(data.agr, taxa = "Ephemeroptera")
#' abuTax(data.agr, taxa = c("Plecopotera", "Ephemeroptera"))

abuTax <- function(x, taxa = NULL, rel = FALSE){
  
  # check if the object d is of class "biomonitoR"
  
  if (class(d) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }
  
  # stop if user does not provide a taxon name
  if(is.null(taxa) == T || taxa == ""){
    stop("Please provide a taxon name")
  }


  # extract taxonomic information from the element Tree in aggregatoR output
  df <- x[["Tree"]][,1:10]
  
  # Position of taxon in the df data.frame
  taxind <- data.frame(row = numeric(), col = numeric())
  for(i in 1:length(taxa)){
    temp <- which(df == taxa[i], arr.ind=T)
    taxind <- rbind(temp, taxind)  
  }
  
  
  taxcol <- unique(taxind[,"col"])
  taxgroup <- names(df)[taxcol]
  taxgroup <- unique(taxgroup)
  
  # check if the data.frame contains the user provided taxon name

  for(i in 1:length(taxa)){
    ctrl <- which(df == taxa[i], arr.ind=T)
    if(nrow(ctrl) == 0){
      stop("Please provide a valid taxon name")
    }
  }
  

  taxsub <- x[taxgroup]

  for(i in 1:length(taxsub)){
    colnames(taxsub[[i]])[1] <- "Taxon"
  }

  taxsub <- do.call( "rbind" , taxsub )
  rownames( taxsub ) <- NULL
  
  abucum <- apply(taxsub[which(taxsub$Taxon %in% taxa), - 1], 2 , sum)
  if(rel == TRUE){
    abucum <- abucum/abu(x)
  }
  abuperc.v <- as.vector(t(abucum))
  names(abuperc.v) <- names(abucum)
  return(round(abuperc.v, 3))
}
