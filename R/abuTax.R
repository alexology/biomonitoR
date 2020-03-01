#' abuTax
#'
#' This function calculates the absolute or relative abundance of a Taxon or of a set Taxa.
#' @param x results of function aggregatoR.
#' @param taxa a Taxon or a vector of taxa.
#' @param rel if TRUE calculates relative abundance. Default to FALSE.
#' @keywords aggregatoR
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' abuTax(data.agR, taxa = "Ephemeroptera")
#' abuTax(data.agR, taxa = c("Setodes", "Orthocladiinae"), rel = TRUE)

abuTax <- function (x, taxa = NULL, rel = FALSE) 
{
  # classCheck(x)
  if (is.null(taxa) == T || (any(taxa == "") & length(taxa) == 
                             1)) {
    stop("Please provide a taxon name")
  }
  taxa <- as.character( taxa )
  df <- x[["Tree"]][, 1:10 ]
  df.vec <- unique( as.character( unlist( df ) ) )
  df.vec <- df.vec[ ! df.vec  %in% ""  ]
  x.name <- deparse( substitute( x ) )
  
  # find taxa not present in the dataset provided by the user
  taxa.sub <- taxa[ ! taxa %in% df.vec ]
  if( length( taxa.sub ) > 0 ){ 
    print( paste( "The following taxon were not find in the ", x.name ," database and has been excluded: ", taxa.sub , sep = "" ) )
    taxa <- taxa[ taxa %in% df.vec ]
    if( length( taxa ) == 0 ){
      stop( "None of the taxa provided were find in the ", x.name ," database" )
    }
  }
  
  
  taxind <- data.frame(row = numeric(), col = numeric())
  for (i in 1:length(taxa)) {
    temp <- which(df == taxa[i], arr.ind = T)
    taxind <- rbind(temp, taxind)
  }
  taxcol <- unique(taxind[, "col"])
  taxgroup <- names(df)[taxcol]
  taxgroup <- unique(taxgroup)
  for (i in 1:length(taxa)) {
    ctrl <- which(df == taxa[i], arr.ind = T)
    if (nrow(ctrl) == 0) {
      stop("Please provide a valid taxon name. Names provided are probably not present in your database")
    }
  }
  taxsub <- x[taxgroup]
  for (i in 1:length(taxsub)) {
    colnames(taxsub[[i]])[1] <- "Taxon"
  }
  taxsub <- do.call("rbind", taxsub)
  rownames(taxsub) <- NULL
  abucum <- apply(taxsub[which(taxsub$Taxon %in% taxa), -1], 
                  2, sum)
  if (rel == TRUE) {
    abucum <- abucum/abu(x)
  }
  abuperc.v <- as.vector(t(abucum))
  names(abuperc.v) <- names(abucum)
  return(abuperc.v)
}
