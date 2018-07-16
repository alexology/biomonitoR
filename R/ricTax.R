#' ricTax
#'
#' This function calculates the absolute richness of a Taxon or of a set Taxa at a user provided taxonomic level.
#' @param x results of function aggregatoR.
#' @param taxa a Taxon or a vector of taxa.
#' @param taxLev taxonomic level at which the richness has to be calculated. It could be also a vector of taxnomici levels.
#' @keywords aggregatoR
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' ricTax(data.agR, taxa = "Ephemeroptera", taxLev = "Family")

ricTax <-  function(x , taxa = NULL, taxLev = NULL){

  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }

  # stop if user does not provide a taxon name
  if(is.null(taxa) == TRUE || taxa == "" || is.null(taxLev) == TRUE || taxLev == ""){
    stop("Please provide a taxon name and/or a taxonomic level")
  }

  # Allow the user to provide a single taxonomic level for all the slected taxa
  if(length(taxa) > 1 & length(taxLev) == 1){
    taxLev <- rep(taxLev, length(taxa))
  }

  # Allow the user to provide a single taxonomic level for all the slected taxa
  if(length(taxa) > 1 & length(taxLev) > 1 & length(taxa) != length(taxLev)){
    stop("taxLev must be of the same length of taxa")
  }

  # extract taxonomic information from the element Tree in aggregatoR output
  df <- x[["Tree"]][,1:10]
  df.names <- names(df)



  for(i in 1:length(taxa)){
    ctrl <- which(df == taxa[i], arr.ind = TRUE)
    if(nrow(ctrl) == 0){
      stop("Please provide a valid taxon name. Names provided can also be absent in your database.")
    }
  }

  # Position of taxon in the df data.frame
  taxind <- rep(0, ncol(x[["Tree"]][ , -c(1:11) ]))

  for(i in 1:length(taxa)){
    temp <- which(df == taxa[i], arr.ind = TRUE)
    tax.lev <- taxLev[i]
    # check if taxa taxonomic level is lower than taxLev taxonomic level
    if(unique(temp[ , "col"]) > which(df.names == tax.lev)){
      stop("Taxonomic level of taxa cannot be lower than taxonomic level of taxLev")
    }

    tax.sel <- as.character(df[temp[, "row"] , tax.lev])
    df.sel <- x[[tax.lev]]
    df.sel <- df.sel[which(df.sel[,1] %in% tax.sel),]
    if("unassigned" %in% df.sel[,1]){
      z <- which(df.sel == "unassigned")
      df.sel <- df.sel[-z,] # remove unassigned row from the species count
    }
    # Column 1 represents Taxa names and it need to be excluded from the calculations

    ntax <- apply(df.sel[, -1, drop = FALSE], 2, FUN = function(x) { length( x[x>0] ) } )
    taxind <- taxind + ntax
  }

 return(taxind)
}
