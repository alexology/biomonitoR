#' remove_taxa
#'
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("experimental") }
#'
#' Remove a taxon or a set of taxa and create a new `data.frame`.
#'
#'
#'
#' @param x result of the function aggregatoR.
#' @param taxa a taxon or a vector of taxa.
#' @details This function does not check for parent-child pairs but it should handle them pretty well.
#' @keywords aggregatoR
#' @export
#' @seealso \code{\link{aggregatoR}} \code{\link{ambiguousSolver}}
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' remove_taxa(data_agr, taxa = "Ephemeroptera")
remove_taxa <- function(x, taxa = NULL) {

  # check if the object d is of class "biomonitoR"
  classCheck(x)

  # check if the taxa argument is empty or contains NULL strings
  if (is.null(taxa) == TRUE || (any(taxa == "") & length(taxa) == 1)) {
    stop("Please provide at least taxon name")
  }

  # transform the taxa argument to character
  taxa <- as.character(taxa)

  # get the list of all the taxa present in the x object including the information stored
  # in the taxonomic tree
  df <- x[["Tree"]][, 1:10]
  df.vec <- unique(as.character(unlist(df)))
  df.vec <- df.vec[!df.vec %in% ""]

  # get the name of the x object
  x.name <- deparse(substitute(x))

  # find any taxa not present in the dataset provided by the user. If no taxa belongs to the list of all the taxa present
  # in the x object the function will be stopped, otherwise it will provide the list of missing taxa
  taxa.sub <- taxa[!taxa %in% df.vec]
  if (length(taxa.sub) > 0) {
    print(paste("The following taxon were not find in the ", x.name, " database and has been excluded: ", taxa.sub, sep = ""))
    taxa <- taxa[taxa %in% df.vec]
    if (length(taxa) == 0) {
      stop("None of the taxa provided were found in the ", x.name, " database")
    }
  }

  # list the row and column numbers at which the desired taxa are stored
  taxind <- data.frame(row = numeric(), col = numeric())
  for (i in 1:length(taxa)) {
    temp <- which(df == taxa[i], arr.ind = T)
    taxind <- rbind(temp, taxind)
  }



  df_sub <- x[["Tree"]][-unique(taxind$row), -c(1:10)]
  df_sub
}
