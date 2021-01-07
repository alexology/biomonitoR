#' get_taxa_richness
#'
#' @description
#'
#' This function calculates the absolute richness of a taxon or of a set taxa at the taxonomic level provided by the user.
#'
#'
#' @param x Result of the function `aggregate_taxa()`.
#' @param taxa a taxon or a vector of taxa.
#' @param tax_lev taxonomic level at which the richness has to be calculated. It could be also a vector of taxonomic levels.
#' @keywords aggregate_taxa
#' @export
#' @seealso \code{\link{aggregate_taxa}}
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' get_taxa_richness(data_agr, taxa = "Ephemeroptera", tax_lev = "Family")
get_taxa_richness <- function(x, taxa = NULL, tax_lev = NULL) {

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  # stop if user does not provide a taxon name
  if (is.null(taxa) == TRUE | is.null(tax_lev) == TRUE) {
    stop("Please provide a taxon name and/or a taxonomic level")
  }

  # Allow the user to provide a single taxonomic level for all the selected taxa
  if (length(taxa) > 1 & length(tax_lev) == 1) {
    tax_lev <- rep(tax_lev, length(taxa))
  }

  # Allow the user to provide a single taxonomic level for all the selcted taxa
  if (length(taxa) > 1 & length(tax_lev) > 1 & length(taxa) != length(tax_lev)) {
    stop("tax_lev must be of the same length of taxa")
  }

  # remove leadine and trailing spaces, moreover the first letter is set to capital
  # while the others to lowercase
  taxa <- trimws(taxa)
  taxa <- unlist(lapply(as.list(taxa), capWords))
  tax_lev <- trimws(tax_lev)
  tax_lev <- unlist(lapply(as.list(tax_lev), capWords))

  # extract taxonomic information from the element Tree in the aggregatoR output
  DF <- x[["Tree"]][, 1:10]
  df.names <- names(DF)


  # check if the taxa provided by the user are in the taxonomic tree
  for (i in 1:length(taxa)) {
    ctrl <- which(DF == taxa[i], arr.ind = TRUE)
    if (nrow(ctrl) == 0) {
      stop("Please provide a valid taxon name. Names provided can also be absent in your database.")
    }
  }

  # Position of taxon in the df data.frame. 1:11 referes to the columns occupied by the tree
  taxind <- rep(0, ncol(x[["Tree"]][, -c(1:11)]))

  # the calculation is made for each taxa throught a for loop. Maybe better solution exists. An apply?
  for (i in 1:length(taxa)) {
    temp <- which(DF == taxa[i], arr.ind = TRUE)
    tax.lev <- tax_lev[i]

    # check if the taxa taxonomic level is lower than tax_lev taxonomic level
    if (unique(temp[, "col"]) > which(df.names == tax.lev)) {
      stop("Taxonomic level of taxa cannot be lower than taxonomic level of tax_lev")
    }

    tax.sel <- as.character(DF[temp[, "row"], tax.lev])
    df.sel <- x[[tax.lev]]
    df.sel <- df.sel[which(df.sel[, 1] %in% tax.sel), ]
    if ("unassigned" %in% df.sel[, 1]) {
      z <- which(df.sel == "unassigned")
      df.sel <- df.sel[-z, ] # remove unassigned row from the species count
    }

    # Column 1 represents Taxa names and must be excluded from the calculations

    ntax <- apply(df.sel[, -1, drop = FALSE], 2, FUN = function(x) {
      length(x[x > 0])
    })
    taxind <- taxind + ntax
  }

  taxind
}
