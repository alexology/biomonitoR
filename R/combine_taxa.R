#' combine_taxa
#'
#' This function returns all the combinations of n taxa at the desired taxonomic resolution.
#' @param x Result of the function `aggregate_taxa()`.
#' @param ntaxa Rumber of taxa.
#' @param tax_lev Taxonomic level on which the calculation has to be made.
#' @details This function is intended to help the user to identify the best subset of taxa that correlates with environmental variables.
#' Metrics based on a specific group of taxa (EPT, 1-GOLD, etc) are currently used in biomonitoring. They rely on available information on the responses
#' of the selected taxa to target strssors. The relationship between an environmental variable and a subset of taxa could exist and not detected due to the scarce
#' availaibility of autoecological data. With more than 4 taxa the calculations should become infeasible, expecially when the number
#' of taxa in the user dataset is high.
#' @keywords aggregatoR
#' @importFrom utils combn
#' @export
#' @seealso \code{\link{aggregate_taxa}}
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' combine_taxa(data_agr, tax_lev = "Family")
combine_taxa <- function(x, ntaxa = 2, tax_lev = "Taxa") {

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  # get the data.frame at the desired taxonomic level
  DF <- x[[tax_lev]]

  if (inherits(x, "bin")) {
    DF <- to_bin(DF)
  }

  # remove unassigned row from the species count if present
  if ("unassigned" %in% DF[, 1]) {
    z <- which(DF[, 1] == "unassigned")
    DF <- DF[-z, ] # remove unassigned row from the species count
  }

  # list all the combination of x taxa taken n at time
  cbn <- combn(nrow(DF), ntaxa)
  DF.sub <- lapply(seq(ncol(cbn)), function(x) DF[cbn[, x], ])

  # sum the abundances of the n-th combination of DF.sub and change the label
  DF.agg <- lapply(DF.sub, FUN = agg_fun)

  # create a data.frame to store the results
  return(do.call(rbind, DF.agg))
}
