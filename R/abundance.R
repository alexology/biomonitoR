#' abundance
#'
#' @description Calculates taxa abundance at the desired taxonomic level.
#' Unassigned taxa can be included or excluded from the calculation.
#'
#' @param x Result of `aggregate_taxa()`.
#' @param tax_lev Taxonomic level at which the calculation has to be performed.
#' @param unassigned Does unassigned taxa need to be taken into account? If yes set unassigned to `TRUE`, otherwise `FALSE`.
#'
#' @keywords abundance
#'
#' @details If unassigned is set to `TRUE` abundances are equal among taxonomic levels.
#'
#' @export
#'
#' @seealso [aggregate_taxa] [get_taxa_abundance]
#'
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' abundance(data_agr, tax_lev = "Family")
#' abundance(data_agr, tax_lev = "Family", unassigned = TRUE)
abundance <- function(x, tax_lev = "Taxa", unassigned = FALSE) {

  # check if the object d is of class "biomonitoR"
  classCheck(x)

  # get the desired taxomic level
  tax <- x[[tax_lev]]

  if (inherits(x, "bin")) {
    tax <- to_bin(tax)
  }

  if (!unassigned) {
    if ("unassigned" %in% tax[, 1]) {
      z <- which(tax[, 1] == "unassigned") # find the row corresponding to the unassigned
      tax <- tax[-z, ] # remove unassigned row from the species count
    }
  }

  # calculate abundance from the desired taxonomic level. The -1 is to remove the Taxa column
  res <- apply(tax[, -1, drop = FALSE], 2, FUN = sum)

  res
}
