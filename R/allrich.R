#' @name allrich
#'
#' @title Function for calculating richness at different taxonomic levels
#'
#' @description
#' This function allows the calculation of richness at different taxonomic levels and of the total abundance.
#' A convenient function for calculating richness at all the taxonomic levels simultaneously is provided.
#'
#' @aliases  richness
#'
#' @param x Result of `aggregate_taxa()`.
#' @param tax_lev Taxonomic level at which the calculation has to be performed.
#'
#' @details `richness()` returns a warning meassage if unassigned taxa are detected at the desired taxonomic level. Unassigned taxa can be present
#' because some taxa are entered at higher taxonomic levels than those at which the calculation is performed.
#' For instance, this problem appears when the calculation is performed at genus level and some taxa has been entered at family level.
#' Unassigned taxa can be present also because of missing levels in the reference database at the desired taxonomic level.
#' This problem is present, for example, for Mollusca in the default database for macroinvertebrates.
#' Mollusca orders are not present because there is no consensus among taxonomist (see the [www.molluscabase.org](https://www.molluscabase.org/index.php) website).
#'
#' @keywords allrich, richness
#'
#' @export
#'
#' @export richness
#'
#' @seealso [aggregate_taxa]
#'
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' richness(data_agr, tax_lev = "Family")
allrich <- function(x) {

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  # extract the names of all the taxonomic levels on which perform the calculation
  n.tree <- names(x[["Tree"]][, sapply(x[["Tree"]], is.character)])

  # initialize the data.frame to store the results
  res <- data.frame(matrix(NA, ncol = length(n.tree), nrow = sum(!names(x[["Tree"]]) %in% n.tree)))

  # a for loop to calculate richness at different taxonomic levels
  # while filling the res data.frame
  for (i in 1:length(n.tree)) {
    res[, i] <- richness(x, tax_lev = n.tree[i])
  }

  # assign column and row names and return the results
  colnames(res) <- n.tree
  rownames(res) <- colnames(x[["Tree"]])[!colnames(x[["Tree"]]) %in% n.tree]
  res
}
