#' ept
#'
#' This function calculates the richness of Ephemeroptera, Plecotera and Trichoptera (EPT) taxa at the desired taxonomic level.
#' @param x Result of `aggregate_taxa()`.
#' @param tax_lev The taxonomic level for calculating EPT richness.
#' @keywords ept
#' @details `tax_lev` must be finer than Order (e.g. Species, Genus, Family ).
#' The `ept()` function works if Ephemeroptera, Trichoptera and Plecoptera are set as orders in the reference database. Otherwise try with `aggregate_richness()`.
#' @importFrom stats aggregate
#' @export
#' @seealso [aggregate_taxa] [get_taxa_richness]
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' ept(data_agr)
#'
#' # ept should return the same results as the function get_taxa_richness
#' ept(data_agr, tax_lev = "Family")
#' get_taxa_richness(data_agr,
#'   taxa = c("Ephemeroptera", "Plecoptera", "Trichoptera"), tax_lev = "Family"
#' )
ept <- function(x, tax_lev = "Taxa") {

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  # get the tree from the aggregatoR object
  x_ept <- x[["Tree"]]

  # list the taxonomic levels in the Tree
  tx <- c("Phylum", "Class", "Subclass", "Order", "Family", "Subfamily", "Tribus", "Genus", "Species", "Subspecies", "Taxa")

  # stop the calculation if tax_lev is coarser than Family
  if (which(tax_lev == tx) <= 4) (stop("Taxonomic level cannot be equal or higher than Order"))

  # get sample names
  stz <- x_ept[!(names(x_ept) %in% tx)]
  stz_n <- names(stz)

  # subset the EPT orders from the Tree
  ept_taxa <- x_ept[which(x_ept$Order == "Plecoptera" | x_ept$Order == "Ephemeroptera" | x_ept$Order == "Trichoptera"), , drop = F]

  # if EPT are absent from the database returns 0 otherwise return the results
  if (nrow(ept_taxa) == 0) {
    res_null <- rep(0, length(stz_n))
    names(res_null) <- stz_n
    return(res_null)
  } else {
    ept_temp <- ept_taxa[, c(tax_lev, stz_n)]
    colnames(ept_temp)[1] <- "selection"
    levels(ept_temp$selection)[levels(ept_temp$selection) == ""] <- "unassigned"
    ept_temp.agg <- aggregate(. ~ selection, ept_temp, FUN = sum)
    if ("unassigned" %in% ept_temp.agg[, 1]) {
      z <- which(ept_temp.agg[, 1] == "unassigned")
      ept_temp.agg <- ept_temp.agg[-z, ] # remove unassigned row from the species count
    }
    temp <- colSums(ept_temp.agg[, -1, drop = F] > 0)
    temp
  }
}
