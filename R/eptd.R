#' eptd
#'
#' This function calculates the log10(Sel_EPTD + 1) index, where EPTD stands for Ephemeroptera, Plecoptera, Trichoptera and Diptera.
#' @param x Result of `aggregate_taxa()`.
#' @param base Base of the logarithm. Default to 10.
#' @param eptd_families A vector of selected EPTD taxa at family level. Optional, default to `NULL`.
#' @param traceB If set to `TRUE` a list as specified below will be returned.
#' @keywords eptd
#' @details log10(Sel_EPTD + 1) is the base-10 logarithm of the abundance of the selected EPTD families plus 1. Selected EPTD families are Heptageniidae, Ephemeridae, Leptophlebiidae, Brachycentridae, Goeridae, Polycentropodidae, Limnephilidae, Odontoceridae, Dolichopodidae, Stratiomyidae, Dixidae, Empididae, Athericidae and Nemouridae.
#' This metric is part of the italian STAR_ICMi index, where it is supposed to be relate to habitat integrity. To accomodate for other EPTD taxa,
#' users can provide their own list of taxa as a character vector. \cr
#' The function `eptd()` will search for the selected EPTD within Ephemeroptera, Plecoptera, Trichoptera and Diptera set at Order level.
#' If store Ephemeroptera, Plecoptera, Trichoptera and Diptera are stored at taxonomic level other than Order, this function will
#' not work properly, see `get_taxa_abundance()` instead.
#'
#' @return If `traceB` is set to `TRUE` a list with the following elements will be returned:
#' \itemize{
#'  \item `results` Results of the log10(Sel_EPTD + 1) index.
#'  \item `taxa_df` The `data.frame` used for the calculation containing the abundance of the EPTD taxa.
#' }
#'
#' @export
#' @seealso [aggregate_taxa] [get_taxa_abundance]
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' eptd(data_agr)
eptd <- function(x, base = 10, eptd_families = NULL, traceB = FALSE) {

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  if (is.null(eptd_families)) {
    eptd_fam <- c("Heptageniidae", "Ephemeridae", "Leptophlebiidae", "Brachycentridae", "Goeridae", "Polycentropodidae", "Limnephilidae", "Odontoceridae", "Dolichopodidae", "Stratiomyidae", "Dixidae", "Empididae", "Athericidae", "Nemouridae")
  } else {
    eptd_fam <- trimws(sapply(as.character(eptd_families), capWords, USE.NAMES = FALSE))
  }

  x_fam <- x[["Family"]]

  if (inherits(x, "bin")) {
    x_fam <- to_bin(x_fam)
  }

  x_eptd <- x_fam[which(x_fam$Family %in% eptd_fam), , drop = F]
  temp <- log(apply(x_eptd[, -1, drop = FALSE], 2, sum) + 1, base = base)

  if (!traceB) {
    temp
  } else {
    rownames(x_eptd) <- NULL
    list(results = temp, taxa_df = x_eptd)
  }
}
