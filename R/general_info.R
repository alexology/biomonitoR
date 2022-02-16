#' general_info
#'
#' Calculates the overall richness at different taxonomic level.
#' @param x Result of `aggregate_taxa()`.
#' @param taxalist If `TRUE` returns the list of taxa for each taxonomic level.
#' @keywords general_info
#' @export
#' @seealso [as_biomonitor]
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' general_info(data_agr)
general_info <- function(x, taxalist = FALSE) {

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  x <- x[["Tree"]]
  tx <- c("Phylum", "Class", "Subclass", "Order", "Family", "Subfamily", "Tribus", "Genus", "Species", "Subspecies", "Taxa")
  stz1 <- x[!(names(x) %in% tx)]
  stz <- apply(stz1, 1, sum)
  stz[stz > 0] <- 1
  stz <- as.logical(stz)
  treeo <- x[stz, names(x) %in% tx]
  tlist <- apply(treeo, 2, unique)
  if (taxalist) {
    lapply(tlist, function(x) x[x != ""])
  }
  else {
    if (any(stz1 > 1)) {
      c(unlist(lapply(tlist, function(x) (sum(x != "")))), Abundance = sum(stz1))
    } else {
      unlist(lapply(tlist, function(x) (sum(x != ""))))
    }
  }
}
