#' totInfo
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("deprecated") }
#'
#' Function to calculate the overall richness at different taxonomic level.
#' @param x results of the function [aggregatoR]
#' @param taxalist if `TRUE` returns the list of taxa for each taxonomic level.
#' @keywords totInfo
#' @export
#' @seealso \code{\link{asBiomonitor}}



totInfo <- function(x, taxalist = FALSE) {
  .Deprecated("general_info", package = "biomonitoR")

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
