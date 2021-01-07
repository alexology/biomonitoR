#' aggregatoR - deprecated
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("deprecated") }
#'
#' This function prepares data for further calculations.
#'
#' @param x results of function asBiomonitor
#' @param FUN the function to be applied for aggregating to higher taxonomic levels.
#' Must be sum for both abundances and presence-absence data.
#' Default to `sum`. Change if you know what you are doing.
#' @keywords aggregatoR
#' @importFrom stats aggregate
#' @export
#' @seealso \code{\link{aggregate_taxa}}


aggregatoR <- function(x, FUN = sum) {
  .Deprecated("aggregate_taxa")

  # check if the object x is of class "biomonitoR"
  if (!inherits(x, "asb")) stop("x is not an object created with asBiomonitor")

  class_x <- class(x)
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  # The following part is to aggregate the dataset provided by the user at different taxonomic levels
  tx <- c("Phylum", "Class", "Subclass", "Order", "Family", "Subfamily", "Tribus", "Genus", "Species", "Subspecies", "Taxa")
  stz <- x[!(names(x) %in% tx)]
  stz_n <- names(stz)

  agg_list <- list()

  for (i in seq_along(tx)) {
    temp.agg <- aggregate(stz, by = list(x[, i]), FUN)
    temp.agg[temp.agg$Group.1 == "", "Group.1"] <- "unassigned"
    names(temp.agg) <- c(tx[i], stz_n)
    agg_list[[i]] <- temp.agg
  }

  names(agg_list) <- tx

  agg_list$Tree <- x

  class_x[class_x %in% "asb"] <- "biomonitoR"
  class(agg_list) <- class_x
  agg_list
}
