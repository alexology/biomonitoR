#' aggregate_taxa
#'
#' @description
#' This function prepares data for further calculations.
#'
#' @param x results of function asBiomonitor
#' @param FUN the function to be applied for aggregating to higher taxonomic levels.
#' Must be sum for both abundances and presence-absence data.
#' Default to `sum`.
#'
#' @keywords aggregatoR
#'
#' @importFrom stats aggregate
#'
#' @export
#'
#' @seealso \code{\link{asBiomonitor}}
#'
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#'
#' # example for macrophytes
#' data(oglio)
#'
#' oglio_asb <- as_biomonitor(oglio, group = "mf")
#' oglio_agg <- aggregate_taxa(oglio_asb)
#' richness(oglio_agg, taxLev = "Species")
#' richness(oglio_agg, taxLev = "Genus")
#' richness(oglio_agg, taxLev = "Family")


aggregate_taxa <- function(x, FUN = sum) {

  # check if the object x is of class "biomonitoR"
  if (!inherits(x, "asb")) stop("x is not an object created with asBiomonitor")

  class_x <- class(x)
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  # The following part is to aggregate the dataset provided by the user at different taxonomic levels
  tx <- c("Phylum", "Class", "Subclass", "Order", "Family", "Subfamily", "Tribus", "Genus", "Species", "Subspecies", "Taxa")
  stz <- x[!(names(x) %in% tx)]
  stz_n <- names(stz)

  agg_list <- vector(mode = "list", length = 12)
  names(agg_list) <- c(tx, "Tree")

  for (i in seq_along(tx)) {
    temp.agg <- aggregate(stz, by = list(x[, i]), FUN)
    temp.agg[temp.agg$Group.1 == "", "Group.1"] <- "unassigned"
    names(temp.agg) <- c(tx[i], stz_n)
    agg_list[[i]] <- temp.agg
  }

  agg_list$Tree <- x

  class_x[class_x %in% "asb"] <- "biomonitoR"
  class(agg_list) <- class_x
  agg_list
}
